"""Strict inventory for registration-only post-Xenium IF processing.

The biological section registry intentionally covers only sections whose antibody roles were
recorded by the maintainer.  Registration needs fewer facts: an image, its exact Xenium serial,
and the pancreas ROI.  Those facts live in ``restain_registration_inputs.csv``; calibration and
channel names come from the qptiff itself.  Keeping the two contracts separate prevents a batch
registration from inventing INS/GCG identities for technical Opal channel names.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from phenocycler.cohort import ensure_eligible_donors, filter_eligible_donors


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
REPO = HERE.parents[3]

IMAGES = Path.home() / "IO60panc2nd" / "PostXeniumIHC" / "Images"
GEOJSON_DIR = Path.home() / "IO60panc2nd" / "PostXeniumIHC" / "geojson"
REGISTRATION_INPUTS = PROJECT / "data" / "integration" / "restain_registration_inputs.csv"
XENIUM_PATHS = PROJECT / "data" / "integration" / "xenium_paths.csv"
XENIUM_ROOT = Path("/data/xenium")


class RegistrationInputError(RuntimeError):
    """A registration input is absent, ambiguous, or contradicts its manifest."""


@dataclass(frozen=True)
class RegistrationInput:
    donor_id: str
    serial: str
    roi: str
    qptiff: Path
    bundle: Path
    geojson: Path | None
    pixel_um: float
    dapi_page: int
    channel_names: tuple[str, ...]
    pyramid_shapes: tuple[tuple[int, int], ...]
    provenance_source: str


def _channel_name(description: str) -> str:
    match = re.search(r"<Name>([^<]*)</Name>", description or "")
    return "" if match is None else match.group(1).strip()


def qptiff_registration_facts(path: Path) -> dict:
    """Read the image facts registration is allowed to use from qptiff XML.

    ``PixelSizeMicrons`` is exact acquisition metadata.  The TIFF resolution rational is only
    approximate (0.25047 or 0.50139 in this cohort) and must not choose a pyramid level.  Every
    raw channel name is retained verbatim; only the literal DAPI channel is interpreted.
    """
    import tifffile

    path = Path(path)
    if not path.is_file():
        raise RegistrationInputError(f"qptiff is absent: {path}")
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = tuple((int(lv.shape[-2]), int(lv.shape[-1])) for lv in levels)
        names = tuple(_channel_name(page.description or "") for page in levels[0].pages)
        pixel_values: list[float] = []
        for page in levels[0].pages:
            matches = re.findall(r"<PixelSizeMicrons>\s*([^<]+?)\s*</PixelSizeMicrons>",
                                 page.description or "")
            pixel_values.extend(float(value) for value in matches)
        dtype = np.dtype(series.dtype)

    if dtype != np.dtype("uint8"):
        raise RegistrationInputError(f"{path.name}: registration expects uint8, found {dtype}")
    if not names or any(not name for name in names):
        raise RegistrationInputError(f"{path.name}: every qptiff channel must have an XML name")
    dapi = [i for i, name in enumerate(names) if name.strip().casefold() == "dapi"]
    if len(dapi) != 1:
        raise RegistrationInputError(
            f"{path.name}: expected exactly one literal DAPI channel, found {len(dapi)} in {names}")
    unique_px = sorted({float(value) for value in pixel_values})
    if len(unique_px) != 1 or not np.isfinite(unique_px[0]) or unique_px[0] <= 0:
        raise RegistrationInputError(
            f"{path.name}: expected one positive PixelSizeMicrons value, found {unique_px}")
    for level, (height, width) in enumerate(shapes):
        factor = 2 ** level
        expected = (shapes[0][0] // factor, shapes[0][1] // factor)
        if (height, width) != expected:
            raise RegistrationInputError(
                f"{path.name}: pyramid level {level} is {(height, width)}, expected floor-halved "
                f"{expected} from level 0")
    return {
        "pixel_um": unique_px[0],
        "dapi_page": dapi[0],
        "channel_names": names,
        "pyramid_shapes": shapes,
        "dtype": str(dtype),
    }


def _serial_match(source_path: str, serial: str) -> bool:
    return re.search(rf"(?:^|__){re.escape(serial)}(?:__|/|$)", source_path) is not None


def load_registration_inputs(
    inventory_csv: Path = REGISTRATION_INPUTS,
    xenium_paths_csv: Path = XENIUM_PATHS,
    images_dir: Path = IMAGES,
    geojson_dir: Path = GEOJSON_DIR,
    xenium_root: Path = XENIUM_ROOT,
    *,
    roi: str = "panc",
) -> dict[str, RegistrationInput]:
    """Load and cross-check every registration input, keyed by donor ID."""
    inventory_csv, xenium_paths_csv = Path(inventory_csv), Path(xenium_paths_csv)
    images_dir, geojson_dir, xenium_root = Path(images_dir), Path(geojson_dir), Path(xenium_root)
    inv = pd.read_csv(inventory_csv, comment="#", dtype=str)
    required = ("image_file", "donor_id", "roi", "xenium_serial", "provenance_source")
    missing = [column for column in required if column not in inv.columns]
    if missing:
        raise RegistrationInputError(f"{inventory_csv} is missing columns {missing}")
    for column in required:
        if inv[column].isna().any() or inv[column].str.strip().eq("").any():
            raise RegistrationInputError(f"{inventory_csv} has a blank {column!r}")
    inv = inv[inv.roi == roi].copy()
    for column in ("image_file", "donor_id", "xenium_serial"):
        dup = sorted(inv.loc[inv[column].duplicated(keep=False), column].unique())
        if dup:
            raise RegistrationInputError(f"{inventory_csv} has duplicate {column}: {dup}")

    xen = pd.read_csv(xenium_paths_csv, dtype=str)
    xen = xen[(xen.roi == roi) & xen.source_path.notna() & xen.source_path.ne("NA")]
    by_donor: dict[str, RegistrationInput] = {}
    for row in inv.itertuples(index=False):
        serial, donor = str(row.xenium_serial), str(row.donor_id)
        hits = xen[xen.source_path.map(lambda value: _serial_match(str(value), serial))]
        if len(hits) != 1:
            raise RegistrationInputError(
                f"serial {serial} has {len(hits)} {roi} rows in {xenium_paths_csv}; expected one")
        hit = hits.iloc[0]
        if str(hit.donor_id) != donor:
            raise RegistrationInputError(
                f"serial {serial}: registration inventory says donor {donor}, but "
                f"{xenium_paths_csv.name} says {hit.donor_id}")
        qptiff = images_dir / str(row.image_file)
        facts = qptiff_registration_facts(qptiff)
        geojsons = sorted(geojson_dir.glob(f"*{serial}*.geojson*")) if geojson_dir.exists() else []
        if len(geojsons) > 1:
            raise RegistrationInputError(f"serial {serial} has multiple GeoJSON exports: {geojsons}")
        bundle = xenium_root / str(hit.source_path).rstrip("/")
        by_donor[donor] = RegistrationInput(
            donor_id=donor,
            serial=serial,
            roi=roi,
            qptiff=qptiff,
            bundle=bundle,
            geojson=geojsons[0] if geojsons else None,
            pixel_um=float(facts["pixel_um"]),
            dapi_page=int(facts["dapi_page"]),
            channel_names=tuple(facts["channel_names"]),
            pyramid_shapes=tuple(facts["pyramid_shapes"]),
            provenance_source=str(row.provenance_source),
        )

    local = {path.name for path in images_dir.glob("*.qptiff")}
    recorded = {record.qptiff.name for record in by_donor.values()}
    extra = sorted(local - recorded)
    absent = sorted(recorded - local)
    if extra or absent:
        raise RegistrationInputError(
            f"registration inventory does not equal the local qptiff set; unrecorded={extra}, "
            f"absent={absent}")
    return by_donor


def select_registration_inputs(
    records: Mapping[str, RegistrationInput],
    donors: Sequence[str] | None,
    *,
    context: str = "post-Xenium IF registration",
) -> list[RegistrationInput]:
    """Select donors deterministically while enforcing the central cohort exclusions."""
    if donors:
        requested = list(ensure_eligible_donors((str(d) for d in donors), context=context))
        unknown = sorted(set(requested) - set(records))
        if unknown:
            raise RegistrationInputError(f"{context} requested donors absent from inventory: {unknown}")
        wanted = set(requested)
    else:
        wanted = set(filter_eligible_donors(records))
    return sorted((record for donor, record in records.items() if donor in wanted),
                  key=lambda record: (int(record.donor_id), record.serial))


def require_registration_files(record: RegistrationInput) -> None:
    """Fail before computation unless the three Xenium inputs are complete."""
    required = (
        record.bundle / "cells.parquet",
        record.bundle / "nucleus_boundaries.parquet",
        record.bundle / "morphology_focus" / "morphology_focus_0000.ome.tif",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RegistrationInputError(
            f"donor {record.donor_id} ({record.serial}) is missing registration files: {missing}")
