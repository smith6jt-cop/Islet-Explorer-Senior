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

from phenocycler.cohort import DONOR_EXCLUSIONS, donor_of, excluded_on

#: What a registration reads: the Xenium DAPI and the post-Xenium IF DAPI scanned on the same
#: slide. No PhenoCycler data enters it, so an exclusion covering only a donor's PhenoCycler data
#: -- 6457's (maintainer, 2026-10-22: "the xenium for 6457 is not excluded") -- does not refuse it.
REGISTRATION_PLATFORMS = ("xenium", "if")


def registration_excluded_on(donor: str) -> list[str]:
    """The platforms this registration reads that the donor's central exclusion covers."""
    return [p for p in REGISTRATION_PLATFORMS if excluded_on(donor, p)]


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
    pixel_size_source: str = "qptiff XML PixelSizeMicrons"


def _channel_name(description: str) -> str:
    match = re.search(r"<Name>([^<]*)</Name>", description or "")
    return "" if match is None else match.group(1).strip()


def qptiff_registration_facts(path: Path, *, allow_dtypes: tuple[str, ...] = ("uint8",),
                              pixel_um_override: float | None = None) -> dict:
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

    # uint8 is a property of the post-Xenium IF cohort, not a requirement of registration:
    # the correlation is computed on mean-subtracted, variance-normalised tiles. A
    # PhenoCycler Fusion export is uint16. Widen only when a caller declares it, so the
    # existing cohort keeps its exact guard.
    if str(dtype) not in allow_dtypes:
        raise RegistrationInputError(
            f"{path.name}: registration expects {allow_dtypes}, found {dtype}")
    if not names or any(not name for name in names):
        raise RegistrationInputError(f"{path.name}: every qptiff channel must have an XML name")
    dapi = [i for i, name in enumerate(names) if name.strip().casefold() == "dapi"]
    if len(dapi) != 1:
        raise RegistrationInputError(
            f"{path.name}: expected exactly one literal DAPI channel, found {len(dapi)} in {names}")
    unique_px = sorted({float(value) for value in pixel_values})
    if pixel_um_override is not None:
        # Some exports carry no PixelSizeMicrons at all. The docstring's warning is about the
        # TIFF resolution rational being approximate -- it is not an argument against a value
        # measured from the data. For 6533 the scale is fitted from the cell GeoJSON against
        # the measurement table's own micron centroids (0.507800 um/px, 24 nm median
        # residual), which is stronger evidence than the tag would have been.
        if not np.isfinite(pixel_um_override) or pixel_um_override <= 0:
            raise RegistrationInputError(f"{path.name}: pixel_um override must be positive")
        if unique_px and not np.isclose(unique_px[0], pixel_um_override, rtol=1e-3):
            raise RegistrationInputError(
                f"{path.name}: declared pixel_um {pixel_um_override} contradicts the XML "
                f"value {unique_px[0]}")
        unique_px = [float(pixel_um_override)]
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
        "pixel_size_source": ("registration inventory (fitted, cross-checked)"
                              if pixel_um_override is not None
                              else "qptiff XML PixelSizeMicrons"),
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
        declared_px = getattr(row, "pixel_um", None)
        declared_dt = getattr(row, "allow_dtype", None)
        facts = qptiff_registration_facts(
            qptiff,
            allow_dtypes=(("uint8",) if declared_dt in (None, "")
                          else tuple(str(declared_dt).split("|"))),
            pixel_um_override=(None if declared_px in (None, "")
                               else float(declared_px)))
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
            pixel_size_source=str(facts["pixel_size_source"]),
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
        requested = [str(d) for d in donors]
        blocked = sorted(d for d in requested if registration_excluded_on(d))
        if blocked:
            details = "; ".join(
                f"{d} (on {', '.join(registration_excluded_on(d))}): "
                f"{DONOR_EXCLUSIONS[donor_of(d)]}" for d in blocked)
            raise ValueError(f"{context} includes excluded donor(s): {details}")
        unknown = sorted(set(requested) - set(records))
        if unknown:
            raise RegistrationInputError(f"{context} requested donors absent from inventory: {unknown}")
        wanted = set(requested)
    else:
        wanted = {d for d in records if not registration_excluded_on(d)}
    return sorted((record for donor, record in records.items() if donor in wanted),
                  key=lambda record: (int(record.donor_id), record.serial))


#: Xenium writes the DAPI morphology plane under two different names depending on the
#: analysis version: ``morphology_focus_0000.ome.tif`` on the cohort this lane was built for,
#: and ``ch0000_dapi.ome.tif`` on newer runs (0080870 is one). Channel 0 is DAPI either way,
#: so resolve by trying both rather than pinning the older spelling.
MORPHOLOGY_DAPI_NAMES = ("morphology_focus_0000.ome.tif", "ch0000_dapi.ome.tif")


def morphology_dapi_path(bundle) -> Path:
    """The bundle's DAPI morphology file, whichever naming this Xenium version used."""
    focus = Path(bundle) / "morphology_focus"
    for name in MORPHOLOGY_DAPI_NAMES:
        candidate = focus / name
        if candidate.is_file():
            return candidate
    return focus / MORPHOLOGY_DAPI_NAMES[0]


def require_registration_files(record: RegistrationInput) -> None:
    """Fail before computation unless the three Xenium inputs are complete."""
    required = (
        record.bundle / "cells.parquet",
        record.bundle / "nucleus_boundaries.parquet",
        morphology_dapi_path(record.bundle),
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RegistrationInputError(
            f"donor {record.donor_id} ({record.serial}) is missing registration files: {missing}")
