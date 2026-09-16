#!/usr/bin/env python3
"""Build the post-Xenium IF layer described in ISLET_LAYERS_SPEC.md §1.4.

``warped_if_um.npy`` contains Xenium centroids warped INTO the IF frame,
not IF points mapped into Xenium. These coordinates are already micrometres:
they must not be rescaled or transformed again. The declared pixel size is read
only from ``restain_sections.csv``; this builder never opens a qptiff.

The output lives entirely under ``$ISLET_LAYERS_DIR/restain``. Its small
``availability.parquet`` records every donor, counts before/after the nucleus
filter, the CSV calibration, and the producer's hormone validation. The
``cells/<donor>.parquet`` files contain the IF-frame overlay. A missing value
always has an explicit availability state and reason beside it.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))
from phenocycler.cohort import filter_eligible_donors  # noqa: E402

MARKERS = ("INS", "GCG")
STATES = ("measured", "blank", "unavailable", "withheld", "failed")


def norm_donor(value: str) -> str:
    """The shared spelling prescribed by spec §0.3: '0112' -> '112'."""
    return str(value).strip().lstrip("0") or "0"


def read_sections(path: Path) -> pd.DataFrame:
    """Read the maintainer's declared calibration, without inspecting images."""
    sections = pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)
    required = {"donor_id", "roi", "image_file", "xenium_serial", "pixel_size_um"}
    if not required.issubset(sections):
        raise ValueError(f"{path}: missing section columns {sorted(required - set(sections))}")
    if (sections[list(required)].apply(lambda col: col.str.strip()) == "").any().any():
        raise ValueError(f"{path}: blank required section value")
    sections["donor_id"] = sections["donor_id"].map(norm_donor)
    sections["pixel_size_um"] = pd.to_numeric(sections["pixel_size_um"], errors="raise")
    if not (np.isfinite(sections.pixel_size_um) & (sections.pixel_size_um > 0)).all():
        raise ValueError(f"{path}: pixel_size_um must be finite and positive")
    if sections.duplicated(["donor_id", "roi"]).any():
        raise ValueError(f"{path}: duplicate donor/roi calibration")
    return sections


def read_source(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    """Read file columns, never infer conflicting Hive partition columns."""
    return ds.dataset(path, format="parquet", partitioning=None).to_table(
        columns=columns
    ).to_pandas()


def source_inventory(root: Path) -> dict[str, list[Path]]:
    """Discover pancreas partitions from their recorded columns, not names."""
    inventory: dict[str, list[Path]] = {}
    paths = sorted(root.rglob("*.parquet"))
    if not paths:
        raise FileNotFoundError(f"No Xenium cell parquet files under {root}")
    for path in paths:
        identities = read_source(path, ["donor_id", "roi"]).drop_duplicates()
        identities["donor_id"] = identities.donor_id.map(norm_donor)
        identities = identities.drop_duplicates()
        if len(identities) != 1:
            raise ValueError(f"{path}: expected one recorded donor/roi per source partition")
        row = identities.iloc[0]
        if row.roi == "panc" and filter_eligible_donors([row.donor_id]):
            inventory.setdefault(row.donor_id, []).append(path)
    return inventory


def warp_inventory(root: Path) -> dict[str, Path]:
    """Use the export's explicit donor and ROI; filenames are not donor metadata."""
    found: dict[str, Path] = {}
    for path in sorted(root.glob("*/export.json")):
        record = json.loads(path.read_text())
        donor = norm_donor(record["donor_id"])
        if record["roi"] != "panc" or not filter_eligible_donors([donor]):
            continue
        if donor in found:
            raise ValueError(f"Duplicate pancreas warp exports for donor {donor}")
        found[donor] = path.parent
    return found


def validate_warp_count(warped: np.ndarray, n_xenium_cells: int, donor: str) -> None:
    """Assert the full source count BEFORE filtering nuclei or writing any file."""
    n_xenium_warped_to_if = len(warped)
    if n_xenium_warped_to_if != n_xenium_cells:
        raise AssertionError(
            f"{donor}: n_xenium_warped_to_if={n_xenium_warped_to_if} "
            f"!= n_xenium_cells={n_xenium_cells}"
        )
    if warped.shape != (n_xenium_cells, 2):
        raise ValueError(f"{donor}: expected an (n_xenium_cells, 2) warp, got {warped.shape}")
    if not np.isfinite(warped).all():
        raise ValueError(f"{donor}: nonfinite coordinates in the IF warp")


def hormone_records(path: Path) -> dict[str, dict]:
    if not path.is_file():
        return {}
    records: dict[str, dict] = {}
    for row in json.loads(path.read_text()):
        donor = norm_donor(row["donor_id"])
        if row["roi"] != "panc":
            continue
        if donor in records:
            raise ValueError(f"{path}: duplicate pancreas hormone validation for {donor}")
        records[donor] = row["hormone"]
    return records


def hormone_availability(record: dict, marker: str) -> dict:
    """Preserve the recorded release decision; never fit or invent positivity."""
    evidence = record.get(marker, {})
    auroc = evidence.get("auroc")
    bar = record.get("preregistered_auroc_bar")
    if auroc is None or bar is None or not np.isfinite([auroc, bar]).all():
        return {"state": "unavailable", "reason": "no_hormone_validation",
                "auroc": None, "bar": None}
    if auroc < bar or evidence.get("clears_preregistered_bar") is not True:
        return {"state": "withheld", "reason": "auroc_below_preregistered_bar",
                "auroc": float(auroc), "bar": float(bar)}
    # A validation result does not itself supply any per-cell released calls.
    return {"state": "unavailable", "reason": "no_released_hormone_calls",
            "auroc": float(auroc), "bar": float(bar)}


def _schema() -> pa.Schema:
    fields = [(name, pa.string()) for name in
              ("donor_id", "state", "reason", "source", "frame", "image_file",
               "object_id_state", "object_id_reason")]
    fields += [("pixel_size_um", pa.float64())]
    fields += [(name, pa.int64()) for name in
               ("n", "n_xenium_cells", "n_xenium_warped_to_if", "n_without_nucleus")]
    for marker in MARKERS:
        fields += [(f"n_measured_if_{marker}", pa.int64())]
        prefix = f"hormone_positivity_{marker}_"
        fields += [(prefix + "state", pa.string()), (prefix + "reason", pa.string()),
                   (prefix + "auroc", pa.float64()), (prefix + "bar", pa.float64())]
    return pa.schema(fields)


def _status(donor: str, source: Path) -> dict:
    row = dict.fromkeys(_schema().names)
    row.update(donor_id=donor, state="unavailable", reason="no_restain_pairing",
               source=str(source), frame=f"if_um:{donor}:panc",
               object_id_state="unavailable", object_id_reason="no_qupath_object_id")
    for marker in MARKERS:
        row[f"hormone_positivity_{marker}_state"] = "unavailable"
        row[f"hormone_positivity_{marker}_reason"] = "no_restain_pairing"
    return row


def inspect_donor(paths: list[Path], warp_dir: Path, donor: str, section: pd.Series) -> tuple:
    """Validate source identity and full warp cardinality without aligning records."""
    frames = []
    for path in paths:
        names = set(pq.ParquetFile(path).schema_arrow.names)
        required = {"donor_id", "roi", "cell_id", *(f"feat__if_{m}" for m in MARKERS)}
        if not required.issubset(names):
            raise ValueError(f"{path}: missing IF source columns {sorted(required - names)}")
        columns = sorted(required | ({"object_id"} if "object_id" in names else set()))
        frames.append(read_source(path, columns))
    cells = pd.concat(frames, ignore_index=True)
    cells["donor_id"] = cells.donor_id.map(norm_donor)
    if not cells.donor_id.eq(donor).all() or not cells.roi.eq("panc").all():
        raise ValueError(f"{donor}: source donor/roi disagrees with requested IF section")
    if cells.cell_id.isna().any() or cells.cell_id.duplicated().any() or cells.cell_id.eq("").any():
        raise ValueError(f"{donor}: missing or duplicate Xenium cell_id")
    export = json.loads((warp_dir / "export.json").read_text())
    if (norm_donor(export["donor_id"]), export["roi"], str(export["serial"])) != (
        donor, "panc", section.xenium_serial
    ):
        raise ValueError(f"{donor}: warp export donor/roi/serial disagrees with sections CSV")
    warped = np.load(warp_dir / "warped_if_um.npy", allow_pickle=False)
    validate_warp_count(warped, len(cells), donor)
    return cells, warped


def prepare_restain_cells(cells: pd.DataFrame, warped: np.ndarray,
                          warp_dir: Path, status: dict) -> tuple:
    """Refuse assembly until the contract's UUID join can actually be performed.

    The measured Xenium inputs have only native cell_id values. Neither a
    fabricated object_id nor positional attachment is a UUID join. Keep this
    failure ahead of every output write until the missing crosswalk or an
    explicit exception to the identifier contract is supplied.
    """
    donor = status["donor_id"]
    if "object_id" not in cells:
        raise ValueError(
            f"{donor}: unavailable(no_qupath_object_id): source has cell_id but no "
            "QuPath object_id; UUID-only assembly cannot join warped_if_um.npy. "
            "A producer-recorded UUID crosswalk or an explicit native-cell_id "
            "exception to BRIEF.md is required. No artifacts were written."
        )
    raise ValueError(
        f"{donor}: unavailable(no_uuid_keyed_warp): {warp_dir / 'warped_if_um.npy'} "
        "has no recorded object_id keys. Refusing positional or centroid matching."
    )


def build_restain_layer(*, cells_xen_root: Path, restain_root: Path,
                        sections_csv: Path, xenium_paths_csv: Path,
                        layers_dir: Path) -> pd.DataFrame:
    """Validate every paired donor before creating any output, then write one root."""
    sections = read_sections(sections_csv)
    sections = sections[sections.roi == "panc"].set_index("donor_id")
    sources = source_inventory(cells_xen_root)
    warps = warp_inventory(restain_root / "valis")
    validation = hormone_records(restain_root / "hormone_panc_refined.json")
    known = pd.read_csv(xenium_paths_csv, dtype=str, comment="#", keep_default_na=False)
    donors = sorted(set(filter_eligible_donors(
        set(known.donor_id.map(norm_donor)) | set(sections.index) | set(sources)
    )))
    statuses, overlays = [], {}
    for donor in donors:
        row = _status(donor, cells_xen_root)
        if donor in sections.index:
            section = sections.loc[donor]
            row.update(pixel_size_um=float(section.pixel_size_um), image_file=section.image_file)
        paths = sources.get(donor, [])
        paired = bool(paths) and all(
            all(f"feat__if_{m}" in pq.ParquetFile(path).schema_arrow.names for m in MARKERS)
            for path in paths
        )
        if not paired:
            statuses.append(row)
            continue
        if donor not in sections.index:
            row.update(reason="no_restain_section_calibration")
            statuses.append(row)
            continue
        if donor not in warps:
            row.update(reason="no_xenium_to_if_warp")
            statuses.append(row)
            continue
        section = sections.loc[donor]
        try:
            cells, warped = inspect_donor(paths, warps[donor], donor, section)
        except FileNotFoundError as exc:
            row.update(reason=f"missing_restain_input:{Path(exc.filename).name}")
            statuses.append(row)
            continue
        row.update(n_xenium_cells=len(cells), n_xenium_warped_to_if=len(warped),
                   source=";".join(str(p) for p in paths))
        for marker in MARKERS:
            evidence = hormone_availability(validation.get(donor, {}), marker)
            for field, value in evidence.items():
                row[f"hormone_positivity_{marker}_{field}"] = value
        # Assembly is kept separate so all shape/identity assertions precede writes.
        overlay, details = prepare_restain_cells(cells, warped, warps[donor], row)
        row.update(details)
        if overlay is not None:
            overlays[donor] = overlay
        statuses.append(row)
    status_table = pa.Table.from_pylist(statuses, schema=_schema())
    destination = Path(layers_dir) / "restain"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "cells").mkdir(exist_ok=True)
    for donor, table in overlays.items():
        pq.write_table(table, destination / "cells" / f"{donor}.parquet", compression="zstd")
    # Publish availability last. Stale donor cell files are never selected without
    # a corresponding measured row in this build's availability table.
    pq.write_table(status_table, destination / "availability.parquet", compression="zstd")
    return status_table.to_pandas()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    core_data = REPO / "core/phenocycler/data"
    parser.add_argument("--cells-xen-root", type=Path,
                        default=core_data / "integration/cells_xen")
    parser.add_argument("--restain-root", type=Path,
                        default=REPO / "data/development/restain")
    parser.add_argument("--sections-csv", type=Path,
                        default=core_data / "integration/restain_sections.csv")
    parser.add_argument("--xenium-paths-csv", type=Path,
                        default=core_data / "integration/xenium_paths.csv")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))
    args = parser.parse_args(argv)
    result = build_restain_layer(**vars(args))
    print(result[["donor_id", "state", "reason", "n", "n_xenium_cells",
                  "n_xenium_warped_to_if", "pixel_size_um"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
