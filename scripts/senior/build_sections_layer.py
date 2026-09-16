#!/usr/bin/env python3
"""Build sections/data.parquet from the two maintainer-recorded inventories.

Each Xenium inventory row names a section, not a bundle-selection rule. Regions
sharing a serial retain different frames; different serials for one donor/ROI
retain different rows. ``image_file`` preserves the recorded Xenium source path
(a bundle directory), or the recorded IF image filename. The explicit NA path
remains a row with null image/serial metadata; availability records its reason.

Only restain_sections.csv supplies pixel_size_um. Xenium has no scale in these
inventories, so its scale stays null. No image is opened by either builder.

The helpers here are shared with build_availability.py so both builders use the
exact donor normalisation in docs/ISLET_LAYERS_SPEC.md, section 0.3.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cohort import filter_eligible_donors
from phenocycler.integration.manifest import parse_bundle_path

SECTIONS_SCHEMA = pa.schema([
    ("donor_id", pa.string()),
    ("roi", pa.string()),
    ("modality", pa.string()),
    ("image_file", pa.string()),
    ("xenium_serial", pa.string()),
    ("pixel_size_um", pa.float64()),
    ("frame", pa.string()),
])


def norm_donor(value: str) -> str:
    """'0112' -> '112', using the shared layer contract without numeric coercion."""
    return str(value).strip().lstrip("0") or "0"


def read_inventory(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    """Read text identifiers, refusing blanks instead of inferring missing facts."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(line for line in handle if not line.lstrip().startswith("#"))
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: required columns are {required}")
        rows = []
        for row in reader:
            if None in row:
                raise ValueError(f"{path}: too many CSV fields")
            row = {key: (value or "").strip() for key, value in row.items()}
            for key in required:
                if not row[key] or (row[key].upper() == "NA" and key != "source_path"):
                    raise ValueError(f"{path}: blank or NA required field {key}")
            if not row["donor_id"].isdigit():
                raise ValueError(f"{path}: donor_id must be a string of digits")
            row["donor_id"] = norm_donor(row["donor_id"])
            if filter_eligible_donors([row["donor_id"]]):
                rows.append(row)
    return rows


def section_inventories(inventory_dir: Path) -> tuple[list[dict], list[dict]]:
    restain = read_inventory(inventory_dir / "restain_sections.csv", (
        "donor_id", "roi", "image_file", "xenium_serial", "pixel_size_um",
    ))
    xenium = read_inventory(inventory_dir / "xenium_paths.csv", (
        "donor_id", "roi", "source_path",
    ))
    return restain, xenium


def build_sections(inventory_dir: Path) -> pa.Table:
    restain, xenium = section_inventories(Path(inventory_dir))
    rows = []
    for row in xenium:
        source = row["source_path"]
        serial = None
        if source.upper() != "NA":
            serial, _ = parse_bundle_path(source)
            if not serial:
                raise ValueError(f"xenium_paths.csv: unrecognised bundle path {source!r}")
        rows.append({
            "donor_id": row["donor_id"], "roi": row["roi"], "modality": "xenium",
            "image_file": source if serial else None, "xenium_serial": serial,
            "pixel_size_um": None,
            "frame": f"xenium_um:{row['donor_id']}:{row['roi']}",
        })
    recorded_serials = {(r["donor_id"], r["roi"], r["xenium_serial"]) for r in rows}
    for row in restain:
        scale = float(row["pixel_size_um"])
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError(f"restain_sections.csv: invalid pixel_size_um for {row['image_file']}")
        if (row["donor_id"], row["roi"], row["xenium_serial"]) not in recorded_serials:
            raise ValueError(f"restain_sections.csv: section absent from xenium_paths.csv: {row}")
        rows.append({
            "donor_id": row["donor_id"], "roi": row["roi"], "modality": "restain_if",
            "image_file": row["image_file"], "xenium_serial": row["xenium_serial"],
            "pixel_size_um": scale, "frame": f"if_um:{row['donor_id']}:{row['roi']}",
        })
    unique = {}
    for row in rows:
        key = tuple(row[k] for k in ("donor_id", "roi", "modality", "xenium_serial"))
        if key in unique and unique[key] != row:
            raise ValueError(f"Conflicting records for section {key}")
        unique[key] = row
    if not unique:
        raise ValueError("The section inventories contain no eligible sections")
    return pa.Table.from_pylist(
        sorted(unique.values(), key=lambda r: tuple(str(r[k] or "") for k in (
            "donor_id", "roi", "modality", "xenium_serial"))),
        schema=SECTIONS_SCHEMA,
    )


def write_layer(table: pa.Table, layers_dir: Path, name: str) -> Path:
    """Replace this builder's single artifact only after the parquet write succeeds."""
    target = Path(layers_dir) / name / "data.parquet"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".parquet.tmp")
    try:
        pq.write_table(table, temporary, compression="zstd")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--repo-root", type=Path, default=REPO)
    parser.add_argument("--inventory-dir", type=Path,
                        help="Default: <repo-root>/core/phenocycler/data/integration")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args(argv)
    inventory = args.inventory_dir or args.repo_root / "core/phenocycler/data/integration"
    table = build_sections(inventory)
    path = write_layer(table, args.layers_dir, "sections")
    print(json.dumps({"artifact": str(path), "sections": table.num_rows,
                      "donors": len(set(table["donor_id"].to_pylist()))}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
