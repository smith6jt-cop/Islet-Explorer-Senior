#!/usr/bin/env python3
"""
3A (Islet Explorer Senior): stream the QuPath per-cell export into a tidy,
donor-partitioned Parquet dataset via DuckDB.

Input  : ~/IO60panc2nd/Cellmeasurements.csv  (~23.3M cells, 1,228 cols, ~138 GB)
Output : <out>/cells/donor_id=<id>/*.parquet  — ~65 cols:
         donor_id, image, object_id, X_centroid, Y_centroid, parent_raw,
         classification, cell_region (core|peri|tissue|other), islet_num,
         cell_area, nucleus_area, and 59 whole-cell marker means
         (`Cell: <marker>: Mean` -> sanitized marker name).

Phenotyping is deferred to a separate step that runs on this Parquet (the raw
QuPath intensities need their own gating thresholds, distinct from the scVI
scale used by Junior's phenotype_rules.csv).

Region scheme mirrors Junior: Parent `Islet_N` = core, `Islet_N_exp20um` = peri,
`Annotation (Tissue)` = tissue.

Usage:
    conda activate scvi-env
    python scripts/senior/build_cells_parquet.py            # full run
    python scripts/senior/build_cells_parquet.py --limit 200000 --out /tmp/senior_test
"""
from __future__ import annotations
import argparse, csv, re, sys
from pathlib import Path
import duckdb

CSV_DEFAULT = Path.home() / "IO60panc2nd" / "Cellmeasurements.csv"
OUT_DEFAULT = Path(__file__).resolve().parents[2] / "data"


def header_cols(csv_path: Path) -> list[str]:
    with open(csv_path, newline="") as f:
        return next(csv.reader(f))


def marker_means(cols: list[str]) -> dict[str, str]:
    """Map 'Cell: <marker>: Mean' -> <marker>."""
    out = {}
    for c in cols:
        m = re.match(r"^Cell: (.+): Mean$", c)
        if m:
            out[c] = m.group(1)
    return out


def sanitize(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def build_sql(csv_path: Path, out_dir: Path, limit: int | None,
              min_cell_area: float) -> tuple[str, dict]:
    cols = header_cols(csv_path)
    mm = marker_means(cols)
    if not mm:
        raise SystemExit("No 'Cell: <marker>: Mean' columns found — check the CSV header.")
    selects = [
        "regexp_extract(\"Image\", '[0-9]+') AS donor_id",  # first digit-run: '6374...'->6374, 'HDL115pancLN...'->115
        '"Image" AS image',
        '"Object ID" AS object_id',
        "TRY_CAST(\"Centroid X µm\" AS DOUBLE) AS X_centroid",
        "TRY_CAST(\"Centroid Y µm\" AS DOUBLE) AS Y_centroid",
        '"Parent" AS parent_raw',
        '"Classification" AS classification',
        ("CASE WHEN \"Parent\" LIKE '%exp20um' THEN 'peri' "
         "WHEN regexp_matches(\"Parent\", '^Islet_[0-9]+$') THEN 'core' "
         "WHEN \"Parent\" = 'Annotation (Tissue)' OR \"Parent\" LIKE 'Root object%' THEN 'tissue' "
         "ELSE 'other' END AS cell_region"),
        "regexp_extract(\"Parent\", 'Islet_([0-9]+)', 1) AS islet_num",
        "TRY_CAST(\"Cell: Area µm^2\" AS DOUBLE) AS cell_area",
        "TRY_CAST(\"Nucleus: Area µm^2\" AS DOUBLE) AS nucleus_area",
        # extra morphology + signed-distance columns (needed by fixed_islet_aggregation / 3B)
        "TRY_CAST(\"Cell: Length µm\" AS DOUBLE) AS cell_length",
        "TRY_CAST(\"Cell: Circularity\" AS DOUBLE) AS cell_circularity",
        "TRY_CAST(\"Cell: Solidity\" AS DOUBLE) AS cell_solidity",
        "TRY_CAST(\"Cell: Max diameter µm\" AS DOUBLE) AS cell_maxdiam",
        "TRY_CAST(\"Cell: Min diameter µm\" AS DOUBLE) AS cell_mindiam",
        "TRY_CAST(\"Nucleus: Length µm\" AS DOUBLE) AS nucleus_length",
        "TRY_CAST(\"Nucleus: Circularity\" AS DOUBLE) AS nucleus_circularity",
        "TRY_CAST(\"Nucleus: Solidity\" AS DOUBLE) AS nucleus_solidity",
        "TRY_CAST(\"Nucleus: Max diameter µm\" AS DOUBLE) AS nucleus_maxdiam",
        "TRY_CAST(\"Nucleus: Min diameter µm\" AS DOUBLE) AS nucleus_mindiam",
        "TRY_CAST(\"Nucleus/Cell area ratio\" AS DOUBLE) AS nc_area_ratio",
        "TRY_CAST(\"Signed distance to annotation Islet µm\" AS DOUBLE) AS dist_islet",
        "TRY_CAST(\"Signed distance to annotation Tissue µm\" AS DOUBLE) AS dist_tissue",
    ]
    for orig, marker in mm.items():
        selects.append(f'TRY_CAST("{orig}" AS DOUBLE) AS "{sanitize(marker)}"')
    sel = ",\n  ".join(selects)
    src = (f"read_csv('{csv_path}', header=true, all_varchar=true, "
           f"ignore_errors=true)")
    inner = f"SELECT\n  {sel}\nFROM {src}"
    if min_cell_area and min_cell_area > 0:
        # Drop degenerate sub-cell segmentation objects: QuPath emits NaN
        # intensities for these tiny artifacts. Real cells have area
        # p1=12.3, median=57.9 µm²; the NaN population is all <1.9 µm².
        area = "TRY_CAST(\"Cell: Area µm^2\" AS DOUBLE)"
        inner += f"\nWHERE isfinite({area}) AND {area} >= {min_cell_area}"
    if limit:
        inner += f"\nLIMIT {limit}"
    cells_dir = out_dir / "cells"
    sql = (f"COPY (\n{inner}\n) TO '{cells_dir}' "
           f"(FORMAT PARQUET, PARTITION_BY (donor_id), "
           f"OVERWRITE_OR_IGNORE, COMPRESSION ZSTD)")
    return sql, mm


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=CSV_DEFAULT)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--min-cell-area", type=float, default=5.0,
                    help="Drop cells below this Cell:Area µm^2 (segmentation debris; 0 to disable)")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    sql, mm = build_sql(a.csv, a.out, a.limit, a.min_cell_area)
    print(f"[info] markers={len(mm)} out={a.out}/cells limit={a.limit} min_cell_area={a.min_cell_area}")
    a.out.mkdir(parents=True, exist_ok=True)  # DuckDB COPY PARTITION_BY needs the parent to exist
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={a.threads}")
    con.execute("SET enable_progress_bar=false")
    con.execute(sql)
    print("[done] wrote", a.out / "cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
