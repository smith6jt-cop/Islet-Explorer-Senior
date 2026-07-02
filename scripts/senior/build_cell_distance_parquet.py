#!/usr/bin/env python3
"""
Build the per-cell "cell_distance" parquet that powers the Spatial-tab Card C
right-hand distance plots (per-cell signed distance to nearest islet, selectable
by phenotype OR RESTORE-positive marker, split by disease status; plus a
marker-composition-vs-distance bar view).

Three facts, three sources, joined here once (offline) so the app never joins:
  * signed distance      -> data/cells/donor_id=<case>/*.parquet   (`dist_islet`,
        negative inside / positive outside; + `object_id`, `cell_region`)
  * app phenotype label  -> data/parquet/tissue/case_id=<case>/*.parquet
        (`phenotype`; same cells, same ROW ORDER as the cells parquet — verified)
  * RESTORE positivity   -> data/restore_gated_redsea/donor_id=<case>/*.parquet
        (`{m}_pos` booleans for 10 lineage markers, gated on REDSEA-corrected MFI;
        joined by `object_id` — NOT row order, gated has a few extra edge cells)

`_pos` supersedes the earlier raw-intensity percentile heuristic: it is a
calibrated per-image threshold (RESTORE, mutually-exclusive negative controls).

Output : data/parquet/cell_distance/case_id=<case>/part-0.parquet
         columns: donor_status, cell_region, phenotype, dist_islet,
                  + {INS,GCG,SST,CD3e,CD20,CD163,CD31,SMA,Vimentin,Pan_Cytokeratin}_pos
         (case_id comes from the hive path)

Run    : conda activate scvi-env && python scripts/senior/build_cell_distance_parquet.py
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
CELLS_ROOT = DATA / "cells"                          # donor_id=<case>/*.parquet (dist + object_id)
TISSUE_ROOT = DATA / "parquet" / "tissue"            # case_id=<case>/*.parquet (phenotype)
GATED_ROOT = DATA / "restore_gated_redsea"           # donor_id=<case>/*.parquet ({m}_pos)
OUT_ROOT = DATA / "parquet" / "cell_distance"        # case_id=<case>/part-0.parquet (this script)
STATUS_H5AD = DATA / "islets_core_senior.h5ad"       # imageid -> donor_status map (tiny, per-islet)

# The 10 markers RESTORE gates (chosen=True in restore_thresholds_redsea.csv).
POS_MARKERS = ["INS", "GCG", "SST", "CD3e", "CD20", "CD163", "CD31", "SMA", "Vimentin", "Pan_Cytokeratin"]
POS_COLS = [f"{m}_pos" for m in POS_MARKERS]

COORD_TOL_UM = 1.0        # per-row |Δ| tolerance between cells/tissue coordinates
MIN_AGREE_FRAC = 0.999    # abort a donor if fewer rows than this agree within COORD_TOL_UM


def read_obs_series(f: h5py.File, key: str) -> np.ndarray:
    node = f["obs"][key]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        cats = node["categories"][:]
        cats = np.array([c.decode() if isinstance(c, (bytes, bytearray)) else c for c in cats])
        codes = node["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], None)
    arr = node[:]
    if arr.dtype.kind in ("S", "O"):
        return np.array([a.decode() if isinstance(a, (bytes, bytearray)) else a for a in arr])
    return arr


def build_status_map() -> dict[str, str]:
    if not STATUS_H5AD.exists():
        sys.exit(f"[cell_distance] status source missing: {STATUS_H5AD}")
    with h5py.File(STATUS_H5AD, "r") as f:
        keys = set(f["obs"].keys())
        img_key = "imageid" if "imageid" in keys else next(k for k in keys if k.lower() == "imageid")
        st_key = "donor_status" if "donor_status" in keys else "Donor Status"
        imageid = read_obs_series(f, img_key).astype(str)
        status = read_obs_series(f, st_key).astype(str)
    m: dict[str, str] = {}
    for i, s in zip(imageid, status):
        m.setdefault(str(i), str(s))
    return m


def read_parquet_ordered(part_dir: Path, columns: list[str]) -> pd.DataFrame:
    files = sorted(part_dir.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no parquet files in {part_dir}")
    frames = [pd.read_parquet(fp, columns=columns) for fp in files]
    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0].reset_index(drop=True)


def process_donor(case: str, status: str) -> tuple[int, int, pd.DataFrame]:
    cells = read_parquet_ordered(
        CELLS_ROOT / f"donor_id={case}",
        ["object_id", "X_centroid", "Y_centroid", "cell_region", "dist_islet"],
    )
    tissue = read_parquet_ordered(
        TISSUE_ROOT / f"case_id={case}", ["X_centroid", "Y_centroid", "phenotype"]
    )
    if len(cells) != len(tissue):
        sys.exit(f"[cell_distance][{case}] ROW COUNT MISMATCH cells={len(cells):,} tissue={len(tissue):,}")

    # Positional (row-order) join of phenotype is only valid if the two parquets
    # describe the same cell at each ordinal position — verify via coordinates.
    dx = np.abs(cells["X_centroid"].to_numpy() - tissue["X_centroid"].to_numpy())
    dy = np.abs(cells["Y_centroid"].to_numpy() - tissue["Y_centroid"].to_numpy())
    agree = float(np.mean((dx < COORD_TOL_UM) & (dy < COORD_TOL_UM)))
    if agree < MIN_AGREE_FRAC:
        sys.exit(f"[cell_distance][{case}] COORD AGREEMENT {agree:.4f} < {MIN_AGREE_FRAC} — cells/tissue no longer row-aligned; aborting")
    cells["phenotype"] = tissue["phenotype"].to_numpy()   # positional attach (order preserved)

    # RESTORE _pos by object_id (unique key → merge preserves cells order).
    gdir = GATED_ROOT / f"donor_id={case}"
    if gdir.exists():
        gated = read_parquet_ordered(gdir, ["object_id"] + POS_COLS)
        merged = cells.merge(gated, on="object_id", how="left")
        n_pos_matched = int(merged[POS_COLS[0]].notna().sum())
    else:
        merged = cells.copy()
        for c in POS_COLS:
            merged[c] = False
        n_pos_matched = 0

    out = pd.DataFrame({
        "donor_status": pd.Categorical([status] * len(merged)),
        "cell_region": merged["cell_region"].astype("category"),
        "phenotype": merged["phenotype"].astype("category"),
        "dist_islet": merged["dist_islet"].astype("float32"),
    })
    for c in POS_COLS:
        out[c] = merged[c].fillna(False).astype(bool).to_numpy()
    return int(round(agree * len(cells))), n_pos_matched, out


def main() -> None:
    ap = argparse.ArgumentParser(description="Build data/parquet/cell_distance/ (Spatial Card C).")
    ap.add_argument("--donor", action="append", help="limit to these case ids (repeatable)")
    args = ap.parse_args()

    if not CELLS_ROOT.exists() or not TISSUE_ROOT.exists():
        sys.exit(f"[cell_distance] need both {CELLS_ROOT} and {TISSUE_ROOT}")
    if not GATED_ROOT.exists():
        print(f"[cell_distance] WARNING: {GATED_ROOT} missing — all _pos will be FALSE (marker modes empty)")

    status_map = build_status_map()
    cells_cases = {p.name.split("=", 1)[1] for p in CELLS_ROOT.glob("donor_id=*")}
    tissue_cases = {p.name.split("=", 1)[1] for p in TISSUE_ROOT.glob("case_id=*")}
    cases = sorted(cells_cases & tissue_cases, key=lambda x: (len(x), x))
    if args.donor:
        cases = [c for c in cases if c in set(args.donor)]
    if not cases:
        sys.exit("[cell_distance] no donors to process")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    total = 0
    print(f"[cell_distance] {len(cases)} donors | markers: {', '.join(POS_MARKERS)}")
    for case in cases:
        status = status_map.get(case)
        if status is None:
            sys.exit(f"[cell_distance][{case}] no donor_status in {STATUS_H5AD.name}; aborting")
        n_agree, n_pos, df = process_donor(case, status)
        ddir = OUT_ROOT / f"case_id={case}"
        if ddir.exists():
            shutil.rmtree(ddir)
        ddir.mkdir(parents=True)
        df.to_parquet(ddir / "part-0.parquet", index=False)
        total += len(df)
        ins = int(df["INS_pos"].sum()); cd3 = int(df["CD3e_pos"].sum())
        print(f"[{case}] {len(df):>9,} cells | {status:>4} | coord {100*n_agree/len(df):5.1f}% | "
              f"_pos matched {100*n_pos/len(df):5.1f}% | INS+ {ins:>7,}  CD3e+ {cd3:>6,}")
    print(f"[cell_distance] DONE — {total:,} cells across {len(cases)} donors -> {OUT_ROOT}")


if __name__ == "__main__":
    main()
