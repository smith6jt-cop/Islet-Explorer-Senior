#!/usr/bin/env python3
"""
Joined per-cell protein + transcript table for a same-slide pair.

One row per Xenium cell. Protein is transferred from PhenoCycler cells by AREA-WEIGHTED
POLYGON OVERLAP in the registered frame: each PhenoCycler cell contributes to each Xenium
cell in proportion to the area they share. That neither forces a bijection nor throws away
extent, which matters because the two segmentations disagree (1,028,140 vs 717,461 cells).

``overlap_fraction`` is the per-row honesty column: a cell assembled from a sliver of one
protein cell is not the same evidence as one covered by a well-matched partner.

Genes go to a sparse sidecar rather than 5,001 parquet columns -- dense would be ~14 GB and
every contract validation would walk all of it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from phenocycler.integration.contract import (  # noqa: E402
    CONTRACT_PROFILE_GEOMETRY, coerce_cell_table, feature_columns, write_cell_table,
)
from phenocycler.integration.transform import Transform  # noqa: E402


def xenium_polygons(bundle: Path) -> tuple[np.ndarray, np.ndarray]:
    """Xenium cell outlines (already microns), one polygon per cell_id."""
    b = pd.read_parquet(bundle / "cell_boundaries.parquet")
    b = b.sort_values("cell_id", kind="stable")
    ids, starts = np.unique(b["cell_id"].to_numpy(), return_index=True)
    order = np.argsort(starts)
    ids, starts = ids[order], starts[order]
    coords = b[["vertex_x", "vertex_y"]].to_numpy(dtype=float)
    idx = np.repeat(np.arange(len(starts)), np.diff(np.append(starts, len(b))))
    polys = shapely.polygons(shapely.linearrings(coords, indices=idx))
    return ids, polys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno", required=True, type=Path)
    ap.add_argument("--pheno-polygons", required=True, type=Path)
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--transform-dir", required=True, type=Path)
    ap.add_argument("--transform-stem", default="refined_field2")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--serial", required=True)
    ap.add_argument("--tissue", default="pancreas")
    args = ap.parse_args(argv)

    t0 = time.time()
    tf = Transform.load(args.transform_dir, stem=args.transform_stem)

    xids, xpolys = xenium_polygons(args.bundle)
    print(f"[joined] xenium polygons {len(xids):,} ({time.time()-t0:.0f}s)", flush=True)
    # Warp Xenium FORWARD into the PhenoCycler frame. Transform.invert() drops the
    # displacement field, so the inverse is a different transform, not just a slower one.
    xw = shapely.transform(xpolys, lambda a: tf.apply(a))
    xarea = shapely.area(xw)

    pp = pd.read_parquet(args.pheno_polygons)
    ppolys = shapely.from_wkb(pp["wkb"].to_numpy())
    pids = pp["cell_id"].to_numpy()
    print(f"[joined] phenocycler polygons {len(pids):,} ({time.time()-t0:.0f}s)", flush=True)

    # Warping through a displacement field can produce self-intersecting rings, and GEOS
    # raises a TopologyException rather than returning a wrong area. Repair only the ones
    # that are actually invalid -- make_valid on all 1.7M is needless work.
    for name, arr in (("xenium", xw), ("phenocycler", ppolys)):
        bad = ~shapely.is_valid(arr)
        if bad.any():
            arr[bad] = shapely.make_valid(arr[bad])
            print(f"[joined] repaired {int(bad.sum()):,} invalid {name} polygons", flush=True)

    tree = shapely.STRtree(ppolys)
    xi, pi = tree.query(xw, predicate="intersects")
    print(f"[joined] candidate pairs {len(xi):,} ({time.time()-t0:.0f}s)", flush=True)
    inter = shapely.area(shapely.intersection(xw[xi], ppolys[pi]))
    keep = inter > 0
    xi, pi, inter = xi[keep], pi[keep], inter[keep]
    print(f"[joined] overlapping pairs {len(xi):,} ({time.time()-t0:.0f}s)", flush=True)

    pheno = pd.read_parquet(args.pheno)
    pos = pd.Series(np.arange(len(pheno)), index=pheno["cell_id"].to_numpy())
    prow = pos.reindex(pids[pi]).to_numpy()
    ok = np.isfinite(prow)
    xi, pi, inter, prow = xi[ok], pi[ok], inter[ok], prow[ok].astype(int)

    n = len(xids)
    wsum = np.bincount(xi, weights=inter, minlength=n)
    feats = feature_columns(pheno)
    out = {"cell_id": xids}
    vals = pheno[feats].to_numpy(dtype=np.float32)
    for j, c in enumerate(feats):
        num = np.bincount(xi, weights=inter * np.nan_to_num(vals[prow, j]), minlength=n)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[c] = np.where(wsum > 0, num / np.maximum(wsum, 1e-12), np.nan).astype(np.float32)
    print(f"[joined] aggregated {len(feats)} features ({time.time()-t0:.0f}s)", flush=True)

    cent = shapely.centroid(xpolys)
    df = pd.DataFrame(out)
    df["x_um"] = shapely.get_x(cent)
    df["y_um"] = shapely.get_y(cent)
    df["area_um2"] = shapely.area(xpolys)
    df["n_pheno_cells"] = np.bincount(xi, minlength=n)
    df["overlap_fraction"] = np.where(xarea > 0, wsum / np.maximum(xarea, 1e-12), 0.0)
    df["xenium_sample_key"] = f"{args.serial}__{args.roi}"

    df = coerce_cell_table(df, modality="xenium", donor_id=args.donor,
                           profile=CONTRACT_PROFILE_GEOMETRY, roi=args.roi, tissue=args.tissue)
    args.out.mkdir(parents=True, exist_ok=True)
    path = write_cell_table(df, args.out / "joined_cells.parquet")

    meta = {"n_xenium_cells": int(n), "n_pheno_cells": int(len(pids)),
            "n_overlapping_pairs": int(len(xi)),
            "transform_stem": args.transform_stem,
            "n_protein_features": len(feats),
            "cells_with_protein": int((df["n_pheno_cells"] > 0).sum()),
            "overlap_fraction_median": float(np.median(df["overlap_fraction"])),
            "n_pheno_cells_median": float(np.median(df["n_pheno_cells"])),
            "seconds": round(time.time() - t0, 1)}
    (args.out / "build_joined.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[joined] wrote {path}")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
