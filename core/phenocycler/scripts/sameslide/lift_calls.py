#!/usr/bin/env python3
"""Carry the PhenoCycler cell-type call onto Xenium cells, with its own uncertainty.

``build_joined.py`` transfers protein INTENSITIES by area-weighted mean, which is the right
operator for a continuous quantity. A cell-type call is categorical and a weighted mean of
labels is meaningless, so it is carried by LARGEST SHARED AREA instead -- the Xenium cell
takes the call of whichever PhenoCycler cell it overlaps most.

Two columns say how much that is worth, and Stage 3 needs both because a prior is only as
good as its certainty:

* ``call_share``  -- the winning PhenoCycler cell's fraction of this Xenium cell's overlapped
  area. Low means the Xenium cell was assembled from several protein cells and the winner
  merely came first.
* ``call_purity`` -- the fraction of overlapped area whose PhenoCycler cells carry the SAME
  call. High purity with low share is a Xenium cell sitting inside a uniform patch, which is
  a perfectly good prior; low purity is a cell straddling a boundary, which is not.

The two segmentations disagree by design (1,028,140 PhenoCycler cells against 717,461 Xenium),
so neither number is a defect to be tuned away. They are the honest record of a many-to-one
transfer.

Geometry is identical to ``build_joined.py`` -- the same transform, the same forward warp of
Xenium polygons into the PhenoCycler frame, the same intersection areas -- so a cell's protein
features and its protein call rest on exactly the same overlap.

    python scripts/sameslide/lift_calls.py \
        --calls  ../../data/development/celltype/sameslide_protein/donor_id=6533/data.parquet \
        --pheno-polygons <dir>/pheno_cell_polygons.parquet \
        --bundle <xenium bundle> --transform-dir <valis dir> --out <dir>
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_joined import xenium_polygons                      # noqa: E402
from phenocycler.integration.transform import Transform       # noqa: E402

#: Columns carried across, when the call table has them.
CARRY = ("call", "reason", "cell_region", "structure", "islet_endocrine", "islet_num")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--pheno-polygons", required=True, type=Path)
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--transform-dir", required=True, type=Path)
    ap.add_argument("--transform-stem", default="refined_field2")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    t0 = time.time()
    tf = Transform.load(args.transform_dir, stem=args.transform_stem)
    xids, xpolys = xenium_polygons(args.bundle)
    xw = shapely.transform(xpolys, lambda a: tf.apply(a))
    xarea = shapely.area(xw)
    print(f"[lift] xenium polygons {len(xids):,} ({time.time()-t0:.0f}s)", flush=True)

    pp = pd.read_parquet(args.pheno_polygons)
    ppolys = shapely.from_wkb(pp["wkb"].to_numpy())
    pids = pp["cell_id"].astype(str).to_numpy()
    for name, arr in (("xenium", xw), ("phenocycler", ppolys)):
        bad = ~shapely.is_valid(arr)
        if bad.any():
            arr[bad] = shapely.make_valid(arr[bad])
            print(f"[lift] repaired {int(bad.sum()):,} invalid {name} polygons", flush=True)

    xi, pi = shapely.STRtree(ppolys).query(xw, predicate="intersects")
    inter = shapely.area(shapely.intersection(xw[xi], ppolys[pi]))
    keep = inter > 0
    xi, pi, inter = xi[keep], pi[keep], inter[keep]
    print(f"[lift] overlapping pairs {len(xi):,} ({time.time()-t0:.0f}s)", flush=True)

    calls = pd.read_parquet(args.calls)
    key = "object_id" if "object_id" in calls.columns else "cell_id"
    pos = pd.Series(np.arange(len(calls)), index=calls[key].astype(str).to_numpy())
    prow = pos.reindex(pids[pi]).to_numpy()
    ok = np.isfinite(prow)
    dropped = int((~ok).sum())
    xi, inter, prow = xi[ok], inter[ok], prow[ok].astype(int)
    print(f"[lift] {len(xi):,} pairs resolve to a called PhenoCycler cell "
          f"({dropped:,} dropped: not in the call table, e.g. QC-ineligible)", flush=True)

    n = len(xids)
    wsum = np.bincount(xi, weights=inter, minlength=n)

    # Winner = largest shared area. np.maximum.at over a (cell, area) reduction rather than a
    # sort: 2.3M pairs sorted by two keys is slower and needs no ordering guarantee here.
    best = np.zeros(n)
    np.maximum.at(best, xi, inter)
    is_best = inter >= best[xi]                      # ties keep the first, deterministically
    first = np.full(n, -1, dtype=np.int64)
    sel = np.flatnonzero(is_best)
    first[xi[sel[::-1]]] = sel[::-1]
    has = first >= 0

    out = pd.DataFrame({"cell_id": xids})
    win = np.where(has, first, 0)
    for c in CARRY:
        if c not in calls.columns:
            continue
        v = calls[c].to_numpy()[prow[win]]
        out[c] = np.where(has, v, None if v.dtype == object else np.nan)
    out["call_share"] = np.where(has & (wsum > 0), inter[win] / np.maximum(wsum, 1e-12), np.nan)

    # Purity: the share of overlapped area agreeing with the winner's call.
    codes, _ = pd.factorize(calls["call"].astype(str))
    agree = codes[prow] == codes[prow[win]][xi]
    num = np.bincount(xi[agree], weights=inter[agree], minlength=n)
    out["call_purity"] = np.where(wsum > 0, num / np.maximum(wsum, 1e-12), np.nan)

    out["n_pheno_cells"] = np.bincount(xi, minlength=n)
    out["overlap_fraction"] = np.where(xarea > 0, wsum / np.maximum(xarea, 1e-12), 0.0)

    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "protein_calls.parquet"
    out.to_parquet(path, index=False)
    meta = {
        "n_xenium_cells": int(n), "n_with_call": int(has.sum()),
        "n_pairs": int(len(xi)), "n_pairs_dropped_uncalled": dropped,
        "transform_stem": args.transform_stem,
        "call_share_median": float(np.nanmedian(out["call_share"])),
        "call_purity_median": float(np.nanmedian(out["call_purity"])),
        "call_purity_ge_0.8": float(np.nanmean(out["call_purity"] >= 0.8)),
        "seconds": round(time.time() - t0, 1),
    }
    (args.out / "lift_calls.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    print(f"[lift] wrote {path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
