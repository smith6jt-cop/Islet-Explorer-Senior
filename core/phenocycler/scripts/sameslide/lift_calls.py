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

from build_joined import xenium_polygons, overlap_edges, require_unique_ids  # noqa: E402
from phenocycler.integration.transform import Transform       # noqa: E402

#: Columns carried across, when the call table has them.
CARRY = ("call", "reason", "cell_region", "structure", "islet_endocrine", "islet_num")


def transfer_calls(xids, edges, calls, xarea):
    """Transfer using the same eligible UUID edges as the continuous measurements."""
    key = "object_id" if "object_id" in calls else "cell_id"
    require_unique_ids(calls, key)
    calls = calls.copy()
    calls[key] = calls[key].astype(str)
    e = edges[edges.object_id.isin(calls[key])].copy()
    if not np.isfinite(e.overlap_area_um2).all() or (e.overlap_area_um2 <= 0).any():
        raise ValueError("overlap areas must be finite and positive")
    if e.duplicated(["cell_id", "object_id"]).any():
        raise ValueError("duplicate overlap edges")
    xi = pd.Index(xids).get_indexer(e.cell_id)
    if (xi < 0).any():
        raise ValueError("overlap edges contain unknown Xenium cells")
    area = np.asarray(xarea, float)
    if area.shape != (len(xids),) or not np.isfinite(area).all() or (area <= 0).any():
        raise ValueError("registered Xenium areas must be finite and positive")
    total = np.bincount(xi, weights=e.overlap_area_um2, minlength=len(xids))
    e = e.merge(calls, left_on="object_id", right_on=key, validate="many_to_one",
                suffixes=("", "_protein"))
    # Resolve equal areas by immutable UUID, not input row order.
    winners = e.sort_values(["overlap_area_um2", "object_id"], ascending=[False, True],
                            kind="stable").drop_duplicates("cell_id").set_index("cell_id")
    winner = winners.reindex(xids)
    out = pd.DataFrame({"cell_id": xids})
    for col in (*CARRY, "status"):
        if col in winner:
            out[col] = winner[col].to_numpy()
    if "status" not in out:
        out["status"] = np.where(winner["call"].notna(), "awarded", "unavailable")
    out["winning_object_id"] = winner.object_id.to_numpy()
    out["call_share"] = np.divide(winner.overlap_area_um2.to_numpy(float), total,
                                 out=np.full(len(xids), np.nan), where=total > 0)
    winning_call = pd.Series(winner["call"].to_numpy(), index=xids)
    agrees = e["call"].notna() & e["call"].eq(e.cell_id.map(winning_call))
    numerator = np.bincount(pd.Index(xids).get_indexer(e.loc[agrees, "cell_id"]),
                            weights=e.loc[agrees, "overlap_area_um2"], minlength=len(xids))
    out["call_purity"] = np.divide(numerator, total, out=np.full(len(xids), np.nan),
                                  where=total > 0)
    out["n_pheno_cells"] = np.bincount(xi, minlength=len(xids))
    out["overlap_fraction"] = total / area
    out.loc[winner["call"].isna().to_numpy(), "status"] = "unavailable"
    return out


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
    print(f"[lift] xenium polygons {len(xids):,} ({time.time()-t0:.0f}s)", flush=True)

    pp = pd.read_parquet(args.pheno_polygons)
    edges, xarea = overlap_edges(xids, xpolys, pp, tf)
    calls = pd.read_parquet(args.calls)
    if "qc_analysis_eligible" in calls:
        calls = calls[calls.qc_analysis_eligible]
    out = transfer_calls(xids, edges, calls, xarea)
    n = len(xids)
    pairs = int(out.n_pheno_cells.sum())

    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "protein_calls.parquet"
    out.to_parquet(path, index=False)
    meta = {
        "n_xenium_cells": int(n), "n_with_call": int(out.winning_object_id.notna().sum()),
        "n_pairs": pairs, "n_pairs_dropped_uncalled": len(edges) - pairs,
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
