#!/usr/bin/env python
"""Step 3 driver: name the SOM clusters, behind the gate and the matched null.

    python scripts/kintsugi/run_step3.py --clusters 6551_armA_n25
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
sys.path.insert(0, str(HERE))

import annotate                                     # noqa: E402
import features                                     # noqa: E402
import panel                                        # noqa: E402


OUT = panel.REPO / "data" / "development" / "kintsugi"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True,
                    help="stem under data/development/kintsugi/clusters, e.g. 6551_armA_n25")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-cells", type=int, default=None,
                    help="subsample for the null gate; recorded in the output when used")
    args = ap.parse_args()

    clusters_path = OUT / "clusters" / f"{args.clusters}.parquet"
    assigned = pq.read_table(clusters_path).to_pandas()
    arm = args.clusters.split("_arm")[1][0]
    markers = features.feature_markers(arm)
    donors = sorted(assigned["donor_id"].unique())

    block, _ = features.load_block(donors, markers, arm)
    block = block.merge(assigned[["donor_id", "object_id", "cluster"]],
                        on=["donor_id", "object_id"], how="inner")
    if len(block) != len(assigned):
        raise RuntimeError(f"join lost cells: {len(assigned)} assigned, {len(block)} joined")

    used = block
    if args.max_cells and len(block) > args.max_cells:
        # Stratified by cluster so the null gate sees every cluster in proportion, with a
        # floor so a small cluster is not sampled below the point where it can be annotated.
        rng = np.random.default_rng(args.seed)
        keep = []
        for cluster, rows in block.groupby("cluster").indices.items():
            take = max(annotate.MIN_CLUSTER_CELLS,
                       int(round(args.max_cells * len(rows) / len(block))))
            keep.append(rows if len(rows) <= take
                        else rng.choice(rows, size=take, replace=False))
        used = block.iloc[np.sort(np.concatenate(keep))]
        print(f"  null gate on {len(used):,} of {len(block):,} cells "
              f"(stratified by cluster)", flush=True)

    t0 = time.time()
    calls, provenance = annotate.annotate(
        used[markers].to_numpy(dtype=np.float64), markers,
        used["cluster"].to_numpy(), used["donor_id"].to_numpy(), seed=args.seed)
    provenance["seconds"] = round(time.time() - t0, 1)
    provenance["subsampled_for_null"] = bool(args.max_cells and len(used) < len(block))
    provenance["n_cells_used"] = int(len(used))
    provenance["n_cells_total"] = int(len(block))

    out_dir = OUT / "calls"
    out_dir.mkdir(parents=True, exist_ok=True)
    calls.to_csv(out_dir / f"{args.clusters}_clusters.csv", index=False)

    per_cell = block[["donor_id", "object_id", "cell_region", "islet_num", "cluster"]].merge(
        calls[["cluster", "call", "reason", "n_candidates"]], on="cluster", how="left")
    per_cell.to_parquet(out_dir / f"{args.clusters}_cells.parquet", index=False)

    (out_dir / f"{args.clusters}_manifest.json").write_text(
        json.dumps(provenance, indent=2, default=str))

    named = calls[~calls["call"].isin([annotate.AMBIGUOUS, annotate.UNAVAILABLE])]
    print(f"\n  {len(named)} of {len(calls)} clusters named, "
          f"{provenance['seconds']}s", flush=True)
    with pd.option_context("display.width", 200, "display.max_colwidth", 34):
        print(calls[["cluster", "n_cells", "call", "reason", "n_candidates"]].to_string(
            index=False))
    share = per_cell["call"].value_counts(normalize=True)
    print("\n  cells by call:")
    for label, frac in share.items():
        print(f"    {label:<38} {frac:6.2%}")
    print(f"\nwrote {out_dir}/{args.clusters}_*.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
