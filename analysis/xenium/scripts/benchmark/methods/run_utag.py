#!/usr/bin/env python3
"""UTAG runner (ElementoLab/utag) — spatial-domain track. Runs in the `bench_py`
env (UTAG pins scanpy<1.12).

UTAG = message-passing over the spatial adjacency graph (A_norm @ X) then leiden.
We run a resolution grid in ONE utag() call (message-passing happens once; leiden
per resolution) and keep the partition whose domain count is closest to
--target-domains, so the head-to-head ARI is granularity-matched. With
--fixed-resolution, a single resolution is used (clean timing for the sweep).

Usage:
  conda activate bench_py
  python scripts/benchmark/methods/run_utag.py --sample 0041435 --variant crop25k
  python scripts/benchmark/methods/run_utag.py --sample 0041435 --variant crop200k --fixed-resolution 0.5
"""
from __future__ import annotations

import argparse

import _method_common as mc
import bench_common as bc
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--variant", default=bc.HEADLINE_VARIANT)
    ap.add_argument("--max-dist", type=float, default=30.0,
                    help="micron adjacency threshold (Xenium cell spacing ~10-15um)")
    ap.add_argument("--resolutions", type=float, nargs="+",
                    default=[0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5])
    ap.add_argument("--fixed-resolution", type=float, default=None)
    ap.add_argument("--target-domains", type=int, default=bc.N_DOMAINS_TARGET)
    ap.add_argument("--parallel", action="store_true",
                    help="enable UTAG multiprocessing (off by default for clean peak-RSS)")
    args = ap.parse_args()

    from utag import utag

    adata = mc.load_input(args.sample, args.variant)
    adata.obsm["spatial"] = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    resolutions = ([args.fixed_resolution] if args.fixed_resolution is not None
                   else list(args.resolutions))
    print(f"  [utag] {adata.n_obs:,} cells, max_dist={args.max_dist}, "
          f"resolutions={resolutions}, parallel={args.parallel}")

    out = utag(
        adata,
        slide_key=None,
        max_dist=args.max_dist,
        normalization_mode="l1_norm",
        apply_clustering=True,
        clustering_method=["leiden"],
        resolutions=resolutions,
        parallel=args.parallel,
    )

    sets = {c: out.obs[c].astype(str).to_numpy()
            for c in out.obs.columns if str(c).startswith("UTAG Label")}
    if not sets:
        raise RuntimeError(f"no 'UTAG Label*' columns produced; got {list(out.obs.columns)}")
    key, lab, n = mc.pick_closest_label_set(sets, args.target_domains)
    print(f"  [utag] selected '{key}' -> {n} domains")
    mc.write_labels(args.sample, args.variant, "utag", mc.cell_ids(adata), lab)


if __name__ == "__main__":
    main()
