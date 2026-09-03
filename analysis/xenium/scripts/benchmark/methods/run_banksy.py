#!/usr/bin/env python3
"""BANKSY runner (prabhakarlab/Banksy_py, pip `pybanksy`) — spatial-domain track.
Runs in the `bench_py` env. REQUIRES LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6
(matplotlib in bench_py needs a newer CXXABI than the system libstdc++) — the
driver / profile_wrap sets this; see env_bench_py.sh.

Pipeline (domain mode, lambda=0.8 + AGF max_m=1):
  initialize_banksy -> generate_banksy_matrix -> pca_umap(20) -> run_Leiden_partition
We sweep a resolution grid in ONE partition call and keep the row whose
num_labels is closest to --target-domains (granularity-matched head-to-head);
--fixed-resolution forces a single resolution (clean timing for the sweep).

Usage (driver sets env):
  python scripts/benchmark/methods/run_banksy.py --sample 0041435 --variant crop25k
"""
from __future__ import annotations

import argparse
import contextlib
import io

import _method_common as mc
import bench_common as bc
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--variant", default=bc.HEADLINE_VARIANT)
    ap.add_argument("--lam", type=float, default=0.8,
                    help="0.8 = spatial-domain segmentation; 0.2 = cell-typing")
    ap.add_argument("--num-neighbours", type=int, default=18)
    ap.add_argument("--pca-dims", type=int, default=20)
    ap.add_argument("--resolutions", type=float, nargs="+",
                    default=[0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.6])
    ap.add_argument("--fixed-resolution", type=float, default=None)
    ap.add_argument("--target-domains", type=int, default=bc.N_DOMAINS_TARGET)
    ap.add_argument("--quiet", action="store_true", default=True)
    args = ap.parse_args()

    from banksy.cluster_methods import run_Leiden_partition
    from banksy.embed_banksy import generate_banksy_matrix
    from banksy.initialize_banksy import initialize_banksy
    from banksy.run_banksy import pca_umap

    adata = mc.load_input(args.sample, args.variant)
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    adata.obs["x"], adata.obs["y"] = xy[:, 0], xy[:, 1]
    resolutions = ([args.fixed_resolution] if args.fixed_resolution is not None
                   else list(args.resolutions))
    print(f"  [banksy] {adata.n_obs:,} cells, lambda={args.lam}, "
          f"k_geom={args.num_neighbours}, pca={args.pca_dims}, res={resolutions}")

    sink = io.StringIO()
    cm_ctx = contextlib.redirect_stdout(sink) if args.quiet else contextlib.nullcontext()
    with cm_ctx:
        bd = initialize_banksy(
            adata, ("x", "y", "spatial"), num_neighbours=args.num_neighbours,
            nbr_weight_decay="scaled_gaussian", max_m=1,
            plt_edge_hist=False, plt_nbr_weights=False,
            plt_agf_angles=False, plt_theta=False,
        )
        bd, _ = generate_banksy_matrix(adata, bd, lambda_list=[args.lam],
                                       max_m=1, verbose=False)
        pca_umap(bd, pca_dims=[args.pca_dims], add_umap=False, plt_remaining_var=False)
        df, _ = run_Leiden_partition(
            bd, resolutions=resolutions, num_nn=50, num_iterations=-1,
            partition_seed=bc.SEED, match_labels=False,
        )

    # pick the partition closest to the target domain count
    nl = df["num_labels"].to_numpy()
    pick = int(np.argmin(np.abs(nl - args.target_domains)))
    row = df.iloc[pick]
    labels = np.asarray(row["labels"].dense).ravel()
    print(f"  [banksy] selected res={row['resolution']} -> {row['num_labels']} domains")
    mc.write_labels(args.sample, args.variant, "banksy", mc.cell_ids(adata), labels)


if __name__ == "__main__":
    main()
