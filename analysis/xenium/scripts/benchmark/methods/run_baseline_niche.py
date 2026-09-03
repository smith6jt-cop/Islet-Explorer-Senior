#!/usr/bin/env python3
"""Baseline domain method (5th method) — the stage-03 "spatial niche", recomputed
on the SAME benchmark subsample so it is comparable to BANKSY/UTAG/KODAMA.

Mechanism (ported from notebooks/03_spatial_analysis.ipynb cell 12):
  1. spatial k-NN (k=50)
  2. per-cell celltype-composition vector = mean one-hot celltype over neighbours
  3. KMeans(n_clusters=12)
  4. 2-pass k=30 majority-vote spatial smoothing

Runs in `xenium_analysis` (sklearn). Seed pinned to bench_common.SEED (=0;
stage-03 used 42 — we pin 0 across ALL benchmark methods for reproducibility).

Usage:
  python scripts/benchmark/methods/run_baseline_niche.py --sample 0041435 --variant crop25k
"""
from __future__ import annotations

import argparse

import _method_common as mc
import bench_common as bc
import numpy as np
import pandas as pd


def _knn_idx(xy: np.ndarray, k: int):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k, algorithm="auto").fit(xy)
    return nn.kneighbors(xy, return_distance=False)


def _majority(labels: np.ndarray, idx: np.ndarray) -> np.ndarray:
    n, k = idx.shape
    K = int(labels.max()) + 1
    counts = np.zeros((n, K), dtype=np.int32)
    neigh = labels[idx]
    rows = np.arange(n)
    for j in range(k):
        counts[rows, neigh[:, j]] += 1
    return counts.argmax(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--variant", default=bc.HEADLINE_VARIANT)
    ap.add_argument("--celltype-col", default="celltype")
    ap.add_argument("--k-comp", type=int, default=50)
    ap.add_argument("--n-clusters", type=int, default=bc.N_DOMAINS_TARGET)
    ap.add_argument("--k-smooth", type=int, default=30)
    ap.add_argument("--smooth-passes", type=int, default=2)
    args = ap.parse_args()

    from sklearn.cluster import KMeans

    adata = mc.load_input(args.sample, args.variant)
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    n = xy.shape[0]
    ct = adata.obs[args.celltype_col].astype(str).to_numpy()
    cats = pd.Categorical(ct)
    print(f"  [baseline] {n:,} cells, {len(cats.categories)} celltypes, "
          f"k_comp={args.k_comp}, n_clusters={args.n_clusters}")

    # composition = mean one-hot celltype over k spatial neighbours
    idx_c = _knn_idx(xy, args.k_comp)
    onehot = np.zeros((n, len(cats.categories)), dtype=np.float32)
    onehot[np.arange(n), cats.codes] = 1.0
    comp = onehot[idx_c].mean(axis=1)

    lab = KMeans(n_clusters=args.n_clusters, n_init=10, max_iter=300,
                 random_state=bc.SEED).fit_predict(comp)

    # 2-pass majority-vote spatial smoothing
    idx_s = _knn_idx(xy, args.k_smooth)
    for _ in range(args.smooth_passes):
        lab = _majority(lab, idx_s)

    mc.write_labels(args.sample, args.variant, "baseline", mc.cell_ids(adata), lab)


if __name__ == "__main__":
    main()
