#!/home/smith6jt/miniconda3/envs/scvi-env/bin/python
"""Group the per-cell KRONOS2 embeddings exactly as CORAL's Tutorial 4 documents it.

Verbatim from ``tutorials/4-Patch-Clustering.ipynb`` (the only label-free route CORAL documents;
Tutorial 5 phenotyping is a supervised classifier on ground-truth cell labels, which this project
does not have and may not use -- QuPath classes are audit context only)::

    LEIDEN_RESOLUTION = 0.5  # higher -> more, smaller clusters
    SEED = 0                 # pin every random step (neighbours, Leiden, UMAP)

    adata = ad.AnnData(X)                                    # X = (n, 768) float32, raw embeddings
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X", random_state=SEED)   # no intermediate PCA
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, random_state=SEED, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=SEED)

Run on the KRONOS2 embeddings AND on the mean-marker baseline the tutorial pairs with them, so the
two are read side by side.  Each cluster is then characterised by its mean-marker profile -- the
"classic per-cell expression vector, the standard input to cell phenotyping" in CORAL's own words.
No label is assigned here: the maintainer reads the profiles.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import pandas as pd

RESOLUTION = 0.5
SEED = 0
N_NEIGHBORS = 15


def documented_clustering(X: np.ndarray, *, resolution: float = RESOLUTION, seed: int = SEED) -> tuple:
    import anndata as ad
    import scanpy as sc
    adata = ad.AnnData(np.asarray(X, dtype=np.float32))
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X", random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=seed)
    return adata.obs["leiden"].astype(int).to_numpy(), np.asarray(adata.obsm["X_umap"], dtype=np.float32)


def profiles(clusters: np.ndarray, M: np.ndarray, markers: list) -> pd.DataFrame:
    """Mean-marker profile per cluster: the raw mean and its z-score ACROSS clusters (display only)."""
    df = pd.DataFrame(M, columns=markers)
    df["cluster"] = clusters
    mean = df.groupby("cluster").mean()
    z = (mean - mean.mean(axis=0)) / mean.std(axis=0).replace(0, np.nan)
    mean.columns = [f"mean_{c}" for c in mean.columns]
    z.columns = [f"z_{c}" for c in z.columns]
    out = mean.join(z)
    out.insert(0, "n_cells", pd.Series(clusters).value_counts().sort_index())
    return out.reset_index()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True, help="the worker's features.npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--resolution", type=float, default=RESOLUTION)
    ap.add_argument("--seed", type=int, default=SEED)
    a = ap.parse_args(argv)

    z = np.load(a.features, allow_pickle=True)
    markers = [str(m) for m in z["markers"]]
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rec = {"resolution": a.resolution, "seed": a.seed, "n_neighbors": N_NEIGHBORS,
           "recipe": "sc.pp.neighbors(n_neighbors=15, use_rep='X') -> sc.tl.leiden(flavor='igraph', n_iterations=2, "
                     "directed=False) -> sc.tl.umap; CORAL Tutorial 4 verbatim",
           "n_cells": int(z["labels"].shape[0]), "spaces": {}}
    for space in ("kronos2", "mean_marker"):
        t0 = time.time()
        cl, um = documented_clustering(z[space], resolution=a.resolution, seed=a.seed)
        prof = profiles(cl, z["mean_marker"], markers)
        pd.DataFrame({"label": z["labels"], "region_idx": z["region_idx"], "cluster": cl,
                      "umap_1": um[:, 0], "umap_2": um[:, 1]}).to_parquet(out / f"clusters__{space}.parquet", index=False)
        prof.to_csv(out / f"profiles__{space}.csv", index=False)
        rec["spaces"][space] = {"n_clusters": int(cl.max() + 1), "dim": int(z[space].shape[1]),
                                "sizes": pd.Series(cl).value_counts().sort_index().to_dict(),
                                "seconds": round(time.time() - t0, 1)}
        print(f"[doc-cluster] {space}: {int(cl.max() + 1)} clusters at resolution {a.resolution} "
              f"({time.time() - t0:.0f}s)", flush=True)
    # descriptive: how much the two groupings agree (never a criterion)
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    ck = pd.read_parquet(out / "clusters__kronos2.parquet")["cluster"].to_numpy()
    cm = pd.read_parquet(out / "clusters__mean_marker.parquet")["cluster"].to_numpy()
    rec["kronos2_vs_mean_marker"] = {"ari": float(adjusted_rand_score(ck, cm)),
                                     "ami": float(adjusted_mutual_info_score(ck, cm))}
    (out / "cluster_record.json").write_text(json.dumps(rec, indent=1, default=str))
    print(json.dumps({k: v for k, v in rec.items() if k != "spaces"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
