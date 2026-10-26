"""Shared pLN partition helpers for protein and Xenium typing.

The functions were moved verbatim from cluster_pln.py. Both historical lanes
import them here, so their seed and clustering methods have one code owner.
"""

from __future__ import annotations

import hashlib

import numpy as np


def stable_seed(*parts: object) -> int:
    """A seed derived from strings that is the same in every interpreter.

    `hash('6461')` is salted per process unless PYTHONHASHSEED is set, so seeding with
    `SEED + hash(donor)` made every run of this script a different run wearing the same
    command.  BLAKE2b over the joined parts has no such salt.  The value is reduced into the
    32-bit range because scanpy, scikit-learn and igraph all take `random_state` as an int.
    """
    digest = hashlib.blake2b("|".join(str(p) for p in parts).encode("utf-8"),
                             digest_size=8).digest()
    return int.from_bytes(digest, "big") % (2 ** 31 - 1)

def leiden_on_subsample(
    features: np.ndarray,
    *,
    n_pcs: int,
    subsample: int,
    resolution: float,
    seed: int,
) -> tuple[np.ndarray, dict]:
    """Leiden on a random subsample, then k-NN assignment of every remaining cell.

    Leiden needs a neighbour graph, and a graph over five million cells is not affordable per
    donor.  The partition is therefore learned on a random subsample and carried to the rest by
    a 15-nearest-neighbour vote in the same PCA space.  The cost of that shortcut is measured
    rather than assumed: the subsample is split 80/20, the classifier is refit on the 80 and
    its agreement with Leiden's own labels on the held-out 20 is reported as
    `holdout_agreement`.  A low number there means the assignment step, not the partition, is
    what a reader should distrust.
    """
    import anndata as ad
    import scanpy as sc
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KNeighborsClassifier

    rng = np.random.default_rng(seed)
    n_cells = features.shape[0]
    n_pcs = int(min(n_pcs, features.shape[1] - 1, n_cells - 1))
    pca = PCA(n_components=n_pcs, random_state=seed)
    coords = pca.fit_transform(features).astype(np.float32)

    take = min(subsample, n_cells)
    idx = rng.choice(n_cells, size=take, replace=False) if take < n_cells else np.arange(n_cells)
    idx = np.sort(idx)

    sub = ad.AnnData(X=coords[idx].copy())
    sub.obsm["X_pca"] = coords[idx].copy()
    sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_pca", random_state=seed)
    sc.tl.leiden(
        sub,
        resolution=resolution,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=seed,
    )
    sub_labels = sub.obs["leiden"].to_numpy().astype(str)

    # measured cost of the assignment shortcut
    perm = rng.permutation(take)
    cut = int(take * 0.8)
    tr, te = perm[:cut], perm[cut:]
    probe = KNeighborsClassifier(n_neighbors=15, n_jobs=-1)
    probe.fit(coords[idx][tr], sub_labels[tr])
    holdout = float((probe.predict(coords[idx][te]) == sub_labels[te]).mean())

    labels = np.empty(n_cells, dtype=object)
    labels[idx] = sub_labels
    rest = np.setdiff1d(np.arange(n_cells), idx, assume_unique=True)
    if rest.size:
        clf = KNeighborsClassifier(n_neighbors=15, n_jobs=-1)
        clf.fit(coords[idx], sub_labels)
        step = 200_000
        for start in range(0, rest.size, step):
            chunk = rest[start:start + step]
            labels[chunk] = clf.predict(coords[chunk])

    meta = {
        "n_cells": int(n_cells),
        "n_pcs": n_pcs,
        "pca_variance_explained": float(pca.explained_variance_ratio_.sum()),
        "subsample_used_for_graph": int(take),
        "leiden_resolution": resolution,
        "n_clusters": int(len(set(sub_labels))),
        "knn_assignment_k": 15,
        "knn_holdout_agreement": holdout,
        "cells_assigned_by_knn": int(rest.size),
        "seed": int(seed),
    }
    return labels.astype(str), meta
