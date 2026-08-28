"""Domain-track metrics for the spatial-method benchmark.

All metrics are REPORTED, never gated (Tier-C). Each carries a bias note in
METRIC_BIAS so the user interprets fairly. Two families:
  * reference-free spatial coherence : CHAOS, PAS, Moran's I of one-hot domains
  * reference-free expression separation : silhouette / Davies-Bouldin /
    Calinski-Harabasz on a NEUTRAL embedding (X_scvi — used by none of the
    methods, so it does not favour any one).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

# Higher-is-better (True) / lower-is-better (False) + one-line bias caveat.
METRIC_DIRECTION = {
    "n_domains": None, "size_gini": None,
    "chaos": False, "pas": False, "moran_domains": True,
    "silhouette_scvi": True, "davies_bouldin_scvi": False,
    "calinski_harabasz_scvi": True,
}
METRIC_BIAS = {
    "n_domains": "descriptive only — reported, never gated",
    "size_gini": "descriptive: 0=equal-sized domains, ->1=one dominant domain",
    "chaos": "favours over-smoothed/few-domain solutions; read with n_domains",
    "pas": "favours over-smoothed solutions (UTAG/baseline smooth); read with n_domains",
    "moran_domains": "favours spatially-contiguous domains; smoothing inflates it",
    "silhouette_scvi": "favours partitions aligned with scVI structure",
    "davies_bouldin_scvi": "favours scVI-aligned, compact clusters",
    "calinski_harabasz_scvi": "grows with cluster count and scVI alignment",
}


def spatial_knn_idx(xy: np.ndarray, k: int) -> np.ndarray:
    """Indices of the k nearest spatial neighbours (self excluded)."""
    nn = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    idx = nn.kneighbors(xy, return_distance=False)
    return idx[:, 1:]


def _codes(labels) -> np.ndarray:
    return pd.factorize(np.asarray(labels))[0]


def chaos(labels, xy: np.ndarray) -> float:
    """Mean distance from each cell to its nearest SAME-domain neighbour
    (spatial CHAOS; lower = more spatially coherent)."""
    labels = np.asarray(labels)
    tot, n = 0.0, 0
    for d in pd.unique(labels):
        m = labels == d
        if m.sum() < 2:
            continue
        nn = NearestNeighbors(n_neighbors=2).fit(xy[m])
        dist, _ = nn.kneighbors(xy[m])
        tot += float(dist[:, 1].sum())
        n += int(m.sum())
    return tot / max(n, 1)


def pas(labels, knn_idx: np.ndarray) -> float:
    """Percentage of Abnormal Spots: fraction of cells whose label differs from
    the majority label among spatial neighbours (lower = smoother)."""
    codes = _codes(labels)
    K = int(codes.max()) + 1
    neigh = codes[knn_idx]
    n, k = neigh.shape
    counts = np.zeros((n, K), dtype=np.int32)
    rows = np.arange(n)
    for j in range(k):
        counts[rows, neigh[:, j]] += 1
    maj = counts.argmax(axis=1)
    return float(np.mean(codes != maj))


def moran_onehot(labels, knn_idx: np.ndarray) -> float:
    """Mean Moran's I across one-hot domain indicators on the kNN graph
    (higher = more spatial autocorrelation of domain membership)."""
    codes = _codes(labels)
    n = len(codes)
    K = int(codes.max()) + 1
    rows = np.repeat(np.arange(n), knn_idx.shape[1])
    cols = knn_idx.ravel()
    W = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    S0 = W.sum()
    Is = []
    for d in range(K):
        z = (codes == d).astype(float)
        z -= z.mean()
        den = float((z * z).sum())
        if den > 0:
            num = float(z @ (W @ z))
            Is.append((n / S0) * (num / den))
    return float(np.mean(Is)) if Is else float("nan")


def expression_separation(labels, X: np.ndarray, seed: int = 0,
                          sample_size: int = 10000) -> dict:
    """silhouette / Davies-Bouldin / Calinski-Harabasz on a neutral embedding."""
    from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
    labels = np.asarray(labels)
    out = {"silhouette_scvi": float("nan"), "davies_bouldin_scvi": float("nan"),
           "calinski_harabasz_scvi": float("nan")}
    if X is None or len(pd.unique(labels)) < 2:
        return out
    ss = None if X.shape[0] <= sample_size else sample_size
    try:
        out["silhouette_scvi"] = float(silhouette_score(
            X, labels, sample_size=ss, random_state=seed))
        out["davies_bouldin_scvi"] = float(davies_bouldin_score(X, labels))
        out["calinski_harabasz_scvi"] = float(calinski_harabasz_score(X, labels))
    except Exception as e:  # noqa: BLE001
        print(f"    [warn] expression separation failed: {e}")
    return out


def _gini(sizes: np.ndarray) -> float:
    s = np.sort(np.asarray(sizes, dtype=float))
    n = len(s)
    if n == 0 or s.sum() == 0:
        return float("nan")
    cum = np.cumsum(s)
    return float((n + 1 - 2 * (cum / cum[-1]).sum()) / n)


def compute_domain_metrics(labels, xy: np.ndarray, X_embed: np.ndarray,
                           knn_idx: np.ndarray, seed: int = 0) -> dict:
    labels = np.asarray(labels).astype(str)
    sizes = pd.Series(labels).value_counts().to_numpy()
    out = {
        "n_domains": int(len(sizes)),
        "size_gini": _gini(sizes),
        "chaos": chaos(labels, xy),
        "pas": pas(labels, knn_idx),
        "moran_domains": moran_onehot(labels, knn_idx),
    }
    out.update(expression_separation(labels, X_embed, seed=seed))
    return out
