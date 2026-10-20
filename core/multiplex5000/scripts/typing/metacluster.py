#!/usr/bin/env python3
"""Merge fine clusters down to a number a person can actually read.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  This module belongs to the Nimbus/Leiden chain
(``substrate.py`` -> ``positivity.py`` -> ``call.py`` -> ``validate.py`` -> ``export_qupath.py``,
with ``metacluster.py``) that the typing README used to describe.  It did not produce the calls on
disk.  The current typing is the subpanel classification in ``subpanel.py`` / ``build_cells.py``;
see ``README.md``, "Retired, and not to be revived as fallbacks".  What follows is kept as a record
of what was measured and why, not as a lane to fall back to.

Two stages, for two different reasons that pull in opposite directions.

*Fine first, because rare cell types vanish otherwise.*  In the validated test window the delta
cells are 44 of 11,659 -- 0.4%.  At a few dozen clusters a population that size is absorbed into a
larger one and disappears, and scattered rare populations are exactly what a study of the immune
environment around islets is about.

*Then merge, because a few dozen is where judgement happens.*  A cluster count in the hundreds or
thousands costs three things: the naming family grows as ``n_clusters x n_lineages`` so the bar
must rise to match (at 1500 clusters that is 13,500 tests, and ~555 clusters would clear some
lineage by chance); each lineage fragments across up to 165 clusters so a "cluster" stops being a
biological unit; and the marker profile table stops being readable, which matters most of all,
because the reviewer looking at these calls in QuPath is the only real check on whether they are
right.  Thirty-five profiles can be read and disagreed with.  Fifteen hundred cannot.

The merge is KINTSUGI's mechanism, unchanged: correlation distance between cluster profiles, Ward
linkage, ``fcluster`` at ``maxclust``.  It works on the profile matrix -- one row per cluster --
so a 44-cell cluster weighs exactly as much as a 40,000-cell one, which is what carries a rare
population through the merge instead of dissolving it.

``n_final`` is **not** a guessed constant.  It is chosen from a stability sweep and the whole sweep
is reported, because a cluster count picked to land near an expected number of cell types is
frequency reasoning.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist


def profiles(Z: np.ndarray, labels: np.ndarray, k: int) -> np.ndarray:
    """Mean marker vector per fine cluster."""
    P = np.full((k, Z.shape[1]), np.nan)
    for c in range(k):
        sel = labels == c
        if sel.any():
            P[c] = np.nanmean(Z[sel], axis=0)
    return P


def merge(P: np.ndarray, n_final: int) -> np.ndarray:
    """Fine cluster -> merged cluster id, 0-based and gapless."""
    ok = np.isfinite(P).all(axis=1)
    if n_final < 2 or n_final > int(ok.sum()):
        raise ValueError(f"n_final={n_final} outside [2, {int(ok.sum())}]")
    Z = linkage(pdist(P[ok], metric="correlation"), method="ward")
    a = fcluster(Z, n_final, criterion="maxclust")
    out = np.full(len(P), -1, dtype=np.int32)
    out[np.where(ok)[0]] = a - 1
    return out


def stability(P: np.ndarray, n_final: int, draws: int, seed: int) -> float:
    """How reproducible the merge is when the profiles are perturbed by their own noise.

    A merge that reorganises under a jitter smaller than the spread between profiles is not a
    structure, whatever it looks like.  The jitter is scaled to each marker's own variation across
    clusters, so it is set by the data rather than chosen.
    """
    from sklearn.metrics import adjusted_rand_score
    ok = np.isfinite(P).all(axis=1)
    base = merge(P, n_final)[ok]
    rng = np.random.default_rng(seed)
    sd = np.nanstd(P[ok], axis=0)
    scores = []
    for _ in range(draws):
        jitter = P.copy()
        jitter[ok] += rng.normal(scale=0.25 * sd, size=(int(ok.sum()), P.shape[1]))
        scores.append(adjusted_rand_score(base, merge(jitter, n_final)[ok]))
    return float(np.mean(scores))


def sweep(P: np.ndarray, lo: int, hi: int, draws: int = 25, seed: int = 0) -> pd.DataFrame:
    """Stability across candidate ``n_final``.  Report all of it; choose within a stated range."""
    rows = []
    for n in range(lo, hi + 1):
        try:
            s = stability(P, n, draws, seed)
        except ValueError:
            continue
        m = merge(P, n)
        sizes = np.bincount(m[m >= 0], minlength=n)
        rows.append({"n_final": n, "stability_ari": round(s, 4),
                     "smallest_merged_cluster": int(sizes.min()),
                     "singletons": int((sizes == 1).sum())})
    return pd.DataFrame(rows)


def choose(sweep_df: pd.DataFrame, lo: int, hi: int, tol: float = 1e-9) -> int:
    """The most stable count inside the reviewable band, ties broken toward MORE clusters.

    The direction of the tie-break is load-bearing and the obvious choice is wrong.  Measured on a
    synthetic set with five common groups and one rare one, ``n_final`` 5 and 6 were equally
    stable (ARI 1.0000 both); breaking toward the smaller count merged the rare group into a
    common one and dropped its purity from 1.00 to 0.11, while the larger count kept it separate
    at no cost to stability.  Within the reviewable band the readability constraint is already
    satisfied, so the only thing a coarser tie buys is lost structure.
    """
    band = sweep_df[(sweep_df.n_final >= lo) & (sweep_df.n_final <= hi)]
    if band.empty:
        raise ValueError(f"no candidate n_final in [{lo}, {hi}]")
    best = band["stability_ari"].max()
    return int(band[band.stability_ari >= best - tol]["n_final"].max())


def rare_population_survives(fine: np.ndarray, merged_of_fine: np.ndarray,
                             is_rare: np.ndarray) -> dict:
    """Did a known rare population keep a cluster of its own through the merge?

    Reported, never optimised against: the point is to see whether the merge dissolved something,
    not to pick the merge that flatters this number.
    """
    merged = merged_of_fine[fine]
    if not is_rare.any():
        return {"n_rare": 0}
    counts = pd.Series(merged[is_rare]).value_counts()
    top = int(counts.index[0])
    in_top = int(counts.iloc[0])
    size_of_top = int((merged == top).sum())
    return {"n_rare": int(is_rare.sum()),
            "best_cluster": top,
            "completeness": round(in_top / int(is_rare.sum()), 4),
            "purity": round(in_top / max(size_of_top, 1), 4),
            "clusters_holding_them": int(len(counts))}
