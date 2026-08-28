"""Cross-method agreement + stability for the domain track.

Reuses the exact sklearn metrics the repo already uses in
scripts/validate_hierarchical_pilot.py (ARI + AMI). AMI is the k-robust
primary; ARI is reported alongside but is granularity-sensitive (penalises
methods that return more/fewer domains), so always read it with n_domains.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
)


def _ami(a, b):
    return adjusted_mutual_info_score(a, b, average_method="arithmetic")


def pairwise_matrices(label_dict: dict) -> dict:
    """{method: 1d-labels (aligned)} -> {'ari','ami','nmi': DataFrame}."""
    methods = list(label_dict)
    mats = {m: pd.DataFrame(index=methods, columns=methods, dtype=float)
            for m in ("ari", "ami", "nmi")}
    for i, a in enumerate(methods):
        for b in methods[i:]:
            x, y = np.asarray(label_dict[a]), np.asarray(label_dict[b])
            vals = {"ari": adjusted_rand_score(x, y), "ami": _ami(x, y),
                    "nmi": normalized_mutual_info_score(x, y)}
            for m, v in vals.items():
                mats[m].loc[a, b] = mats[m].loc[b, a] = v
    return mats


def mean_offdiag_agreement(mat: pd.DataFrame) -> pd.Series:
    """Per-method mean agreement with the OTHER methods (consensus-ness)."""
    m = mat.copy().astype(float)
    np.fill_diagonal(m.values, np.nan)
    return m.mean(axis=1, skipna=True)


def stability(seed_label_sets: dict) -> dict:
    """{method: {seed: labels}} -> {method: {'mean_ari','sd_ari','n_seeds'}}.

    Pairwise ARI across the available seeds for each method (same cells, so
    directly comparable). Deterministic methods -> ARI 1.0.
    """
    out = {}
    for method, by_seed in seed_label_sets.items():
        seeds = sorted(by_seed)
        if len(seeds) < 2:
            out[method] = {"mean_ari": float("nan"), "sd_ari": float("nan"),
                           "n_seeds": len(seeds)}
            continue
        aris = []
        for i in range(len(seeds)):
            for j in range(i + 1, len(seeds)):
                aris.append(adjusted_rand_score(
                    np.asarray(by_seed[seeds[i]]), np.asarray(by_seed[seeds[j]])))
        out[method] = {"mean_ari": float(np.mean(aris)),
                       "sd_ari": float(np.std(aris)), "n_seeds": len(seeds)}
    return out
