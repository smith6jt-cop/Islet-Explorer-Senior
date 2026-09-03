"""Biological-recovery (SOFT-reference, DESCRIPTIVE) track.

Describes WHAT biology each method's domains capture, using the pipeline's
celltype labels as a soft reference. This is intentionally NOT a quality verdict:
no proportion-plausibility judgement, no over/under-call language (Tier-C;
memories no-proportion-based-quality-judgments / categorize-by-expressed-genes).
The user reads the enrichment heatmaps and decides.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score


def contingency(domains, celltypes) -> pd.DataFrame:
    return pd.crosstab(pd.Series(np.asarray(domains).astype(str), name="domain"),
                       pd.Series(np.asarray(celltypes).astype(str), name="celltype"))


def enrichment(domains, celltypes):
    """Return (counts, log2_obs_over_exp, pearson_residual) domain x celltype."""
    ct = contingency(domains, celltypes)
    obs = ct.values.astype(float)
    row = obs.sum(1, keepdims=True)
    col = obs.sum(0, keepdims=True)
    tot = obs.sum()
    exp = row @ col / max(tot, 1e-9)
    with np.errstate(divide="ignore", invalid="ignore"):
        log2oe = np.log2((obs + 1e-9) / (exp + 1e-9))
        pearson = (obs - exp) / np.sqrt(exp + 1e-9)
    idx, cols = ct.index, ct.columns
    return (ct, pd.DataFrame(log2oe, index=idx, columns=cols),
            pd.DataFrame(pearson, index=idx, columns=cols))


def label_nmi(domains, celltypes) -> dict:
    d, c = np.asarray(domains), np.asarray(celltypes)
    return {"nmi": float(normalized_mutual_info_score(d, c)),
            "ami": float(adjusted_mutual_info_score(d, c, average_method="arithmetic"))}


def top_celltypes_per_domain(domains, celltypes, top: int = 3) -> dict:
    ct = contingency(domains, celltypes)
    comp = ct.div(ct.sum(axis=1), axis=0)
    out = {}
    for d in comp.index:
        s = comp.loc[d].sort_values(ascending=False).head(top)
        out[str(d)] = [(str(k), float(v)) for k, v in s.items()]
    return out
