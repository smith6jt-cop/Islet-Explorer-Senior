"""Stage 3 — detection conditioned on what the cell is, not on one threshold for everything.

A uniform per-cell rule was built, measured against a ground-truth negative set, and rejected
as the decision rule (`detect.py`). Every mode traded background for sensitivity almost
proportionally, because no per-cell threshold can separate a real single count from a
background one when 81.5% of nonzero entries ARE single counts. The information a uniform
rule discards is prior knowledge about which genes a given cell can plausibly express.

    P(real | k, g, c)  proportional to   pi(g|k_c) * P(k | expressed)
                       versus        (1 - pi(g|k_c)) * P(k | lambda_c)

``lambda_c`` is Stage 1's per-cell gDNA background rate; ``pi(g|k_c)`` is the probability that
gene g is genuinely expressed by a cell of class ``k_c``. The threshold then varies by
(gene, cell) instead of being one number: an exocrine transcript in a beta cell needs far more
evidence than an endocrine transcript in the same cell, which is exactly what a uniform rule
throws away.

**Why this is not circular.** The classes come from the PROTEIN arm -- 70-plex on the same
physical section, registered DAPI-to-DAPI at a measured 0.20 um median residual -- so the
partition of cells is external to the transcript counts it is used to interpret. Estimating
pi from RNA *within* protein-defined classes therefore borrows no information from the RNA
about which class a cell is in.

**What it cannot do**, stated because the estimator is empirical Bayes and the limitation is
structural: pi is estimated from the same cells it is applied to, so a gene whose background
exactly mimics its expression in a class is indistinguishable there. The prior sharpens a
decision; it does not manufacture evidence.

**One context variable, not several.** Class is the protein call alone. Islet membership is
already inside it -- the endocrine calls exist only inside islets by construction of the
cascade -- and stacking a second prior on top would double-count the same fact.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import poisson

#: A class must have at least this many cells before it gets its own prior. Below it, the
#: pooled prior is used and the cell is recorded as having had no usable context -- a small
#: class would otherwise estimate pi from a handful of cells and state it with confidence.
MIN_CLASS_CELLS = 200

#: Purity below this and the protein call is not a usable context: the Xenium cell straddles
#: several PhenoCycler cells that disagree, so its "class" is an artefact of the overlap.
#: These cells fall back to the pooled prior rather than being dropped -- they still have
#: counts, and a pooled prior is exactly the uniform assumption, which is the honest default.
MIN_CALL_PURITY = 0.8

#: Floor on the expressed-state mean, so a class where a gene is background-only still has a
#: defined alternative hypothesis rather than a degenerate one.
MU_FLOOR = 0.05

POOLED = "__pooled__"


class ContextError(ValueError):
    pass


def class_labels(calls: pd.DataFrame, *, min_purity: float = MIN_CALL_PURITY,
                 min_cells: int = MIN_CLASS_CELLS) -> pd.Series:
    """Protein call per cell, blanked to POOLED where it is not usable as context.

    Two ways a call fails to be context: the transfer was impure (the Xenium cell overlaps
    PhenoCycler cells that disagree), or the class is too small to estimate a prior from.
    Both fall back to POOLED rather than being dropped.
    """
    lab = calls["call"].astype(str).where(
        calls["call_purity"].fillna(0.0) >= min_purity, POOLED)
    counts = lab.value_counts()
    small = set(counts[counts < min_cells].index) - {POOLED}
    return lab.where(~lab.isin(small), POOLED)


def background_detection_rate(lam: np.ndarray) -> float:
    """Mean P(at least one background count) over a set of cells."""
    return float(np.mean(1.0 - np.exp(-np.asarray(lam, dtype=float))))


def class_tables(counts: sp.csr_matrix, lam: np.ndarray, labels: pd.Series
                 ) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Per (class, gene): the expressing prior pi, and the expressed-state mean mu.

    ``pi`` is the background-CORRECTED detection rate:

        r_obs = pi * P(detect | expressed) + (1 - pi) * r_bg

    with ``P(detect | expressed)`` taken as 1, which makes pi a LOWER bound -- a gene that is
    expressed but often undetected reads as less prevalent than it is. Conservative in the
    right direction: it understates the prior rather than inventing one.

    ``mu`` is the background-corrected mean count, floored so the alternative hypothesis is
    always defined.
    """
    X = sp.csr_matrix(counts)
    if X.shape[0] != len(lam) or X.shape[0] != len(labels):
        raise ContextError(f"counts has {X.shape[0]} rows, lambda {len(lam)}, "
                           f"labels {len(labels)}")
    lam = np.asarray(lam, dtype=float)
    det = (X > 0)
    pis, mus, meta = {}, {}, {}
    for cls, idx in labels.groupby(labels).groups.items():
        rows = labels.index.get_indexer(idx)
        r_obs = np.asarray(det[rows].mean(axis=0)).ravel()
        r_bg = background_detection_rate(lam[rows])
        pi = np.clip((r_obs - r_bg) / max(1.0 - r_bg, 1e-9), 0.0, 1.0)
        mean_obs = np.asarray(X[rows].mean(axis=0)).ravel()
        mu = np.maximum(mean_obs - lam[rows].mean(), MU_FLOOR)
        pis[cls], mus[cls] = pi, mu
        meta[cls] = {"n_cells": int(len(rows)), "r_bg": r_bg,
                     "pi_median": float(np.median(pi)),
                     "genes_pi_zero": int((pi <= 0).sum())}
    return pd.DataFrame(pis).T, pd.DataFrame(mus).T, meta


def detected(counts: sp.csr_matrix, lam: np.ndarray, labels: pd.Series,
             pi: pd.DataFrame, mu: pd.DataFrame, *, odds: float = 1.0) -> sp.csr_matrix:
    """Boolean CSR: is each observed count better explained by expression than by background?

    Evaluated only where ``k > 0`` -- a zero count is not a detection under any rule, and
    restricting to the nonzeros keeps this to the ~215M entries that exist rather than the
    3.6 billion that do not.
    """
    X = sp.csr_matrix(counts).tocoo()
    lam = np.asarray(lam, dtype=float)
    cls_of_row = labels.to_numpy()
    order = {c: i for i, c in enumerate(pi.index)}
    ci = np.array([order[c] for c in cls_of_row], dtype=np.int64)

    k = X.data.astype(np.float64)
    p_i, g_i = ci[X.row], X.col
    prior = pi.to_numpy()[p_i, g_i]
    mu_e = mu.to_numpy()[p_i, g_i]
    lam_c = lam[X.row]

    # Poisson on both sides: background at lambda_c, expression at mu. Both are counts of
    # the same kind, so the likelihood ratio is the only thing that differs between them.
    ll_expressed = poisson.pmf(k, mu_e) * prior
    ll_background = poisson.pmf(k, lam_c) * (1.0 - prior)
    keep = ll_expressed > odds * ll_background
    return sp.csr_matrix(
        (np.ones(int(keep.sum()), bool), (X.row[keep], X.col[keep])), shape=X.shape)
