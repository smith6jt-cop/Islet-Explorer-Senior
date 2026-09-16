#!/usr/bin/env python
"""One detection predicate, aware of each cell's genomic-DNA background.

Four incompatible definitions of "detected" exist across the three typing lanes:
`hierarchical_phenotype` uses raw counts >= 1 (and >= 2 in one place) plus a lognorm family
log1p(1)/log1p(3)/log1p(5); `nullgate` has no per-cell notion at all; `celltype` uses a
scaled value >= 0.10 of the way above a marker's background mode. This module supplies the
single predicate the unified pipeline uses, and every mode is evaluated against a ground
truth rather than argued for.

**The ground truth.** Both Pheno-Xeno donors are female, so the six Y probes have no template
of any kind -- every count they register is background. `TBL1Y`, `USP9Y` and `EIF1AY` read
0.21-0.35% at raw >= 1, which is the non-gDNA floor (optical + decoding). `DDX3Y` and `KDM5D`
read higher (1.47-2.08%) through cross-hybridisation with their highly-expressed X paralogs
`DDX3X` / `KDM5C`, so they measure a different thing and are reported separately.

**Why `raw >= 1` fails here.** At lambda ~ 0.0252 counts per gene per cell, P(>=1 count from
gDNA alone) is 2.49% of all gene-cell pairs -- about half of the median gene's 5.02% detection
rate. Worse, the failure is size-dependent: satisfying `ANCHOR_MIN_GENES = 2` by background
alone has probability 0.016 for the 8-gene Exocrine core and 0.834 for the 129-gene Immune_T
core, a 50x differential that systematically favours large cores. That is a mechanism for the
independently measured failure in `nullgate`, where the small cores came out below their null.

**The modes.**

``raw``      counts >= floor. The current behaviour, kept so the regression is measurable.
``count``    the same with a higher floor; ``floor=2`` drops background detection 80x.
``poisson``  per cell, the smallest k with ``P(K >= k | lambda_c) < alpha``. This is the only
             mode that adapts: a small clean cell (lambda at p5 = 0.0038) keeps k=1 because
             one count is already unlikely there, while a large background-heavy cell (p95 =
             0.0655) is required to show more. Sensitivity is spent where it buys something.
``lognorm``  value >= thr on the normalised layer. Included because it is measurably
             background-resistant on its own -- on the zero-template Y probes, raw >= 1 reads
             1.88% where lognorm >= log1p(1) reads 0.55%, a 3.4x reduction, since
             `normalize_total` scales a stray count below 1 in any cell of above-median depth
             and gDNA scales with cell size.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import poisson

#: Per-cell false-positive rate a Poisson call is allowed. One gene, one cell.
DEFAULT_ALPHA = 0.01
MODES = ("raw", "count", "poisson", "lognorm")


class DetectError(ValueError):
    pass


def poisson_min_counts(lam: np.ndarray, alpha: float = DEFAULT_ALPHA,
                       k_max: int = 8) -> np.ndarray:
    """Smallest k per cell with ``P(K >= k | lambda) < alpha``.

    Vectorised over cells and capped: lambda here is ~0.025, so k is 1 or 2 for essentially
    every cell, and the cap only guards a pathological background.
    """
    lam = np.asarray(lam, dtype=float)
    need = np.full(lam.shape, k_max, dtype=np.int16)
    for k in range(k_max, 0, -1):
        # P(K >= k) = sf(k-1)
        need = np.where(poisson.sf(k - 1, lam) < alpha, k, need)
    return np.maximum(need, 1).astype(np.int16)


def detected(matrix, *, mode: str = "poisson", gdna_lambda: np.ndarray | None = None,
             floor: int = 1, alpha: float = DEFAULT_ALPHA, thr: float | None = None):
    """Boolean sparse matrix of per-(cell, gene) detections.

    ``matrix`` is counts for every mode except ``lognorm``, where it is the normalised layer.
    Returns a CSR boolean matrix of the same shape, so callers can sum over gene sets.
    """
    if mode not in MODES:
        raise DetectError(f"mode must be one of {MODES}, got {mode!r}")
    X = sp.csr_matrix(matrix)

    if mode == "raw":
        return (X >= max(int(floor), 1)).astype(bool)
    if mode == "count":
        if floor < 2:
            raise DetectError("mode='count' exists to raise the floor; use mode='raw' for 1")
        return (X >= int(floor)).astype(bool)
    if mode == "lognorm":
        if thr is None:
            raise DetectError("mode='lognorm' needs thr, e.g. np.log1p(1.0)")
        return (X >= float(thr)).astype(bool)

    # poisson: a per-CELL threshold, broadcast down the row
    if gdna_lambda is None:
        raise DetectError("mode='poisson' needs gdna_lambda (obs['gdna_lambda'])")
    need = poisson_min_counts(gdna_lambda, alpha)
    if len(need) != X.shape[0]:
        raise DetectError(f"gdna_lambda has {len(need)} entries for {X.shape[0]} cells")
    out = X.tocoo()
    keep = out.data >= need[out.row]
    return sp.csr_matrix((np.ones(keep.sum(), bool), (out.row[keep], out.col[keep])),
                         shape=X.shape)


def detection_rate(det: sp.csr_matrix, cols: np.ndarray | None = None) -> np.ndarray:
    """Per-gene fraction of cells detecting it."""
    d = det if cols is None else det[:, cols]
    return np.asarray(d.mean(axis=0)).ravel()


def core_satisfied(det: sp.csr_matrix, cols, min_genes: int = 2) -> np.ndarray:
    """Per cell: does it detect >= min_genes of this gene set?"""
    return np.asarray(det[:, cols].sum(axis=1)).ravel() >= int(min_genes)
