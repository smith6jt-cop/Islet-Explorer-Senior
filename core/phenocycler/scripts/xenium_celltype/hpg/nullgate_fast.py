"""The same gate as :mod:`nullgate`, reduced to one sparse pass per grouping.

Not a different statistic and not an approximation. ``nullgate._contrast`` is

    mean over cells in the group of (mean over the core's genes)
      minus the same outside the group

and the two means commute, so the contrast of a gene set is exactly the mean, over that
set's genes, of a per-gene quantity that does not depend on which set is being scored:

    delta[g] = mean of gene g inside the group  -  mean of gene g outside it
    contrast(core) = delta[core].mean()

Verified against the reference implementation to 1.3e-18 on a random sparse matrix.

Why it matters is arithmetic, not elegance. The reference evaluates one column-slice-and-mean
over the whole cell x gene matrix per (group, lineage, draw). The production run wants
k up to a few hundred groups, 8 lineages and 200 null draws over 1.1M cells: ~320,000 passes
over a 1.1M x 5,101 matrix, which is not a slow analysis, it is an analysis that does not
finish. Here the groups partition the cells, so every per-gene group sum comes from ONE
sparse product ``G @ X`` with a k x n indicator, and each of the 320,000 draws becomes a mean
over ~100 entries of a length-5,101 vector.

The consequence is scientific, not just operational: the draw count and the cluster
granularity stop being things to economise on. 200 draws puts the null's 95th percentile on
10 draws a side instead of 50 draws' 2.5, and k can go fine enough that a rare immune
population is its own group instead of a minority inside a large one.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# The reference lives one directory up. Importing it rather than restating its constants is
# what keeps MIN_GROUP_CELLS, the reasons and the GroupCall shape from drifting apart between
# the two implementations while the tests still pass on each of them separately.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from nullgate import (DEFAULT_DRAWS, MIN_GROUP_CELLS, GroupCall, NullGateError,
                      build_null_pools, decide_group, detection_deciles, match_strata,
                      resolve_cores, standardise_genes)


def group_gene_deltas(X, groups: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-group, per-gene in-minus-out mean, from a single sparse product.

    Returns ``(group_labels, sizes, delta)`` with ``delta`` shaped (n_groups, n_genes).
    """
    X = sp.csr_matrix(X)
    n_cells, n_genes = X.shape
    if groups.shape[0] != n_cells:
        raise NullGateError(f"groups has {groups.shape[0]} entries for {n_cells} cells")
    labels, codes = np.unique(groups, return_inverse=True)
    k = len(labels)
    sizes = np.bincount(codes, minlength=k).astype(np.int64)
    if (sizes == n_cells).any():
        raise NullGateError("a group holds every cell; the contrast has no outside")

    G = sp.csr_matrix((np.ones(n_cells), (codes, np.arange(n_cells))), shape=(k, n_cells))
    sums = np.asarray((G @ X).todense()) if sp.issparse(G @ X) else np.asarray(G @ X)
    total = np.asarray(X.sum(axis=0)).ravel()

    inside = sums / sizes[:, None]
    outside = (total[None, :] - sums) / (n_cells - sizes)[:, None]
    return labels, sizes, inside - outside


def annotate_groups(
    X,
    genes,
    groups: np.ndarray,
    core_genes,
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 0,
    min_cells: int = MIN_GROUP_CELLS,
) -> pd.DataFrame:
    """Drop-in replacement for :func:`nullgate.annotate_groups` with identical semantics.

    Every scientific decision is the reference's, imported rather than restated:
    ``standardise_genes`` so no single gene carries a core, ``match_strata`` so the null is
    matched on detection AND magnitude, ``build_null_pools`` so a gene never stands in for
    itself, and ``decide_group`` for the award -- clearing the null, reading higher inside
    the group than outside, and being separable from the runner-up. The draws are the same
    draws, because the RNG is consumed in the same order. Only the arithmetic differs.

    **This file drifted once and it matters how.** It was written before the reference gained
    the negative-contrast and separability corrections, kept ``real > p95`` with no paired
    null, and therefore awarded lineages on groups the reference refuses -- in the production
    HiPerGator path, on every donor, silently. ``tests/test_nullgate_fast.py`` caught it and
    was red. Nothing here may restate a rule the reference owns.
    """
    X = standardise_genes(X)
    dec = match_strata(X)
    coarse = detection_deciles(X)
    rng = np.random.default_rng(seed)

    resolved = resolve_cores(genes, core_genes)
    pools = build_null_pools(genes, resolved, dec, coarse)

    labels, sizes, delta = group_gene_deltas(X, groups)
    order = {lab: i for i, lab in enumerate(labels)}

    out: list[GroupCall] = []
    for g in pd.unique(groups):
        i = order[g]
        n = int(sizes[i])
        if n < min_cells:
            out.append(GroupCall(g, n, "", "too_few_cells", np.nan, np.nan, np.nan, np.nan, ""))
            continue
        d = delta[i]
        stats, nulls = {}, {}
        for lin, cols in resolved.items():
            cols = np.asarray(cols)
            real = float(d[cols].mean())
            # Drawn inside the group loop, lineage by lineage, one index per core gene per
            # draw -- the reference's exact order of rng.choice calls. Matching the RNG
            # consumption is what lets tests/test_nullgate_fast.py assert the two
            # implementations agree frame for frame rather than merely agree in distribution.
            picks = np.stack([[rng.choice(pools[lin][c]) for c in cols]
                              for _ in range(n_draws)])
            null = d[picks].mean(axis=1)
            sd = float(null.std(ddof=1))
            p95 = float(np.percentile(null, 95))
            z = (real - float(null.mean())) / sd if sd > 0 else np.nan
            stats[lin] = (real, float(null.mean()), p95, z)
            nulls[lin] = null

        out.append(decide_group(g, n, stats, nulls))
    return pd.DataFrame([vars(c) for c in out])
