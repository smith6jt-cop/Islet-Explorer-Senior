"""Cluster -> lineage annotation, gated on a size- and detection-matched null.

This module exists because the obvious thing does not work. A per-cell argmax over
``identity_core_scores`` was measured on donor 6539 (50,000 cells) against random gene sets
matched on core size *and* per-gene detection decile, 10 draws:

    lineage     genes   REAL%   null mean   null p5-p95
    Endocrine      18    2.83        5.84   4.49- 7.66   <- BELOW its null
    Exocrine        8   75.46       61.18  38.68-77.76
    Immune_B       58    1.91        6.74   2.65-12.87   <- BELOW its null
    Immune_T      129    1.08        5.09   0.75-18.50
    Myeloid       111    2.78        3.08   1.53- 4.56
    Neural         17    1.23        1.42   0.95- 1.92
    Stromal        48    7.02        4.36   2.47- 6.37
    Vascular       87    7.68       12.30   5.51-34.87

Five of eight sit inside the matched null and two do worse than random: at ~96% per-cell
zeros the argmax reports which core happens to hold a high expresser, which is a property a
random size-matched set has too. So the null is not a post-hoc check here. It is the gate.

Two consequences are built in rather than documented:

* **Annotation happens on groups, not cells.** Pooling hundreds of cells is what buys the
  signal-to-noise a single 316-transcript cell does not have. The caller supplies the
  grouping (a KODAMA embedding clustered, say); this module never invents one.
* **A lineage is only awarded when the real core beats its own matched null.** A group whose
  best core is inside the null gets ``Unavailable`` with a recorded reason, never the nearest
  lineage. That is the difference between a measurement and an argmax.
* **The statistic is a contrast, not a level** — how much higher the core reads inside the
  group than outside. Scoring the level would hand the award to whichever core is brightest;
  a core detected in every cell reads high everywhere while distinguishing nothing, and its
  contrast is ~0 in every group. A unit test pins this.

The matching is on the detection decile because that is the confound: a core of ubiquitous
genes scores high everywhere regardless of identity. Matching size alone would leave that in.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

#: Draws of the matched null per lineage. 50 puts the 5th/95th percentile on ~2.5 draws a
#: side, which is the coarsest that still resolves "inside" from "outside" honestly.
DEFAULT_DRAWS = 50

#: A group smaller than this is not annotated at all. Pooling is the entire mechanism, and a
#: handful of cells does not pool. Reported as ``too_few_cells``, never as a weak call.
MIN_GROUP_CELLS = 50

#: Fewest genes a matched stratum must offer, after the identity cores are removed from it,
#: before the match is widened to the detection decile alone. Matching on detection AND
#: magnitude makes strata narrow, and a narrow stratum on a small panel can leave a pool of
#: one -- which is not a null.
MIN_NULL_CANDIDATES = 8


class NullGateError(RuntimeError):
    """The gate cannot be evaluated — never downgraded to a pass."""


@dataclass(frozen=True)
class GroupCall:
    group: object
    n_cells: int
    lineage: str          # "" when nothing cleared its null
    reason: str           # "" only for an awarded lineage
    real: float           # winning core's in-group minus out-of-group mean
    null_mean: float
    null_p95: float
    z: float              # (real - null_mean) / null_sd, for ranking only
    runners_up: str       # other lineages that also cleared, pipe-joined


def detection_deciles(X, n_bins: int = 10) -> np.ndarray:
    """Per-gene detection rate binned into equal-count deciles.

    Ranked before binning so that ties (very common at this sparsity — thousands of genes
    share a detection rate) cannot collapse several deciles into one.
    """
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    ranks = pd.Series(det).rank(method="first")
    return pd.qcut(ranks, n_bins, labels=False).to_numpy()


def _gene_sd(X) -> np.ndarray:
    """Per-gene standard deviation, computed without densifying a sparse matrix."""
    n = X.shape[0]
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = X.multiply(X) if hasattr(X, "multiply") else np.square(X)
    mean_sq = np.asarray(sq.mean(axis=0)).ravel()
    return np.sqrt(np.maximum(mean_sq - mean ** 2, 0.0)) * np.sqrt(n / max(n - 1, 1))


def standardise_genes(X):
    """Scale each gene to unit SD, so no single gene can carry a core.

    Measured on donor 6539: CUZD1 alone supplies 57% of the eight-gene Exocrine core's mean
    and XBP1 39% of the 58-gene Immune_B core's. A core's contrast was therefore largely one
    gene's contrast, and no random set of the same size shares that structure — so the null
    could not see it.

    Only the SCALE is applied, never the centring, and that is exact rather than a
    concession to sparsity: the statistic is a difference of means, ``inside - outside``,
    so a per-gene offset cancels between the two terms. Dividing by the SD alone therefore
    gives precisely the z-scored contrast while leaving a sparse matrix sparse.
    """
    sd = _gene_sd(X)
    inv = np.where(sd > 0, 1.0 / np.where(sd > 0, sd, 1.0), 0.0)
    if hasattr(X, "multiply"):
        import scipy.sparse as sp
        return X @ sp.diags(inv)
    return np.asarray(X) * inv[None, :]


def match_strata(X, n_bins: int = 10) -> np.ndarray:
    """Gene strata matched on BOTH detection and magnitude, not detection alone.

    Matching on detection alone leaves the confound it was meant to remove. Measured on
    donor 6539, the top detection decile spans 0.120–0.874 detection and 0.075–2.069 mean
    lognorm — a 28x range in magnitude — and the genes that actually carry the cores sit at
    its very top: XBP1 at the 97.8th percentile of its own decile, CUZD1 at the 99.4th. A
    "matched" draw therefore replaced each core gene with one ~8–13x weaker, and the null was
    systematically too easy for exactly the cores under test.
    """
    det = np.asarray((X > 0).mean(axis=0)).ravel()
    mag = np.asarray(X.mean(axis=0)).ravel()
    d = pd.qcut(pd.Series(det).rank(method="first"), n_bins, labels=False).to_numpy()
    m = pd.qcut(pd.Series(mag).rank(method="first"), n_bins, labels=False).to_numpy()
    return d * n_bins + m


def _contrast(X, cols: Sequence[int], in_group: np.ndarray) -> float:
    """How much higher this core reads INSIDE the group than outside it.

    KNOWN LIMIT, measured 2026-08-27 -- this statistic is blind to the MAJORITY compartment,
    and that is structural rather than a tuning problem. The comparison is against the
    complement, so when a compartment is most of the tissue, "outside" is also mostly that
    compartment and the difference collapses however unambiguous the identity is.

    On donor 6539 at k=30 the nine largest Unavailable groups (2,189-3,193 cells each) carry
    acinar (CUZD1/CTRL) z of 1.51-1.82 -- the highest in the donor bar one -- against 0.43
    for the awarded groups and 0.99 cohort-wide. Group 29 at acinar z=1.878 was awarded
    Exocrine while group 2 at 1.823 was refused: near-identical identity, opposite outcome,
    decided by which side of the mostly-acinar remainder each happened to fall.

    So this gate answers "is this group ENRICHED for lineage L relative to the rest of the
    tissue?", which is a real question, and NOT "what is this group?". For minority
    compartments -- immune, vascular, neural, endocrine, which is what the islet work needs
    -- the two coincide. For the dominant compartment they do not, and no composition
    denominator computed from these calls is meaningful until that is addressed.
    A contrast, not a level. Scoring the level would award any core that happens to be
    brightly expressed -- a core detected in every cell reads high in every group while
    distinguishing none of them, which is the opposite of identity. The difference of means
    is ~0 for such a core in every group, so it clears no null anywhere.
    """
    if len(cols) == 0:
        raise NullGateError("a core with no gene on this panel reached the gate")
    per_cell = np.asarray(X[:, list(cols)].mean(axis=1)).ravel()
    inside = per_cell[in_group]
    outside = per_cell[~in_group]
    if inside.size == 0 or outside.size == 0:
        raise NullGateError("contrast needs cells on both sides of the group boundary")
    return float(inside.mean() - outside.mean())


def annotate_groups(
    X,
    genes: Sequence[str],
    groups: np.ndarray,
    core_genes: Mapping[str, Sequence[str]],
    *,
    n_draws: int = DEFAULT_DRAWS,
    seed: int = 0,
    min_cells: int = MIN_GROUP_CELLS,
) -> pd.DataFrame:
    """Award each group a lineage only where the real core beats its matched null.

    ``X`` is cells x genes (sparse or dense) on whatever scale the caller chose; the gate is
    a comparison of the real core against random cores on that same scale, so it inherits the
    caller's normalisation rather than imposing one.
    """
    genes = list(map(str, genes))
    index = {g: i for i, g in enumerate(genes)}
    X = standardise_genes(X)                       # no gene may carry a core; see above
    dec = match_strata(X)                          # matched on detection AND magnitude
    coarse = detection_deciles(X)                  # the fallback when a stratum runs thin
    rng = np.random.default_rng(seed)

    n_cells = X.shape[0]
    resolved = {}
    for lin, gs in core_genes.items():
        cols = [index[g] for g in map(str, gs) if g in index]
        if not cols:
            raise NullGateError(
                f"lineage {lin!r} has no identity-core gene on this panel; it cannot be "
                f"scored and must be declared absent rather than silently skipped"
            )
        resolved[lin] = cols

    # A gene must never stand in for itself. Drawing the replacement from the same stratum
    # can otherwise return the core's own gene -- certainly on a small panel, where a 10x10
    # joint stratification leaves only a handful of genes per stratum -- and the "null" then
    # contains the very signal it exists to be compared against.
    #
    # The exclusion is per core, not global. A gene belonging to some OTHER core is a
    # perfectly good matched draw; excluding every core's genes from every pool is what
    # empties a stratum that a bright core happens to monopolise. Where the fine stratum runs
    # thin it widens to the detection decile rather than drawing from a pool of one.
    def _pool(col: int, own: frozenset) -> np.ndarray:
        fine = np.array([i for i in np.where(dec == dec[col])[0] if i not in own], dtype=int)
        if fine.size >= MIN_NULL_CANDIDATES:
            return fine
        wide = np.array([i for i in np.where(coarse == coarse[col])[0] if i not in own],
                        dtype=int)
        if wide.size == 0:
            raise NullGateError(
                f"no gene outside this core shares {genes[col]!r}'s detection decile, so "
                "its null cannot be matched"
            )
        return wide

    pools = {
        lin: {c: _pool(c, frozenset(cols)) for c in cols}
        for lin, cols in resolved.items()
    }

    out: list[GroupCall] = []
    for g in pd.unique(groups):
        rows = np.where(groups == g)[0]
        if len(rows) < min_cells:
            out.append(GroupCall(g, len(rows), "", "too_few_cells",
                                 np.nan, np.nan, np.nan, np.nan, ""))
            continue

        in_group = np.zeros(n_cells, dtype=bool)
        in_group[rows] = True
        cleared, stats, nulls = [], {}, {}
        for lin, cols in resolved.items():
            real = _contrast(X, cols, in_group)
            null = np.array([
                _contrast(X, [rng.choice(pools[lin][c]) for c in cols], in_group)
                for _ in range(n_draws)
            ])
            sd = float(null.std(ddof=1))
            p95 = float(np.percentile(null, 95))
            z = (real - float(null.mean())) / sd if sd > 0 else np.nan
            stats[lin] = (real, float(null.mean()), p95, z)
            nulls[lin] = null
            # A core must beat its null AND read higher inside the group than outside it.
            # Without the second clause a group reading LOWER on a core than the rest of the
            # tissue was awarded it whenever the random cores read lower still: measured,
            # nonspatial k=30 group 11 took Vascular at real=-0.0134 and spatial k=30
            # group 12 took Myeloid at real=-0.0045 across 833 cells.
            if real > max(p95, 0.0):
                cleared.append(lin)

        if not cleared:
            best = max(stats, key=lambda k: (stats[k][3] if np.isfinite(stats[k][3]) else -np.inf))
            r, nm, p95, z = stats[best]
            out.append(GroupCall(g, len(rows), "", "no_core_beats_its_matched_null",
                                 r, nm, p95, z, ""))
            continue

        cleared.sort(key=lambda k: stats[k][3], reverse=True)
        win = cleared[0]
        r, nm, p95, z = stats[win]

        # Clearing is not winning. Two cores can both beat their nulls while being
        # indistinguishable from each other, and ranking them by z then hands the group
        # entirely to one: measured on donor 6539, Exocrine at z=3.597 and Immune_B at
        # z=3.539 -- a gap of 0.06 -- decided a 100%-vs-0% assignment, and four groups
        # swapped between k=20 and k=30 on gaps that size. The paired null says whether the
        # gap means anything: for each draw, the two cores' random contrasts differ by some
        # amount, and the real gap has to beat the 95th percentile of that.
        if len(cleared) > 1:
            runner = cleared[1]
            paired = nulls[win] - nulls[runner]
            gap = r - stats[runner][0]
            if gap <= float(np.percentile(paired, 95)):
                out.append(GroupCall(
                    g, len(rows), "", "cores_cleared_but_are_not_separable",
                    r, nm, p95, z, "|".join(cleared),
                ))
                continue

        out.append(GroupCall(
            g, len(rows), win,
            "" if len(cleared) == 1 else "multiple_cores_cleared",
            r, nm, p95, z, "|".join(cleared[1:]),
        ))

    return pd.DataFrame([vars(c) for c in out])
