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


def _contrast(X, cols: Sequence[int], in_group: np.ndarray) -> float:
    """How much higher this core reads INSIDE the group than outside it.

    A contrast, not a level. Scoring the level would award any core that happens to be
    brightly expressed — a core detected in every cell reads high in every group while
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
    dec = detection_deciles(X)
    by_dec = {d: np.where(dec == d)[0] for d in np.unique(dec)}
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

    out: list[GroupCall] = []
    for g in pd.unique(groups):
        rows = np.where(groups == g)[0]
        if len(rows) < min_cells:
            out.append(GroupCall(g, len(rows), "", "too_few_cells",
                                 np.nan, np.nan, np.nan, np.nan, ""))
            continue

        in_group = np.zeros(n_cells, dtype=bool)
        in_group[rows] = True
        cleared, stats = [], {}
        for lin, cols in resolved.items():
            real = _contrast(X, cols, in_group)
            null = np.array([
                _contrast(X, [rng.choice(by_dec[dec[c]]) for c in cols], in_group)
                for _ in range(n_draws)
            ])
            sd = float(null.std(ddof=1))
            p95 = float(np.percentile(null, 95))
            z = (real - float(null.mean())) / sd if sd > 0 else np.nan
            stats[lin] = (real, float(null.mean()), p95, z)
            if real > p95:
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
        out.append(GroupCall(
            g, len(rows), win,
            "" if len(cleared) == 1 else "multiple_cores_cleared",
            r, nm, p95, z, "|".join(cleared[1:]),
        ))

    return pd.DataFrame([vars(c) for c in out])
