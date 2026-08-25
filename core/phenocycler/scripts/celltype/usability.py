"""Per-donor, per-marker statistics that do not assume a scale.

Every statistic here is computed WITHIN a donor and is invariant to any monotone
per-donor rescaling of a marker, because the thing we do not yet have is a common scale
across donors.  That rules out means, variances and anything compared to a fixed number,
and it rules in ranks, tails and spatial structure.

Four families:

* **shape** -- quantiles, the fraction reading exactly zero, and ``top1_signal_share``: the
  share of a marker's total signal carried by its brightest 1% of cells.  DAPI, which is on
  every nucleus, sits at ~1.9% and is the reference for "flat".
* **tail overlap** -- how much more often than chance a marker's brightest cells are also
  its declared partner's brightest cells, reported beside the same statistic for declared
  NON-partners.  An enrichment with no control beside it is not a measurement.
* **spatial coherence** -- Moran's I over the cell centroids.  Real staining is spatially
  clustered because cells of a type sit together; noise is not.  Computed twice: on ranks
  (overall structure) and on the top-1% indicator (structure of the bright tail, which is
  the only part a sparse marker has).  This is what separates "flat because it is on
  everything" from "flat because it is on nothing" -- ``top1_signal_share`` cannot.
* **acquisition-group structure** -- whether a marker's level steps with the acquisition
  group (calibration epoch x panel variant) rather than with donors.  A rigid shift models
  a per-donor scale difference; it does not model a different antibody clone, so a marker
  that steps with the group must be flagged whatever its shift statistic looks like.
"""
from __future__ import annotations

import sys
import pathlib

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "dct"))
from gates import Verdict, require_evidence          # noqa: E402  (reused, never copied)

__all__ = ["shape_stats", "tail_overlap", "neighbour_graph", "morans_i", "morans_null",
           "spatial_stats", "acquisition_group_ratio", "judge_marker", "MIN_CELLS", "BARS"]

#: Below this many usable cells a donor-marker is unjudged, never passed.
MIN_CELLS = 10_000

#: Pre-registered bars.  Fixed before the first full run, and stated as properties of the
#: statistic rather than as counts over an enumerated set, so they do not change meaning if
#: the marker or partner list grows.
BARS = {
    "partner_min": 3.0,        # a partner tail must beat chance by this much ...
    "control_max": 2.0,        # ... and its controls must stay below this
    "partner_over_control": 3.0,   # and the partner must beat its worst control by this
    "morans_null_draws": 10,   # permutation draws per probe marker, for the spatial null
    "morans_null_probes": 3,   # probe markers spanning the concentration range
    "morans_null_q": 0.99,
    "acq_ratio_flag": 2.0,     # above this the level steps with the acquisition group
    "zero_frac_max": 0.50,     # more than half the cells reading exactly 0 is not a measurement
    "tail_q": 0.001,           # the tail the partner statistic is computed on
    "morans_k": 12,
}


# --------------------------------------------------------------------------- shape
def shape_stats(x: np.ndarray) -> dict:
    """Quantiles, zero/null fractions and the brightest-1% signal share."""
    x = np.asarray(x, dtype=float)
    n_total = len(x)
    ok = np.isfinite(x)
    v = x[ok]
    if len(v) == 0:
        return {"n": 0, "null_frac": 1.0, "zero_frac": np.nan, "top1_signal_share": np.nan}
    qs = np.percentile(v, [10, 25, 50, 75, 90, 99, 99.9])
    s = np.sort(v)[::-1]
    total = s.sum()
    k = max(1, len(s) // 100)
    return {
        "n": int(len(v)),
        "null_frac": float(1.0 - len(v) / max(n_total, 1)),
        "zero_frac": float((v == 0).mean()),
        **{f"p{q}": float(a) for q, a in zip((10, 25, 50, 75, 90, 99, 99.9), qs)},
        "top1_signal_share": float(s[:k].sum() / total) if total > 0 else np.nan,
    }


# ---------------------------------------------------------------------- tail overlap
def tail_overlap(a: np.ndarray, b: np.ndarray, q: float = 0.001) -> float:
    """Fold-enrichment over chance that the top ``q`` of ``a`` is also the top ``q`` of ``b``.

    1.0 means independent.  Scale-free in both arguments: it depends only on which cells
    rank highest, so no normalisation is presupposed.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < MIN_CELLS:
        return np.nan
    a, b = a[ok], b[ok]
    n = len(a)
    k = max(20, int(n * q))
    ta = np.argpartition(-a, k)[:k]
    tb = np.argpartition(-b, k)[:k]
    return float(len(np.intersect1d(ta, tb, assume_unique=False)) / (k * (k / n)))


# ------------------------------------------------------------------ spatial coherence
def neighbour_graph(xy: np.ndarray, k: int = 12):
    """Indices of each cell's ``k`` nearest neighbours, excluding itself."""
    tree = cKDTree(np.asarray(xy, dtype=float))
    _, idx = tree.query(xy, k=k + 1, workers=-1)
    return idx[:, 1:]


def morans_i(values: np.ndarray, nbr: np.ndarray) -> float:
    """Moran's I with row-standardised k-nearest-neighbour weights.

    With row-standardised weights the normaliser cancels, leaving
    ``I = sum(z * mean_neighbour(z)) / sum(z**2)``.  0 is no spatial structure;
    the expectation under a random permutation is ``-1/(n-1)``, i.e. ~0 at this scale.
    """
    z = np.asarray(values, dtype=float)
    z = z - z.mean()
    denom = float(z @ z)
    if denom <= 0:
        return np.nan
    return float(z @ z[nbr].mean(axis=1) / denom)


def morans_null(probe_values: list, nbr: np.ndarray, *, draws: int = 10,
                q: float = 0.99, seed: int = 0) -> float:
    """What Moran's I looks like when there is nothing to see.

    Permuting a marker's values across cells destroys spatial structure while keeping the
    value distribution and the neighbour graph exactly as they are, so the spread of |I| over
    permutations is the "no structure" scale for THIS donor at THIS cell count.  Several
    probe markers are used because the null depends mildly on how heavy-tailed the values
    are; the bar is the worst of them.
    """
    rng = np.random.default_rng(seed)
    worst = 0.0
    for v in probe_values:
        v = np.asarray(v, dtype=float)
        draws_i = [abs(morans_i(rng.permutation(v), nbr)) for _ in range(draws)]
        draws_i = [d for d in draws_i if np.isfinite(d)]
        if draws_i:
            worst = max(worst, float(np.quantile(draws_i, q)))
    return worst


def spatial_stats(x: np.ndarray, nbr: np.ndarray) -> dict:
    """Moran's I on ranks (overall) and on the top-1% indicator (the bright tail)."""
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    if ok.sum() < MIN_CELLS or ok.sum() != len(x):
        # nbr indexes the full array; a partial mask would point at the wrong cells
        return {"morans_rank": np.nan, "morans_top1": np.nan}
    r = rankdata(x) / len(x)
    thr = np.quantile(x, 0.99)
    top = (x >= thr).astype(float)
    return {"morans_rank": morans_i(r, nbr), "morans_top1": morans_i(top, nbr)}


# ------------------------------------------------------- acquisition-group structure
def acquisition_group_ratio(level_by_donor: pd.Series, group_by_donor: pd.Series) -> dict:
    """Between-group spread over the worst within-group spread, on donor-level medians.

    A marker whose level is set by the donor gives a ratio near or below 1.  A marker whose
    level steps with the acquisition group -- a different antibody clone or conjugate lot in
    one panel build -- gives a ratio well above it.
    """
    df = pd.DataFrame({"level": level_by_donor, "group": group_by_donor}).dropna()
    df = df[df["level"] > 0]
    if df["group"].nunique() < 2:
        return {"acq_between": np.nan, "acq_within": np.nan, "acq_ratio": np.nan}
    g = df.groupby("group")["level"]
    between = float(g.median().max() / max(g.median().min(), 1e-12))
    within = float(g.apply(lambda s: s.max() / max(s.min(), 1e-12)).max())
    return {"acq_between": between, "acq_within": within,
            "acq_ratio": float(between / max(within, 1e-12))}


# ------------------------------------------------------------------------- the verdict
def judge_marker(row: dict, bars: dict | None = None) -> Verdict:
    """usable / degraded / unusable for one (donor, marker), fail-closed.

    Absent, non-finite or under-powered evidence is a FAILURE, never a pass -- the rule
    ``scripts/dct/gates.py`` exists to enforce.  ``Verdict.rows`` carries one
    ``(check, value, ok, why)`` tuple per check so the reason survives into the output.
    """
    B = dict(BARS if bars is None else bars)
    checks, fatal, soft = [], [], []

    why = require_evidence(int(row.get("n", 0) or 0), MIN_CELLS, "cells")
    checks.append(("n_cells", row.get("n"), why is None, why or ""))
    if why:
        fatal.append(why)

    z = row.get("zero_frac", np.nan)
    ok = np.isfinite(z) and z <= B["zero_frac_max"]
    checks.append(("zero_frac", z, ok, "" if ok else f"{z:.3f} of cells read exactly 0"))
    if not ok:
        fatal.append("too many cells read exactly 0")

    # Spatial structure, judged against this donor's own permutation null rather than a
    # fixed number.  The two statistics answer different questions and a marker only has to
    # pass one: ``morans_rank`` sees structure anywhere in the distribution, ``morans_top1``
    # sees it in the bright tail.  A marker on ~90% of cells (E-cadherin, Keratin_5) has a
    # scattered top 1% and low ``morans_top1`` while being perfectly real -- that is the
    # "flat because it is on everything" case, and failing it here would be wrong.
    null = row.get("morans_null", np.nan)
    mr, mt = row.get("morans_rank", np.nan), row.get("morans_top1", np.nan)
    structured = np.isfinite(null) and ((np.isfinite(mr) and mr > null)
                                        or (np.isfinite(mt) and mt > null))
    checks.append(("spatial_vs_null", max([v for v in (mr, mt) if np.isfinite(v)] or [np.nan]),
                   structured, "" if structured else f"neither Moran's I exceeds the null {null:.4f}"))

    # Partner evidence is only judged where a partner was declared.  A marker with none is
    # not thereby degraded -- it is recorded as unjudged on that axis.
    p, c = row.get("partner_best", np.nan), row.get("control_worst", np.nan)
    has_partner = np.isfinite(p) and np.isfinite(c)
    partner_ok = False
    if has_partner:
        partner_ok = bool(p >= B["partner_min"] and c <= B["control_max"]
                          and p >= B["partner_over_control"] * c)
        checks.append(("partner_vs_control", p / max(c, 1e-9), partner_ok,
                       "" if partner_ok else f"partner {p:.1f}x vs control {c:.1f}x"))
        # Name the clause that fired.  "does not beat its controls" can mean three different
        # things and they call for different follow-up: a weak partner is a marker problem, a
        # hot control is a specificity problem, and a thin ratio is neither on its own.
        if p < c:
            soft.append(f"bright tail overlaps {row.get('control_worst_with', 'a non-partner')} "
                        f"({c:.1f}x) more than its partner {row.get('partner_best_with', '')} ({p:.1f}x)")
        elif not partner_ok:
            if p < B["partner_min"]:
                soft.append(f"partner tail is only {p:.1f}x chance "
                            f"(bar {B['partner_min']:.0f}x)")
            elif c > B["control_max"]:
                soft.append(f"non-partner {row.get('control_worst_with', '')} sits at {c:.1f}x "
                            f"(bar {B['control_max']:.0f}x)")
            else:
                soft.append(f"partner {p:.1f}x is under {B['partner_over_control']:.0f}x "
                            f"its control {c:.1f}x")
    else:
        checks.append(("partner_vs_control", np.nan, None, "no declared partner"))

    # Unusable only when NOTHING shows signal: no spatial structure above the null and no
    # partner tail either.  Either one alone is evidence the marker is measuring something.
    if not structured and not partner_ok:
        fatal.append("no spatial structure above the null and no partner tail")

    ar = row.get("acq_ratio", np.nan)
    ok = np.isfinite(ar) and ar < B["acq_ratio_flag"]
    checks.append(("acq_group", ar, ok,
                   "" if ok else "level steps with the acquisition group, not with donors"))
    if not ok:
        soft.append("reagent_change_between_acquisition_groups")

    if fatal:
        return Verdict(False, checks, "unusable: " + "; ".join(fatal))
    if soft:
        return Verdict(False, checks, "degraded: " + "; ".join(soft))
    return Verdict(True, checks, "usable")
