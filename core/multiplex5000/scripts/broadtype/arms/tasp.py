"""Arm H ``tasp`` -- tail-anchored seed propagation over a within-donor rank kNN graph.

WHY THIS ARM EXISTS.  Every other table arm carries a model of what a "positive" cell looks
like (a signature, a mixture, a trained network).  This one carries none.  The only thing it
asserts is that the cells at the extreme top of a declared marker's distribution -- the Step-1
tail the whole lane already stands on (``gate.tail_mask``, k = max(20, 0.001 n)) -- are of that
marker's lineage, and it lets those few cells label their neighbours in marker space.  There
is no threshold on any marker, no prior on any label's share, no class weighting and no
prevalence term: a label's final count is whatever the neighbourhood structure of THIS donor's
cells hands it.  The seed count is fixed by the tail definition; the called share is not.

Method (frozen 2026-09-05; any change is a new arm name)
  features     F = ``arm_base.within_donor_rank(ctx.X)`` -- per-marker rank in [0, 1] over the
               loaded cells, so no marker's dynamic range dominates the euclidean distance.
  seeds        tail_L = union over m in M_L ∩ ctx.markers of ``gate.tail_mask(ctx.X[m])`` on the
               ORIGINAL aligned values.  A cell seeds L iff it is in tail_L and in NO other
               label's tail; a cell in two labels' tails is not a seed and is propagated like
               any other cell.  A label with no seeds is absent from the output; a label with
               fewer than 20 seeds is recorded (``seeds_below_min``) and still propagates.
  graph        kNN (k = 15, euclidean) on F with pynndescent, ``random_state = ctx.seed``,
               single-threaded so the graph is reproducible for the seed.
  propagation  synchronous majority vote.  Every unlabelled cell takes the majority label of
               its LABELLED neighbours iff at least 3 neighbours are labelled and the majority
               holds >= 50 % of them while strictly beating the runner-up.  Repeat until no
               cell changes or 20 iterations.  A final pass applies the same rule once more,
               records the vote fraction, and classifies the residue:
                 Ambiguous(low_margin)               labelled neighbours, no label holds the margin
                 Ambiguous(few_labelled_neighbours)  only 1-2 labelled neighbours
                 Ambiguous(no_labelled_neighbour)    none of the k neighbours is labelled
  score        the vote fraction at assignment; seeds keep their label with score 1.0.
  NaN          never in a tail (``tail_mask``).  For the kNN only, a missing rank is imputed to
               the column's median rank over measured cells and the count is recorded.  A cell
               with NO finite identity marker is ``Unavailable(no_identity_measurement)`` and is
               left out of the graph (see Deviations).

Records (``params["run"]`` in the JSONL record, and ``ctx.extra["tasp"]`` for a driver):
n_seeds per label, labels below 20 seeds, labels absent, cells in >= 2 labels' tails, the
per-marker tail size, iterations and cells changed per iteration, convergence, the number and
fraction of cells assigned in the final pass, the Ambiguous reason counts, imputation counts,
and the pynndescent / numba versions.

Deviations from the planner's text (recorded here and in the report)
  * A cell whose EVERY identity marker is NaN is ``Unavailable(no_identity_measurement)`` rather
    than imputed to the all-median point.  Imputing one missing marker of a measured cell is a
    recorded distortion of a real cell; imputing all of them invents a cell at the centroid of
    feature space and hands it whatever label lives there, which root ``CLAUDE.md`` forbids
    ("a missing measurement is unavailable").  On donor 6436 the NaN set is the SAME cells for
    every identity marker (12.9 % of the uniform table), i.e. cells with no measurement at all.
  * ``Ambiguous(few_labelled_neighbours)`` is split out of ``low_margin`` for the 1-2 labelled
    neighbour case so the reason names the actual cause; the label is Ambiguous either way.
  * A tie for first place is not a majority: it is ``Ambiguous(low_margin)`` with the tied
    fraction as score, and the count is recorded as ``n_final_tie``.

Pins: pynndescent 0.6.0 (numba 0.67.0) from the repo ``.venv``; the installed versions are read
at import and written into every record.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/broadtype/arms/tasp.py --donor 6436 --scope frame
"""
from __future__ import annotations

import importlib.metadata
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import gate  # noqa: E402
import records  # noqa: E402
from criteria import CRITERIA  # noqa: E402

ARM = "tasp"
STOCHASTIC = True   # honours ctx.seed; the screen runs seeds 1..5 as replicates
METHOD_VERSION = "tasp-2026-09-05"

# Reason vocabulary.  Every Ambiguous / Unavailable row carries exactly one of these.
R_LOW_MARGIN = "low_margin"
R_FEW_LABELLED = "few_labelled_neighbours"
R_NO_LABELLED = "no_labelled_neighbour"
R_NO_MEASUREMENT = "no_identity_measurement"
R_NO_SEEDS = "no_label_has_seeds"
R_MULTI = "multi_lineage"
REASONS = (R_LOW_MARGIN, R_FEW_LABELLED, R_NO_LABELLED, R_NO_MEASUREMENT, R_NO_SEEDS, R_MULTI)


def _dist_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "missing"


#: The frozen configuration.  No key is a share, a prior or a floor on any label's frequency.
DEFAULT_PARAMS = {
    "method_version": METHOD_VERSION,
    "features": ("gate.standardise_within_donor(ctx.X): each declared marker centred and scaled by its own "
                 "spread over the loaded cells. Ranks were tried first and rejected before any arm ran: a rank "
                 "flattens a lineage's marker shift to a fraction of a rank unit, so kNN neighbourhoods were "
                 "dominated by the other markers and votes crossed lineages (planted-block accuracy 0.52)"),
    "propagation": "one kNN graph over every declared marker; iterative majority vote from the seeds until no cell changes",
    "margin_frac": 0.25,
    "tail_q": CRITERIA["tail_q"],
    "tail_min_k": CRITERIA["tail_min_k"],
    "seed_rule": "in the tail of >= 1 declared marker of exactly one label and in no other label's tail",
    "knn": "pynndescent.NNDescent",
    "k_neighbours": 15,
    "metric": "euclidean",
    "n_jobs": 1,
    "min_labelled_neighbours": 3,
    "majority_frac": 0.5,
    "max_iterations": 500,   # until no cell changes; each iteration touches only the unlabelled cells
    "seeds_min_for_note": 20,
    "seed_min_same_label_neighbours": 1,   # a seed with no same-label seed among its k neighbours is isolated and dropped
    "nan_policy": ("never in a tail; column-median rank imputed for the kNN only; "
                   "a cell with no finite identity marker -> Unavailable(no_identity_measurement)"),
    "weighting": "none: every labelled neighbour is one vote; no class weight and no term for any label's abundance",
    "pins": {"pynndescent": _dist_version("pynndescent"), "numba": _dist_version("numba")},
}
#: ``arm_base.cli`` hands this dict to the record AFTER ``run`` returns, so ``run`` stores its
#: diagnostics under ``RUN_PARAMS["run"]`` and they land in ``records.jsonl`` without the CLI
#: needing an arm-specific hook.  ``ctx.extra["tasp"]`` carries the same dict for a driver.
RUN_PARAMS = dict(DEFAULT_PARAMS)


# --------------------------------------------------------------------------- seeds
def seed_table(X: pd.DataFrame, label_markers: dict) -> tuple:
    """Seed code per cell (-1 = not a seed) from the Step-1 tails of the ORIGINAL values.

    ``label_markers`` maps label -> the declared markers present in ``X`` (declaration order,
    which fixes the integer code of each label).  A cell in the tails of two labels is not a
    seed for either: the method has no rule for choosing between them and must not invent one.
    """
    labels = list(label_markers)
    n = len(X)
    tails = np.zeros((n, len(labels)), dtype=bool)
    tail_k: dict = {}
    for j, L in enumerate(labels):
        for m in label_markers[L]:
            t = gate.tail_mask(X[m].to_numpy(float))
            tails[:, j] |= t
            tail_k[m] = int(t.sum())
    n_in = tails.sum(axis=1)
    seed = np.full(n, -1, dtype=np.int64)
    only = n_in == 1
    seed[only] = tails[only].argmax(axis=1)
    n_seeds = {L: int((seed == j).sum()) for j, L in enumerate(labels)}
    diag = {"n_seeds": n_seeds, "n_multi_tail": int((n_in >= 2).sum()), "tail_k": tail_k,
            "n_seed_cells": int(only.sum())}
    return seed, diag


# --------------------------------------------------------------------------- features
def feature_matrix(X: pd.DataFrame) -> tuple:
    """Rank features for the kNN.  Returns (F over measured cells, measured mask, diagnostics).

    A measured cell has at least one finite identity marker; its missing ranks are imputed to
    the column median rank over measured cells, for the graph only (the tails never see an
    imputed value).  A marker with no finite value anywhere is dropped from the features and
    recorded -- it cannot carry a tail either, so it contributes nothing to any call.
    """
    R = gate.standardise_within_donor(X.to_numpy(float)) if X.shape[1] else np.zeros((len(X), 0))
    R = np.where(np.isfinite(X.to_numpy(float)), R, np.nan)   # standardise_within_donor zero-fills NaN; keep them NaN here
    finite = np.isfinite(R)
    have_col = finite.any(axis=0)
    dropped = [c for c, ok in zip(X.columns, have_col) if not ok]
    R, finite = R[:, have_col], finite[:, have_col]
    measured = finite.any(axis=1) if R.shape[1] else np.zeros(len(X), dtype=bool)
    if R.shape[1] and measured.any():
        med = np.nanmedian(R[measured], axis=0)   # standardised scale: the median is ~0
    else:
        med = np.zeros(R.shape[1])
    impute = measured[:, None] & ~finite
    F = np.where(impute, med[None, :], R)[measured]
    diag = {"n_cells_measured": int(measured.sum()), "n_cells_unmeasured": int((~measured).sum()),
            "n_cells_partially_imputed": int(impute.any(axis=1).sum()), "n_values_imputed": int(impute.sum()),
            "markers_without_measurement": dropped, "n_features": int(R.shape[1])}
    return F, measured, diag


# --------------------------------------------------------------------------- graph
def knn_graph(F: np.ndarray, k: int, seed: int, *, n_jobs: int = 1) -> np.ndarray:
    """(n, k) neighbour indices without the cell itself; -1 marks a slot pynndescent left empty.

    pynndescent returns the query cell as its own first neighbour on the training data, so k+1
    are requested and the self entry is removed per row (by identity, not by position, because
    the search is approximate).  A row without a self entry loses its farthest neighbour
    instead so every row keeps exactly k.
    """
    import pynndescent
    n = F.shape[0]
    k_eff = int(min(k, n - 1))
    index = pynndescent.NNDescent(np.ascontiguousarray(F, dtype=np.float32), n_neighbors=k_eff + 1,
                                  metric="euclidean", random_state=seed, n_jobs=n_jobs)
    nbr = np.asarray(index.neighbor_graph[0], dtype=np.int64)
    is_self = nbr == np.arange(n)[:, None]
    keep = ~is_self
    keep[~is_self.any(axis=1), -1] = False
    return nbr[keep].reshape(n, k_eff)


# --------------------------------------------------------------------------- propagation
def _votes(nbr_labels: np.ndarray, n_labels: int) -> tuple:
    """Per row: number of labelled neighbours, the leading label, its count, the runner-up count."""
    mask = nbr_labels >= 0
    m = nbr_labels.shape[0]
    rows = np.broadcast_to(np.arange(m)[:, None], nbr_labels.shape)[mask]
    counts = np.bincount(rows * n_labels + nbr_labels[mask], minlength=m * n_labels).reshape(m, n_labels)
    n_lab = mask.sum(axis=1)
    top = counts.argmax(axis=1)
    top_count = counts.max(axis=1)
    second = np.sort(counts, axis=1)[:, -2] if n_labels >= 2 else np.zeros(m, dtype=np.int64)
    return n_lab, top, top_count, second


def _neighbour_labels(nbr: np.ndarray, label: np.ndarray, rows: np.ndarray) -> np.ndarray:
    nb = nbr[rows]
    out = np.full(nb.shape, -1, dtype=np.int64)
    valid = nb >= 0
    out[valid] = label[nb[valid]]
    return out


def propagate(nbr: np.ndarray, seed: np.ndarray, n_labels: int, *, min_labelled: int = DEFAULT_PARAMS["min_labelled_neighbours"],
              majority: float = DEFAULT_PARAMS["majority_frac"], max_iter: int = DEFAULT_PARAMS["max_iterations"]) -> tuple:
    """Iterative majority vote from the seeds; returns (label code, score, reason, diagnostics).

    Updates are synchronous (every unlabelled cell votes on the labels as they stood at the
    start of the iteration), so the result does not depend on cell order.  ``label`` is -1 for
    a cell that never met the rule; its ``reason`` says why.
    """
    n = len(seed)
    label = seed.astype(np.int64).copy()
    score = np.where(label >= 0, 1.0, np.nan)
    reason = np.full(n, "", dtype=object)
    changed: list = []
    converged = False
    for _ in range(max_iter):
        unl = np.flatnonzero(label < 0)
        if unl.size == 0:
            converged = True
            break
        n_lab, top, top_count, second = _votes(_neighbour_labels(nbr, label, unl), n_labels)
        ok = (n_lab >= min_labelled) & (top_count > second) & (top_count >= majority * n_lab)
        if not ok.any():
            converged = True
            break
        label[unl[ok]] = top[ok]
        score[unl[ok]] = top_count[ok] / n_lab[ok]
        changed.append(int(ok.sum()))
    # final pass: the same rule once more, then classify the residue with its margin recorded
    unl = np.flatnonzero(label < 0)
    n_final, n_tie = 0, 0
    counts = {R_LOW_MARGIN: 0, R_FEW_LABELLED: 0, R_NO_LABELLED: 0}
    if unl.size:
        n_lab, top, top_count, second = _votes(_neighbour_labels(nbr, label, unl), n_labels)
        ok = (n_lab >= min_labelled) & (top_count > second) & (top_count >= majority * n_lab)
        label[unl[ok]] = top[ok]
        score[unl[ok]] = top_count[ok] / n_lab[ok]
        n_final = int(ok.sum())
        rest = ~ok
        frac = np.where(n_lab > 0, top_count / np.maximum(n_lab, 1), np.nan)
        none = rest & (n_lab == 0)
        few = rest & (n_lab > 0) & (n_lab < min_labelled)
        low = rest & (n_lab >= min_labelled)
        for sel, r in ((none, R_NO_LABELLED), (few, R_FEW_LABELLED), (low, R_LOW_MARGIN)):
            reason[unl[sel]] = r
            score[unl[sel]] = frac[sel]
            counts[r] = int(sel.sum())
        n_tie = int((low & (top_count == second)).sum())
    diag = {"n_iterations": len(changed), "changed_per_iteration": changed, "converged": converged,
            "n_final_pass_assigned": n_final, "frac_final_pass_assigned": float(n_final / n) if n else 0.0,
            "n_final_tie": n_tie, "ambiguous_reasons": counts}
    return label, score, reason, diag


# --------------------------------------------------------------------------- the arm
def coherent_seeds(seed: np.ndarray, nbr: np.ndarray, min_same: int, labels: list) -> tuple:
    """Drop isolated seeds: a seed with fewer than ``min_same`` same-label seeds among its neighbours.

    A lineage's bright tail is a population -- its cells co-express the lineage's markers and sit
    together in feature space -- while the tail of a marker that carries no lineage structure is a
    random scatter of cells whose neighbourhoods hold no other seed of that label.  Requiring one
    same-label seed neighbour removes the scatter without consulting any label's abundance
    (measured 2026-09-05 on planted blocks: random-tail seeds fell from ~20 to ~0 per label while
    planted-block seeds were kept).  Dropped seeds are counted per label in the diagnostics.
    """
    seed = seed.astype(np.int64).copy()
    is_seed = seed >= 0
    nb = nbr[is_seed]
    nb_lab = np.where(nb >= 0, seed[np.maximum(nb, 0)], -1)
    same = (nb_lab == seed[is_seed][:, None]).sum(axis=1)
    keep = same >= min_same
    dropped = seed[is_seed][~keep]
    out = seed.copy()
    out[np.flatnonzero(is_seed)[~keep]] = -1
    diag = {"n_seeds_isolated_dropped": {L: int((dropped == j).sum()) for j, L in enumerate(labels)},
            "n_seeds_kept": {L: int((out == j).sum()) for j, L in enumerate(labels)}}
    return out, diag


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    """Calls for every ``ctx.ids`` row from ``ctx.X`` alone (declared identity markers).

    One kNN graph over every declared marker on the standardised scale; the seeds are the Step-1
    tails; labels grow by synchronous majority vote until no cell changes.  Nothing here weights
    a class or consults any label's abundance: the seed count is fixed by the tail definition and
    the called shares are whatever the graph structure yields.
    """
    dec = declarations.load()
    label_markers = {L: [m for m in spec["markers"] if m in ctx.markers] for L, spec in dec["labels"].items()}
    absent = [L for L, ms in label_markers.items() if not ms]
    present = {L: ms for L, ms in label_markers.items() if ms}
    labels = list(present)
    n = len(ctx.ids)

    out_label = np.full(n, records.UNAVAILABLE, dtype=object)
    out_score = np.full(n, np.nan)
    out_reason = np.full(n, "", dtype=object)

    F_all, measured, fdiag = feature_matrix(ctx.X)
    kept_cols = [c for c in ctx.X.columns if c not in set(fdiag["markers_without_measurement"])]
    out_reason[~measured] = R_NO_MEASUREMENT
    seed_all, sdiag = seed_table(ctx.X, present) if present else (np.full(n, -1, dtype=np.int64),
                                                                    {"n_seeds": {}, "n_multi_tail": 0, "tail_k": {}, "n_seed_cells": 0})
    diag = {"labels": labels, "labels_absent": absent, "label_markers": present, **sdiag, **fdiag,
            "seeds_below_min": [L for L, c in sdiag["n_seeds"].items() if c < DEFAULT_PARAMS["seeds_min_for_note"]],
            "k_neighbours_used": 0, "per_label": {}}
    idx = np.flatnonzero(measured)
    if idx.size >= 2 and sdiag["n_seed_cells"] > 0:
        k = DEFAULT_PARAMS["k_neighbours"]
        nbr = knn_graph(F_all, k, ctx.seed, n_jobs=DEFAULT_PARAMS["n_jobs"])
        diag["k_neighbours_used"] = int(nbr.shape[1])
        seed_m, cdiag = coherent_seeds(seed_all[measured], nbr, DEFAULT_PARAMS["seed_min_same_label_neighbours"], labels)
        diag.update(cdiag)
        code, score, reason, pdiag = propagate(nbr, seed_m, len(labels))
        diag.update(pdiag)
        diag["seed_index"] = [int(i) for i in idx[seed_m >= 0]]      # frozen-order indices of the seeds that propagated
        called = code >= 0
        lab_arr = np.asarray(labels, dtype=object)
        out_label[idx[called]] = lab_arr[code[called]]
        out_label[idx[~called]] = records.AMBIGUOUS
        out_score[idx] = score
        out_reason[idx] = reason
    elif idx.size:
        out_reason[idx] = R_NO_SEEDS
        diag.update({"n_iterations": 0, "n_final_pass_assigned": 0, "frac_final_pass_assigned": 0.0,
                     "ambiguous_reasons": {R_LOW_MARGIN: 0, R_FEW_LABELLED: 0, R_NO_LABELLED: 0, R_MULTI: 0}})
    diag["n_called_per_label"] = {L: int((out_label == L).sum()) for L in labels}
    diag["n_ambiguous"] = int((out_label == records.AMBIGUOUS).sum())
    diag["n_unavailable"] = int((out_label == records.UNAVAILABLE).sum())
    ctx.extra["tasp"] = diag
    RUN_PARAMS["run"] = diag
    return pd.DataFrame({"object_id": ctx.ids.astype(str), "label": out_label, "score": out_score, "reason": out_reason})


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=sys.executable, default_params=RUN_PARAMS))
