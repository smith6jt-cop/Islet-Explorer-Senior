"""Step 3: name each cluster, behind a gate and a matched null.

**Why the cluster and not the cell.**  KINTSUGI's lane is cluster-then-annotate, and so is
the Xenium lane -- for a reason recorded in ``scripts/xenium_celltype/nullgate.py``: a
per-cell argmax over identity cores was measured against size- and detection-matched random
cores and five of eight lineages sat inside their null, two doing worse than random.  This
lane inherits that finding rather than re-testing it.

**Two stages, and they answer different questions.**

1. *The gate sets the candidates.*  A cluster's median profile either shows what a reference
   demands or it does not.  ``rules.requirements`` supplies exactly that -- ``all_of``
   markers and ``any_of`` groups, read off the declaring row so an OR-group already satisfied
   by an ancestor's commitment is dropped rather than turned into a hard requirement.  A
   marker counts as elevated when it sits ``MIN_ELEVATED`` of the way from its own measured
   background to its robust top, on the UniFORM-aligned scale.  That is a statement about
   background, not about prevalence -- a fixed quantile would assume how common each cell
   type is, which is the trap ``CLAUDE.md`` forbids and the one that made the DeepCellTypes
   decile contrast unreadable.
2. *The null decides the award.*  Clearing a gate is not the same as being that thing.
   ``nullgate.annotate_groups`` is imported and called, never copied, so this lane inherits
   its four corrections at once: no award on a negative contrast, per-marker unit scaling so
   one marker cannot carry a core, nulls matched on detection **and** magnitude with no
   self-substitution, and a winner that must beat the runner-up by more than the paired
   null difference or be refused as ``cores_cleared_but_are_not_separable``.

**The limit this inherits, stated where the number is used.**  The null-gate contrast is
``inside - outside``, so it asks "is this cluster enriched for L relative to the rest of the
tissue?" and not "what is this cluster?".  For minority compartments -- immune, vascular,
neural, endocrine, which is what the islet question needs -- the two coincide.  For the
dominant compartment they do not: pancreas is ~85% acinar and this 59-plex carries no acinar
enzyme marker, which is already the leading explanation for the cosine lane's 169,734 no-rule
cells.  **No composition denominator over the exocrine compartment is meaningful from these
calls**, and :func:`annotate` records that in the output rather than leaving it to a reader.

**No Unknown bucket.**  A cluster that clears nothing is ``Ambiguous`` with a reason, and a
cluster whose evidence is absent is ``Unavailable`` with a reason.  Both are recoverable
statements; "Unknown" is not.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
sys.path.insert(0, str(HERE.parent / "xenium_celltype"))
import panel                                        # noqa: E402
import rules                                        # noqa: E402
import nullgate                                     # noqa: E402
import score                                        # noqa: E402

#: How far above a marker's own measured background a cluster's median must sit to count as
#: elevated.  One number for the whole panel, stated and swept -- the same choice
#: ``celltype/score.py`` makes, and for the same reason.  It is a background statement, not a
#: prevalence one.
MIN_ELEVATED = 0.10

#: A cluster smaller than this is not annotated.  ``nullgate`` applies the same floor.
MIN_CLUSTER_CELLS = nullgate.MIN_GROUP_CELLS

AMBIGUOUS = "Ambiguous"
UNAVAILABLE = "Unavailable"


def background_mode(values: np.ndarray) -> float:
    """A marker's own background mode on the aligned scale.

    Mirrors ``celltype/run_step3.build_scale``, which is the floor the cosine lane's
    ``REQUIRE_MIN`` is measured from.  It is imported by value rather than by import because
    that module pulls in matplotlib and shares a filename with this lane's own driver.

    Using ``min()`` instead -- which an earlier draft of this file did -- makes the floor the
    REDSEA clamp at zero rather than the background, and then 56.8% of all (cell, marker)
    values sit above ``REQUIRE_MIN``: measured on donor 6551, every cluster showed 28-43 of
    55 markers, so every reference's negative constraints were violated and 94% of cells came
    back Ambiguous.  The floor is the whole gate.
    """
    v = np.asarray(values, dtype="float64")
    v = v[np.isfinite(v) & (v > 0)]
    if len(v) < 1000:
        return 0.0
    counts, edges = np.histogram(np.log1p(v), bins=512)
    return float(np.expm1(edges[1:][int(np.argmax(counts[1:])) + 1]))


def marker_scale(values: np.ndarray, markers: list[str]) -> pd.DataFrame:
    """[background mode, robust top] per marker -- ``score.marker_scale``'s own frame."""
    by_marker = {m: values[:, j] for j, m in enumerate(markers)}
    floors = {m: background_mode(v) for m, v in by_marker.items()}
    return score.marker_scale(by_marker, floors)


def elevated_markers(profile_row: pd.Series, scale: pd.DataFrame) -> set[str]:
    """Which markers this cluster's median profile shows, on the [0, 1] background scale."""
    lo, span = scale["floor"], scale["span"]
    scaled = ((profile_row - lo) / span).clip(0.0, 1.0)
    return set(scaled.index[scaled >= MIN_ELEVATED])


def candidate_references(shown: set[str], reqs: dict, flat: dict,
                         available: set[str]) -> tuple[list[str], dict]:
    """References whose requirements this cluster meets, and why the others were refused.

    A reference naming a marker this donor never acquired is **unavailable**, not failed:
    the panel is three variants, so evaluability is a per-donor property.
    """
    passed, refused = [], {}
    for phenotype, (all_of, any_groups) in reqs.items():
        named = set(all_of) | {m for g in any_groups for m in g} \
            | set(flat.get(phenotype, {}))
        missing_from_panel = named - available
        if missing_from_panel:
            refused[phenotype] = f"unavailable:{','.join(sorted(missing_from_panel))}"
            continue
        if not all(m in shown for m in all_of):
            refused[phenotype] = "all_of_not_met"
            continue
        if not all(any(m in shown for m in group) for group in any_groups):
            refused[phenotype] = "any_of_not_met"
            continue
        negatives = [m for m, t in flat.get(phenotype, {}).items() if rules.is_negative(str(t))]
        if any(m in shown for m in negatives):
            refused[phenotype] = "negative_constraint_violated"
            continue
        passed.append(phenotype)
    return passed, refused


def ancestry(table: pd.DataFrame) -> dict[str, list[str]]:
    """phenotype -> its chain of ancestors, nearest first, excluding ``all``."""
    parent = dict(zip(table["phenotype"], table["parent"]))
    out = {}
    for phenotype in table["phenotype"]:
        chain, cur = [], parent.get(phenotype)
        while cur and cur != "all":
            chain.append(cur)
            cur = parent.get(cur)
            if cur in chain:
                raise ValueError(f"cycle in the rule hierarchy at {phenotype}")
        out[phenotype] = chain
    return out


def resolve(passed: list[str], anc: dict[str, list[str]]) -> tuple[str, str]:
    """The most specific statement a set of satisfied references supports.

    The gate returns nested chains, because ``rules.flatten`` composes each phenotype with its
    ancestors: a cluster that satisfies ``Mac-M2-like`` necessarily satisfies
    ``Imm-macrophage`` and ``Immune`` too.  Those are not competing candidates and treating
    them as such is what made the matched null return
    ``cores_cleared_but_are_not_separable`` -- correctly, since their cores share most of
    their markers by construction.

    So the hierarchy resolves them, not a statistic.  Take the maximal candidates (those no
    other candidate descends from); if there is one, that is the call.  If there are several
    -- four T-cell subsets under ``Imm-T-cell``, say -- the call is their lowest common
    ancestor **within the candidate set**: the deepest thing the evidence does support,
    rather than an arbitrary pick among siblings it cannot separate.  Candidates spanning
    different top-level branches have no common ancestor below ``all`` and are genuinely
    ambiguous.
    """
    if not passed:
        return "", "no_reference_requirements_met"
    candidates = set(passed)
    maximal = [p for p in passed
               if not any(p in anc.get(q, []) for q in candidates if q != p)]
    if len(maximal) == 1:
        return maximal[0], ""
    common = [a for a in anc[maximal[0]] if all(a in anc[m] for m in maximal[1:])]
    common = [a for a in common if a in candidates]
    if common:
        return common[0], "resolved_to_common_ancestor"
    return "", "candidates_span_multiple_branches"


def core_markers(reqs: dict, flat: dict) -> dict[str, list[str]]:
    """The positive evidence each reference rests on -- what the null is matched against.

    Negative constraints are deliberately excluded: the null asks whether a cluster is
    *enriched* for this reference's markers, and a marker the reference forbids carries no
    enrichment to test.
    """
    out = {}
    for phenotype, (all_of, any_groups) in reqs.items():
        markers = sorted(set(all_of) | {m for g in any_groups for m in g})
        if markers:
            out[phenotype] = markers
    return out


def annotate(values: np.ndarray, markers: list[str], clusters: np.ndarray,
             donor_ids: np.ndarray, *, seed: int = 0,
             rules_csv: pathlib.Path | None = None) -> tuple[pd.DataFrame, dict]:
    """Name every cluster, or refuse it with a reason.  Returns (per-cluster, provenance)."""
    table = rules.load_rules(rules_csv)
    flat = rules.flatten(table)
    reqs = rules.requirements(table)
    available = set(markers)

    frame = pd.DataFrame(values, columns=markers)
    frame["cluster"] = clusters
    profile = frame.groupby("cluster")[markers].median()
    sizes = frame.groupby("cluster").size()

    scale = marker_scale(values, markers)

    anc = ancestry(table)

    # Stage 1 -- the gate proposes, the hierarchy resolves.
    proposals = {}
    for cluster in profile.index:
        n = int(sizes.loc[cluster])
        if n < MIN_CLUSTER_CELLS:
            proposals[cluster] = (n, set(), "", "too_few_cells", {})
            continue
        shown = elevated_markers(profile.loc[cluster], scale)
        passed, _refused = candidate_references(shown, reqs, flat, available)
        chosen, why = resolve(passed, anc)
        proposals[cluster] = (n, shown, chosen, why, passed)

    # Stage 2 -- the null confirms the resolved call, one reference at a time.  Scored
    # against ALL 37 at once it becomes a competition between nested references, which is a
    # question the hierarchy has already answered; asked one at a time it answers the
    # question that is actually open -- is this cluster enriched for THIS reference?
    all_cores = core_markers(reqs, flat)
    all_cores = {p: [m for m in ms if m in available] for p, ms in all_cores.items()}
    null_stats: dict[str, pd.DataFrame] = {}
    for reference in sorted({c for _, _, c, _, _ in proposals.values() if c}):
        core = all_cores.get(reference)
        if not core:
            continue
        null_stats[reference] = nullgate.annotate_groups(
            values, markers, clusters, {reference: core}, seed=seed).set_index("group")

    rows = []
    for cluster, (n, shown, chosen, why, passed) in proposals.items():
        row = {"cluster": cluster, "n_cells": n,
               "n_candidates": len(passed),
               "candidates": ",".join(sorted(passed)) if passed else "",
               "markers_shown": ",".join(sorted(shown))}
        if why == "too_few_cells":
            rows.append({**row, "call": UNAVAILABLE, "reason": "too_few_cells"})
            continue
        if not chosen:
            rows.append({**row, "call": AMBIGUOUS, "reason": why})
            continue
        stat = null_stats.get(chosen)
        nrow = stat.loc[cluster].to_dict() if stat is not None and cluster in stat.index else {}
        confirmed = bool(nrow.get("lineage"))
        rows.append({**row,
                     "call": chosen if confirmed else AMBIGUOUS,
                     "reason": why if confirmed else (
                         nrow.get("reason") or "not_confirmed_by_matched_null"),
                     "null_z": nrow.get("z", np.nan),
                     "null_real": nrow.get("real", np.nan),
                     "null_p95": nrow.get("null_p95", np.nan)})

    provenance = {
        "min_elevated": MIN_ELEVATED,
        "min_cluster_cells": MIN_CLUSTER_CELLS,
        "scale_quantile": score.SCALE_Q,
        "n_references": len(reqs),
        "n_references_with_core": len(all_cores),
        "rules_csv": str(rules_csv or rules.RULES_CSV),
        "majority_compartment_limit": (
            "The null-gate contrast is inside-minus-outside, so it measures enrichment "
            "relative to the rest of the tissue, not identity. Pancreas is ~85% acinar and "
            "this panel carries no acinar enzyme marker, so no composition denominator over "
            "the exocrine compartment is meaningful from these calls."),
    }
    return pd.DataFrame(rows), provenance
