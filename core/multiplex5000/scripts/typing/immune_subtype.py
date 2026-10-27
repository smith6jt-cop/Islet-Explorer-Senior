#!/usr/bin/env python3
"""Phase 2 -- name the immune cells, inside the class they are already in.

The cells that carry the single label ``Immune`` and nothing finer are the population the study
question is about.  This resolves them against the established rule table
(``phenotype_workflow.csv``, the phenotyping guide) and **re-decides nothing outside them**: a cell
that is not ``Immune`` is not read, not scored and not written.  ``cell_class`` is never modified.

The hierarchy is the table's own, three levels deep:

    Immune -> Imm-T-cell | Imm-B-cell | Imm-macrophage | Imm-dendritic-cell
              Imm-NK-cell | Imm-neutrophil | Imm-plasma-cell
      Imm-T-cell     -> T-regulatory | T-exhausted | T-progenitor-memory
                        T-CD8-cytotoxic | T-CD8 | T-CD4-helper
      Imm-macrophage -> Mac-M1-like | Mac-M2-like

**CD4 and CD8 are read only after the T-cell parent is established.** The recorded
rules still name the more specific T-cell leaves. If none fires, the higher
background-to-bright-end scaled CD8 or CD4 names the coarse subtype without
another 0.10 cutoff (owner, 2026-10-27). The comparison uses values before
clipping. Missing measurements stay unavailable; exact ties keep the T-cell
parent and require review of the total marker profile.

**Panel variants are availability, not absence.** A rule whose markers a donor did not ACQUIRE is
``unavailable`` for that donor and recorded as such -- never scored as a negative, which would
quietly call every cell "not a plasma cell" on exactly the donors that cannot be asked.  Which
markers a donor acquired is measured from that donor's own aligned partition, per-row-group null
counts included: a column exists for every marker in the cohort union, and a column that is
entirely null is a marker the slide never carried, not a marker that read zero.

**A cell whose value is missing is undecided, not negative.** ``value >= 0.10`` is False for NaN,
so a sparse missing measurement used to slide silently into "this rule does not fire".  Every rule
naming a marker that is NaN *on that cell* is undecided for that cell, the cell's call is withheld
with that reason, and the count is reported.  **"Undecidable" is its own outcome, distinct from
"no rule fired"**: ``immune_class_undecided`` / ``immune_subtype_undecided`` carry the first and
the summary counts ``(uncalled: no rule fired)`` separately, so the two can never be added
together into one "unnamed" number.

**NO APP ARTIFACT READS ``typing/immune/``.** This tree is a second pass beside the typing, not
part of it.  The four app builders -- ``build_app_data.py``, ``build_workbook.py``,
``neighborhood.py``, ``build_parquet_views.py`` -- contain no reference to it, and their donor
coverage comes from each builder's own ``COHORT.resolve``.  Nothing here is on the app's path,
and ``cell_class`` is never modified.

**Multiple rules can fire.** That is not resolved by a silent precedence: every match is counted,
the most specific rule wins (most required markers, ties broken by name for determinism), and both
the match count and the losing rules are written per cell so the choice is auditable.

**Nothing is written as 0 where it was not decidable.** A subtype whose rule is unavailable on a
donor is reported EMPTY with the reason, the way the pLN lane reports NK_ILC.  A zero is a
measurement only where the rule could actually be evaluated.

**But EMPTY is not the end of the question.** A rule the cohort cannot ask everywhere is also
evaluated in its REDUCED form -- the same rule with the markers no donor of this run is missing,
and only those, taken out -- and that count is reported on EVERY donor, including the ones that
do carry the marker.  It is a different question with a different name, never written into
``immune_class``, and it is cross-donor comparable precisely because it asks the same thing of
everyone.  `Imm-plasma-cell` is `CD38 + CD79a` with CD20 negative; on the seven slides without
CD38 what remains is CD79a-positive and CD20-negative, which is the B-lineage marker plasma cells
keep with the one they lose -- worth reporting as itself, and not worth calling a plasma cell.

WHAT TRANSFERS FROM THE pLN LANE (``build_pln_immune_pheno.py``), and what does not.

  Transfers: (1) EMPTY-with-a-reason instead of 0; (2) the per-cell offset -- each cell's own
  median rank over the whole panel -- because a dim cell reads low on *everything*, so "did not
  clear the level on any naming marker" and "is not that cell type" are different statements and
  the offset is what tells them apart; (3) the discipline that a ``SUPPORTING_ONLY`` marker may
  not carry a call alone, so what a FOXP3 or CD56 rule names has to be said in the words of the
  measurement.

  Does NOT transfer: the clustering that substitutes for context, and the within-section
  standout statistic.  A node has no partition of the section to draw a reference population
  from; a pancreas does -- the QuPath objects -- and this lane's levels are fitted against it
  already.  The offset is carried here as a DIAGNOSTIC on the calls, not as the gate: the gate
  stays the rule table at the level the coarse immune definitions were evaluated at.

    python core/multiplex5000/scripts/typing/immune_subtype.py --include-marginal
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import rankdata

from multiplex5000.tissues.t_cell import (
    CD4,
    CD8,
    REVIEW_TIE,
    UNAVAILABLE,
    compare_cd8_cd4,
    unbounded_scale,
)

HERE = __import__("pathlib").Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "celltype"))
import cohort as COHORT
import lineages as LIN
import rules as RULES
import subpanel_io as S

OUT = S.OUT / "immune"

#: The same level the coarse immune definitions were evaluated at (RECONSTRUCTION.md).
REQUIRE_MIN = 0.10

LEVEL2_PARENT = "Immune"
LEVEL3_PARENTS = ("Imm-T-cell", "Imm-macrophage")

#: A term is REPORTED SEPARATELY, with an alternative count, when it clears the level in more
#: than this fraction of the cells its own rule is read on.  That is a within-donor observation
#: and NOT a verdict: CD4 above the level in most T cells is biology, and `(CD68 or Iba1)` above
#: the level in most macrophages is what made them macrophages.  It only says the rule's other
#: terms did the separating, so the count with that term dropped is worth having beside the count
#: with it kept.
MAJORITY_TERM = 0.5

#: The VERDICT is cross-donor and is a property of the TERM, not of any donor.  A term that
#: clears the level in most of one donor's cells and almost none of another's is not measuring
#: the same thing on the two, so any count whose rule uses it is not comparable across them.
#: The level was fitted as a p99 of a reference population -- built to admit about 1% of it --
#: so `most` against `almost none` is stated as this pair rather than as a tuned cut.
#:
#: Measured 2026-09-17 over the twenty-donor build, and this is not hypothetical.  CD57
#: corroborates `Imm-NK-cell` through `(CD57 or Granzyme_B)`.  It clears the level in 86.03-95.95%
#: of immune cells on four donors and in 0.013-1.543% on the other sixteen -- and on fifteen of
#: those sixteen in 0.013-0.245%, the sole exception reaching 1.543%.  The raw CD57 median over
#: the same immune cells is 358-533 on the four against 0.21-9.67 on the sixteen; over each
#: donor's whole aligned partition it is 496-625 against 0.08-10.49, so the split is the channel
#: and not the immune selection.  Those four are exactly two acq_run values, which is
#: the finer batch variable correction is licensed on, never the Fusion version.  Their
#: `Imm-NK-cell` count is a measurement of that channel on that run.  The donors are not named
#: here: `term_comparability.csv` and `coverage.json` carry them, derived from the run.
NOT_COMPARABLE_HIGH = 0.5
NOT_COMPARABLE_LOW = 0.05

#: Permutation draws for the eta-squared null on a term's cross-donor positive fraction.
PERMUTATIONS = 10_000

#: How many cells to rank at once when building the per-cell offset.  The offset is a median over
#: the panel *per cell*, so the panel has to be resident for the rows being reduced; a donor's
#: whole aligned matrix is not.
OFFSET_BLOCK = 200_000


def immune_rules() -> tuple[pd.DataFrame, dict]:
    table = RULES.load_rules()
    reqs = S.shared_negatives(RULES.requirements(table, with_negatives=True))
    return table, reqs


def children_of(table: pd.DataFrame, parent: str) -> list[str]:
    return [p for p in table.loc[table["parent"] == parent, "phenotype"].astype(str).unique()]


def rule_markers(requirement) -> list[str]:
    """Every marker one rule names, positive, or-group and negative alike."""
    all_of, any_groups, none_of = requirement
    return sorted(set(all_of) | set(none_of) | {m for g in any_groups for m in g})


def supporting_only_note(name: str, requirement) -> str:
    """What a rule resting on a ``SUPPORTING_ONLY`` marker actually names, in its own terms.

    ``lineages.SUPPORTING_ONLY`` is the maintainer's statement that a marker may not carry a call
    alone -- FOXP3 and CD56 are both on it.  The rule table licenses them anyway, but only in
    company, so the honest description of the label is the company: the terms the rule requires
    ALONGSIDE the supporting-only marker.  Derived from the rule, so it cannot drift from it.

    LOAD-BEARING is not the same as MENTIONED.  A supporting-only marker inside an OR group that
    has other members is not what the rule rests on -- any other member satisfies the group -- so
    only ``all_of`` membership and a group whose SOLE member it is count.  Used negatively it gets
    its own sentence, because an exclusion inherits the marker's uncertainty just as a requirement
    does.
    """
    all_of, any_groups, none_of = requirement
    required = sorted(set(all_of) & LIN.SUPPORTING_ONLY)
    sole = sorted({g[0] for g in any_groups if len(g) == 1 and g[0] in LIN.SUPPORTING_ONLY})
    excluded = sorted(set(none_of) & LIN.SUPPORTING_ONLY)
    load_bearing = sorted(set(required) | set(sole))
    if not load_bearing and not excluded:
        return ""
    parts = []
    if load_bearing:
        corroborating = sorted(set(all_of) - set(load_bearing))
        groups = ["(" + " or ".join(sorted(g)) + ")" for g in any_groups
                  if not (len(g) == 1 and g[0] in load_bearing)]
        company = " and ".join([x for x in [", ".join(corroborating)] if x] + groups)
        parts.append(
            f"{name} requires {', '.join(load_bearing)}, which lineages.SUPPORTING_ONLY says may "
            f"not carry a call alone. It names the cells that clear the level on "
            f"{company or 'no other required marker'} and ALSO on {', '.join(load_bearing)}. Read "
            f"it as \"the cells of that population reading highest on {'/'.join(load_bearing)}\", "
            f"not as a count of {name} by an independent marker, because this panel carries none.")
    if excluded:
        parts.append(
            f"{name} excludes cells above the level on {', '.join(excluded)}, a supporting-only "
            "marker used negatively; the exclusion inherits the same uncertainty as a requirement "
            "on it would.")
    return " ".join(parts)


def term_fractions(scaled: pd.DataFrame, requirement, available: set[str]) -> dict:
    """marker -> fraction of the evaluated cells that clear the level, for every term of a rule."""
    markers = [m for m in rule_markers(requirement) if m in available]
    if not markers or not len(scaled):
        return {}
    block = scaled[markers].to_numpy()
    with np.errstate(invalid="ignore"):
        positive = block >= REQUIRE_MIN
    finite = np.isfinite(block)
    out = {}
    for j, m in enumerate(markers):
        n = int(finite[:, j].sum())
        out[m] = float(positive[finite[:, j], j].mean()) if n else float("nan")
    return out


def majority_terms(requirement, fractions: dict, parent_requirement=None) -> list[dict]:
    """The terms of one rule that clear the level in most of the cells the rule is read on.

    Reported for every role, because the consequence differs: a REQUIRED term that admits
    everyone silently removes itself from the rule; an OR-GROUP member that admits everyone
    satisfies the group for every cell, so the group's other members never speak; a NEGATIVE term
    that admits everyone excludes nearly every cell.  Only the OR-group case has a well-defined
    alternative -- drop that member and let the group's other members decide -- so only that case
    carries an alternative count.

    A term the PARENT rule already required is not reported: `Mac-M1-like` is read only on cells
    `Imm-macrophage` already named, and `(CD68 or Iba1)` is above the level in most of them
    BECAUSE that is what named them.  Reporting it would be reporting the selection.
    """
    all_of, any_groups, none_of = requirement
    inherited = set()
    if parent_requirement is not None:
        inherited = set(rule_markers(parent_requirement))
    fractions = {m: f for m, f in fractions.items() if m not in inherited}
    flagged = []
    for m in sorted(set(all_of)):
        f = fractions.get(m)
        if f is not None and np.isfinite(f) and f > MAJORITY_TERM:
            flagged.append({"marker": m, "role": "required", "positive_fraction": round(f, 4),
                            "has_alternative": False})
    for group in any_groups:
        members = [m for m in group if m in fractions]
        hot = [m for m in members
               if np.isfinite(fractions[m]) and fractions[m] > MAJORITY_TERM]
        for m in sorted(hot):
            flagged.append({"marker": m, "role": "or_group",
                            "positive_fraction": round(fractions[m], 4),
                            "group": sorted(members),
                            "has_alternative": len(members) - len(hot) > 0})
    for m in sorted(set(none_of)):
        f = fractions.get(m)
        if f is not None and np.isfinite(f) and f > MAJORITY_TERM:
            flagged.append({"marker": m, "role": "excluded", "positive_fraction": round(f, 4),
                            "has_alternative": False})
    return flagged


def restrict_required(requirement, drop: set[str]):
    """The same rule with ``drop`` removed from every position it holds -- a DIFFERENT question.

    This is what remains of a rule when a marker the cohort does not share is taken out of it.
    It is never written into ``immune_class``: `CD79a positive, CD20 negative` is not
    `CD38 positive, CD79a positive, CD20 negative`, and reporting one as the other would let a
    donor that cannot be asked look like a donor that answered.  It is reported beside the EMPTY,
    under its own name, on every donor -- including the ones that do carry the marker, which is
    what makes the two columns comparable to each other.
    """
    all_of, any_groups, none_of = requirement
    groups = [[m for m in g if m not in drop] for g in any_groups]
    return ([m for m in all_of if m not in drop],
            [g for g in groups if g],
            [m for m in none_of if m not in drop])


def restrict(requirement, drop: set[str]):
    """The same rule with ``drop`` removed from any OR group that keeps at least one member.

    Nothing else is touched: a required or a negative term is the rule's meaning, and dropping
    one would answer a different question rather than the same question better.
    """
    all_of, any_groups, none_of = requirement
    groups = []
    for group in any_groups:
        kept = [m for m in group if m not in drop]
        groups.append(kept if kept else list(group))
    return (all_of, groups, none_of)


def acquired_markers(donor: str, wanted: set[str]) -> tuple[list[str], list[str], list[str]]:
    """(acquired, column-absent, column-present-but-entirely-null) for one donor's own partition.

    The aligned schema is the cohort union, so a column's existence is not evidence the marker was
    measured.  Per-row-group null counts settle it exactly for the all-null case, which is the one
    that matters, and cost no data read -- the same test ``panel.donor_markers`` makes, applied to
    the partition this lane actually reads.  A donor with no partition at all is credited with
    NOTHING, because "not evidence a marker is absent" is not evidence it is present.
    """
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    if not path.exists():
        return [], sorted(wanted), []
    md = pq.ParquetFile(path).metadata
    index = {md.schema.column(i).name: i for i in range(md.num_columns)}
    got, absent, empty = [], [], []
    for marker in sorted(wanted):
        i = index.get(marker)
        if i is None:
            absent.append(marker)
            continue
        nulls = sum(md.row_group(g).column(i).statistics.null_count
                    for g in range(md.num_row_groups)
                    if md.row_group(g).column(i).statistics is not None)
        (empty if nulls >= md.num_rows else got).append(marker)
    return got, absent, empty


def evaluate(scaled: pd.DataFrame, requirement, available: set[str]):
    """(fires, undecided, unavailable) for one rule over every row.

    ``requirement`` is (all_of, [any_of groups], none_of) as the table declares it.
    A rule naming a marker this donor does not carry is unavailable for the whole donor, not false.
    A rule naming a marker that is NaN on one cell is undecided for that cell, not false.
    """
    all_of, any_groups, none_of = requirement
    needed = set(rule_markers(requirement))
    if not needed <= available:
        return None, None, sorted(needed - available)

    named = sorted(needed)
    block = scaled[named].to_numpy()
    undecided = ~np.isfinite(block).all(axis=1)

    fires = np.ones(len(scaled), dtype=bool)
    for marker in all_of:
        fires &= scaled[marker].to_numpy() >= REQUIRE_MIN
    for group in any_groups:
        if not group:
            continue
        hit = np.zeros(len(scaled), dtype=bool)
        for marker in group:
            hit |= scaled[marker].to_numpy() >= REQUIRE_MIN
        fires &= hit
    for marker in none_of:
        fires &= scaled[marker].to_numpy() < REQUIRE_MIN
    return fires & ~undecided, undecided, []


def resolve(frame: pd.DataFrame, candidates: list[str], reqs, scaled, available, specificity):
    """Assign the most specific firing rule; record the count, the runners-up and the undecided."""
    hits, undecided, skipped = {}, np.zeros(len(frame), dtype=bool), {}
    for name in candidates:
        fires, undec, missing = evaluate(scaled, reqs[name], available)
        if fires is None:
            skipped[name] = missing
            continue
        hits[name] = fires
        undecided |= undec
    blank = (pd.Series([None] * len(frame), index=frame.index),
             pd.Series(0, index=frame.index), pd.Series([""] * len(frame), index=frame.index),
             pd.Series(undecided, index=frame.index), skipped)
    if not hits:
        return blank

    order = sorted(hits, key=lambda n: (-specificity[n], n))
    stacked = np.vstack([hits[n] for n in order])
    n_match = stacked.sum(axis=0)
    first = np.argmax(stacked, axis=0)
    call = np.where(n_match > 0, np.array(order)[first], None)
    others = [
        ";".join(n for k, n in enumerate(order) if stacked[k, i] and k != first[i])
        if n_match[i] > 1 else ""
        for i in range(len(frame))
    ]
    # A cell that a rule DID name is decided, whatever some other rule could not read on it.
    undecided &= ~(n_match > 0)
    return (pd.Series(call, index=frame.index), pd.Series(n_match, index=frame.index),
            pd.Series(others, index=frame.index), pd.Series(undecided, index=frame.index), skipped)


def per_cell_offset(values: np.ndarray) -> np.ndarray:
    """Each cell's own median place on this donor's panel, in [0, 1].

    The pLN lane's term, and its reason: a dim cell reads low on every channel, so a cell that
    clears no rule may be a cell type this panel cannot name OR simply a dim cell.  Ranks are
    taken per marker over this donor's analysis-eligible cells, so the statistic is invariant to
    multiplying any channel by any positive constant; the per-cell median over the panel is then
    that cell's brightness, in the units of the panel rather than of a channel.

    Diagnostic only.  Nothing here decides a call.
    """
    n, k = values.shape
    if n == 0 or k == 0:
        return np.full(n, np.nan, dtype=np.float32)
    unit = np.empty((n, k), dtype=np.float32)
    for j in range(k):
        column = values[:, j]
        finite = np.isfinite(column)
        unit[:, j] = np.nan
        if finite.sum() > 1:
            unit[finite, j] = (rankdata(column[finite]) / finite.sum()).astype(np.float32)
    out = np.full(n, np.nan, dtype=np.float32)
    # A cell with no finite value anywhere on the panel has no offset. That is a measurement --
    # it is counted and reported -- not a warning to be stepped over, and `nanmedian` would emit
    # `RuntimeWarning: All-NaN slice` for it and return NaN anyway.
    any_finite = np.isfinite(unit).any(axis=1)
    for a in range(0, n, OFFSET_BLOCK):
        block = unit[a:a + OFFSET_BLOCK]
        rows = any_finite[a:a + OFFSET_BLOCK]
        if rows.any():
            out[a:a + OFFSET_BLOCK][rows] = np.nanmedian(block[rows], axis=1)
    return out


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """Between-group share of the total variance -- the effect size, not a test."""
    values = np.asarray(values, dtype=float)
    total = ((values - values.mean()) ** 2).sum()
    if total <= 0:
        return float("nan")
    between = 0.0
    for g in np.unique(groups):
        take = groups == g
        between += take.sum() * (values[take].mean() - values.mean()) ** 2
    return float(between / total)


def eta_squared_with_null(values, groups, draws: int = PERMUTATIONS, seed: int = 0):
    """(eta2, permutation p, null median) -- the null is the same labels, shuffled.

    Acquisition batch is aliased with disease status at the coarse level, so a batch effect
    reported as a bare eta-squared means nothing: with 20 donors over many runs a large share of
    the variance falls between groups BY CONSTRUCTION.  The null says how large it falls there
    when the labels carry no information at all, and the p is the share of shuffles that reach
    the observed value.
    """
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    observed = eta_squared(values, groups)
    if not np.isfinite(observed):
        return observed, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    null = np.empty(draws)
    for i in range(draws):
        null[i] = eta_squared(values, rng.permutation(groups))
    return observed, float((null >= observed).mean()), float(np.median(null))


def _count_or_empty(counts: pd.Series, label: str, unavailable: dict[str, list[str]],
                    parent_present: bool, parent: str,
                    flagged: dict[str, list[dict]]) -> tuple[object, str, str]:
    """(count, status, reason) for one label on one donor -- an integer, or EMPTY with a reason."""
    if label in unavailable:
        return "", "EMPTY_rule_unavailable", (
            f"the rule names {', '.join(unavailable[label])}, which this donor did not acquire; "
            "an unasked question is not a negative answer")
    if not parent_present:
        return "", "EMPTY_parent_absent", (
            f"no cell of this donor was named {parent}, so this subtype was never evaluated; "
            "the rule is evaluable but was not reached")
    if label in flagged:
        terms = "; ".join(f"{t['marker']} ({t['role']}) clears the level in "
                          f"{100 * t['positive_fraction']:.1f}% of the cells this rule was read on"
                          for t in flagged[label])
        return int(counts.get(label, 0)), "measured_majority_term", (
            f"{terms}. The rule's other terms did the separating; the count with that term "
            "dropped is in count_without_majority_term. This is an observation about this "
            "donor, not a verdict about the count.")
    return int(counts.get(label, 0)), "measured", ""


def apply_t_cd8_cd4_comparison(out: pd.DataFrame, values: pd.DataFrame,
                              markers: list[str], scale: pd.DataFrame) -> np.ndarray:
    """Name otherwise untyped T cells by the higher measured marker, or flag a tie.

    The parent and any more specific rule calls have already been decided. This
    comparison never treats a missing marker as zero and never clips before
    comparing. It mutates only the subtype and its audit fields in ``out``.
    """
    t_untyped = ((out["immune_class"] == "Imm-T-cell")
                 & out["immune_subtype"].isna()).to_numpy(dtype=bool)
    compared = {}
    for marker in ("CD8", "CD4"):
        compared[marker] = (unbounded_scale(
            values[marker].to_numpy(dtype=float),
            floor=float(scale.loc[marker, "floor"]),
            span=float(scale.loc[marker, "span"]))
            if marker in markers else np.full(len(values), np.nan))
    comparison = compare_cd8_cd4(compared["CD8"], compared["CD4"], t_untyped)
    out["immune_subtype_resolved_by"] = np.where(out["immune_subtype"].notna(), "rule", "")
    out["immune_subtype_review_reason"] = ""
    for state, subtype in ((CD8, "T-CD8"), (CD4, "T-CD4-helper")):
        selected = comparison == state
        out.loc[selected, "immune_subtype"] = subtype
        out.loc[selected, "immune_subtype_resolved_by"] = "higher_scaled_CD8_CD4"
    out.loc[comparison == REVIEW_TIE, "immune_subtype_review_reason"] = (
        "exact_CD8_CD4_tie_review_total_marker_profile")
    out.loc[comparison == UNAVAILABLE, "immune_subtype_review_reason"] = (
        "CD8_or_CD4_measurement_unavailable")
    return comparison


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--all", action="store_true", help="every manifest donor, not just clean")
    args = parser.parse_args()

    marginal_left_out: dict[str, str] = {}
    if args.all:
        manifest = json.loads(S.QUPATH_MANIFEST.read_text())
        sys.path.insert(0, str(S.REPO / "core/multiplex5000"))
        from phenocycler import cohort as pc
        donors = sorted({r["donor_id"] for r in manifest["images"]} - set(pc.DONOR_EXCLUSIONS))
    else:
        donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
        # A run that silently omits donors the cohort module documents is how this lane came to
        # cover ten of twenty for three days without anything saying so.  Name them here.
        marginal_left_out = {d: why for d, why in COHORT.MARGINAL.items() if d not in donors}

    table, reqs = immune_rules()
    level2 = children_of(table, LEVEL2_PARENT)
    level3 = {p: children_of(table, p) for p in LEVEL3_PARENTS}
    specificity = {n: len(reqs[n][0]) + sum(len(g) for g in reqs[n][1]) + len(reqs[n][2])
                   for n in reqs}
    every_rule = level2 + [k for kids in level3.values() for k in kids]
    #: child rule -> its parent's requirement, so a term the parent already required is not
    #: reported as "above the level in most cells" when the parent is what selected those cells.
    parent_req = {k: reqs[parent] for parent, kids in level3.items() for k in kids}
    needed_markers = sorted({m for n in every_rule for m in rule_markers(reqs[n])})
    print(f"level 2 ({len(level2)}): {', '.join(level2)}")
    for parent, kids in level3.items():
        print(f"  {parent} ({len(kids)}): {', '.join(kids)}")
    print(f"\n{len(every_rule)} rules over {len(needed_markers)} markers; "
          f"{len(donors)} donors requested")
    if marginal_left_out:
        print("NOT requested, and the cohort module says why:")
        for d, why in sorted(marginal_left_out.items()):
            print(f"  {d}: {why}  (--include-marginal to build it)")
    notes = {n: supporting_only_note(n, reqs[n]) for n in every_rule}
    for n, note in notes.items():
        if note:
            print(f"\nSUPPORTING_ONLY -- {note}")
    print()

    # Which of the rules' markers is some donor of THIS run missing?  Metadata only, no data
    # read.  A rule naming one of them cannot be asked of the whole cohort, so its reduced form
    # is evaluated on everyone -- see the module docstring.
    scarce: set[str] = set()
    for donor in donors:
        have, _, _ = acquired_markers(donor, set(needed_markers))
        scarce |= set(needed_markers) - set(have)
    reduced = {n: restrict_required(reqs[n], scarce) for n in every_rule
               if set(rule_markers(reqs[n])) & scarce}
    if reduced:
        print(f"markers no donor of this run is missing: none of {sorted(scarce)}; "
              f"reduced forms evaluated on every donor for "
              f"{', '.join(sorted(reduced))}\n")

    scale = S.load_scale()
    OUT.mkdir(parents=True, exist_ok=True)
    summary, skipped_report, coverage, rows = [], {}, {}, []

    for donor in donors:
        started = time.time()
        cells = pq.read_table(S.OUT / f"cells/donor_id={donor}/data.parquet",
                              columns=["object_id", "cell_class", "acq_run",
                                       "stain_group"]).to_pandas()
        batch = {v: sorted(set(cells[v].dropna().astype(str))) for v in ("acq_run", "stain_group")}
        immune = cells[cells["cell_class"] == "Immune"]
        if immune.empty:
            coverage[donor] = {"state": "empty", "n_cells": int(len(cells)), "n_immune": 0,
                               "reason": "no cell of this donor carries the class Immune"}
            print(f"  {donor}: no Immune cells"); continue

        panel_all = sorted(set(pq.ParquetFile(
            S.ALIGNED / f"donor_id={donor}/data.parquet").schema_arrow.names) & set(scale.index))
        markers, absent, empty_cols = acquired_markers(donor, set(panel_all) - {"DAPI"})
        if not markers:
            coverage[donor] = {"state": "unavailable", "n_cells": int(len(cells)),
                               "n_immune": int(len(immune)),
                               "reason": "this donor has no aligned expression partition, so it "
                                         "is credited with no acquired marker at all"}
            print(f"  {donor}: UNAVAILABLE -- no aligned partition"); continue

        values = pq.read_table(S.ALIGNED / f"donor_id={donor}/data.parquet",
                               columns=["object_id"] + markers).to_pandas()
        values = S.mask_fragment_hormones(values, donor)
        # The offset is a property of this donor's analysis-eligible cells, which is what
        # `cells/` already is -- a cell with no nucleus is not in it and does not enter the ranks.
        eligible = cells[["object_id"]].merge(values, on="object_id", how="left",
                                              validate="one_to_one")
        offset = per_cell_offset(eligible[markers].to_numpy())
        n_no_measurement = int(np.isnan(offset).sum())
        offset = pd.Series(offset, index=eligible["object_id"].to_numpy())

        values = immune.merge(values, on="object_id", how="left", validate="one_to_one")
        scaled = pd.DataFrame(
            S.apply_scale(values[markers].to_numpy(), markers, scale),
            columns=markers, index=values.index)
        # apply_scale clips, and np.clip propagates NaN -- checked, not assumed, because a NaN
        # silently becoming 0.0 would turn every missing measurement into a negative answer.
        raw_missing = ~np.isfinite(values[markers].to_numpy())
        assert (~np.isfinite(scaled.to_numpy()) == raw_missing).all(), \
            "apply_scale did not preserve missingness -- refusing to score"

        call2, n2, other2, undec2, skip2 = resolve(
            values, level2, reqs, scaled, set(markers), specificity)
        out = pd.DataFrame({
            "donor_id": donor, "object_id": values["object_id"],
            "immune_class": call2, "immune_class_matches": n2, "immune_class_also": other2,
            "immune_class_undecided": undec2,
            "immune_subtype": None, "immune_subtype_matches": 0, "immune_subtype_also": "",
            "immune_subtype_undecided": False,
            "panel_offset": offset.reindex(values["object_id"].to_numpy()).to_numpy(),
        })
        for parent, kids in level3.items():
            mask = (out["immune_class"] == parent).to_numpy()
            if not mask.any():
                continue
            sub, n3, other3, undec3, skip3 = resolve(values[mask], kids, reqs,
                                                     scaled[mask], set(markers), specificity)
            out.loc[mask, "immune_subtype"] = sub.to_numpy()
            out.loc[mask, "immune_subtype_matches"] = n3.to_numpy()
            out.loc[mask, "immune_subtype_also"] = other3.to_numpy()
            out.loc[mask, "immune_subtype_undecided"] = undec3.to_numpy()
            skip2.update(skip3)

        comparison = apply_t_cd8_cd4_comparison(out, values, markers, scale)

        # ---- Do the rule's own terms separate the cells it is read on, ON THIS DONOR?
        evaluated = {n: np.ones(len(out), dtype=bool) for n in level2}
        for parent, kids in level3.items():
            for k in kids:
                evaluated[k] = (out["immune_class"] == parent).to_numpy()
        fractions, flagged, repairable = {}, {}, set()
        for name in every_rule:
            if name in skip2:
                continue
            rows_here = evaluated[name]
            if not rows_here.any():
                continue
            fractions[name] = term_fractions(scaled[rows_here], reqs[name], set(markers))
            hits = majority_terms(reqs[name], fractions[name], parent_req.get(name))
            if hits:
                flagged[name] = hits
                repairable |= {t["marker"] for t in hits if t["has_alternative"]}

        restricted_counts: dict[str, int] = {}
        if repairable:
            tight = {n: restrict(r, repairable) for n, r in reqs.items()}
            # Specificity is recomputed on the restricted rules: "most specific wins" has to be
            # decided on the rules actually being compared, or a rule that lost a term keeps the
            # precedence it had with it.
            tight_specificity = {n: len(r[0]) + sum(len(g) for g in r[1]) + len(r[2])
                                 for n, r in tight.items()}
            r2, _, _, _, _ = resolve(values, level2, tight, scaled, set(markers),
                                     tight_specificity)
            restricted_counts = {k: int(v) for k, v in r2.value_counts().items()}
            for parent, kids in level3.items():
                mask2 = (r2 == parent).to_numpy()
                if not mask2.any():
                    continue
                r3, _, _, _, _ = resolve(values[mask2], kids, tight, scaled[mask2],
                                         set(markers), tight_specificity)
                restricted_counts.update({k: int(v) for k, v in r3.value_counts().items()})
            print(f"    {donor}: terms above the level in most of the cells their rule is "
                  f"read on -- {sorted(repairable)}; re-resolved without them for comparison",
                  flush=True)

        skipped_report[donor] = skip2
        target = OUT / f"donor_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        out.to_parquet(target / "data.parquet", index=False)
        tied = comparison == REVIEW_TIE
        if tied.any():
            review = out.loc[tied, ["donor_id", "object_id", "immune_class",
                                    "immune_subtype_review_reason"]].copy()
            for marker in ("CD8", "CD4"):
                review[f"{marker}_scaled_before_clip"] = unbounded_scale(
                    values.loc[tied, marker].to_numpy(dtype=float),
                    floor=float(scale.loc[marker, "floor"]),
                    span=float(scale.loc[marker, "span"]))
            review["full_profile_source"] = str(
                S.ALIGNED / f"donor_id={donor}/data.parquet")
            review.to_csv(target / "t_cd8_cd4_profile_review.csv", index=False)

        class_counts = out["immune_class"].value_counts()
        sub_counts = out["immune_subtype"].value_counts()
        named = int(out["immune_class"].notna().sum())
        undecided = int(out["immune_class_undecided"].sum())
        unnamed = len(out) - named - undecided
        for label in level2:
            count, status, why = _count_or_empty(class_counts, label, skip2, True,
                                                 LEVEL2_PARENT, flagged)
            rows.append({"donor_id": donor, "level": "immune_class", "label": label,
                         "count": count, "status": status, "reason": why,
                         "count_without_majority_term":
                             restricted_counts.get(label, "") if repairable else "",
                         "supporting_only_note": notes[label]})
        for parent, kids in level3.items():
            present = bool(class_counts.get(parent, 0))
            for label in kids:
                count, status, why = _count_or_empty(sub_counts, label, skip2, present,
                                                     parent, flagged)
                rows.append({"donor_id": donor, "level": "immune_subtype", "label": label,
                             "count": count, "status": status, "reason": why,
                             "count_without_majority_term":
                                 restricted_counts.get(label, "") if repairable else "",
                             "supporting_only_note": notes[label]})
        for name, requirement in reduced.items():
            rows_here = evaluated.get(name)
            if rows_here is None or not rows_here.any():
                continue
            fires, _, missing = evaluate(scaled[rows_here], requirement, set(markers))
            dropped = sorted(set(rule_markers(reqs[name])) - set(rule_markers(requirement)))
            rows.append({
                "donor_id": donor, "level": "reduced_rule",
                "label": f"{name} [without {', '.join(dropped)}]",
                "count": "" if fires is None else int(fires.sum()),
                "status": "EMPTY_rule_unavailable" if fires is None else "measured_reduced_rule",
                "reason": (f"the reduced form still names {', '.join(missing)}, which this donor "
                           "did not acquire") if fires is None else
                          (f"{name} with {', '.join(dropped)} taken out, evaluated on every donor "
                           "of this run so the column is comparable across them. It is a "
                           "DIFFERENT question from " + name + " and is never written into "
                           "immune_class; on a donor that did acquire " + ", ".join(dropped) +
                           " it is the same cells asked a weaker question."),
                "count_without_majority_term": "",
                "supporting_only_note": notes.get(name, "")})
        rows.append({"donor_id": donor, "level": "immune_class", "label": "(uncalled: no rule fired)",
                     "count": unnamed, "status": "measured",
                     "reason": "every evaluable rule was read on this cell and none fired",
                     "count_without_majority_term": "",
                     "supporting_only_note": ""})
        rows.append({"donor_id": donor, "level": "immune_class", "label": "(uncalled: undecidable)",
                     "count": undecided, "status": "measured",
                     "reason": "a marker some rule names is missing on this cell, so the rule "
                               "could not be read; a missing measurement is not a negative",
                     "count_without_majority_term": "",
                     "supporting_only_note": ""})

        called_offset = out.loc[out["immune_class"].notna(), "panel_offset"]
        blank_offset = out.loc[out["immune_class"].isna(), "panel_offset"]
        coverage[donor] = {
            "state": "built",
            "cohort_note": COHORT.MARGINAL.get(donor, "clean"),
            "acq_run": batch["acq_run"], "stain_group": batch["stain_group"],
            "n_cells": int(len(cells)), "n_immune": int(len(out)),
            "n_markers_acquired": len(markers),
            "n_eligible_cells_with_no_finite_measurement": n_no_measurement,
            "rule_markers_absent": sorted(set(needed_markers) - set(markers)),
            "panel_columns_entirely_null": empty_cols,
            "unavailable_rules": {k: v for k, v in skip2.items()},
            "named": named, "uncalled_no_rule_fired": unnamed, "uncalled_undecidable": undecided,
            "panel_offset_median_named": float(np.nanmedian(called_offset)) if named else None,
            "panel_offset_median_uncalled": float(np.nanmedian(blank_offset)) if unnamed else None,
            "immune_class": {k: int(v) for k, v in class_counts.items()},
            "immune_subtype": {k: int(v) for k, v in sub_counts.items()},
            "t_subtype_by_higher_marker": int((out["immune_subtype_resolved_by"]
                                                == "higher_scaled_CD8_CD4").sum()),
            "t_cd8_cd4_exact_ties_for_profile_review": int(tied.sum()),
            "term_positive_fraction": {n: {m: round(f, 5) for m, f in v.items()}
                                       for n, v in fractions.items()},
            "majority_terms": flagged,
            "counts_without_majority_terms": restricted_counts,
        }
        print(f"  {donor}: {len(out):>9,} Immune | named {100*named/len(out):5.1f}% | "
              f"undecidable {undecided:>6,} | "
              f"multi-match {100*int((out['immune_class_matches'] > 1).sum())/len(out):4.1f}% | "
              f"offset named {np.nanmedian(called_offset):.3f} vs uncalled "
              f"{np.nanmedian(blank_offset) if unnamed else float('nan'):.3f} | "
              f"unavailable rules: {sorted(skip2) or 'none'} | {time.time()-started:.0f}s",
              flush=True)
        summary.append(out)

    if not summary:
        return 1
    everything = pd.concat(summary, ignore_index=True)
    built = sorted(d for d, c in coverage.items() if c["state"] == "built")
    print(f"\n{len(everything):,} Immune cells over {len(built)} donors")
    print("\nimmune_class:")
    # The two ways a cell goes unnamed are different statements and are never pooled into one
    # "(uncalled)" row: one rule was read and did not fire, the other could not be read at all.
    labelled = everything["immune_class"].fillna(
        pd.Series(np.where(everything["immune_class_undecided"].to_numpy(),
                           "(uncalled: undecidable)", "(uncalled: no rule fired)"),
                  index=everything.index))
    print(labelled.value_counts().to_string())
    sub = everything[everything["immune_subtype"].notna()]
    if len(sub):
        print("\nimmune_subtype:")
        print(sub["immune_subtype"].value_counts().to_string())

    per_donor = pd.DataFrame(rows)

    # ---- THE CROSS-DONOR VERDICT.  A term that clears the level in most of one donor's cells
    # and almost none of another's is not measuring the same thing on the two, whatever it does
    # within either.  This is the only place a statement about more than one donor is made, and
    # it comes with eta-squared against a permutation null on the finer batch variables, never
    # on the Fusion version.
    term_rows, not_comparable = [], {}
    fraction_of = {d: coverage[d].get("term_positive_fraction", {}) for d in built}
    for rule in every_rule:
        markers_here = sorted({m for d in built for m in fraction_of[d].get(rule, {})})
        for marker in markers_here:
            have = [(d, fraction_of[d][rule][marker]) for d in built
                    if marker in fraction_of[d].get(rule, {})
                    and np.isfinite(fraction_of[d][rule][marker])]
            if len(have) < 2:
                continue
            ids = [d for d, _ in have]
            fracs = np.array([f for _, f in have])
            row = {"rule": rule, "marker": marker, "n_donors": len(have),
                   "min_positive_fraction": round(float(fracs.min()), 5),
                   "max_positive_fraction": round(float(fracs.max()), 5),
                   "donor_at_min": ids[int(fracs.argmin())],
                   "donor_at_max": ids[int(fracs.argmax())]}
            for variable in ("acq_run", "stain_group"):
                labels = np.array([";".join(coverage[d][variable]) for d in ids])
                e2, pval, null = eta_squared_with_null(fracs, labels)
                row[f"eta2_{variable}"] = round(e2, 4) if np.isfinite(e2) else ""
                row[f"p_perm_{variable}"] = round(pval, 4) if np.isfinite(pval) else ""
                row[f"eta2_null_median_{variable}"] = round(null, 4) if np.isfinite(null) else ""
            row["cross_donor_comparable"] = not (
                fracs.max() > NOT_COMPARABLE_HIGH and fracs.min() < NOT_COMPARABLE_LOW)
            term_rows.append(row)
            if not row["cross_donor_comparable"]:
                high = sorted(d for d, f in have if f > NOT_COMPARABLE_HIGH)
                not_comparable.setdefault(rule, []).append({"marker": marker, "donors": high,
                                                            **{k: row[k] for k in row
                                                               if k.startswith(("min_", "max_",
                                                                                "eta2", "p_perm"))}})
    terms = pd.DataFrame(term_rows)
    if len(terms):
        terms.to_csv(OUT / "term_comparability.csv", index=False)

    for rule, entries in not_comparable.items():
        for entry in entries:
            hit = ((per_donor["label"] == rule) & (per_donor["level"] != "reduced_rule")
                   & per_donor["donor_id"].isin(entry["donors"])
                   & ~per_donor["status"].str.startswith("EMPTY"))
            why = (f"{entry['marker']} is a term of this rule and clears the level in "
                   f"{100 * entry['max_positive_fraction']:.1f}% of immune cells on its highest "
                   f"donor against {100 * entry['min_positive_fraction']:.3f}% on its lowest, so "
                   f"it is not measuring the same thing on the two. eta2 by acq_run "
                   f"{entry['eta2_acq_run']} (permutation p {entry['p_perm_acq_run']}), by "
                   f"stain_group {entry['eta2_stain_group']} (p {entry['p_perm_stain_group']}). "
                   f"This count is a measurement of this donor and is NOT comparable with a "
                   f"donor where that term is at background.")
            per_donor.loc[hit, "status"] = "measured_not_cross_donor_comparable"
            per_donor.loc[hit, "reason"] = why

    per_donor.to_csv(OUT / "per_donor_counts.csv", index=False)
    (OUT / "unavailable_rules.json").write_text(json.dumps(skipped_report, indent=2))
    (OUT / "coverage.json").write_text(json.dumps(
        {"donors_built": built,
         "donors_requested": list(donors),
         "donors_not_requested": marginal_left_out,
         "require_min": REQUIRE_MIN,
         "majority_term": MAJORITY_TERM,
         "not_comparable_criterion": {"high": NOT_COMPARABLE_HIGH, "low": NOT_COMPARABLE_LOW,
                                      "permutations": PERMUTATIONS},
         "terms_not_cross_donor_comparable": not_comparable,
         "supporting_only_notes": {k: v for k, v in notes.items() if v},
         "per_donor": coverage}, indent=2))

    hot = per_donor[per_donor["status"] == "measured_not_cross_donor_comparable"]
    if len(hot):
        print(f"\n{hot['donor_id'].nunique()} donor(s) x {hot['label'].nunique()} label(s) are "
              f"NOT cross-donor comparable, because a term of the rule is at background on part "
              f"of the cohort and above the level in most cells on the rest:")
        for rule, entries in not_comparable.items():
            for entry in entries:
                print(f"  {rule} via {entry['marker']}: {', '.join(entry['donors'])} "
                      f"(max {entry['max_positive_fraction']:.3f} vs min "
                      f"{entry['min_positive_fraction']:.5f}; eta2 acq_run "
                      f"{entry['eta2_acq_run']}, permutation p {entry['p_perm_acq_run']})")
    empties = per_donor[per_donor["status"].str.startswith("EMPTY")]
    print(f"\n{len(empties)} (donor, label) pairs written EMPTY rather than 0:")
    if len(empties):
        print(empties.groupby(["label", "status"]).size().to_string())
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
