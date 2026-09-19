#!/usr/bin/env python3
"""The refinement pass: a class is kept only where its own defining marker is present.

WRITTEN AND NOT WIRED (status, 2026-09-18).  This module runs and nothing consumes what it
writes.  ``refined_class`` has ZERO consumers; no cell in ``typing/cells/`` carries
``Unresolved (no lineage evidence)`` or any ``refined_*`` column; ``build_cells.py --refine
<dir>`` is the join that would put the pass in front of the app and it has not been run.  What
follows describes a PROPOSAL, and the open decisions in ``REFINEMENT.md`` section 10 are open.

THE BIGGEST OF THEM, because it qualifies every count this module reports: ``bound_X`` is the
5th percentile of a class's own concordant members, so **a class with no low tail at all still
loses about 5 % of it**.  Measured per class in ``REFINEMENT.md`` section 6: Lymphatic, Nerve,
Immune and Delta are demoted far above that floor -- the rule finding cells without their own
evidence -- while Exocrine's excess is only 0.2-0.5 percentage points, so essentially the whole
Exocrine demotion is the percentile relabelling a fixed share of the class by construction.
That is an OPEN QUESTION for the maintainer, not a decision this module makes; it implements
the recovered rule and reports the floor in every log.

RECOVERED, NOT REDESIGNED.  This module is the 2026-09-11 refinement pass made maintainable.
The rule, the constants and the routing are the code as it ran; what is new here is the scale
it runs at (the cohort, not a 320,000-cell subsample), the per-donor / per-stratum fitting the
maintainer ruled on 2026-09-18, and the guard that stops a caller iterating it.  Every
departure from the recovered run is named in ``REFINEMENT.md`` beside this file, with numbers.

WHAT IT DOES

    concordance(cell) = share of its K=30 nearest neighbours IN MARKER SPACE with the same label

Scatter is DETECTED in an embedding -- that is where it is visible -- but every cell is
CHARACTERISED and acted on in marker space, so the rule is a statement about measurements and
a tidier picture is a consequence, never the target.  This module builds no embedding at all.

    bound_X = 5th percentile of max(OWN markers of X) over the cells of X that are concordant
              (concordance >= 0.40), fitted within one donor and one islet size stratum
    a cell of X whose max(OWN markers of X) falls below bound_X is DEMOTED

Region membership -- a centroid inside a nerve, lymphatic or vessel object -- and the islet
prior are evidence of WHERE a cell is.  They are not evidence that the cell IS that thing when
its own marker is absent.  The bound is taken from the class's own clear members, never from a
constant, and a demoted cell falls back into the complement it came out of:

    demoted Exocrine          -> "Unresolved (no lineage evidence)"   (no complement beneath it)
    demoted anything else     -> "Immune" if it satisfies an immune leaf combination
                              -> "Exocrine" otherwise

``immune_leaf`` is what makes that routing correct.  Without it every demoted region cell went
to Exocrine and none could fall back to Immune; adding it moved 5,131 cells on the recovered
subsample.

IT IS APPLIED ONCE.  Step 2 was tested on the recovered run and rejected with its reason
recorded: re-applying the rule to an already-cleaned population is a ratchet, not a new cause.
The cells a second pass removes carry 3+ SD of their own defining marker and lie above the
Exocrine 95th percentile in 91-100 % of cases, for a 0.12 percentage-point gain.  Three
hypotheses were falsified before it was rejected.  :func:`refine_donor` therefore refuses an
input that has already been refined, and there is no step argument to raise.

WHAT IS FITTED PER DONOR AND PER STRATUM, AND WHY

``bound_X`` is a reference level -- a value derived from a reference population -- so the
standing ruling applies to it: per donor, and since 2026-09-18 per islet size stratum as well.
Per donor because a level fitted across donors splits by acquisition batch, and batch is
aliased with disease status at the coarse level, so a pooled level would carry a disease
signal by construction.  Per stratum because the same 2026-09-18 measurement that produced the
ruling shows the islet-core population is not one population across sizes.  A stratum exists
only where it is defined -- for islet-core cells -- and every cell outside an islet core sits
in one unstratified pool, because "the size of the islet this cell is not in" is not a number.

A (donor, stratum, class) group with fewer than :data:`MIN_CLEAR` concordant members gets NO
bound.  It is not borrowed from another stratum, not replaced by the donor's, and not
defaulted: demotion requires positive evidence that a cell falls below its own class's clear
members, and an unfittable group supplies none.  Those cells keep their class and the reason
is recorded, so an under-powered group is visible as one.

NEIGHBOURS ARE DRAWN WITHIN ONE DONOR

The recovered run pooled 20 donors into one 320,000-cell graph.  This module builds one graph
per donor, for three reasons, and the difference it makes is measured in ``REFINEMENT.md``
rather than hidden:

  * A neighbour drawn from another donor's staining is the same hazard that makes reference
    levels per-donor.  In the pooled graph a cell's concordance is decided in part by cells
    from another acquisition run, and acquisition batch is aliased with disease status.
  * Availability is a per-donor fact.  A within-donor graph uses exactly the markers that
    donor carries and imputes nothing.  The pooled matrix had to substitute a value for a
    marker a donor never acquired; a missing measurement is unavailable, not a number.
  * 30-NN over 25,230,109 x 54 in one graph is an order of magnitude beyond what was run.
    Twenty graphs of ~1.26M are tractable and independent.

The REPRESENTATION is still cohort-wide: one cofactor and one standardisation per marker,
from a seeded equal-per-donor sample.  That is deliberate and is the opposite decision to the
bound, for the reason ``run_step3.build_scale`` already gives -- a per-donor representation
would re-impose the distribution matching UniFORM's step 2 avoided.  The representation says
what a marker's units are; the bound says how much of it this donor's clear members show.

    python core/phenocycler/scripts/typing/refine.py --out <dir> [--donors 6450 6539]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "celltype"))

import subpanel_io as S                                      # noqa: E402
from islet_identity import STRUCTURAL                        # noqa: E402
import panel                                                 # noqa: E402
import rules as R                                            # noqa: E402
import score as SC                                           # noqa: E402
import run_step3 as S3                                       # noqa: E402

#: The share of its K nearest neighbours a cell needs to share its label with before it counts
#: as clear.  Pre-registered on the recovered run; not a tunable.
LOW = 0.40
K = 30
BOUND_PCT = 5.0
SEED = 0

#: Fewest concordant members a (donor, stratum, class) group must hold before a bound is
#: fitted from them.  The recovered run's own guard.  It is NOT ``subpanel_io``'s 500 -- that
#: guard belongs to the endocrine level, which is a per-cell decision made against the level
#: itself, whereas this one only decides whether the 5th percentile of a group is worth
#: taking.  The two are kept apart on purpose and both are reported, so a group that clears 50
#: but not 500 is visible in the log rather than assumed adequate.
MIN_CLEAR = 50

#: The label a demoted Exocrine cell takes.  It has no complement beneath it -- the whole point
#: of the demotion is that it carries no lineage evidence -- so it is recorded as unresolved
#: with that reason attached, not reassigned to a class by assertion.
UNRESOLVED = "Unresolved (no lineage evidence)"

#: Each class's OWN defining markers, verbatim from the recovered run.  A class is kept only
#: where the MAXIMUM over these reaches the class's own bound; it is a "does this cell carry
#: any of its own evidence" test, not a signature match.
OWN_MARKERS: dict[str, tuple[str, ...]] = {
    "Nerve": ("B3TUBB", "CD56"),
    "Lymphatic": ("Podoplanin",),
    "Blood vessel": ("CD31", "Caveolin", "CD34"),
    "Beta": ("INS",),
    "Alpha": ("GCG",),
    "Delta": ("SST",),
    "Immune": ("CD3e", "CD68", "Iba1", "CD20", "CD79a", "MPO", "CD66", "CD11c", "CD209"),
    "Exocrine": ("Ker8_18", "E_cadherin", "Pan_Cytokeratin"),
}

#: Classes the pass does not touch, and why.  ``Excluded`` is the
#: ``exocrine_inside_islet_border`` set the maintainer ruled out of analysis: an excluded cell
#: is excluded from every calculation, so it is neither demoted nor allowed to be somebody
#: else's neighbour.  ``UNRESOLVED`` cannot appear in an input that has not already been
#: refined, and its presence is how :func:`refine_donor` detects a second pass.
UNTOUCHED = ("Excluded",)

#: Markers dropped by name before the matrix is built.  CD38 and b_Catenin1 are panel variants
#: -- not acquired on all twenty donors -- and a marker only some donors carry cannot enter a
#: representation shared by all of them.  Everything else dropped is DERIVED from
#: ``shifts.json`` provenance, never restated here.
DROP_BY_NAME = ("CD38", "b_Catenin1")
DROP_PROVENANCE = "refused_reagent_change"

SHIFTS = S.REPO / "data/development/celltype/shifts.json"
CELLS = S.OUT / "cells"


class RefinementAlreadyApplied(RuntimeError):
    """Raised when a caller tries to refine an input the pass has already been applied to.

    The refinement is a single pass by measurement, not by convention -- see the module
    docstring and ``REFINEMENT.md`` section "Step 2, tested and rejected".  There is no step
    argument to raise and no loop to enter; a caller that wants a second pass has to delete
    this guard, and then has to answer the measurement that put it here.
    """


# --------------------------------------------------------------------------- inputs


def derive_markers() -> list[str]:
    """The markers the pass measures in, derived from the rule table -- never a list.

    Every marker that carries a non-zero entry in the flattened rule matrix, minus the panel
    variants and minus anything UniFORM refused to shift.  On 2026-09-18 this derives 54
    markers and drops three: CD38, b_Catenin1 (by name) and CD57 (by provenance).
    """
    flat = R.flatten(R.load_rules())
    Rm, _ = R.reference_matrix(flat, panel.markers(include_variant=True), S3.ENCODING,
                               strict=False)
    Rm = Rm.loc[(Rm != 0).any(axis=1)]
    shifts = json.loads(SHIFTS.read_text())["shifts"]
    drop = set(DROP_BY_NAME) | {m for m, r in shifts.items()
                                if r.get("provenance") == DROP_PROVENANCE}
    return [m for m in Rm.index if m not in drop]


def cohort_donors() -> list[str]:
    """The donors the pass runs over: those the current typing has written cells for.

    Read off the stage root rather than ``panel.cohort()``, which carries 22 donors -- 6533
    and 112 are not in this typing and a marker they do or do not share cannot enter a
    representation the other twenty are measured in.
    """
    return sorted(p.name.split("=", 1)[1] for p in CELLS.glob("donor_id=*")
                  if (p / "data.parquet").exists())


def load_cells(donor: str) -> pd.DataFrame:
    """One donor's current typing, with the islet size stratum attached.

    Islet size is counted exactly as ``islet_identity.load`` counts it -- over the cells the
    analysis ADMITS: islet core, non-structural, not Excluded -- so it is the size of the thing
    being looked at rather than QuPath's detection count.  A stratum is attached only to
    islet-core cells; every other cell carries ``""``, which is not a fifth stratum but the
    absence of one.
    """
    df = pq.read_table(CELLS / f"donor_id={donor}/data.parquet",
                       columns=["object_id", "cell_class", "cell_region_20um", "islet_num",
                                "qc_analysis_eligible"]).to_pandas()
    core = ((df["cell_region_20um"] == "core") & (~df["cell_class"].isin(STRUCTURAL))
            & (df["cell_class"] != "Excluded"))
    size = df.loc[core].groupby("islet_num")["object_id"].size()
    n = df["islet_num"].map(size)
    n = n.where(df["cell_region_20um"] == "core", other=np.nan)
    df["islet_cells"] = n.to_numpy(dtype=float)
    df["size_stratum"] = S.size_stratum(df["islet_cells"].to_numpy())
    return df


def donor_values(donor: str, markers: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """The donor's UniFORM-aligned values for ``markers``, and which of them it carries.

    Read through the celltype lane's own reader so availability cannot diverge from what the
    scoring that produced these classes saw.  ``ratio=False``: the recovered run measured the
    raw aligned values, and the DAPI ratio is a per-reference decision that has no meaning for
    a "does this cell carry its own marker at all" test.
    """
    df, have, _ = S3._read_donor(donor, markers, ["object_id", "qc_analysis_eligible"], False)
    return df, have


# ------------------------------------------------------------------- representation


def representation_params(donors: list[str], markers: list[str], n_per_donor: int = 200_000,
                          seed: int = SEED) -> pd.DataFrame:
    """Cohort-wide cofactor, centre and spread per marker -- the units the pass measures in.

    Verbatim from the recovered run: an arcsinh on each marker's own measured background
    (median / 5), then standardised per marker, so a bound reads in SD of that marker.  Equal
    cells per donor, seeded, so no donor's size decides the units.

    Cohort-wide on purpose, and the opposite decision to the bound: the alignment has already
    put the donors on one scale, and a per-donor representation would re-impose the
    distribution matching that step avoided (``run_step3.build_scale`` says the same of the
    marker scale).  The units are cohort-wide; the LEVEL a donor's own clear members show is
    not.
    """
    rng = np.random.default_rng(seed)
    pool: dict[str, list[np.ndarray]] = {m: [] for m in markers}
    for d in donors:
        t, have = donor_values(d, markers)
        t = t[t["qc_analysis_eligible"]]
        if len(t) > n_per_donor:
            t = t.iloc[rng.choice(len(t), n_per_donor, replace=False)]
        for m in have:
            pool[m].append(t[m].to_numpy("float64"))
    rows = {}
    for m, chunks in pool.items():
        if not chunks:
            continue                                    # unavailable, not zero: no row at all
        v = np.concatenate(chunks)
        cof = np.nanmedian(v) / 5.0
        if not np.isfinite(cof) or cof <= 0:
            cof = 1.0
        z = np.arcsinh(v / cof)
        rows[m] = {"cofactor": float(cof), "centre": float(np.nanmean(z)),
                   "spread": float(np.nanstd(z)), "n": int(np.isfinite(v).sum())}
    out = pd.DataFrame(rows).T
    out.index.name = "marker"
    return out


def to_Z(values: np.ndarray, markers: list[str], params: pd.DataFrame) -> np.ndarray:
    """(cells x markers) -> standardised arcsinh, in SD of each marker.

    Refuses a marker it has no parameters for rather than filling one in, and refuses a cell
    carrying a non-finite value rather than substituting a number for it -- the caller is
    expected to have removed those cells and to have said how many.
    """
    missing = [m for m in markers if m not in params.index]
    if missing:
        raise KeyError(f"no representation parameters for {missing} -- unavailable, not zero")
    cof = params.loc[markers, "cofactor"].to_numpy(float)
    centre = params.loc[markers, "centre"].to_numpy(float)
    spread = params.loc[markers, "spread"].to_numpy(float)
    spread = np.where(np.isfinite(spread) & (spread > 0), spread, 1.0)
    Z = (np.arcsinh(np.asarray(values, dtype="float64") / cof) - centre) / spread
    if not np.isfinite(Z).all():
        raise ValueError("non-finite value in the representation -- a missing measurement is "
                         "unavailable and its cell must be dropped, not filled in")
    return Z.astype("float32")


# ------------------------------------------------------------------------ the rule


def neighbours(Z: np.ndarray, k: int = K, seed: int = SEED) -> np.ndarray:
    """Index of each cell's ``k`` nearest neighbours in marker space, self dropped."""
    from pynndescent import NNDescent
    idx = NNDescent(Z, n_neighbors=k + 1, metric="euclidean", random_state=seed, n_jobs=-1)
    nbr, _ = idx.neighbor_graph
    return nbr[:, 1:]


def concordance(nbr: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Share of each cell's neighbours carrying the same label.

    Recomputed on the NEW labels after the demotion, over the SAME neighbour index: the
    representation does not change when a label does, and NNDescent on the same matrix with
    the same seed returns the same graph, so rebuilding it would cost a second search to
    reproduce the first one exactly.
    """
    return (labels[nbr] == labels[:, None]).mean(axis=1)


def own_marker_max(Z: np.ndarray, markers: list[str], cls: str) -> np.ndarray | None:
    """max(OWN markers of ``cls``) per cell, or None when the donor carries none of them."""
    cols = [markers.index(m) for m in OWN_MARKERS[cls] if m in markers]
    if not cols:
        return None
    return Z[:, cols].max(axis=1)


def fit_bounds(Z: np.ndarray, markers: list[str], labels: np.ndarray, strata: np.ndarray,
               conc: np.ndarray, low: float = LOW, pct: float = BOUND_PCT,
               min_clear: int = MIN_CLEAR) -> dict:
    """(class, stratum) -> the level that class's own clear members show, in that stratum.

    Returns an entry for every (class, stratum) group that exists, with ``bound`` None where
    the group holds too few concordant members to fit one.  None means UNAVAILABLE: the caller
    must not demote in that group, must not borrow another group's bound, and must say so.
    """
    out: dict = {}
    for cls in OWN_MARKERS:
        v = own_marker_max(Z, markers, cls)
        if v is None:
            continue
        in_class = labels == cls
        for st in pd.unique(strata[in_class]) if in_class.any() else []:
            m = in_class & (strata == st)
            clear = m & (conc >= low)
            bound = (float(np.percentile(v[clear], pct))
                     if int(clear.sum()) >= min_clear else None)
            out[(cls, str(st))] = {"n": int(m.sum()), "n_clear": int(clear.sum()),
                                   "bound_sd": None if bound is None else round(bound, 4)}
    return out


def _demote(Z: np.ndarray, markers: list[str], labels: np.ndarray, strata: np.ndarray,
            bounds: dict) -> tuple[np.ndarray, np.ndarray]:
    """Which cells fall below their own (class, stratum) bound, and why not where they do not.

    Private, and it takes the labels it is given rather than reading them: the public entry
    point applies it exactly once.
    """
    demoted = np.zeros(len(labels), bool)
    why = np.full(len(labels), "", dtype=object)
    for cls in OWN_MARKERS:
        v = own_marker_max(Z, markers, cls)
        in_class = labels == cls
        if v is None:
            why[in_class] = "no own marker on this donor"
            continue
        for st in (pd.unique(strata[in_class]) if in_class.any() else []):
            m = in_class & (strata == st)
            rec = bounds.get((cls, str(st)))
            if rec is None or rec["bound_sd"] is None:
                why[m] = f"bound unavailable ({rec['n_clear'] if rec else 0} clear < {MIN_CLEAR})"
                continue
            below = m & (v < rec["bound_sd"])
            demoted |= below
            why[below] = f"below {cls} bound {rec['bound_sd']:+.3f} SD in stratum {st or 'none'}"
    return demoted, why


def apply_once(Z: np.ndarray, markers: list[str], labels: np.ndarray, strata: np.ndarray,
               immune_leaf: np.ndarray, nbr: np.ndarray | None = None) -> dict:
    """THE pass, applied ONCE to a matrix that is already in the representation.

    The whole rule in one call, so a caller cannot assemble half of it or run it twice by
    assembling it twice: fit the bounds on the labels it is given, demote against them, route
    the demoted cells, and recompute concordance on the result.  It takes labels as an
    argument and returns new ones; it never reads or writes a table, which is what makes it
    testable on a synthetic matrix.

    Returns ``{labels, demoted, note, bounds, concordance_before, concordance, neighbours}``.
    """
    if len(labels) != len(Z) or len(strata) != len(Z) or len(immune_leaf) != len(Z):
        raise ValueError("labels, strata and immune_leaf must be one per row of Z")
    if (labels == UNRESOLVED).any():
        raise RefinementAlreadyApplied(
            f"{UNRESOLVED!r} is produced only by this pass -- refusing to apply it twice")
    nbr = neighbours(Z) if nbr is None else nbr
    before = concordance(nbr, labels)
    bounds = fit_bounds(Z, markers, labels, strata, before)
    demoted, note = _demote(Z, markers, labels, strata, bounds)
    new = route(labels, demoted, immune_leaf)
    return {"labels": new, "demoted": demoted, "note": note, "bounds": bounds,
            "concordance_before": before, "concordance": concordance(nbr, new),
            "neighbours": nbr}


def route(labels: np.ndarray, demoted: np.ndarray, immune_leaf: np.ndarray) -> np.ndarray:
    """Where a demoted cell goes: back into the complement it came out of.

    A demoted Exocrine cell has no complement beneath it and becomes ``UNRESOLVED``.  Anything
    else loses only the override and is then Immune if it satisfies an immune leaf
    combination, else Exocrine.  It is not relabelled Exocrine by assertion.
    """
    new = labels.copy()
    exo = demoted & (labels == "Exocrine")
    oth = demoted & (labels != "Exocrine")
    new[exo] = UNRESOLVED
    new[oth & immune_leaf] = "Immune"
    new[oth & ~immune_leaf] = "Exocrine"
    return new


# ---------------------------------------------------------------------- immune_leaf


def immune_leaves() -> tuple[dict, list[str]]:
    """The rule table's immune leaf phenotypes, and what each of them requires.

    ``Imm-`` / ``T-`` / ``Mac-`` prefixed rows that survive into ``requirements`` -- read off
    the table, so a table that gains a leaf gains it here.
    """
    raw = R.load_rules()
    reqs = R.requirements(raw, with_negatives=True)
    phenotype = raw.columns[1]
    leaves = [p for p in raw[phenotype].astype(str)
              if p.startswith(("Imm-", "T-", "Mac-")) and p in reqs]
    return reqs, leaves


def immune_leaf(donor: str, markers: list[str] | None = None,
                require_min: float = S3.REQUIRE_MIN) -> pd.DataFrame:
    """object_id -> does this cell satisfy AT LEAST ONE immune leaf combination.

    Per donor, over ``qc_analysis_eligible`` cells, joined on ``object_id`` and never on
    position or row order.  ``require_min`` is ``run_step3.REQUIRE_MIN`` -- 0.10, pre-registered
    and swept once, imported rather than restated so it cannot drift from the value the
    classification's own immune tier uses.

    The scale is the frozen cohort scale the classification was made on
    (``subpanel_io.load_scale``), not a recomputed one: ``build_scale`` draws a seeded sample,
    so recomputing it would drift from the values these classes were called against.  That is
    a deliberate departure from the recovered run, which recomputed; the difference it makes is
    reported in ``REFINEMENT.md``.
    """
    reqs, leaves = immune_leaves()
    markers = markers or derive_markers()
    t, have = donor_values(donor, markers)
    t = t[t["qc_analysis_eligible"]]
    Xs = S.apply_scale(t[have].to_numpy("float32"), have, S.load_scale())
    ok = SC.meets_requirements(Xs, have, reqs, leaves, require_min).any(axis=1)
    return pd.DataFrame({"object_id": t["object_id"].to_numpy(), "immune_leaf": ok})


# -------------------------------------------------------------------------- the pass


def refine_donor(donor: str, params: pd.DataFrame, markers: list[str] | None = None,
                 cells: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict]:
    """Apply the pass ONCE to one donor.  Returns (per-cell frame, log).

    The frame is a sidecar keyed on ``object_id``, every row of the donor's typing, in the
    typing's own order.  It does not rewrite the typing: a cell the pass could not measure
    keeps its class and carries the reason it was not measured, which a table that had been
    overwritten in place could no longer say.
    """
    markers = markers or derive_markers()
    cells = load_cells(donor) if cells is None else cells
    _refuse_second_pass(cells)

    values, have = donor_values(donor, markers)
    elig = values.loc[values["qc_analysis_eligible"], "object_id"]
    frame = cells.merge(values[["object_id"] + have], on="object_id", how="left",
                        validate="one_to_one")

    # Who the pass may act on. A cell excluded from the analysis is excluded from this
    # calculation too -- it is neither demoted nor allowed to be somebody else's neighbour.
    measurable = (frame["object_id"].isin(set(elig)) & frame["qc_analysis_eligible"]
                  & ~frame["cell_class"].isin(UNTOUCHED)
                  & np.isfinite(frame[have].to_numpy("float64")).all(axis=1)).to_numpy()

    reason = np.full(len(frame), "", dtype=object)
    reason[~measurable] = "not measured (excluded, ineligible or no aligned value)"
    for c in UNTOUCHED:
        reason[(frame["cell_class"] == c).to_numpy()] = f"{c}: outside the pass"

    Z = to_Z(frame.loc[measurable, have].to_numpy("float32"), have, params)
    labels = frame.loc[measurable, "cell_class"].to_numpy().astype(object)
    strata = frame.loc[measurable, "size_stratum"].to_numpy().astype(object)

    leaf = immune_leaf(donor, markers)
    flag = (frame[["object_id"]].merge(leaf, on="object_id", how="left")["immune_leaf"]
            .fillna(False).to_numpy(bool))

    t0 = time.time()
    got = apply_once(Z, have, labels, strata, flag[measurable])
    t_graph = time.time() - t0
    new, demoted, why = got["labels"], got["demoted"], got["note"]
    conc_before, conc_after = got["concordance_before"], got["concordance"]
    bounds = got["bounds"]

    out = pd.DataFrame({
        "object_id": frame["object_id"].to_numpy(),
        "donor_id": donor,
        "cell_class": frame["cell_class"].to_numpy(),
        "refined_class": frame["cell_class"].to_numpy().astype(object),
        "size_stratum": frame["size_stratum"].to_numpy(),
        "islet_cells": frame["islet_cells"].to_numpy(),
        "immune_leaf": flag,
        "concordance": np.full(len(frame), np.nan),
        "concordance_before": np.full(len(frame), np.nan),
        "demoted": np.zeros(len(frame), bool),
        "refine_note": reason,
        "refine_pass": np.where(measurable, 1, 0),
    })
    out.loc[measurable, "refined_class"] = new
    out.loc[measurable, "concordance"] = conc_after
    out.loc[measurable, "concordance_before"] = conc_before
    out.loc[measurable, "demoted"] = demoted
    out.loc[measurable, "refine_note"] = why

    log = _log(donor, frame, measurable, have, bounds, labels, new, demoted, flag[measurable],
               conc_before, conc_after, strata, t_graph)
    return out, log


def _refuse_second_pass(cells: pd.DataFrame) -> None:
    """Refuse an input the pass has already been applied to.

    The refinement converged at one pass and the measurement that says so is recorded in
    ``REFINEMENT.md``.  Iterating it to convergence is a ratchet: each pass re-derives the
    bound from a population the previous pass already cleaned, so it keeps finding a 5th
    percentile to cut whether or not anything is wrong with the cells below it.
    """
    applied = [c for c in ("refined_class", "refine_pass") if c in cells.columns]
    if applied:
        raise RefinementAlreadyApplied(
            f"input already carries {applied} -- the refinement is a single pass by "
            "measurement (step 2 removed cells carrying 3+ SD of their own defining marker "
            "for a 0.12 pp gain and was rejected); refine the ORIGINAL typing, not its output")
    if (cells["cell_class"] == UNRESOLVED).any():
        raise RefinementAlreadyApplied(
            f"input carries {UNRESOLVED!r}, which only this pass produces -- refusing to "
            "apply it a second time")


def _log(donor, frame, measurable, have, bounds, labels, new, demoted, flag, before, after,
         strata, t_graph) -> dict:
    rows = []
    for (cls, st), rec in sorted(bounds.items()):
        m = (labels == cls) & (strata == st)
        n_dem = int((demoted & m).sum())
        rows.append({"cell_class": cls, "size_stratum": st, **rec, "n_demoted": n_dem,
                     "pct_demoted": round(100 * n_dem / max(int(m.sum()), 1), 2)})
    exo_dem = int((demoted & (labels == "Exocrine")).sum())
    oth = demoted & (labels != "Exocrine")
    core = (frame.loc[measurable, "size_stratum"].to_numpy() != "")
    endo = np.isin(labels, ["Beta", "Alpha", "Delta"])
    return {
        "donor_id": donor,
        "rule": ("keep a class only if its own defining marker reaches the "
                 f"{BOUND_PCT}th percentile of that class's own concordant members "
                 f"(concordance >= {LOW}) in the same donor and size stratum; demoted cells "
                 "fall back to the complement"),
        "own_markers": {k: list(v) for k, v in OWN_MARKERS.items()},
        "markers_used": list(have),
        "neighbours": {"k": K, "scope": "within donor", "seconds": round(t_graph, 1)},
        "n_cells": int(len(frame)),
        "n_measured": int(measurable.sum()),
        "n_not_measured": int((~measurable).sum()),
        "groups": rows,
        "n_groups_unfittable": sum(1 for r in rows if r["bound_sd"] is None),
        "n_groups_under_500": sum(1 for r in rows
                                  if r["bound_sd"] is not None and r["n_clear"] < 500),
        "n_demoted": int(demoted.sum()),
        "routed": {"Unresolved": exo_dem,
                   "Immune": int((oth & flag).sum()),
                   "Exocrine": int((oth & ~flag).sum())},
        "demoted_islet_core_endocrine_to_exocrine": int((oth & ~flag & core & endo).sum()),
        "scattered_pct_before": round(100 * float((before < LOW).mean()), 2),
        "scattered_pct_after": round(100 * float((after < LOW).mean()), 2),
        "before": [{"cell_class": c, "n": int((labels == c).sum()),
                    "median_concordance": round(float(np.median(before[labels == c])), 3),
                    "pct_scattered": round(100 * float((before[labels == c] < LOW).mean()), 1)}
                   for c in sorted(pd.unique(labels))],
        "after": [{"cell_class": c, "n": int((new == c).sum()),
                   "median_concordance": round(float(np.median(after[new == c])), 3),
                   "pct_scattered": round(100 * float((after[new == c] < LOW).mean()), 1)}
                  for c in sorted(pd.unique(new))],
        "iterated": False,
        "single_pass_reason": ("step 2 was tested and rejected: it removed cells carrying 3+ SD "
                               "of their own defining marker, above the Exocrine 95th "
                               "percentile in 91-100% of cases, for a 0.12 pp gain"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, type=pathlib.Path,
                    help="directory to write donor_id=<d>/data.parquet and refine_log.json")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--params", type=pathlib.Path, default=None,
                    help="cached representation parameters; built and written if absent")
    args = ap.parse_args(argv)

    donors = args.donors or cohort_donors()
    markers = derive_markers()
    print(f"refinement: {len(donors)} donors, {len(markers)} markers, K={K}, "
          f"bound = p{BOUND_PCT:g} of each class's concordant members per donor x stratum")

    if args.params and args.params.exists():
        params = pd.read_parquet(args.params)
        print(f"  representation read from {args.params}")
    else:
        t0 = time.time()
        params = representation_params(cohort_donors(), markers)
        print(f"  representation built over {len(cohort_donors())} donors "
              f"({time.time() - t0:.0f}s)")
        if args.params:
            args.params.parent.mkdir(parents=True, exist_ok=True)
            params.to_parquet(args.params)

    args.out.mkdir(parents=True, exist_ok=True)
    logs = []
    for d in donors:
        t0 = time.time()
        out, log = refine_donor(d, params, markers)
        p = args.out / f"donor_id={d}"
        p.mkdir(parents=True, exist_ok=True)
        out.to_parquet(p / "data.parquet", index=False)
        log["seconds"] = round(time.time() - t0, 1)
        logs.append(log)
        print(f"  {d}: {log['n_measured']:,} measured, {log['n_demoted']:,} demoted "
              f"-> Unresolved {log['routed']['Unresolved']:,} · "
              f"Immune {log['routed']['Immune']:,} · Exocrine {log['routed']['Exocrine']:,}"
              f"   scatter {log['scattered_pct_before']}% -> {log['scattered_pct_after']}%"
              f"   {log['seconds']}s", flush=True)
    (args.out / "refine_log.json").write_text(json.dumps(logs, indent=1))
    print(f"recorded to {args.out / 'refine_log.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
