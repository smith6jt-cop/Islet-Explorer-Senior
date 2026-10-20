#!/usr/bin/env python3
"""The refinement pass: a class is kept only where its own defining marker is present.

RULED IN, WIRED THROUGH THE ``--refine`` JOIN (status, 2026-09-23).  The maintainer ruled that
the pass "needs to be included in the cell typing for phenocycler data" and scoped it the same
day to one route: ``refine.py --out data/development/typing/refine``, then ``build_cells.py
--refine data/development/typing/refine``, which joins ``refined_class`` and this pass's
per-cell columns BESIDE ``cell_class``.  The islet chain's ``cell_class``, ``is_endocrine`` and
``endocrine_route`` are not touched, and this module reads no CD99, elongation or hormone
route: it demotes after ``subpanel.py``'s classes exist, on its own 54-marker neighbour
concordance.  ON DISK the join has not run until the DATA step -- no cell in ``typing/cells/``
carries ``Unresolved (no lineage evidence)`` or a filled ``refined_*`` column yet -- and every
builder still reads ``cell_class``; whether one should read the refined class instead is a
separate open question (``docs/PRODUCERS.md``, "Not decided here").

WHAT QUALIFIES EVERY COUNT THIS MODULE REPORTS: ``bound_X`` is the 5th percentile of a class's
own concordant members, so **a class with no low tail at all still loses about 5 % of it**.
Measured per class in ``REFINEMENT.md`` section 6: Lymphatic, Nerve, Immune and Delta are
demoted far above that floor -- the rule finding cells without their own evidence -- while
Exocrine's excess is only 0.2-0.5 percentage points, so essentially the whole Exocrine demotion
is the percentile relabelling a fixed share of the class by construction.  The ruling included
the pass without revisiting that floor; the module implements the recovered rule unchanged and
reports the floor in every log.

THE ORPHAN IT REPLACES IS CLOSED (2026-09-19).  The recovered run's outcome was written to
``data/development/typing/umap/refinement_log.json``, a file NOTHING in this repository
produces, and until 2026-09-19 ``trajectory_singlecell_features.py`` parameterised its entire
feature set from it -- superseded work used as a live input.  That file is now refused: this
module installs a process-wide guard at import, so any read of it from any code raises
:class:`SupersededRefinementLog` naming its date, the absence of a producer, the rejected step
2 it carries, its subsample counts against the cohort's, and what to read instead.  Nothing is
deleted; the record stays on disk and stops being an input.  :data:`OWN_MARKERS` here is the
committed authority for the defining-marker table, and it is identical to the one that file
carries (checked 2026-09-19, both steps, all eight classes).

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

    python core/multiplex5000/scripts/typing/refine.py --out <dir> [--donors 6450 6539]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
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

SHIFTS = pathlib.Path(os.environ.get(
    "IO60_SHIFTS_PATH", S.REPO / "data/development/celltype/shifts.json"))
#: How the pass reads CD8: with UniFORM's shift, like every other marker and like the primary calls.
#: Ruled 2026-10-19 for all 25 donors ("Do the shift"); the typing compares markers within a cell,
#: so one marker read in other units than the rest would skew those comparisons. ``unshifted``
#: divides the shift back out and exists only to reproduce the 2026-10-08 refinement.
CD8_ALIGNMENT = "aligned"
CELLS = S.OUT / "cells"

# --------------------------------------------------- the 2026-09-11 orphan, made unreadable
#
# THE ARTIFACT :data:`OWN_MARKERS` WAS RECOVERED FROM, AND WHY IT IS NOT AN INPUT.
#
# ``data/development/typing/umap/refinement_log.json`` was written 2026-09-11 13:38 by the
# recovered run, at a time when NO PRODUCER for it existed anywhere in the repository -- that
# is the failure this module was written to end.  It is kept on disk because it is the only
# surviving record of that run and of the step-2 rejection, and this lane deletes no data.  It
# is a RECORD.  It is not a parameter table, and until 2026-09-19 it was being read as one:
# ``trajectory_singlecell_features.py`` parameterised its whole feature set from it.
#
# Three separate reasons it may not be read as an input, each measured:
#
#   1. No producer.  Nothing in the repository writes it, so nothing can rebuild it, and a
#      parameter you cannot rebuild is not a parameter.  ``refine.py --out <dir>`` writes
#      ``refine_log.json``; that is a run's outcome and it has a producer.
#   2. It carries TWO steps and step 2 is ``"accepted": false`` -- rejected on the recovered
#      run, with the reason recorded in the file itself.  A reader that iterates its steps
#      reads the rejected one, and the old reader did exactly that.
#   3. Its counts are the 2026-09-11 320,000-cell subsample (Exocrine n=236,927, Immune
#      n=26,433, Blood vessel n=34,684, Nerve n=2,073).  The cohort is 25,230,109 cells over
#      20 donors in ``typing/cells/``, written 2026-09-15.  Quoting the file's numbers as
#      cohort numbers is off by two orders of magnitude.
#
# The defining-marker table it carries is identical to :data:`OWN_MARKERS` (checked
# 2026-09-19, both steps, all eight classes, marker for marker).  That is why closing this
# read changes no feature and no call -- and it is exactly why the read had to be closed while
# it still changes nothing, rather than after ``OWN_MARKERS`` moves.
#
# TWO GUARDS, NEITHER REDUNDANT.  Since 2026-09-19 ``subpanel_io`` hands out ``OUT`` and
# ``RUN`` as ``superseded.GuardedPath``, so every path this lane BUILDS from them refuses the
# registered artifacts, this file among them.  The hook below is the other half: it catches a
# literal path -- ``Path("/home/.../umap/refinement_log.json").read_text()``, a hard-coded
# string, a path from a config -- which never passes through a ``GuardedPath`` at all.  The
# two cover different routes to the same file and are deliberately kept separate.
SUPERSEDED_REFINEMENT_LOG = S.OUT / "umap/refinement_log.json"

SUPERSEDED_REFINEMENT_LOG_REFUSAL = (
    "Refusing data/development/typing/umap/refinement_log.json. "
    "WRITTEN 2026-09-11 13:38, and it HAS NO PRODUCER in this repository: nothing writes it, "
    "so nothing can rebuild it. It is a record of the recovered run, never an input. "
    "It also carries two steps and step 2 is \"accepted\": false -- rejected on that run -- so "
    "a reader that walks its steps reads a rejected step. Its counts are the 320,000-cell "
    "subsample (Exocrine n=236,927, Immune n=26,433, Blood vessel n=34,684, Nerve n=2,073); "
    "the cohort is 25,230,109 cells over 20 donors in data/development/typing/cells/, written "
    "2026-09-15. "
    "READ INSTEAD: core/multiplex5000/scripts/typing/refine.py -- OWN_MARKERS is the committed "
    "defining-marker table (identical to this file's, checked 2026-09-19) and "
    "UNRESOLVED/UNTOUCHED are the committed class vocabulary; "
    "core/multiplex5000/scripts/typing/REFINEMENT.md for what the pass does, what it measured at "
    "cohort scale and which decisions are still open; and the refine_log.json that "
    "`refine.py --out <dir>` writes for a run's outcome. "
    "The file stays on disk. It is not readable as an input."
)


class SupersededRefinementLog(RuntimeError):
    """Raised on any attempt to read the 2026-09-11 orphan refinement log.

    Carries :data:`SUPERSEDED_REFINEMENT_LOG_REFUSAL`, which names the file, its date, the
    fact that it has no producer, the rejected step it carries, its subsample counts against
    the cohort's, and what to read instead.  A silent fallback to the orphan is the defect
    this class exists to convert into a loud failure.
    """


def _is_superseded_refinement_log(path: object) -> bool:
    """True for the orphan log, however it is spelled -- str, bytes, Path, relative.

    Matched by SHAPE (``<anything>/umap/refinement_log.json``) rather than by one absolute
    path, so a copy of the file refuses too; a relative open is resolved against the working
    directory first, because ``open("refinement_log.json")`` from inside ``umap/`` is the same
    read.  The basename test is the only cost an unrelated ``open`` pays.
    """
    if isinstance(path, int):            # an already-open fd; nothing to match on
        return False
    try:
        text = os.fsdecode(path)
    except TypeError:
        return False
    if not text.replace("\\", "/").endswith("refinement_log.json"):
        return False
    p = pathlib.Path(text)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    return p.parent.name == "umap"


def refuse_superseded_refinement_log(path: object) -> None:
    """Refuse, by name, to read the 2026-09-11 orphan refinement log.

    The named authority for that refusal.  Other lanes import THIS rather than restating the
    path or the reason, so there is one place that knows why the file is closed.
    """
    if _is_superseded_refinement_log(path):
        raise SupersededRefinementLog(SUPERSEDED_REFINEMENT_LOG_REFUSAL)


#: `open` modes that only READ.  The audit hook refuses these and nothing else.
_READ_MODES = frozenset({"r", "rb", "rt", "br", "tr"})


def _refusal_audit_hook(event: str, args: tuple) -> None:
    """Refuse a READ of the orphan.  Never a write, and this is not a detail.

    The refusal's own premise is that the file HAS NO PRODUCER, and the remedy superseded.py
    names is "either commit the producer, or make the read refuse".  A hook that refused every
    mode made the first remedy impossible from any process that imports this module: measured
    2026-09-19, modes 'w', 'a', 'x' and 'r+' all raised and no file was created.  A guard that
    blocks the fix it demands is not a guard.
    """
    if event != "open":                  # one string compare on every other audited event
        return
    mode = args[1] if len(args) > 1 else "r"
    if mode is None:                     # os.open and friends pass no textual mode
        mode = "r"
    if isinstance(mode, str) and mode not in _READ_MODES:
        return                           # a write, an append or an update -- let it through
    refuse_superseded_refinement_log(args[0])


def install_refusal_hook() -> bool:
    """Make the orphan refuse to be read ANYWHERE in this process, not just here.

    A guard only the callers who remember to call it reach is not a closed door.  This adds a
    process-wide audit hook, so ``open()``, ``os.open()`` and ``pathlib.Path.read_text()`` on
    ``umap/refinement_log.json`` raise :class:`SupersededRefinementLog` from whatever code
    attempts them -- including code written later that never heard of this module.

    Installed at import, once, and NOT removable: an escape hatch is how the file would come
    back as an input.  Anything that genuinely needs the bytes (an archive, a forensic read of
    the recovered run) copies them outside Python or from a process that does not import this
    lane; the hook's blast radius is exactly the refinement lane, on purpose.

    Cost, measured 2026-09-19 on this machine: 3.1 % over 20,000 ``open()`` calls (0.187 s ->
    0.193 s), about 0.3 us per audited event.  Bulk reads here go through pyarrow's C++ IO and
    raise no ``open`` event at all.
    """
    global _HOOK_INSTALLED
    if _HOOK_INSTALLED:
        return False
    sys.addaudithook(_refusal_audit_hook)
    _HOOK_INSTALLED = True
    return True


_HOOK_INSTALLED = False
install_refusal_hook()


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
                       columns=["object_id", "cell_class", "cell_region", "islet_num",
                                "qc_analysis_eligible"]).to_pandas()
    core = ((df["cell_region"] == "core") & (~df["cell_class"].isin(STRUCTURAL))
            & (df["cell_class"] != "Excluded"))
    size = df.loc[core].groupby("islet_num")["object_id"].size()
    n = df["islet_num"].map(size)
    n = n.where(df["cell_region"] == "core", other=np.nan)
    df["islet_cells"] = n.to_numpy(dtype=float)
    df["size_stratum"] = S.size_stratum(df["islet_cells"].to_numpy())
    return df


def cd8_multiplier(donor: str, shifts: dict) -> float:
    """Require this exact donor's finite positive CD8 multiplier."""
    try:
        multiplier = float(shifts["shifts"]["CD8"]["scale"][str(donor)])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{donor}: matching CD8 alignment multiplier unavailable") from exc
    if not np.isfinite(multiplier) or multiplier <= 0:
        raise ValueError(f"{donor}: invalid CD8 alignment multiplier {multiplier}")
    return multiplier


def restore_unshifted_cd8(values: pd.DataFrame, donor: str, shifts: dict) -> pd.DataFrame:
    """Reverse only CD8's declared multiplier, preserving keys and missing measurements.

    The 2026-10-08 generation read shifted CD8 for primary calls and unshifted CD8
    throughout refinement; since 2026-10-19 refinement reads the shift too, and this
    survives only to reproduce that run. The matching alignment record is required: a donor
    multiplier from another expression arm is not interchangeable. Other channels
    and the primary aligned artifact remain unchanged.
    """
    if "CD8" not in values:
        return values
    multiplier = cd8_multiplier(donor, shifts)
    result = values.copy()
    result["CD8"] = result["CD8"] / multiplier
    return result


def scale_for_unshifted_cd8(scale: pd.DataFrame, donor: str, shifts: dict) -> pd.DataFrame:
    """Express frozen CD8 floor/span in the same units as the restored measurement.

    The immune fallback retains the frozen reference criterion. Dividing only the
    measurement would compare an unshifted value to a shifted reference. Dividing
    both by the same factor changes units without changing the criterion. The
    neighborhood representation is separately fitted from unshifted measurements.
    """
    if "CD8" not in scale.index:
        return scale
    multiplier = cd8_multiplier(donor, shifts)
    result = scale.copy()
    for field in ("floor", "top", "span"):
        if field in result.columns:
            result.loc["CD8", field] = result.loc["CD8", field] / multiplier
    return result


def donor_values(donor: str, markers: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """The donor's UniFORM-aligned values for ``markers``, and which of them it carries.

    Read through the celltype lane's own reader so availability cannot diverge from what the
    scoring that produced these classes saw.  ``ratio=False``: the recovered run measured the
    raw aligned values, and the DAPI ratio is a per-reference decision that has no meaning for
    a "does this cell carry its own marker at all" test.
    """
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    file = pq.ParquetFile(path)
    have = [m for m in markers if m in file.schema_arrow.names]
    df = file.read(columns=["object_id"] + have).to_pandas()
    if CD8_ALIGNMENT == "unshifted" and "CD8" in have:
        df = restore_unshifted_cd8(df, donor, json.loads(SHIFTS.read_text()))
    elif CD8_ALIGNMENT != "aligned" and CD8_ALIGNMENT != "unshifted":
        raise ValueError(f"unknown refinement CD8 alignment mode {CD8_ALIGNMENT!r}")
    geometry = pq.ParquetFile(S.cell_source(donor)[0]).read(
        columns=["object_id", "qc_analysis_eligible"]).to_pandas()
    df = df.merge(geometry, on="object_id", how="inner", validate="one_to_one")
    df = S.mask_fragment_hormones(df, donor)
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
    reqs = S.shared_negatives(R.requirements(raw, with_negatives=True))
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
    scale = S.load_scale()
    if CD8_ALIGNMENT == "unshifted" and "CD8" in have:
        scale = scale_for_unshifted_cd8(scale, donor, json.loads(SHIFTS.read_text()))
    Xs = S.apply_scale(t[have].to_numpy("float32"), have, scale)
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
    endocrine_uncertainty = ["Endocrine (unresolved subtype)",
                             "Endocrine (ambiguous subtype)"]
    endo = np.isin(labels, ["Beta", "Alpha", "Delta", *endocrine_uncertainty])
    return {
        "donor_id": donor,
        "rule": ("keep a class only if its own defining marker reaches the "
                 f"{BOUND_PCT}th percentile of that class's own concordant members "
                 f"(concordance >= {LOW}) in the same donor and size stratum; demoted cells "
                 "fall back to the complement"),
        "own_markers": {k: list(v) for k, v in OWN_MARKERS.items()},
        "endocrine_subtype_uncertainty": {
            "classes": endocrine_uncertainty,
            "n_before": int(np.isin(labels, endocrine_uncertainty).sum()),
            "n_after": int(np.isin(new, endocrine_uncertainty).sum()),
            "subtype_defining_marker_test_applied": False,
            "concordance_interpretation": (
                "neighbour agreement for a retained broad endocrine class; "
                "does not establish an Alpha, Beta or Delta subtype"),
        },
        "markers_used": list(have),
        "measurement_policy": {
            "cd8_alignment": CD8_ALIGNMENT,
            "aligned_root": str(S.ALIGNED),
            "shifts_path": str(SHIFTS),
            "shifts_sha256": hashlib.sha256(SHIFTS.read_bytes()).hexdigest()
            if SHIFTS.is_file() else None,
            "cd8_multiplier_reversed": float(json.loads(SHIFTS.read_text())
                                             ["shifts"]["CD8"]["scale"][str(donor)])
            if CD8_ALIGNMENT == "unshifted" and "CD8" in have else None,
            "scale_path": str(S.SCALE_CSV),
            "immune_fallback_cd8_scale": "frozen floor/top/span expressed in unshifted donor units"
            if CD8_ALIGNMENT == "unshifted" else "declared aligned scale",
            "primary_calls_unchanged": True,
        },
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
    global CD8_ALIGNMENT, SHIFTS
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, type=pathlib.Path,
                    help="directory to write donor_id=<d>/data.parquet and refine_log.json")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--params", type=pathlib.Path, default=None,
                    help="cached representation parameters; built and written if absent")
    ap.add_argument("--cd8-alignment", choices=["aligned", "unshifted"], default="aligned",
                    help="'aligned' is the method (ruled 2026-10-19: UniFORM's CD8 shift on all 25 "
                         "donors, as in the primary calls); 'unshifted' only reproduces the 2026-10-08 "
                         "refinement")
    ap.add_argument("--shifts", type=pathlib.Path, default=SHIFTS,
                    help="matching UniFORM multiplier record for this aligned input")
    args = ap.parse_args(argv)
    CD8_ALIGNMENT = args.cd8_alignment
    SHIFTS = args.shifts.resolve()
    if CD8_ALIGNMENT == "unshifted" and args.params:
        raise ValueError("unshifted CD8 requires fresh representation parameters; "
                         "omit --params to avoid reusing an aligned-CD8 representation")

    donors = args.donors or cohort_donors()
    representation_donors = cohort_donors()
    markers = derive_markers()
    print(f"refinement: {len(donors)} donors, {len(markers)} markers, K={K}, "
          f"bound = p{BOUND_PCT:g} of each class's concordant members per donor x stratum")

    using_cached_params = bool(args.params and args.params.exists())
    if using_cached_params:
        params = pd.read_parquet(args.params)
        print(f"  representation read from {args.params}")
    else:
        t0 = time.time()
        params = representation_params(representation_donors, markers)
        print(f"  representation built over {len(representation_donors)} donors "
              f"({time.time() - t0:.0f}s)")
        if args.params:
            args.params.parent.mkdir(parents=True, exist_ok=True)
            params.to_parquet(args.params)

    args.out.mkdir(parents=True, exist_ok=True)
    parameter_path = args.out / "representation_params.parquet"
    params.to_parquet(parameter_path)
    (args.out / "representation_provenance.json").write_text(json.dumps({
        "fit_state": "provided_cache" if using_cached_params else "fresh",
        "fit_donors": None if using_cached_params else representation_donors,
        "aligned_root": str(S.ALIGNED),
        "cd8_alignment": CD8_ALIGNMENT,
        "shifts_path": str(SHIFTS),
        "shifts_sha256": hashlib.sha256(SHIFTS.read_bytes()).hexdigest(),
        "scale_path": str(S.SCALE_CSV),
        "scale_sha256": hashlib.sha256(S.SCALE_CSV.read_bytes()).hexdigest(),
        "parameter_path": str(parameter_path),
        "parameter_sha256": hashlib.sha256(parameter_path.read_bytes()).hexdigest(),
        "source_script_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "immune_fallback_cd8_scale": "frozen reference converted by inverse matching donor multiplier"
        if CD8_ALIGNMENT == "unshifted" else "declared aligned scale",
    }, indent=2) + "\n")
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
