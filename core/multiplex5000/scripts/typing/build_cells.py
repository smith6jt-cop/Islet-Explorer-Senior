#!/usr/bin/env python3
"""Assemble the per-cell table every app artifact is aggregated from.

Joins the subpanel calls to the run's own coordinates and QC by ``(donor_id, object_id)`` -- never
row order, never a rounded centroid, never a filename -- and applies the corrections that the
calls on disk do not carry:

  G  no nucleus, no cell.  The calls already cover only ``qc_analysis_eligible`` cells, and a
     detection without a nucleus is never eligible, so this is satisfied by construction rather
     than by a filter.  The dropped count is recorded per donor so the denominator stays auditable.
  Spatial regions retain only true islet-core membership and ordinary tissue. Continuous
     nearest-islet distance stays available as a measurement; no expansion ring or replacement
     distance band is reconstructed here.
  C  which object a cell sits in, recorded.  The classification already uses object membership --
     including inside islets -- but ``context`` overwrites it with ``islet`` for every islet cell,
     so the evidence the call was made on is not in the output.  Kept as its own column.

THIS FILE CARRIES; IT DOES NOT DECIDE.  Every column the producer writes reaches the output --
that is the contract, and it is asserted in ``tests/current/test_build_cells_carry.py``.  A
producer column that this file quietly dropped would be indistinguishable downstream from a
decision that was never made, which is exactly how ``is_endocrine`` came to be a re-derivation.

    is_endocrine IS THE ISLET CHAIN'S DECISION, READ OFF ``resolved_by``.

It used to be ``cell_class in (Alpha, Beta, Delta)`` -- the NAME the hormone argmax wrote, not the
question the chain answered.  Under the CD99 / elongation chain those are two different facts: a
cell can be endocrine because CD99 cleared its donor's and stratum's level and then be NAMED Beta
by the hormone argmax, and a cell in a stratum whose CD99 level could not be fitted is endocrine
on the hormone question alone.  Both are endocrine; they were decided by different evidence, and
``endocrine_route`` is which.  Re-deriving the flag from the class name loses that and would go on
reading the same whatever the chain decided, so the flag is carried and the two are cross-checked:
a cell NAMED endocrine with no route, or carrying a route without the name, stops the build.

    A KNOWN-DEFECTIVE COLUMN IS CARRIED WITH ITS DEFECT NAMED, NOT REPAIRED HERE.

``top_identity_marker`` reads ``FOXP3`` on 21.1 % (115), 21.7 % (6442) and 25.3 % (6482) of the
deployed table's cells, which no Treg nuclear factor can be.  Measured 2026-09-18 over all
eligible cells of those three donors, recomputing the argmax from the aligned values and the
cohort scale (19.6 / 20.2 / 23.9 % over all 34 identity markers -- see the last point):

  * it is NOT a missing-value artefact.  0 cells have no identity marker measured, FOXP3 is
    missing on 0-1 cells per donor, and 0 cells are won by the only marker they carry.
  * it is NOT a tie or a saturation.  3.1-3.4 k cells per donor have any tie at the maximum and
    almost none of them is FOXP3; 0.01-0.03 % of cells saturate FOXP3 at the cohort top.  The
    median FOXP3 winner beats its runner-up by 0.142 scaled units -- a systematic offset, not a
    hairline.
  * it IS the cohort scale.  FOXP3's floor is 533.0 and its top 1108.1, a span of 575.1 -- the
    narrowest dynamic range of any identity marker.  A donor's own FOXP3 runs p50 508-513,
    p75 614-616, p99 831-854, so the "background" floor sits between its median and its third
    quartile and an ORDINARY cell scores 0.14 scaled where the same 1.15x-over-background reads
    0.003 on CD8 (floor 78.8, span 3432.8).  Under a per-marker rank representation, which
    removes exactly that span inflation, FOXP3 tops 4.2-5.6 % instead of 19.6-23.9 %.
  * what it is measuring is nuclear signal, not Treg identity.  On 115, log FOXP3 correlates
    0.788 with DAPI and 0.762 with PCNA against 0.372 with CD4, and FOXP3-topped cells are
    87.6 % Exocrine against a 77.1 % cohort share -- enriched in the bulk tissue, not in immune.
  * and the argmax does not range over what ``IDENTITY_MARKERS`` says it does.  ``classify``
    intersects it with the markers the run happened to LOAD, and CD31, Collagen_IV, Keratin_5,
    SMA, TP63 and Vimentin are in no subpanel and no immune rule, so 28 of the 34 compete.  On
    115 that is the difference between 181,557 FOXP3 wins and 168,772: 12,785 cells whose top
    identity marker is a marker that was never allowed to enter.

The column is READ BY THE PRODUCER -- ``subpanel._structure_calls`` withholds a structure call
when the cell's top identity marker is not one that drew the object -- so it is not cosmetic, and
it is not this file's to fix.  It is carried verbatim, and the per-donor FOXP3 share is written
into the build ledger so the defect is in front of whoever reads the build rather than in a note.

    python core/multiplex5000/scripts/typing/build_cells.py [--donors 6442 ...] [--include-marginal]
    python core/multiplex5000/scripts/typing/build_cells.py --refine <dir> --audit <audit.json>
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import islet_membership as MEMBERSHIP
import subpanel_io as S

DISTANCE = S.REPO / "data/development/kintsugi/distance"
OUT = S.OUT / "cells"

#: The label ``refine.UNRESOLVED`` gives a demoted Exocrine cell.  Restated here rather than
#: imported because importing ``refine`` pulls in the whole celltype scoring stack, and asserted
#: equal to it in the tests so the two cannot drift.  ``REFINEMENT.md``'s open decision 4 is that
#: nothing downstream of ``typing/cells`` knows this value; this is that file being told.
UNRESOLVED = "Unresolved (no lineage evidence)"

#: Display label per call, and the region a cell is reported in.  ``Excluded`` keeps its reason.
DISPLAY = {
    "Alpha": "Alpha_cell", "Beta": "Beta_cell", "Delta": "Delta_cell",
    "Blood vessel": "Blood_Vessel", "Lymphatic": "Lymphatic", "Nerve": "Neural",
    "Immune": "Immune", "Exocrine": "Exocrine", "Excluded": "Excluded",
    "Endocrine (unresolved subtype)": "Endocrine_unresolved",
    "Endocrine (ambiguous subtype)": "Endocrine_ambiguous",
    UNRESOLVED: "Unresolved",
}
ENDOCRINE = ("Alpha", "Beta", "Delta", "Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)")

#: The values ``subpanel.classify`` writes into ``resolved_by`` when the ISLET CHAIN calls a cell
#: endocrine -- one per route, and the route is the evidence:
#:
#:   cd99        CD99 cleared the level fitted from this donor's and this stratum's
#:               hormone-confirmed cells.  The maintainer's key differentiator.
#:   elongation  CD99 did not reach the cell, and its nucleus is as elongated as that stratum's
#:               confirmed endocrine cells.  Large stratum only, by construction upstream.
#:   hormone     the stratum could supply no CD99 level, so the cell kept the hormone question.
#:               Also the only route the pre-chain producer and ``--no-islet-chain`` ever write.
#:
#: ``is_endocrine`` is exactly "one of these three is written on this cell".
ENDOCRINE_ROUTES = ("cd99", "elongation", "hormone", "single_hormone")

#: Everything else ``resolved_by`` can carry.  A value in NEITHER tuple stops the build: a route
#: added upstream must not be able to arrive here and be silently counted as not-endocrine.
NON_ENDOCRINE_RESOLVED = ("structure", "definition", "supporting", "fallthrough",
                         "excluded_single_hormone")

#: The producer's per-cell islet number is kept, under its own name, beside this file's.
#:
#: They are the same call -- ``islet_membership.membership`` over the outlines -- when the
#: producer ran ``--membership geometry``, and they are NOT the same call when it ran
#: ``--membership parent``, which QuPath mislabels on 2.7 % of islet-core cells.  ``islet_num``
#: stays this file's geometric answer, because ``islet_key``, ``combined_islet_id`` and every
#: consumer of them are built on it.  ``typing_islet_num`` is the number the SIZE and so the
#: STRATUM were counted on, which is the only thing that can say whether a stratum was counted
#: on the same islets the cell is reported in.  Dropping either would make that unanswerable.
CARRY_RENAME = {"islet_num": "typing_islet_num"}

#: Columns the producer and this file both write, reconciled rather than suffixed.  Both read
#: ``cell_region`` out of the same export, so a disagreement means the two were built against
#: different exports and the build stops.
RECONCILED = ("cell_region",)

#: This file's own columns.  None of them may arrive from the producer: if one does, the producer
#: has started deciding something this file also decides, and that has to be ruled on, not merged.
DERIVED = ("X_centroid", "Y_centroid", "qc_analysis_eligible", "islet_num",
           "signed_dist_islet_um", "phenotype", "is_endocrine",
           "endocrine_route", "islet_key", "combined_islet_id", "acq_run", "stain_group",
           "islet_cd99_level", "islet_elongation_level", "refined_phenotype")

#: The refinement sidecar's own columns, carried verbatim.  ``refine.py`` writes these keyed on
#: ``object_id``; it deliberately does not rewrite ``typing/cells``, so joining them here is what
#: puts the pass in front of the app instead of in a directory nobody reads.
REFINE_CARRY = ("refined_class", "refine_note", "refine_pass", "demoted",
                "concordance", "concordance_before", "immune_leaf")

#: Sidecar column -> the producer column that already answers the same question.
#:
#: Where the producer answered it too, the two are CHECKED and the sidecar's is dropped only if
#: they AGREE: two columns meant to say the same thing are how a table comes to disagree with
#: itself.  Where they disagree the sidecar's is kept, prefixed, because a disagreement about the
#: stratum a bound was fitted in is exactly what a reader would need.  Where the producer never
#: answered it -- the pre-islet-chain calls carry no ``islet_stratum`` -- the sidecar's column is
#: left entirely alone, because a column dropped for being redundant with nothing is a loss.
REFINE_CHECK = {"cell_class": "cell_class", "size_stratum": "islet_stratum",
                "islet_cells": "islet_cells"}

#: ``subpanel.SHAPE_FEATURE``.  Restated for the same reason as ``UNRESOLVED`` and asserted equal.
SHAPE_FEATURE = "nucleus_elongation"


def _manifest() -> tuple[dict, dict]:
    raw = json.loads(S.QUPATH_MANIFEST.read_text())
    pixel = {r["donor_id"]: r["pixel_size_um_x"] for r in raw["images"]}
    stem = {r["donor_id"]: r["image_id"].split(" - ")[0].replace(".qptiff", "")
            for r in raw["images"]}
    return pixel, stem


def islet_levels(audit: pathlib.Path | str, donor: str) -> dict:
    """The per-stratum levels one donor's islet chain actually fitted, from its audit.

    ``subpanel.py --audit-out`` is where the fitted level lives; it is not a per-cell column, so
    a cell cannot otherwise say what it was measured against.  A stratum whose level is None came
    back UNAVAILABLE upstream and stays unavailable here -- never 0.0, which for a scaled marker
    is a level every measured cell clears.
    """
    entries = json.loads(pathlib.Path(audit).read_text())
    if isinstance(entries, dict):
        entries = [entries]
    for entry in entries:
        if str(entry.get("donor_id")) == str(donor):
            return entry.get("islet_levels") or {}
    raise KeyError(f"{donor}: no entry in {audit} -- the audit is from a different run")


def _carry(calls: pd.DataFrame, donor: str) -> pd.DataFrame:
    """The producer's table, renamed where this file also writes a column of that name."""
    if "resolved_by" not in calls.columns:
        raise AssertionError(
            f"{donor}: the calls carry no `resolved_by`, so the islet chain's decision is not in "
            "them and `is_endocrine` could only be re-derived from the class NAME -- which is the "
            "defect this build was corrected for. Rerun subpanel.py with --emit-detail.")
    taken = set(calls.columns) & set(CARRY_RENAME.values())
    if taken:
        raise AssertionError(
            f"{donor}: the producer already writes {sorted(taken)}, which is the name this file "
            "renames its own column to. Two columns would end up sharing one name.")
    carried = calls.rename(columns=CARRY_RENAME)
    unexpected = set(carried.columns) & (set(DERIVED) | set(REFINE_CARRY))
    if unexpected:
        raise AssertionError(
            f"{donor}: the producer now writes {sorted(unexpected)}, which this file also writes. "
            "Two answers to one question must be ruled on, not merged silently.")
    return carried


def _join_refinement(frame: pd.DataFrame, donor: str,
                     refine: pathlib.Path) -> tuple[pd.DataFrame, dict]:
    """Join the refinement sidecar onto the typed cells.  Returns (frame, what it found).

    The pass is a sidecar on purpose -- a cell it could not measure keeps its class and carries
    the reason -- so this join ADDS ``refined_class`` beside ``cell_class`` and never overwrites
    it.  ``is_endocrine`` is untouched: it is the islet chain's decision, and the refinement's
    verdict on the same cell is ``demoted`` / ``refined_class``.  Whether a demoted endocrine
    cell should still count as endocrine downstream is a maintainer's ruling, not a merge.

    The sidecar restates three things the producer may also have decided -- the class, the
    stratum and the islet size.  Where both answered, the two are CHECKED: identical, the
    sidecar's goes, because two columns meant to say the same thing are how a table comes to
    disagree with itself; different, the sidecar's is KEPT under a prefix, because that
    disagreement is the reader's whole answer.  Where only the sidecar answered -- the
    pre-islet-chain calls carry no stratum -- its column is left alone.  A disagreement about the
    CLASS stops the build: it means the pass ran against a different generation of the typing,
    and this tree stacks generations.
    """
    path = pathlib.Path(refine) / f"donor_id={donor}/data.parquet"
    if not path.exists():
        raise FileNotFoundError(f"{donor}: no refinement partition at {path}")
    side = pq.read_table(path).to_pandas()
    missing = [c for c in REFINE_CARRY if c not in side.columns]
    if missing:
        raise AssertionError(f"{donor}: the refinement sidecar is missing {missing}")

    restated = [c for c in REFINE_CHECK if c in side.columns]
    wanted = ["object_id"] + list(REFINE_CARRY) + restated
    merged = frame.merge(side[wanted], on="object_id", how="left", validate="one_to_one",
                         suffixes=("", "__refine"))
    absent = int(merged["refine_pass"].isna().sum())
    if absent:
        raise AssertionError(
            f"{donor}: {absent:,} typed cells have no row in the refinement sidecar -- the pass "
            "was run against a different generation of the typing")

    found: dict = {"refine_demoted": int(merged["demoted"].fillna(False).to_numpy(bool).sum())}
    for their in restated:
        ours = REFINE_CHECK[their]
        if ours not in frame.columns:
            # The producer never wrote the counterpart, so there is nothing to check against and
            # nothing that would be said twice. The sidecar's column stays: dropping a column
            # whose only justification for being dropped is that something else says the same
            # thing, when nothing else says it, is the leak this whole file exists to prevent.
            continue
        theirs = f"{their}__refine" if their in frame.columns else their
        a, b = merged[theirs], merged[ours]
        numeric = (pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b))
        if numeric:
            differ = int((~np.isclose(a.to_numpy(float), b.to_numpy(float),
                                      equal_nan=True)).sum())
        else:
            differ = int((a.fillna("").astype(str) != b.fillna("").astype(str)).sum())
        found[f"refine_{their}_differs"] = differ
        if their == "cell_class" and differ:
            raise AssertionError(
                f"{donor}: the refinement sidecar disagrees with the typing about `cell_class` on "
                f"{differ:,} cells -- it was run on a different generation of the typing")
        # Identical: one of the two is redundant and the producer's is the authority. Different:
        # the sidecar's is KEPT under its own prefixed name, because a disagreement about the
        # stratum a bound was fitted in is the one thing a reader would need and the one thing
        # dropping it would hide.
        merged = (merged.drop(columns=[theirs]) if differ == 0
                  else merged.rename(columns={theirs: f"refine_{their}"}))
    return merged, found


#: What ``refine_note`` says on a build the refinement was not joined into.
#:
#: NOT APPLIED IS A STATE AND IT BELONGS ON THE ROW.  The pass reads ``typing/cells`` and writes a
#: sidecar, so the FIRST build of a generation cannot have it -- the flag cannot be required.  But
#: a table with no ``refined_class`` column at all is indistinguishable from one where the pass ran
#: and demoted nothing, and it is the table every consumer would then serve.  So the columns are
#: written either way and the absent case says so: ``refine_pass`` 0, and everything the pass would
#: have decided left UNAVAILABLE rather than filled in with the unrefined answer.  Writing
#: ``refined_class = cell_class`` here would assert that the pass looked at the cell and kept it.
NOT_REFINED = "the refinement pass was not applied to this build"


def _no_refinement(frame: pd.DataFrame) -> pd.DataFrame:
    """The refinement columns on a build made without ``--refine``: unavailable, and saying so."""
    n = len(frame)
    frame["refined_class"] = pd.array([None] * n, dtype="string")
    frame["refine_note"] = NOT_REFINED
    frame["refine_pass"] = np.zeros(n, dtype="int64")
    frame["demoted"] = pd.array([None] * n, dtype="boolean")
    frame["immune_leaf"] = pd.array([None] * n, dtype="boolean")
    frame["concordance"] = np.full(n, np.nan)
    frame["concordance_before"] = np.full(n, np.nan)
    frame["refined_phenotype"] = pd.array([None] * n, dtype="string")
    return frame


def build_donor(donor: str, pixel: dict, stem: dict, runs: pd.DataFrame,
                refine: pathlib.Path | None = None,
                audit: pathlib.Path | None = None,
                calls_root: pathlib.Path | None = None) -> pd.DataFrame:
    # The corrected typing, written per donor by `subpanel.py`. It carries `object_context`,
    # `immune_definition`, `resolved_by`, `top_identity_marker` and the `passed_through_*`
    # FRACTIONS -- how much of the cell each object covers, float and NaN-where-unmeasured since
    # 2026-09-19; they were booleans before, and a consumer that tests them for truth will read
    # every non-zero overlap as a passthrough,
    # all of which the old monolithic table collapsed away or never had -- and, since the islet
    # chain, `islet_cells`, `islet_stratum`, `islet_route` and `excluded_reason` on every row.
    # Read with no column list: a column added upstream must arrive here by default, because a
    # list is a place a new column goes to be forgotten.
    root = pathlib.Path(calls_root) if calls_root else S.TYPING
    calls = pq.read_table(root / f"donor_id={donor}/data.parquet").to_pandas()
    carried = _carry(calls, donor)

    # Cell geometry from whichever export is authoritative for this donor -- the revised rebuild
    # for the seven re-exported ones. Reading the run unconditionally here took coordinates and
    # islet numbers from the defective export for exactly those donors.
    source, which = S.cell_source(donor)
    expression = pq.read_table(
        source, columns=["object_id", "X_centroid", "Y_centroid",
                         "cell_region", "qc_analysis_eligible"]).to_pandas()

    if which == "revised":
        geometry = pq.read_table(source, columns=["object_id", "analysis_reason"]).to_pandas()
        distance = pq.read_table(source,
                                 columns=["object_id", "signed_dist_islet_um"]).to_pandas()
    else:
        geometry = pq.read_table(
            S.RUN / f"geometry_qc/donor_id={donor}/data.parquet",
            columns=["object_id", "analysis_reason"]).to_pandas()
        distance = pq.read_table(
            DISTANCE / f"donor_id={donor}/data_0.parquet",
            columns=["object_id", "signed_dist_islet_um"]).to_pandas()

    # Islet membership from the OUTLINES, the same call the typing made. Taking `islet_num` from
    # QuPath's `Parent` here would have re-imported the defect the typing just removed: the cell
    # would carry a corrected class and the wrong islet.
    member = MEMBERSHIP.membership(donor, expression["X_centroid"].to_numpy(),
                                   expression["Y_centroid"].to_numpy(),
                                   pixel[donor], stem[donor])
    expression["islet_num"] = member["islet_num"].to_numpy()

    # --- correction G, as an accounting not a filter -------------------------------------
    no_nucleus = int((geometry["analysis_reason"] == "missing_nucleus_geometry").sum())
    ineligible = int((~expression["qc_analysis_eligible"]).sum())
    if not expression.loc[expression["object_id"].isin(calls["object_id"]),
                          "qc_analysis_eligible"].all():
        raise AssertionError(f"{donor}: a typed cell is not analysis-eligible -- correction G "
                             "cannot be satisfied by construction, so it must become a filter")

    collide = (set(carried.columns) & set(expression.columns)) - {"object_id"}
    if collide - set(RECONCILED):
        raise AssertionError(
            f"{donor}: unreconciled column(s) {sorted(collide - set(RECONCILED))} in both the "
            "calls and the geometry. Name which one wins rather than letting a suffix decide.")
    frame = carried.merge(expression, on="object_id", how="left", validate="one_to_one",
                          suffixes=("__typing", ""))
    for column in RECONCILED:
        theirs = f"{column}__typing"
        if theirs not in frame.columns:
            continue
        differ = int((frame[theirs].fillna("") != frame[column].fillna("")).sum())
        if differ:
            raise AssertionError(
                f"{donor}: the calls and the geometry disagree about `{column}` on {differ:,} "
                f"cells. Both read it from {source}, so they were built against different "
                "exports and one of them is stale.")
        frame = frame.drop(columns=[theirs])

    frame = frame.merge(distance, on="object_id", how="left", validate="one_to_one")
    if frame["X_centroid"].isna().any():
        raise AssertionError(f"{donor}: {int(frame['X_centroid'].isna().sum())} typed cells have "
                             "no coordinate -- the join on object_id did not resolve")
    unmatched = int(frame["signed_dist_islet_um"].isna().sum())
    if unmatched:
        raise AssertionError(f"{donor}: {unmatched} typed cells have no signed distance")

    # Spatial membership has no peri-islet ring. Keep the continuous signed distance and
    # distinguish true islet core from surrounding tissue without binning it.
    d = frame["signed_dist_islet_um"].to_numpy()
    frame["cell_region"] = pd.Series(np.where(d <= 0, "core", "tissue"), index=frame.index).where(np.isfinite(d))

    # --- correction C: the object the call was actually made on ---------------------------
    # Carried through from the typing rather than recomputed. Recomputing it invited the two
    # tables to disagree about the same cell, and the classification is the one that used it.

    frame["phenotype"] = frame["cell_class"].map(DISPLAY)
    if frame["phenotype"].isna().any():
        missing = sorted(frame.loc[frame["phenotype"].isna(), "cell_class"].unique())
        raise KeyError(f"{donor}: no display label for {missing}")

    # --- the islet chain's decision, CARRIED -----------------------------------------------
    resolved = frame["resolved_by"].fillna("").astype(str)
    unknown = sorted(set(resolved.unique()) - set(ENDOCRINE_ROUTES)
                     - set(NON_ENDOCRINE_RESOLVED) - {""})
    if unknown:
        raise AssertionError(
            f"{donor}: `resolved_by` carries {unknown}, which this file has never been told "
            "about. A new route must be declared endocrine or not -- an unknown one would be "
            "counted as not-endocrine by default, which is the leak this check exists for.")
    frame["endocrine_route"] = np.where(resolved.isin(ENDOCRINE_ROUTES), resolved, "")
    frame["is_endocrine"] = frame["endocrine_route"].to_numpy() != ""
    named = frame["cell_class"].isin(ENDOCRINE).to_numpy()
    if int((named & ~frame["is_endocrine"].to_numpy()).sum()):
        n = int((named & ~frame["is_endocrine"].to_numpy()).sum())
        raise AssertionError(
            f"{donor}: {n:,} cells are NAMED {ENDOCRINE} with no endocrine route on them. The "
            "hormone argmax names a cell; it does not decide one, and a name with no decision "
            "behind it is the pre-chain behaviour arriving through the back door.")
    if int((~named & frame["is_endocrine"].to_numpy()).sum()):
        n = int((~named & frame["is_endocrine"].to_numpy()).sum())
        raise AssertionError(
            f"{donor}: {n:,} cells carry an endocrine route but were not named endocrine.")

    # islet_num is a string in the run, empty for a cell in no islet -- not a nullable int.
    islet_num = frame["islet_num"].fillna("").astype(str).str.strip()
    frame["islet_num"] = islet_num.replace("", None)
    frame["islet_key"] = np.where(islet_num != "", "Islet_" + islet_num, None)
    frame["combined_islet_id"] = np.where(islet_num != "", f"{donor}_Islet_" + islet_num, None)

    run = runs.loc[runs["donor_id"] == donor]
    if run.empty:
        raise KeyError(f"{donor}: absent from donor_acquisition_runs.csv")
    frame["acq_run"] = run["acq_run"].iloc[0]
    frame["stain_group"] = run["stain_group"].iloc[0]

    # --- the level each cell was measured against, where the audit says ---------------------
    if audit is not None:
        if "islet_stratum" not in frame.columns:
            raise AssertionError(
                f"{donor}: --audit needs `islet_stratum` on the calls, and these calls predate "
                "the islet chain. Rebuild the typing before attaching its levels.")
        levels = islet_levels(audit, donor)
        # The audit also carries a donor-wide reference level beside the nested size strata.
        # Keep only the stratum mappings; the donor-wide level is not the level an islet cell
        # was measured against.
        for column, key in (("islet_cd99_level", "CD99"),
                            ("islet_elongation_level", SHAPE_FEATURE)):
            by_stratum = {s: v.get(key) for s, v in levels.items()
                          if isinstance(v, dict)}
            frame[column] = pd.to_numeric(frame["islet_stratum"].map(by_stratum),
                                          errors="coerce")

    # --- the refinement pass, where it has been run -----------------------------------------
    found: dict = {}
    if refine is not None:
        frame, found = _join_refinement(frame, donor, pathlib.Path(refine))
        frame["refined_phenotype"] = frame["refined_class"].map(DISPLAY)
        if frame["refined_phenotype"].isna().any():
            missing = sorted(frame.loc[frame["refined_phenotype"].isna(),
                                       "refined_class"].dropna().unique())
            raise KeyError(f"{donor}: no display label for refined class {missing}")
        found["refine_applied"] = True
    else:
        frame = _no_refinement(frame)
        # No `refine_demoted: 0` here on purpose: nothing was demoted because nothing was
        # asked, and a 0 in that column would read as a measured zero.
        found = {"refine_applied": False}

    # Compare the islet identity carried by the typing call with the current QuPath object
    # membership. Continuous nearest-islet distance remains a separate measurement; no distance
    # band is reconstructed here.
    counts = {"islet_num_only_here": 0, "islet_num_only_in_typing": 0, "islet_num_conflicts": 0}
    if "typing_islet_num" in frame.columns:
        theirs = frame["typing_islet_num"].fillna("").astype(str).str.strip()
        ours = frame["islet_num"].fillna("").astype(str).str.strip()
        counts["islet_num_only_here"] = int(((ours != "") & (theirs == "")).sum())
        counts["islet_num_only_in_typing"] = int(((ours == "") & (theirs != "")).sum())
        counts["islet_num_conflicts"] = int(((ours != "") & (theirs != "")
                                             & (ours != theirs)).sum())

    frame.attrs["dropped_no_nucleus"] = no_nucleus
    frame.attrs["dropped_ineligible"] = ineligible
    frame.attrs.update(counts)
    for route in ENDOCRINE_ROUTES:
        frame.attrs[f"endocrine_{route}"] = int((frame["endocrine_route"] == route).sum())
    frame.attrs["top_marker_FOXP3"] = int(
        (frame["top_identity_marker"] == "FOXP3").sum()) if "top_identity_marker" in frame else 0
    frame.attrs.update({k: v for k, v in found.items()})
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--calls", default=None,
                        help="the subpanel calls to carry (default: the typing stage root)")
    parser.add_argument("--out", default=None,
                        help="output directory (default: the typing/cells stage root)")
    parser.add_argument("--refine", default=None,
                        help="a refine.py output directory, joined on object_id so the "
                             "refinement's verdict reaches the app instead of a sidecar")
    parser.add_argument("--audit", default=None,
                        help="subpanel.py --audit-out JSON, for the per-stratum level each "
                             "cell's endocrine question was measured against")
    args = parser.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    pixel, stem = _manifest()
    runs_path = pathlib.Path(os.environ.get(
        "IO60_DONOR_ACQUISITION_RUNS", S.REPO / "data/donor_acquisition_runs.csv"))
    runs = pd.read_csv(runs_path, dtype={"donor_id": str})

    out = pathlib.Path(args.out) if args.out else OUT
    out.mkdir(parents=True, exist_ok=True)
    ledger = []
    for donor in donors:
        frame = build_donor(donor, pixel, stem, runs, refine=args.refine, audit=args.audit,
                            calls_root=args.calls)
        target = out / f"donor_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pandas(frame, preserve_index=False),
                       target / "data.parquet", compression="zstd")
        counts = frame["cell_region"].value_counts()
        row = {
            "donor_id": donor, "typed_cells": len(frame),
            "dropped_no_nucleus": frame.attrs["dropped_no_nucleus"],
            "dropped_ineligible": frame.attrs["dropped_ineligible"],
            "core": int(counts.get("core", 0)),
            "tissue": int(counts.get("tissue", 0)),
            "in_object": int((frame["object_context"] != "tissue").sum()),
            "endocrine": int(frame["is_endocrine"].sum()),
        }
        row.update({k: v for k, v in frame.attrs.items()
                    if k not in ("dropped_no_nucleus", "dropped_ineligible")})
        ledger.append(row)
        print(f"  {donor}: {len(frame):>9,} cells "
              f"(islet core {row['core']:,} / surrounding tissue {row['tissue']:,})  endocrine "
              + " / ".join(f"{r} {frame.attrs[f'endocrine_{r}']:,}" for r in ENDOCRINE_ROUTES),
              flush=True)

    book = pd.DataFrame(ledger)
    book.to_csv(out / "build_ledger.csv", index=False)
    print(f"\n{len(donors)} donors, {book.typed_cells.sum():,} cells -> {out}")
    print(f"excluded and recorded: {book.dropped_no_nucleus.sum():,} without a nucleus, "
          f"{book.dropped_ineligible.sum():,} ineligible overall")
    conflicts = int(book.get("islet_num_conflicts", pd.Series([0])).sum())
    orphans = int(book.get("islet_num_only_in_typing", pd.Series([0])).sum())
    if conflicts or orphans:
        print(f"ISLET NUMBER: {conflicts:,} cells are reported in one islet and had their size "
              f"stratum counted in another, and {orphans:,} are numbered by the typing and not "
              "here. The typing was not built on this file's membership.")
    if "top_marker_FOXP3" in book:
        share = 100.0 * book.top_marker_FOXP3.sum() / max(int(book.typed_cells.sum()), 1)
        print(f"CARRIED WITH A KNOWN DEFECT: top_identity_marker reads FOXP3 on "
              f"{book.top_marker_FOXP3.sum():,} cells ({share:.2f}%). That is the cohort scale "
              "(FOXP3 floor 533.0, top 1108.1), not Treg biology -- see this file's docstring.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
