#!/usr/bin/env python3
"""The compartment gradient: islet core, peri-islet expansion, exocrine tissue, pLN.

Part 3 of the study question -- immune cells near and in islets against those in
the pancreatic lymph node and in the exocrine tissue.

What this builds
----------------
``$ISLET_LAYERS_DIR/compartments`` -- defaulting to ``/data/islet_layers`` and
overridable with ``--out`` -- holds, per donor:

* ``compartment_composition.csv`` -- immune fraction and immune-subtype
  composition for every compartment, one row per donor x compartment x subtype.
  ``n_lineage_unavailable`` on every row is the count of eligible cells in that
  compartment whose lineage was never written (a ``cell_class`` that is null,
  string-shaped missing, the declared absence, or outside the declared
  vocabulary in the pancreas; a null ``is_immune`` in the pLN); they are outside
  ``n_cells`` and ``lineage_unavailable_reason`` says why, by kind and count.
* ``islet_size_strata.csv`` -- the same core-compartment numbers recomputed on
  islets of at least ``--min-islet-cells`` core cells, which is the stratum in
  which an islet-level statement is defensible at all.
* ``islet_size_summary.csv`` -- per donor, the islet count, how many clear
  ``--min-islet-cells``, the median core cells per islet, the one-cell fraction
  and the ``fragmented_islets`` flag, so a pooled per-islet number can be read
  against its own denominator.
* ``batch_eta_squared.csv`` -- eta-squared of every cross-donor quantity on
  ``stain_group``, ``acq_run`` and ``donor_status``, each against its own
  permutation null.
* ``manifest.json`` -- the reasons, the denominators and the withheld donors.

It is NOT written under ``data/development/integration/``.  That is where it
used to land, and no card reads there: ``apps/islet/R/card_s3_compartments.R``
resolves ``$ISLET_LAYERS_DIR``, and a layer nothing can read is the whole
failure.

Where the pancreas subtype comes from
-------------------------------------
The pancreas immune subtype is ``immune_definition`` in
``data/development/typing/cells`` -- the per-marker immune definition, carried
on all 20 donors of the 2026-09-19 build.  ``data/development/typing/immune`` is
NOT read here: that second pass is a finer subtype (rewritten 2026-09-17, 20
donors, 2,934,898 Immune cells), and it never modifies ``cell_class``.
``refined_class`` is not read either: the refinement pass was ruled in on
2026-09-23 and is wired as the ``build_cells.py --refine`` join, beside
``cell_class``; the join has not run, and which class a builder serves is open.

That build carries eight classes plus ``Excluded``, with the 2026-09-18 islet
chain applied to every donor; no cell carries ``Unresolved``.

Three compartments are pancreas and one is not
----------------------------------------------
The pancreas compartments are **spatial**: core, the 20 um peri-islet expansion
and the remaining exocrine tissue are cut out of one section by one signed
distance.  The pLN compartment is **compositional only**.  A lymph node is a
different block from the pancreas beside it -- for most donors a different slide
entirely -- so it has no islet, and a signed distance to one is undefined there
rather than zero.  The pLN column joins to the pancreas columns on ``donor_id``
and ``donor_status`` and on nothing else.  No distance, adjacency,
neighbourhood or KDE crosses the tissue boundary, and nothing is co-plotted in
space across tissues.  ``refuse_cross_tissue`` makes that structural: asking
this module for a pLN distance raises.

  One exception is coming and is not inferred here.  The maintainer is
  annotating the pLN portions of pancreas slides as a distinct dataset, and some
  of those lymphoid structures are anatomically attached to the section, where
  an islet-to-lymphoid distance *is* meaningful.  That classification arrives as
  data.  ``LYMPHOID_ATTACHMENT_COLUMN`` names the column this reader will honour
  when it lands; until then every pancreas/pLN pair is treated as separate.

Compartments come from the signed distance, never from an annotation
--------------------------------------------------------------------
Every expansion is 20 um.  One donor's drawn expansion annotation measures
23.1 um against 18.8-19.3 for every other, and all of them are named
``_exp20um``, so the annotation cannot be trusted to mean what it is called.
The compartment is therefore cut from ``signed_dist_islet_um``:

    core       signed distance <= 0
    peri       0 < signed distance <= 20
    exocrine   signed distance > 20

QuPath's ``Parent`` is never read.  Where ``typing/islet_membership`` exists the
geometric outline membership is read back as an audit and the disagreement is
reported per donor -- it is an audit, not an override.

Two stratifications, both measured, both applied before anything pools
----------------------------------------------------------------------
1. The endocrine shape signal is real only above roughly ten cells per islet, so
   every core-compartment quantity is emitted a second time over islets with at
   least ``--min-islet-cells`` core cells.
2. 6539's islets are mostly one-cell classifier fragments.  It is measured, it
   is reported, and it is flagged ``fragmented_islets`` so that no cohort-pooled
   per-islet number silently carries it.

Neither stratification drops a donor.  A cohort-pooled per-islet number without
both is wrong; a donor withheld for being inconvenient is worse.

Labels
------
PhenoCycler pancreas and PhenoCycler pLN are named by different lanes, so their
labels are harmonised through
``core/phenocycler/data/integration/vocab_crosswalk.csv`` (generated; never
hand-edited) to the common lineages ``T_NK``, ``B_Plasma``, ``Myeloid`` and
``Granulocyte``.  A lineage the crosswalk does not classify is withheld by name
and never defaulted to non-immune.

Negative evidence is evidence, and "unavailable" is reserved for a measurement
that does not exist.  A donor whose pLN export carries no intensity measurement
has nothing to cluster and nothing to name; that is unavailable with a reason,
and it is never zero immune cells.
"""

from __future__ import annotations

import argparse
import json
import sys
import os
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]

from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

sys.path.insert(0, str(REPO / "core" / "multiplex5000" / "scripts" / "typing"))

import build_cells as TYPING_CELLS  # noqa: E402  -- the declared cell_class table

sys.path.insert(0, str(REPO))

from scripts.senior import crosswalk as XW  # noqa: E402
from scripts.senior import xenium_release  # noqa: E402


@lru_cache(maxsize=None)
def xenium_source(tissue: str):
    """The current Xenium release for a tissue, resolved from its own declaration.

    This file used to pin ``release_pln_immune_v2``.  v2's own successor records
    ``supersedes: .../release_pln_immune_v2`` in its ``release.json``, so the pin
    outlived the release it named; the resolver honours that declaration.
    """
    return xenium_release.resolve(tissue)

CELLS = REPO / "data" / "development" / "typing" / "cells"
MEMBERSHIP = REPO / "data" / "development" / "typing" / "islet_membership"
PLN_PHENO = REPO / "data" / "development" / "typing" / "pln" / "pheno"
DONOR_META = REPO / "data" / "development" / "typing" / "app" / "donor_metadata.csv"
ACQ_RUNS = REPO / "data" / "donor_acquisition_runs.csv"
CROSSWALK = REPO / "core" / "phenocycler" / "data" / "integration" / "vocab_crosswalk.csv"
#: The layers root its READER uses. apps/islet/R/card_s3_compartments.R resolves
#: $ISLET_LAYERS_DIR defaulting to /data/islet_layers, and every other layer already sits
#: there. Writing beside the source trees instead put this artifact where no card looks --
#: the same defect the pLN node layer had, and a layer nothing can read is the whole failure.
OUT = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")) / "compartments"

#: The maintainer's pLN-attachment classification is not in the checkout yet.
#: When it arrives it arrives as a column on the pLN section table; this is the
#: name this reader will honour.  Nothing infers it.
LYMPHOID_ATTACHMENT_COLUMN = "attached_to_pancreas_section"

EXPANSION_UM = 20.0
STATUS_ORDER = ["ND", "Aab+", "T1D"]

#: A cell whose lineage was never written is not non-immune.  It leaves the
#: compartment denominator and is counted on every row as
#: ``n_lineage_unavailable`` with one of these reasons.
PANC_LINEAGE_UNAVAILABLE_REASON = (
    "cell_class is null: no class was written for this cell; its lineage is "
    "unavailable, not non-immune, and it is outside n_cells")
PANC_LINEAGE_MISSING_STRING_REASON = (
    "cell_class is a string-shaped missing value ('', 'nan', 'None', 'NA'): no class "
    "was written for this cell; its lineage is unavailable, not non-immune, and it is "
    "outside n_cells")
PANC_LINEAGE_OUTSIDE_VOCABULARY_REASON = (
    "cell_class is outside the declared vocabulary (build_cells.DISPLAY): it can be read "
    "neither as immune nor as any non-immune class; its lineage is unavailable and it is "
    "outside n_cells")
#: The ``cell_class`` vocabulary is the one declared table, ``DISPLAY`` in
#: core/multiplex5000/scripts/typing/build_cells.py, whose KEYS are what ``typing/cells``
#: carries ('Alpha', 'Blood vessel', 'Nerve', ...): build_cells.py:394 maps those keys to
#: the display labels and build_cells.py:396 reports a value the table lacks.  It is
#: derived here, never restated, so the two cannot drift.  ``Immune`` is the immune
#: class; every other declared class is non-immune, except the declared absence below.
PANC_CLASS_VOCABULARY = frozenset(TYPING_CELLS.DISPLAY)
PANC_IMMUNE_CLASS = "Immune"
#: In the table so it can be displayed, and not a lineage: the "no lineage evidence"
#: entry.  A cell carrying it is lineage-unavailable, never non-immune.  No cell
#: carries it today: refine.py writes it into refined_class, beside cell_class, only
#: through the build_cells.py --refine join, which has not run.
PANC_LINEAGE_ABSENT = frozenset({TYPING_CELLS.UNRESOLVED})
PANC_LINEAGE_DECLARED_ABSENCE_REASON = (
    f"cell_class is {TYPING_CELLS.UNRESOLVED!r}: the declared absence of lineage "
    "evidence; its lineage is unavailable, not non-immune, and it is outside n_cells")
#: A class written as one of these strings (after strip and lower-case) is a missing
#: value wearing a string, not a class.
PANC_MISSING_CLASS_STRINGS = frozenset({"", "nan", "none", "na", "null", "<na>"})
#: Why a cell has no lineage, in the order the reasons are reported.
PANC_LINEAGE_KINDS = {
    "null": PANC_LINEAGE_UNAVAILABLE_REASON,
    "missing_string": PANC_LINEAGE_MISSING_STRING_REASON,
    "declared_absence": PANC_LINEAGE_DECLARED_ABSENCE_REASON,
    "outside_vocabulary": PANC_LINEAGE_OUTSIDE_VOCABULARY_REASON,
}
PLN_LINEAGE_UNAVAILABLE_REASON = (
    "is_immune is null: the cell's cluster carries no immune call (the producer's "
    "n_uncalled); its lineage is unavailable, not non-immune, and it is outside n_cells")
#: The row-level eligibility columns the pancreas path filters on.  The pLN path
#: filters on them the same way where its partition carries them.
ELIGIBILITY_COLUMNS = ("qc_analysis_eligible", "excluded_reason")
PANCREAS_COMPARTMENTS = ["core", "peri", "exocrine"]
COMMON_LINEAGES = ["T_NK", "B_Plasma", "Myeloid", "Granulocyte"]

#: PhenoCycler pancreas fine immune calls -> common lineage.  NK is T_NK because
#: the crosswalk collapses T and NK on the Xenium side under one common lineage;
#: the fine call is kept in the subtype column, not thrown away.
PANC_IMMUNE_TO_COMMON = {
    "Imm-T-cell": "T_NK",
    "Imm-NK-cell": "T_NK",
    "Imm-B-cell": "B_Plasma",
    "Imm-plasma-cell": "B_Plasma",
    "Imm-macrophage": "Myeloid",
    "Imm-dendritic-cell": "Myeloid",
    "Imm-neutrophil": "Granulocyte",
}

#: PhenoCycler pLN cluster-named immune classes -> common lineage.
PLN_CLASS_TO_COMMON = {
    "T_cell": "T_NK",
    "NK_ILC": "T_NK",
    "B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}

#: Xenium pLN release labels -> common lineage.
XENIUM_LABEL_TO_COMMON = {
    "Immune_T": "T_NK",
    "Immune_NK_ILC": "T_NK",
    "Immune_B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}


class CrossTissueDistance(RuntimeError):
    """Raised when something asks for a distance that spans a tissue boundary."""


def refuse_cross_tissue(tissue: str, quantity: str) -> None:
    """Structural refusal, not a convention.

    A lymph node has no islet.  A signed distance to one is undefined there, not
    zero, so the refusal is an exception rather than a NaN that a later mean
    could swallow.
    """
    if tissue != "pancreas":
        raise CrossTissueDistance(
            f"{quantity} is a spatial quantity and {tissue!r} is not the pancreas "
            "section it was measured on. The pLN is a different block: it has no "
            "islet, so a signed distance to one is undefined, not zero. Join on "
            "donor_id and donor_status instead."
        )


def read_parquet_retrying(path: Path, columns=None, attempts: int = 5):
    """A concurrent workflow writes the pLN tree while this reads it.

    A partial read is a wrong number, not a small one, so retry and then withhold
    -- never record what was on disk halfway through a write.
    """
    import time

    last = None
    for attempt in range(attempts):
        try:
            return pd.read_parquet(path, columns=columns)
        except Exception as exc:  # pragma: no cover - timing dependent
            last = exc
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"could not read {path} after {attempts} attempts: {last}")


def eligible_donors() -> list[str]:
    donors = sorted(p.name.split("=", 1)[1] for p in CELLS.glob("donor_id=*"))
    return [d for d in donors if d not in DONOR_EXCLUSIONS]


def load_donor_meta() -> pd.DataFrame:
    meta = pd.read_csv(DONOR_META, dtype={"case_id": str})
    meta = meta.rename(columns={"case_id": "donor_id"})
    acq = pd.read_csv(ACQ_RUNS, dtype={"donor_id": str})
    acq = acq[acq["excluded"].astype(str).str.lower() != "true"]
    # One donor can carry several scans; the compartment quantities are per
    # donor, so the batch covariate is the donor's scan when it is unambiguous
    # and unrecorded when it is not.  Unrecorded is not a default value.
    agg = acq.groupby("donor_id").agg(
        acq_run=("acq_run", lambda s: s.iloc[0] if s.nunique() == 1 else pd.NA),
        stain_group=("stain_group", lambda s: s.iloc[0] if s.nunique() == 1 else pd.NA),
    )
    return meta.merge(agg, on="donor_id", how="left")


def load_crosswalk_immune() -> tuple[set[str], dict[str, str]]:
    """Which common lineages the crosswalk records as immune, and in which tissues.

    Read through the one strict reader (scripts/senior/crosswalk.py): ``is_immune`` must say
    exactly ``yes`` -- it was matched case-insensitively here until 2026-09-24 -- and a blank
    scope is ``""``, not the string ``"nan"`` this reader used to write into the manifest.
    """
    lin = XW.lineage_rows(CROSSWALK)
    immune = lin[lin["is_immune"].eq("yes")]
    return set(immune["common"]), dict(zip(immune["common"], immune["immune_tissue_scope"]))


# ---------------------------------------------------------------------------
# pancreas


def compartment_from_signed_distance(dist: pd.Series) -> pd.Series:
    """core / peri / exocrine, cut from the signed distance alone.

    The drawn expansion annotation is not consulted: one donor's is 23.1 um
    while every other is 18.8-19.3 um, and all of them are named _exp20um.
    A cell with no signed distance has an unavailable compartment, never a
    default one.
    """
    out = pd.Series(pd.NA, index=dist.index, dtype="object")
    have = dist.notna()
    out[have & (dist <= 0)] = "core"
    out[have & (dist > 0) & (dist <= EXPANSION_UM)] = "peri"
    out[have & (dist > EXPANSION_UM)] = "exocrine"
    return out


def pancreas_donor(donor: str, min_islet_cells: int) -> dict:
    path = CELLS / f"donor_id={donor}" / "data.parquet"
    cols = [
        "object_id", "cell_class", "immune_definition", "excluded_reason",
        "qc_analysis_eligible", "signed_dist_islet_um", "islet_num",
        "cell_region_20um", "acq_run", "stain_group",
    ]
    df = read_parquet_retrying(path, columns=cols)
    refuse_cross_tissue("pancreas", "signed_dist_islet_um")

    n_rows = len(df)
    # A cell with no nucleus never reached this tree; what is here and still
    # ineligible is excluded by name, and the name is kept.
    keep = df["qc_analysis_eligible"].fillna(False)
    keep &= df["excluded_reason"].fillna("").str.strip() == ""
    dropped = df.loc[~keep, "excluded_reason"].fillna("ineligible").value_counts().to_dict()
    df = df.loc[keep].copy()

    df["compartment"] = compartment_from_signed_distance(df["signed_dist_islet_um"])
    no_distance = int(df["compartment"].isna().sum())

    # A cell whose class was never written has no lineage.  ``== "Immune"`` on a
    # null is False, which would count it non-immune and keep it in the
    # denominator; so is it on a class written as '' / 'nan' / 'None' / 'NA', and
    # on any value the declared vocabulary does not contain.  Each leaves every
    # calculation here and is counted by name and by kind, per compartment, on
    # the rows it would have been in.
    kind = _lineage_unavailable_kind(df["cell_class"])
    lineage_unavailable = kind != ""
    null_lineage = df.loc[lineage_unavailable, ["compartment", "islet_num", "cell_class"]]
    null_lineage = null_lineage.assign(kind=kind[lineage_unavailable])
    df = df.loc[~lineage_unavailable].copy()

    def _lineage_kwargs(comp, islets=None) -> dict:
        sel = null_lineage["compartment"] == comp
        if islets is not None:
            sel &= null_lineage["islet_num"].isin(islets)
        return {"n_lineage_unavailable": int(sel.sum()),
                "lineage_unavailable_reason": _lineage_reason(null_lineage.loc[sel])}

    df["is_immune"] = df["cell_class"] == PANC_IMMUNE_CLASS
    fine = df["immune_definition"].fillna("").str.strip()
    df["fine"] = fine.where(df["is_immune"], "")
    df["common"] = df["fine"].map(PANC_IMMUNE_TO_COMMON)
    unmapped = sorted(set(df.loc[df["is_immune"] & df["common"].isna(), "fine"]) - {""})

    rows = []
    for comp in PANCREAS_COMPARTMENTS:
        sub = df[df["compartment"] == comp]
        rows.extend(_composition_rows(donor, "pancreas", "phenocycler", comp, sub, "all_islets",
                                      **_lineage_kwargs(comp)))

    # Stratum 2: only islets with at least min_islet_cells core cells.  An
    # islet-level statement below that is a statement about classifier
    # fragments, so it is computed separately rather than pooled in.
    core = df[df["compartment"] == "core"]
    sizes = core["islet_num"].dropna().value_counts()
    big = set(sizes[sizes >= min_islet_cells].index)
    core_big = core[core["islet_num"].isin(big)]
    strat_rows = _composition_rows(
        donor, "pancreas", "phenocycler", "core", core_big,
        f"islets_ge_{min_islet_cells}_core_cells", **_lineage_kwargs("core", big),
    )

    # The peri-islet expansion inherits the same stratification: a 20 um
    # expansion drawn around a one-cell fragment is not a peri-islet expansion.
    peri = df[df["compartment"] == "peri"]
    peri_big = peri[peri["islet_num"].isin(big)]
    strat_rows.extend(_composition_rows(
        donor, "pancreas", "phenocycler", "peri", peri_big,
        f"islets_ge_{min_islet_cells}_core_cells", **_lineage_kwargs("peri", big),
    ))

    n_islets = int(len(sizes))
    n_big = int(len(big))
    frag = float((sizes <= 1).sum() / len(sizes)) if len(sizes) else float("nan")
    median = float(sizes.median()) if n_islets else float("nan")
    # The flag is the measured condition, not a judgement: a donor whose median
    # islet is smaller than the threshold at which an islet-level statement is
    # defensible has islets that are mostly classifier fragments.
    flag = "fragmented_islets" if (n_islets and median < min_islet_cells) else ""

    audit = _membership_audit(donor, df)

    return {
        "rows": rows,
        "strat_rows": strat_rows,
        "islets": {
            "donor_id": donor,
            "n_islets": n_islets,
            f"n_islets_ge_{min_islet_cells}": n_big,
            "median_core_cells_per_islet": median,
            "fraction_islets_one_cell": frag,
            "flag": flag,
        },
        "diagnostics": {
            "donor_id": donor,
            "n_rows": n_rows,
            "n_eligible": int(len(df)),
            "dropped": dropped,
            "cells_without_signed_distance": no_distance,
            "n_lineage_unavailable": int(len(null_lineage)),
            "lineage_unavailable_reason": _lineage_reason(null_lineage),
            "lineage_unavailable_by_kind": {
                k: int((null_lineage["kind"] == k).sum()) for k in PANC_LINEAGE_KINDS},
            "cell_class_values_outside_vocabulary": _outside_vocabulary_values(null_lineage),
            "cell_class_vocabulary_from":
                "core/multiplex5000/scripts/typing/build_cells.py: DISPLAY (keys)",
            "immune_calls_not_in_crosswalk": unmapped,
            "acq_run": _one(df["acq_run"]),
            "stain_group": _one(df["stain_group"]),
            "membership_audit": audit,
        },
    }


def _one(series: pd.Series):
    vals = series.dropna().unique()
    return str(vals[0]) if len(vals) == 1 else None


def _lineage_unavailable_kind(cell_class: pd.Series) -> pd.Series:
    """'' where the class is a declared lineage; otherwise why it is not.

    ``null``: never written.  ``missing_string``: written as one of
    ``PANC_MISSING_CLASS_STRINGS``.  ``declared_absence``: the table's own
    "no lineage evidence" entry.  ``outside_vocabulary``: a value
    ``build_cells.DISPLAY`` does not declare.  None of the four is a non-immune
    cell, and a value is compared to the vocabulary exactly as written.
    """
    text = cell_class.astype("string")
    is_null = text.isna().fillna(False).astype(bool).to_numpy()
    missing = (text.str.strip().str.lower().isin(PANC_MISSING_CLASS_STRINGS)
               .fillna(False).astype(bool).to_numpy())
    absent = text.isin(PANC_LINEAGE_ABSENT).fillna(False).astype(bool).to_numpy()
    declared = text.isin(PANC_CLASS_VOCABULARY).fillna(False).astype(bool).to_numpy()
    kind = np.select([is_null, missing, absent, ~declared],
                     ["null", "missing_string", "declared_absence", "outside_vocabulary"],
                     default="")
    return pd.Series(kind, index=cell_class.index)


def _outside_vocabulary_values(null_lineage: pd.DataFrame, limit: int = 10) -> list[str]:
    hit = null_lineage.loc[null_lineage["kind"] == "outside_vocabulary", "cell_class"]
    values = sorted(str(v) for v in hit.unique())
    if len(values) > limit:
        values = values[:limit] + [f"... and {len(values) - limit} more"]
    return values


def _lineage_reason(sel: pd.DataFrame) -> str:
    """One reason per kind present, each with its count; offending values named."""
    parts = []
    for kind, reason in PANC_LINEAGE_KINDS.items():
        n = int((sel["kind"] == kind).sum()) if len(sel) else 0
        if not n:
            continue
        if kind == "outside_vocabulary":
            reason = f"{reason}; values: {_outside_vocabulary_values(sel)}"
        parts.append(f"{reason} [n={n}]")
    return " | ".join(parts)


def _boolean_valued(series: pd.Series) -> bool:
    """True when every non-null value is a bool.

    A numpy bool, a nullable ``boolean``, or an object column holding only
    True / False / null -- which is how an arrow bool with nulls reads back when
    the file carries no pandas dtype metadata.  A string or integer column is
    not boolean-valued however it spells True.
    """
    if pd.api.types.is_bool_dtype(series):
        return True
    if series.dtype != object:
        return False
    values = series.dropna()
    return bool(values.map(lambda v: isinstance(v, (bool, np.bool_))).all())


def _membership_audit(donor: str, df: pd.DataFrame) -> dict:
    """Read the geometric outline membership back, where it exists.

    This is an audit of the distance cut, not an override of it, and it is
    reported rather than silently reconciled.  QuPath's Parent column is not
    read here or anywhere in this module.
    """
    path = MEMBERSHIP / f"donor_id={donor}" / "data.parquet"
    if not path.exists():
        return {"state": "unavailable",
                "reason": "no islet_membership partition for this donor"}
    m = read_parquet_retrying(path, columns=["object_id", "in_islet"])
    j = df[["object_id", "compartment"]].merge(m, on="object_id", how="left")
    covered = j["in_islet"].notna()
    core_cut = j["compartment"] == "core"
    return {
        "state": "measured",
        "n_covered": int(covered.sum()),
        "n_core_by_distance": int(core_cut.sum()),
        "n_core_by_outline": int((j["in_islet"] == True).sum()),  # noqa: E712
        "n_disagree": int((covered & (core_cut != (j["in_islet"] == True))).sum()),  # noqa: E712
    }


def _composition_rows(donor, tissue, modality, compartment, sub, stratum, *,
                      n_lineage_unavailable: int = 0,
                      lineage_unavailable_reason: str = "") -> list[dict]:
    """One denominator per compartment, repeated across subtype rows.

    ``n_cells`` is the eligible denominator for the compartment.  Do not sum it
    across subtype rows -- the subtype counts sum to the immune count, and the
    immune count plus the non-immune count is the denominator.  That internal
    coherence is checked by the test suite.

    ``is_immune`` must be boolean-valued on every row of ``sub``.  A cell whose
    lineage is null has already left the frame in the caller and arrives here as
    ``n_lineage_unavailable`` with its reason; that count is repeated on every
    row of the compartment, beside ``n_cells`` and outside it.  A null that
    reaches this function is refused rather than summed as False.
    """
    n = int(len(sub))
    if n and sub["is_immune"].isna().any():
        raise ValueError(
            f"{donor} {tissue} {compartment}: is_immune is null on "
            f"{int(sub['is_immune'].isna().sum())} rows. A null lineage is unavailable, "
            "not non-immune: remove it and pass it as n_lineage_unavailable.")
    n_immune = int(sub["is_immune"].sum()) if n else 0
    rows = [{
        "donor_id": donor, "tissue": tissue, "modality": modality,
        "compartment": compartment, "stratum": stratum,
        "lineage": "Immune_total", "n_cells": n, "n_lineage": n_immune,
        "n_immune": n_immune,
        "fraction_of_compartment": (n_immune / n) if n else float("nan"),
        "fraction_of_immune": 1.0 if n_immune else float("nan"),
        "state": "measured" if n else "unavailable",
        "reason": "" if n else f"no eligible cell in compartment {compartment}",
        "n_lineage_unavailable": int(n_lineage_unavailable),
        "lineage_unavailable_reason": lineage_unavailable_reason,
    }]
    for lineage in COMMON_LINEAGES:
        k = int((sub["common"] == lineage).sum()) if n else 0
        rows.append({
            "donor_id": donor, "tissue": tissue, "modality": modality,
            "compartment": compartment, "stratum": stratum,
            "lineage": lineage, "n_cells": n, "n_lineage": k, "n_immune": n_immune,
            "fraction_of_compartment": (k / n) if n else float("nan"),
            "fraction_of_immune": (k / n_immune) if n_immune else float("nan"),
            "state": "measured" if n else "unavailable",
            "reason": "" if n else f"no eligible cell in compartment {compartment}",
            "n_lineage_unavailable": int(n_lineage_unavailable),
            "lineage_unavailable_reason": lineage_unavailable_reason,
        })
    return rows


def _pln_unavailable_rows(donor, modality: str, state: str, reason: str) -> list[dict]:
    """The pLN rows for a donor that has no number: every lineage, one reason.

    ``n_cells`` is null, never 0, and the reason travels on the row so the card
    can name it.
    """
    return [{
        "donor_id": donor, "tissue": "pancreatic_lymph_node",
        "modality": modality, "compartment": "pln",
        "stratum": "not_applicable_no_islets", "lineage": lineage,
        "n_cells": pd.NA, "n_lineage": pd.NA, "n_immune": pd.NA,
        "fraction_of_compartment": float("nan"),
        "fraction_of_immune": float("nan"),
        "state": state, "reason": reason,
        "n_lineage_unavailable": pd.NA, "lineage_unavailable_reason": "",
    } for lineage in ["Immune_total"] + COMMON_LINEAGES]


# ---------------------------------------------------------------------------
# pLN -- compositional only


def pln_phenocycler(donor: str) -> tuple[list[dict], dict]:
    base = PLN_PHENO / f"donor_id={donor}"
    status = pd.read_csv(PLN_PHENO / "per_donor_status.csv", dtype={"donor_id": str})
    row = status[status["donor_id"] == donor]
    if not len(row):
        return [], {"donor_id": donor, "state": "unavailable",
                    "reason": "donor has no pLN section in the pheno tree"}
    state = str(row.iloc[0]["state"])
    if state != "typed":
        # The measured reason travels with the row so the card can name it.
        # "the export carries no intensity measurement, so there is nothing to
        # cluster and nothing to name" is a statement about the measurement,
        # and it is never rendered as zero immune cells.
        reason = str(row.iloc[0].get("reason", "") or f"pLN state is {state}")
        rows = _pln_unavailable_rows(
            donor, "phenocycler",
            state if state in {"unavailable", "withheld", "failed"} else "unavailable",
            reason)
        return rows, {"donor_id": donor, "state": state, "reason": reason}

    parts = sorted(base.glob("roi=*/data.parquet"))
    if not parts:
        return [], {"donor_id": donor, "state": "unavailable",
                    "reason": "pLN status says typed but no roi partition is on disk"}

    # Eligibility is filtered on the row where the partition carries it, with
    # the pancreas path's filter verbatim.  Today's partitions carry neither
    # column: core/multiplex5000/scripts/typing/build_pln_immune_pheno.py applies
    # ``qc_analysis_eligible`` before it writes and records n_rows / n_eligible
    # in per_donor_counts.csv and per_donor_status.csv beside the tree.  Which
    # of the two happened is written into the diagnostics, never assumed.
    carried = [c for c in ELIGIBILITY_COLUMNS if c in set(pq.read_schema(parts[0]).names)]
    if len(carried) == 1:
        missing = next(c for c in ELIGIBILITY_COLUMNS if c not in carried)
        reason = (f"partition carries {carried[0]} but not {missing}; eligibility "
                  "cannot be evaluated on the row and is not assumed")
        return (_pln_unavailable_rows(donor, "phenocycler", "unavailable", reason),
                {"donor_id": donor, "state": "unavailable", "reason": reason})
    frames = [read_parquet_retrying(part, columns=["is_immune", "immune_class"] + carried)
              for part in parts]
    df = pd.concat(frames, ignore_index=True)
    n_rows = len(df)
    n_rows_on_disk = n_rows
    if carried:
        # A cell with no nucleus never reached this tree; what is here and still
        # ineligible is excluded by name, and the name is kept.
        keep = df["qc_analysis_eligible"].fillna(False)
        keep &= df["excluded_reason"].fillna("").str.strip() == ""
        dropped = df.loc[~keep, "excluded_reason"].fillna("ineligible").value_counts().to_dict()
        df = df.loc[keep].copy()
        n_eligible = int(len(df))
        eligibility_filter = ("row-level: qc_analysis_eligible and excluded_reason, the "
                              "pancreas path's filter")
        dropped_reason = ""
    else:
        recorded_eligible = pd.to_numeric(row.iloc[0].get("n_eligible", pd.NA), errors="coerce")
        if pd.notna(recorded_eligible) and n_rows > int(recorded_eligible):
            reason = (f"partition holds {n_rows} rows but per_donor_status.csv records "
                      f"n_eligible={int(recorded_eligible)}; the partition carries no "
                      "qc_analysis_eligible column, so the excess cannot be filtered here")
            return (_pln_unavailable_rows(donor, "phenocycler", "unavailable", reason),
                    {"donor_id": donor, "state": "unavailable", "reason": reason})
        eligibility_filter = (
            "upstream: build_pln_immune_pheno.py filtered on qc_analysis_eligible before "
            "writing; the partition carries no qc_analysis_eligible or excluded_reason "
            "column, so the row-level filter was not re-evaluated here")
        recorded = _pln_recorded_counts(donor)
        if recorded is None:
            dropped = None
            if pd.notna(recorded_eligible):
                n_eligible = int(recorded_eligible)
                dropped_reason = ("per_donor_counts.csv does not record this donor; the "
                                  "ineligible count is not reported")
            else:
                # Nobody recorded an eligible count: per_donor_status.csv carries
                # nothing this reader can parse as one and per_donor_counts.csv
                # does not record the donor.  The rows on disk are
                # n_rows_on_disk, a count of what is here; they are never written
                # up as n_eligible, which would be a count nobody recorded.
                n_eligible = None
                raw = row.iloc[0].get("n_eligible", pd.NA)
                spelled = "blank" if pd.isna(raw) else repr(raw)
                dropped_reason = (
                    "the recorded eligible count is absent: per_donor_status.csv carries "
                    f"n_eligible {spelled}, which is not a count, and per_donor_counts.csv "
                    "does not record this donor; n_eligible is not reported and the "
                    "ineligible count is not reported")
        else:
            n_rows, n_eligible = recorded
            dropped = {"ineligible_recorded_by_producer": n_rows - n_eligible}
            dropped_reason = ""

    # The producer's recorded n_eligible is read back beside what is on disk.  A
    # partition SMALLER than that count is measured as it stands -- what is on
    # disk is real -- but the note says so, in numbers, so n_eligible can never
    # sit beside a smaller n_cells silently.  (Larger is refused above.)
    eligible_count_disagreement = ""
    if not carried and n_eligible is not None and n_eligible > n_rows_on_disk:
        eligible_count_disagreement = (
            f"partition holds {n_rows_on_disk} rows but the producer recorded "
            f"n_eligible={n_eligible}; the {n_eligible - n_rows_on_disk} eligible cells not "
            "on disk are outside n_cells, and nothing here can say which they are")

    # is_immune must be boolean-valued before anything reads it.  ``astype(bool)``
    # on a string column makes every non-empty value True and on an integer
    # column every non-zero value True, so ['True', 'NA'] and [1, -1] each came
    # out as two immune cells, measured.  A partition this reader cannot read is
    # refused the way the pancreas path refuses one -- read_parquet_retrying
    # raises and nothing is recorded -- rather than reported as a number.
    if not _boolean_valued(df["is_immune"]):
        raise ValueError(
            f"{donor} pancreatic_lymph_node: is_immune is {df['is_immune'].dtype}, not boolean. "
            "Only True, False and null can be read as an immune call; coercing this column "
            "would count every non-empty or non-zero value as immune. The partition is "
            "unreadable and is refused, as the pancreas path refuses an unreadable partition.")
    # A null is_immune is a cell whose cluster carries no immune call.  Its
    # lineage is unavailable, not non-immune: it leaves the denominator and is
    # counted under its own name on every row.
    lineage_unavailable = df["is_immune"].isna()
    n_lineage_unavailable = int(lineage_unavailable.sum())
    df = df.loc[~lineage_unavailable].copy()
    df["is_immune"] = df["is_immune"].astype(bool)
    cls = df["immune_class"].astype(str).str.strip()
    df["common"] = cls.map(PLN_CLASS_TO_COMMON)
    unmapped = sorted(set(cls[df["is_immune"] & df["common"].isna()]) - {"", "nan"})
    lineage_reason = PLN_LINEAGE_UNAVAILABLE_REASON if n_lineage_unavailable else ""
    rows = _composition_rows(donor, "pancreatic_lymph_node", "phenocycler",
                             "pln", df, "all_islets",
                             n_lineage_unavailable=n_lineage_unavailable,
                             lineage_unavailable_reason=lineage_reason)
    for r in rows:
        r["stratum"] = "not_applicable_no_islets"
        if not len(df) and n_lineage_unavailable:
            r["reason"] = (f"every eligible cell's lineage is unavailable "
                           f"({n_lineage_unavailable} cells): {lineage_reason}")
    return rows, {"donor_id": donor, "state": "measured" if len(df) else "unavailable",
                  "reason": "" if len(df) else rows[0]["reason"],
                  "n_rows": int(n_rows),
                  "n_eligible": None if n_eligible is None else int(n_eligible),
                  "n_rows_on_disk": int(n_rows_on_disk),
                  "eligible_count_disagreement": eligible_count_disagreement,
                  "dropped": dropped, "dropped_reason": dropped_reason,
                  "eligibility_filter": eligibility_filter,
                  "n_cells": int(len(df)),
                  "n_lineage_unavailable": n_lineage_unavailable,
                  "lineage_unavailable_reason": lineage_reason,
                  "classes_not_in_crosswalk": unmapped}


def _pln_recorded_counts(donor: str) -> tuple[int, int] | None:
    """(n_rows, n_eligible) the producer recorded for this donor, or None.

    ``per_donor_counts.csv`` is written by build_pln_immune_pheno.py in the same
    pass as ``per_donor_status.csv``; a tree without it, or a row without both
    counts, records nothing and nothing is invented for it.
    """
    path = PLN_PHENO / "per_donor_counts.csv"
    if not path.exists():
        return None
    counts = pd.read_csv(path, dtype={"donor_id": str})
    hit = counts[counts["donor_id"] == donor]
    if not len(hit) or not {"n_rows", "n_eligible"} <= set(counts.columns):
        return None
    n_rows, n_eligible = hit.iloc[0]["n_rows"], hit.iloc[0]["n_eligible"]
    if pd.isna(n_rows) or pd.isna(n_eligible):
        return None
    return int(n_rows), int(n_eligible)


def _xenium_count(value, column: str) -> tuple[int | None, str]:
    """The boundary for one per_donor_summary.csv count.

    ``(count, "")`` for a non-negative integer -- an integral float passes,
    because one null anywhere in a CSV column makes pandas read the whole column
    as float, and a numeric string passes, because a CSV token has no type --
    ``(None, "")`` for a null, which the caller reports as not reported, and
    ``(None, reason)`` for a value that cannot be a count: a boolean, a
    non-number, a non-finite number, a fraction, or a negative.
    """
    if value is None or pd.isna(value):
        return None, ""
    if isinstance(value, (bool, np.bool_)):
        return None, f"{column}={value!r} is a boolean, not a count"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None, f"{column}={value!r} is not a number"
    if not np.isfinite(number):
        return None, f"{column}={value!r} is not finite"
    if number != int(number):
        return None, f"{column}={value!r} is not an integer"
    if number < 0:
        return None, f"{column}={value!r} is negative"
    return int(number), ""


def _xenium_counts(r, columns, col_map: dict[str, list[list[str]]]) -> tuple[dict, str]:
    """Every count on one per_donor_summary.csv row, through ONE boundary.

    Reads ``cells``, ``immune``, ``no_nucleus_excluded`` and every lineage
    column the release writes.  Each must be a non-negative integer or null,
    and the decomposition must hold: ``no_nucleus_excluded <= cells``, and
    ``immune`` and each lineage ``<= cells - no_nucleus_excluded``.  Returns
    ``(counts, "")`` with None for every absent or null count, or ``(counts,
    reason)`` naming each column and value that fails; a row with a reason is
    unavailable and nothing downstream reads its numbers.  Before this boundary
    ``immune=-1`` measured a fraction of -0.125, ``no_nucleus_excluded=12`` over
    ``cells=10`` measured ``n_cells=-2``, ``immune=9`` over 8 nucleated cells
    measured 1.125, and ``immune=4.7`` truncated to 4.
    """
    lineage_columns = [c for groups in col_map.values() for group in groups for c in group]
    counts: dict[str, int | None] = {}
    failures: list[str] = []
    for column in ["cells", "immune", "no_nucleus_excluded"] + lineage_columns:
        if column not in columns:
            counts[column] = None
            continue
        counts[column], failure = _xenium_count(r[column], column)
        if failure:
            failures.append(failure)
    cells, nne, immune = counts["cells"], counts["no_nucleus_excluded"], counts["immune"]
    nucleated = None
    if cells is not None and nne is not None:
        if nne > cells:
            failures.append(f"no_nucleus_excluded={nne} exceeds cells={cells}")
        else:
            nucleated = cells - nne
    if nucleated is not None:
        if immune is not None and immune > nucleated:
            failures.append(f"immune={immune} exceeds nucleated cells={nucleated} "
                            f"(cells={cells} - no_nucleus_excluded={nne})")
        for lineage, groups in col_map.items():
            have = next((g for g in groups if all(c in columns for c in g)), None)
            if have is None or any(counts[c] is None for c in have):
                continue
            k = sum(counts[c] for c in have)
            if k > nucleated:
                spelled = " + ".join(f"{c}={counts[c]}" for c in have)
                failures.append(f"{lineage} ({spelled}) = {k} exceeds nucleated cells={nucleated}")
    if not failures:
        return counts, ""
    return counts, ("per_donor_summary.csv row cannot be a measurement: "
                    + "; ".join(failures))


def pln_xenium() -> tuple[list[dict], dict]:
    """The RNA view of the same compartment, as its own modality column.

    It is never merged with the protein pLN numbers and never tuned against
    them.  Agreement is informative and disagreement is informative.
    Cell ids do not cross: a QuPath UUID and a 10x cell_id share no namespace,
    so nothing here joins on a cell.
    """
    release = xenium_source("pln")
    if not release.available:
        return [], {"state": "unavailable", "reason": str(release),
                    "roots_searched": [str(r) for r in release.roots],
                    "candidates_rejected": list(release.rejected)}
    path = release.sidecar("per_donor_summary.csv")
    if path is None:
        return [], {"state": "unavailable",
                    "reason": f"{release.release_id} writes no per_donor_summary.csv at "
                              f"{release.root}",
                    "release": release.describe()}
    s = read_csv_retrying(path)
    s["donor_id"] = s["donor_id"].astype(str)
    rows = []
    # Each lineage names ALTERNATIVE column groups, newest spelling first, and the
    # first group whose columns are ALL present wins.  Two releases spell the same
    # quantity differently -- the superseding one writes `T_NK_compartment` and
    # `Myeloid`, its predecessor wrote `s2__Immune_T` + `s2__Immune_NK_ILC` and
    # `s2__Myeloid` -- and reading a spelling that is not there would have printed
    # every lineage as unavailable rather than as the number the release carries.
    # Groups are all-or-nothing so a half-present group can never half-count.
    col_map = {
        "T_NK": [["T_NK_compartment"], ["s2__Immune_T", "s2__Immune_NK_ILC"]],
        "B_Plasma": [["Immune_B_Plasma"], ["s2__Immune_B_Plasma"]],
        "Myeloid": [["Myeloid"], ["s2__Myeloid"]],
        "Granulocyte": [["Granulocyte"], ["s2__Granulocyte"]],
    }
    donors_without_nne = []
    donors_unavailable: dict[str, str] = {}

    for _, r in s.iterrows():
        # Every count on the row passes ONE boundary at the read: a non-negative
        # integer or null, and a decomposition that holds.  A row that fails is
        # unavailable with the column and values named, and none of the
        # arithmetic below is hardened against it site by site.
        counts, refusal = _xenium_counts(r, s.columns, col_map)
        if refusal:
            donors_unavailable[str(r["donor_id"])] = refusal
            rows.extend(_pln_unavailable_rows(r["donor_id"], "xenium", "unavailable", refusal))
            continue
        # The denominator excludes cells with no nucleus, which are not cells we
        # can be sure of, and excludes nothing else by value.  A release that
        # does not report that exclusion has not reported the denominator: an
        # absent count is not zero exclusions, and an unreduced denominator is a
        # wrong number, not a slightly large one.
        nne = counts["no_nucleus_excluded"]
        if nne is None:
            donors_without_nne.append(str(r["donor_id"]))
            donors_unavailable[str(r["donor_id"])] = "no_nucleus_excluded not reported"
            rows.extend(_pln_unavailable_rows(r["donor_id"], "xenium", "unavailable",
                                              "no_nucleus_excluded not reported"))
            continue
        # The numerator and the denominator are read the same way: a count the
        # summary does not carry, or carries null, is not a count.  int() on a
        # null raised mid-build and an absent column raised KeyError; neither is
        # a reported absence, and this row is one.
        cells_reported = counts["cells"]
        immune_reported = counts["immune"]
        if cells_reported is None or immune_reported is None:
            missing = [c for c, v in (("cells", cells_reported), ("immune", immune_reported))
                       if v is None]
            reason = f"{' and '.join(missing)} not reported"
            donors_unavailable[str(r["donor_id"])] = reason
            rows.extend(_pln_unavailable_rows(r["donor_id"], "xenium", "unavailable", reason))
            continue
        n = cells_reported - nne
        n_immune = immune_reported
        rows.append({
            "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
            "modality": "xenium", "compartment": "pln",
            "stratum": "not_applicable_no_islets", "lineage": "Immune_total",
            "n_cells": n, "n_lineage": n_immune, "n_immune": n_immune,
            "fraction_of_compartment": n_immune / n if n else float("nan"),
            "fraction_of_immune": 1.0 if n_immune else float("nan"),
            "state": "measured", "reason": "",
            "n_lineage_unavailable": pd.NA, "lineage_unavailable_reason": "",
        })
        for lineage, groups in col_map.items():
            have = next((g for g in groups if all(c in s.columns for c in g)), None)
            if have is None:
                rows.append({
                    "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
                    "modality": "xenium", "compartment": "pln",
                    "stratum": "not_applicable_no_islets", "lineage": lineage,
                    "n_cells": n, "n_lineage": pd.NA, "n_immune": n_immune,
                    "fraction_of_compartment": float("nan"),
                    "fraction_of_immune": float("nan"),
                    "state": "unavailable",
                    "reason": f"{lineage} has no column in {release.release_id}; that "
                              f"release writes {sorted(s.columns)}",
                    "n_lineage_unavailable": pd.NA, "lineage_unavailable_reason": "",
                })
                continue
            unreported = [c for c in have if counts[c] is None]
            if unreported:
                # The column is there and the value is null: this lineage is not
                # reported for this donor.  ``int(r[c] or 0)`` raised on a NaN
                # here, which is neither a count nor a reported absence.
                rows.append({
                    "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
                    "modality": "xenium", "compartment": "pln",
                    "stratum": "not_applicable_no_islets", "lineage": lineage,
                    "n_cells": n, "n_lineage": pd.NA, "n_immune": n_immune,
                    "fraction_of_compartment": float("nan"),
                    "fraction_of_immune": float("nan"),
                    "state": "unavailable",
                    "reason": f"{lineage} not reported: {', '.join(unreported)} is null in "
                              f"{release.release_id}",
                    "n_lineage_unavailable": pd.NA, "lineage_unavailable_reason": "",
                })
                continue
            k = sum(counts[c] for c in have)
            rows.append({
                "donor_id": r["donor_id"], "tissue": "pancreatic_lymph_node",
                "modality": "xenium", "compartment": "pln",
                "stratum": "not_applicable_no_islets", "lineage": lineage,
                "n_cells": n, "n_lineage": k, "n_immune": n_immune,
                "fraction_of_compartment": k / n if n else float("nan"),
                "fraction_of_immune": k / n_immune if n_immune else float("nan"),
                "state": "measured", "reason": "",
                "n_lineage_unavailable": pd.NA, "lineage_unavailable_reason": "",
            })
    return rows, {"state": "measured", "n_donors": int(s["donor_id"].nunique()),
                  "donors_without_no_nucleus_excluded": donors_without_nne,
                  "donors_unavailable": donors_unavailable,
                  "release": release.describe(),
                  "per_donor_summary": str(path)}


def read_csv_retrying(path: Path, attempts: int = 5) -> pd.DataFrame:
    import time
    last = None
    for attempt in range(attempts):
        try:
            return pd.read_csv(path)
        except Exception as exc:  # pragma: no cover
            last = exc
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"could not read {path}: {last}")


# ---------------------------------------------------------------------------
# batch


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """Between-group share of the total sum of squares. NaN when undefined."""
    ok = ~pd.isna(values) & ~pd.isna(groups)
    v, g = np.asarray(values, float)[ok], np.asarray(groups)[ok]
    if v.size < 3 or len(set(g.tolist())) < 2:
        return float("nan")
    grand = v.mean()
    sst = float(((v - grand) ** 2).sum())
    if sst <= 0:
        return float("nan")
    ssb = 0.0
    for level in set(g.tolist()):
        sel = v[g == level]
        ssb += sel.size * (sel.mean() - grand) ** 2
    return float(ssb / sst)


def eta_with_null(values, groups, n_perm=2000, seed=0) -> dict:
    """Eta-squared against a permutation null.

    With 14-20 donors over 9-11 acquisition levels, a large eta-squared is the
    *expected* value under no effect at all -- the null median says so.  So the
    observed statistic is reported next to that null and its permutation p, and
    neither is allowed to gate anything.  A confounded quantity stays in.
    """
    obs = eta_squared(values, groups)
    if not np.isfinite(obs):
        return {"eta_squared": float("nan"), "null_median": float("nan"),
                "p_permutation": float("nan"), "n_perm": 0}
    rng = np.random.default_rng(seed)
    ok = ~pd.isna(values) & ~pd.isna(groups)
    v, g = np.asarray(values, float)[ok], np.asarray(groups)[ok]
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = eta_squared(rng.permutation(v), g)
    return {
        "eta_squared": obs,
        "null_median": float(np.nanmedian(null)),
        "p_permutation": float((np.nansum(null >= obs) + 1) / (n_perm + 1)),
        "n_perm": n_perm,
        "n_levels": int(len(set(g.tolist()))),
        "n_donors_used": int(v.size),
    }


def batch_table(comp: pd.DataFrame, meta: pd.DataFrame, n_perm: int) -> pd.DataFrame:
    m = meta.set_index("donor_id")
    out = []
    keys = comp[["tissue", "modality", "compartment", "stratum", "lineage"]].drop_duplicates()
    for _, k in keys.iterrows():
        sel = comp
        for col in ["tissue", "modality", "compartment", "stratum", "lineage"]:
            sel = sel[sel[col] == k[col]]
        sel = sel[sel["state"] == "measured"]
        if len(sel) < 3:
            continue
        vals = sel["fraction_of_compartment"].to_numpy(float)
        donors = sel["donor_id"].tolist()
        for grouping in ["stain_group", "acq_run", "donor_status"]:
            g = np.array([m[grouping].get(d, pd.NA) if d in m.index else pd.NA
                          for d in donors], dtype=object)
            stat = eta_with_null(vals, g, n_perm=n_perm)
            out.append({
                "quantity": f"{k['compartment']}_{k['lineage']}_fraction",
                "tissue": k["tissue"], "modality": k["modality"],
                "compartment": k["compartment"], "stratum": k["stratum"],
                "lineage": k["lineage"], "grouping": grouping,
                "n_donors_unrecorded": int(sum(pd.isna(x) for x in g)),
                **stat,
            })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--min-islet-cells", type=int, default=10,
                    help="core cells an islet needs before an islet-level number "
                         "is defensible; the measured threshold is about 10")
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--domain-root", default=str(REPO / "data" / "development" / "typing" / "app"),
                    help="the app domain whose cells these are. It is recorded in the "
                         "manifest so the card can refuse an artifact built for the "
                         "other domain: the two datasets share five donors and no "
                         "vocabulary, and a donor id is not a licence to read.")
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = load_donor_meta()
    immune_lineages, immune_scope = load_crosswalk_immune()
    donors = args.donors or eligible_donors()

    rows, strat, islets, diags, pln_notes = [], [], [], [], []
    for donor in donors:
        res = pancreas_donor(donor, args.min_islet_cells)
        rows.extend(res["rows"])
        strat.extend(res["strat_rows"])
        islets.append(res["islets"])
        diags.append(res["diagnostics"])
        pln_rows, note = pln_phenocycler(donor)
        rows.extend(pln_rows)
        pln_notes.append(note)

    xen_rows, xen_note = pln_xenium()
    cohort = set(donors)
    xen_in = [r for r in xen_rows if r["donor_id"] in cohort]
    xen_out = sorted({r["donor_id"] for r in xen_rows} - cohort)
    rows.extend(xen_in)

    comp = pd.DataFrame(rows)
    comp = comp.merge(meta[["donor_id", "donor_status", "stain_group", "acq_run"]],
                      on="donor_id", how="left")
    comp["donor_status"] = pd.Categorical(comp["donor_status"], STATUS_ORDER, ordered=True)
    comp = comp.sort_values(["tissue", "modality", "compartment", "stratum",
                             "donor_status", "donor_id", "lineage"])
    comp.to_csv(out / "compartment_composition.csv", index=False)

    st = pd.DataFrame(strat)
    st = st.merge(meta[["donor_id", "donor_status"]], on="donor_id", how="left")
    st.to_csv(out / "islet_size_strata.csv", index=False)
    pd.DataFrame(islets).to_csv(out / "islet_size_summary.csv", index=False)

    eta = batch_table(pd.concat([comp, st.assign(tissue="pancreas",
                                                 modality="phenocycler")],
                                ignore_index=True),
                      meta, args.n_perm)
    eta.to_csv(out / "batch_eta_squared.csv", index=False)

    manifest = {
        "generated_by": "scripts/senior/build_compartment_gradient.py",
        "source_domain_root": str(Path(args.domain_root).resolve()),
        "compartments": {
            "core": "signed_dist_islet_um <= 0",
            "peri": f"0 < signed_dist_islet_um <= {EXPANSION_UM}",
            "exocrine": f"signed_dist_islet_um > {EXPANSION_UM}",
            "pln": "compositional only; no distance exists in a lymph node",
        },
        "expansion_from": "signed distance, never the drawn annotation (one donor's "
                          "is 23.1 um against 18.8-19.3 for every other, and all are "
                          "named _exp20um)",
        "parent_column": "never read; membership is geometric",
        "cross_tissue": "pancreas and pLN join on donor_id and donor_status ONLY. "
                        "No distance, adjacency, neighbourhood or KDE crosses the "
                        "boundary. refuse_cross_tissue() makes that structural.",
        "pending_maintainer_classification": {
            "column": LYMPHOID_ATTACHMENT_COLUMN,
            "note": "lymphoid structures anatomically attached to a pancreas section "
                    "will license an islet-to-lymphoid distance. It arrives as data "
                    "and is not inferred here; until it lands every pancreas/pLN pair "
                    "is separate.",
        },
        "stratifications": {
            "min_islet_cells": args.min_islet_cells,
            "why": "the endocrine shape signal is real only above about 10 cells per "
                   "islet; 6539's islets are mostly one-cell classifier fragments",
        },
        "donor_exclusions_from": "phenocycler.cohort.DONOR_EXCLUSIONS",
        "donors_excluded": sorted(DONOR_EXCLUSIONS),
        "donors_built": donors,
        "immune_lineages_from_crosswalk": sorted(immune_lineages),
        "immune_tissue_scope": immune_scope,
        "pln_phenocycler": pln_notes,
        "pln_xenium": xen_note,
        "pln_xenium_donors_outside_cohort": xen_out,
        "pancreas_diagnostics": diags,
        "batch": "eta-squared on stain_group and acq_run, each against a "
                 "permutation null; coarse batch is aliased with disease status so "
                 "it is not corrected on. Nothing is gated on these numbers.",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
