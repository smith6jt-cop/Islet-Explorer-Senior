#!/usr/bin/env python
"""Protein (PhenoCycler) vs RNA (Xenium) lineage-composition concordance.

WHAT THIS MEASURES
------------------
Per donor, per lineage, the share of that lineage within the *comparable* lineage
set, computed independently in each modality, and the difference between them.

These are DIFFERENT SECTIONS OF DIFFERENT BLOCKS.  They are not two views of one
cell, and no cell identifier crosses: a PhenoCycler ``object_id`` is a QuPath UUID
and a Xenium ``cell_id`` is 10x segmentation, unique only within a bundle.  The
join is on ``donor_id`` and nothing else.  Composition against composition is the
only comparison these data license.

Nothing here is tuned to make the modalities agree.  Agreement is informative and
disagreement is informative; the deliverable is which lineages may carry a
cross-modality claim and which may not.

THE DISAGREEMENT THRESHOLD IS FIXED IN THIS FILE BEFORE ANY NUMBER IS COMPUTED.
See ``THRESHOLD_*`` below and the rationale beside them.

WHAT IT READS
-------------
* protein pancreas -- ``data/development/typing/cells``.  The lineage comes from
  ``cell_class`` and, for the immune split, from ``immune_definition``, which the
  crosswalk records as ``pheno_resolution=subclass``: recovered from the
  per-marker gates, not from the broad label.  ``typing/immune`` is NOT read --
  that second pass is a finer subtype that never modifies ``cell_class`` -- and
  neither is ``refined_class``, which no consumer reads at all.
* protein pLN -- ``data/development/typing/pln/pheno``.
* RNA, both tissues -- the current release for that tissue, resolved by
  ``scripts/senior/xenium_release.py`` from the release's own declared metadata.
* ``core/phenocycler/data/integration/vocab_crosswalk.csv`` (generated) for which
  lineages are comparable at all.

Exocrine is ONE label on both sides.  Neither panel carries a marker or a probe
set that separates acinar from ductal, so no such claim exists to compare.

``typing/cells`` here is the 2026-09-15 build: 25,230,109 cells, 20 donors,
eight classes plus ``Excluded``.  No donor has been re-run with the 2026-09-18
islet chain, so nothing read here carries a CD99 question or an ``Unresolved``
class.

EVERY CLASS IS COUNTED, INCLUDING THE ONES THAT CARRY NO SHARE.  A class with no
common lineage used to be dropped by ``value_counts(dropna=True)`` while
``n_cells`` went on counting it: measured 2026-09-18, 316,752 of the 25,230,109
cells -- ``Lymphatic`` 122,884 and ``Excluded`` 193,868 -- were in the total and
in no column.  They now have a named column each, on every donor, and the reader
refuses a donor whose parts do not sum to its ``n_cells``.  Which classes carry a
share is unchanged: see ``PANC_PHENO_CLASS_NOT_COMPARED`` for what each one is
and why it is out, ``Lymphatic`` among them as an open question rather than a
decision taken here.

WHAT IT WRITES
--------------
``--out``, defaulting to ``data/development/integration/concordance``: the
per-donor counts, per-donor-per-lineage differences, lineage summary and batch
eta-squared for each tissue, ``pln_rna_controls.csv``, and ``manifest.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402
from phenocycler.integration import vocab  # noqa: E402

sys.path.insert(0, str(REPO))

from scripts.senior import xenium_release  # noqa: E402

CROSSWALK = REPO / "core/phenocycler/data/integration/vocab_crosswalk.csv"
ACQ_RUNS = REPO / "data/donor_acquisition_runs.csv"
DONOR_META = REPO / "data/development/typing/app/donor_metadata.csv"

PLN_PHENO = REPO / "data/development/typing/pln/pheno"
PANC_PHENO = REPO / "data/development/typing/cells"


@lru_cache(maxsize=None)
def xenium_source(tissue: str):
    """The current Xenium release for a tissue, from its own declared metadata.

    THIS FILE HELD THREE PINS, and the pancreas one was the damaging kind: it
    named ``release_panc_v1`` first and ``release`` as the fallback.
    ``release_panc_v1`` holds only ``_work/6539.*`` -- no ``calls/`` and no
    ``release.json`` -- so the pancreas comparison fell silently through to the
    two-donor 2026-09-07 tree and every lineage read INSUFFICIENT DONORS FOR A
    VERDICT while a 21-donor release sat on disk.  A fallback chain of literals
    cannot tell "the newer one is not here yet" from "the newer one is not a
    release"; the resolver says which, by name.
    """
    return xenium_release.resolve(tissue)


def display_path(path: Path) -> str:
    """Repo-relative where that is meaningful, absolute where it is not.

    The pancreas release is written OUTSIDE the repository, so ``relative_to(REPO)``
    raises on it.  A provenance field that crashes on the current release is worse
    than one that prints an absolute path.
    """
    try:
        return str(Path(path).relative_to(REPO))
    except ValueError:
        return str(path)

OUT = REPO / "data/development/integration/concordance"

# ---------------------------------------------------------------------------
# THRESHOLD, FIXED BEFORE THE NUMBERS
# ---------------------------------------------------------------------------
# A lineage share is a part of a composition over K comparable lineages, so the
# mean share is 100/K percentage points.  For the pLN that is three parts and a
# mean share of 33.3 pp.
#
# These sections are different blocks of the same donor, so *some* disagreement
# is sectioning and not method: which part of a node a section cuts -- cortex,
# paracortex, medulla, capsule -- moves a B-vs-T share, and no modality can be
# blamed for that.  What sectioning cannot do is reorder the composition.
#
#   |median over donors of (protein share - RNA share)| <= 10 pp  -> CONCORDANT
#       10 pp is under a third of the mean share and inside the range regional
#       sampling of one node plausibly produces.
#   10 pp < |median| <= 20 pp                                     -> GREY
#       usable only with the offset stated alongside every number.
#   |median| > 20 pp                                              -> STRUCTURAL
#       20 pp is 60% of the mean share; a shift that large changes which lineage
#       dominates, which sectioning does not do.
#
# Independently of the median, a DISPERSION test: if the inter-quartile range of
# the per-donor difference exceeds 20 pp the disagreement is donor-dependent, and
# a donor-dependent offset cannot be carried by any pooled cross-modality claim
# even when it happens to cancel in the median.
#
# And a RANK test: if the two modalities disagree about the rank order of the
# comparable lineages in more than half the donors, the composition is structurally
# different whatever the per-lineage medians say.
THRESHOLD_CONCORDANT_PP = 10.0
THRESHOLD_STRUCTURAL_PP = 20.0
THRESHOLD_DISPERSION_IQR_PP = 20.0
THRESHOLD_RANK_DISAGREE_FRAC = 0.50

N_PERM = 10000
RNG_SEED = 20260917

# PhenoCycler pLN immune_class -> crosswalk common lineage.
PLN_PHENO_CLASS_TO_COMMON = {
    "T_cell": "T_NK",
    "NK_ILC": "T_NK",
    "B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}

# PhenoCycler pancreas cell_class -> common lineage (non-immune).
PANC_PHENO_CLASS_TO_COMMON = {
    "Beta": "Endocrine",
    "Alpha": "Endocrine",
    "Delta": "Endocrine",
    "Exocrine": "Exocrine",
    "Blood vessel": "Vascular",
    "Nerve": "Neural",
}

# Pancreas `cell_class` values that get NO common lineage, and why each one gets none.
#
# Naming them is the whole point.  `common.value_counts(dropna=True)` drops every row that
# mapped to nothing, so a class that is in neither table left the per-lineage counts while
# `n_cells` went on counting it, and the parts stopped summing to the total they are reported
# against.  Measured 2026-09-18 over `typing/cells` (25,230,109 cells, 20 donors, all
# `qc_analysis_eligible`): 316,752 cells -- Lymphatic 122,884 and Excluded 193,868 -- were
# counted in `n_cells` and appeared in no column at all.  They are now counted by name.
#
# AUTHORITY.  The class vocabulary is `core/phenocycler/scripts/typing/refine.py`
# (`OWN_MARKERS` + `UNTOUCHED` + `UNRESOLVED`), restated here rather than imported for the
# reason `core/phenocycler/scripts/typing/build_cells.py` gives for restating `UNRESOLVED`:
# importing `refine` pulls in the whole celltype scoring stack.  What a non-call IS belongs to
# `core/phenocycler/phenocycler/integration/vocab.py` (`PHENO_UNRESOLVED`, `NON_CALLS`), which
# is imported above; `check_pancreas_class_vocabulary` checks these spellings against it and
# refuses a non-call vocab records and this file does not name.  Together with the
# reconciliation in `panc_protein_counts`, a class added upstream can no longer arrive unnamed:
# it is counted as unnamed and said out loud.
PANC_PHENO_CLASS_NOT_COMPARED = {
    "Excluded": (
        "excluded from every calculation. It is the exocrine-inside-an-islet annotation error "
        "the maintainer ruled out of analysis, not a cell population, so it has no share and "
        "may not enter a comparable total"
    ),
    "Lymphatic": (
        "OPEN, and deliberately NOT folded into Vascular here. The crosswalk names one protein "
        "Vascular label and this file feeds it from `Blood vessel` alone; whether the protein "
        "Lymphatic class belongs in the same quantity as the Xenium broad Vascular class is a "
        "method question about what is being compared, and answering it in this file would move "
        "every pancreas share (`n_comparable` is the denominator of all of them). Counted by "
        "name and left out of the comparison until that is ruled on"
    ),
    "Unresolved (no lineage evidence)": (
        "a recorded NON-CALL, not a lineage: vocab.NON_CALLS reason_id `no_identity_evidence` -- "
        "the producer looked and named no identity, so there is no share to compare. It is also "
        "not reachable in this column: `refine.py` is the only producer, it is not wired, and "
        "`build_cells._join_refinement` adds `refined_class` BESIDE `cell_class` without "
        "overwriting it. So this count is UNAVAILABLE, not zero, until something writes the "
        "class into the column this builder reads"
    ),
}

#: The not-compared labels that NOTHING CAN WRITE into `cell_class` today, so their count is
#: UNAVAILABLE rather than zero.
#:
#: `Unresolved (no lineage evidence)` is produced only by `core/phenocycler/scripts/typing/
#: refine.py`, which is written and not wired, and `build_cells._join_refinement` adds
#: `refined_class` BESIDE `cell_class` without overwriting it.  Printing 0 here would say "this
#: builder looked and found none", when it did not look -- the column cannot carry the class at
#: all.  That is the missing-measurement rule, and this file's own entry for the label states it.
#:
#: This is not a permanent claim.  `not_compared_count` below re-measures on every donor: the
#: moment a producer writes the class into the column, the count stops being unavailable and
#: becomes the real number, with no edit here.
PANC_PHENO_CLASS_NOT_IN_COLUMN = frozenset({"Unresolved (no lineage evidence)"})

#: Classes counted by a branch of their own rather than by the map: `Immune` goes through
#: `immune_definition`, which is where the crosswalk says the protein immune split lives.
PANC_PHENO_CLASS_COUNTED_ELSEWHERE = ("Immune",)


def not_compared_count(cell_class, label: str):
    """The count for one not-compared label, or ``None`` where the column cannot carry it.

    A zero and an unavailable are different answers and this repo does not let one stand for the
    other.  ``0`` means this builder asked and no cell of that class was there.  ``None`` means
    no producer can put the class in ``cell_class``, so the question was never asked.

    The reachability is re-measured, never assumed: a label declared unreachable that turns up in
    the column is reported as the count it actually has, because the declaration is about today's
    producers and a producer can land without this file hearing about it.
    """
    n = int(cell_class.eq(label).sum())
    if n == 0 and label in PANC_PHENO_CLASS_NOT_IN_COLUMN:
        return None
    return n


def not_compared_column(label: str) -> str:
    """The per-donor column a not-compared class is counted in.

    Derived from the label so the column cannot drift from the class it counts, and stripped
    of the parenthetical reason so `Unresolved (no lineage evidence)` is a usable column name.
    """
    return "n_notcompared_" + label.split(" (")[0].replace(" ", "_")


def check_pancreas_class_vocabulary() -> None:
    """The two tables above have to agree with the authority, and cover it.

    Refuses rather than absorbing: a non-call this file does not name would be counted as an
    unnamed class, which is visible but says the wrong thing about it -- vocab already recorded
    the producer's answer, and reporting it as an unreadable label calls that answer a defect.
    """
    misdeclared = sorted(l for l in PANC_PHENO_CLASS_NOT_COMPARED if l in PANC_PHENO_CLASS_TO_COMMON)
    if misdeclared:
        raise RuntimeError(
            f"{misdeclared} are both mapped to a common lineage and listed as not compared. "
            f"One class cannot be in and out of the comparison at the same time."
        )
    declared = {nc.phenocycler for nc in vocab.NON_CALLS if nc.phenocycler}
    unnamed = sorted(declared - set(PANC_PHENO_CLASS_NOT_COMPARED))
    if unnamed:
        raise RuntimeError(
            f"phenocycler.integration.vocab records the protein-side non-call(s) {unnamed}, "
            f"which this file does not name. vocab.NON_CALLS is the authority on what is a "
            f"recorded absence rather than an identity; add it to PANC_PHENO_CLASS_NOT_COMPARED "
            f"with the reason it carries no share."
        )
    named_noncalls = [l for l in PANC_PHENO_CLASS_NOT_COMPARED if vocab.is_noncall(l)]
    if not named_noncalls:
        raise RuntimeError(
            "no label in PANC_PHENO_CLASS_NOT_COMPARED is a non-call vocab recognises, which "
            "means the spelling here and the spelling in vocab have drifted apart into two "
            "classes."
        )
# PhenoCycler pancreas immune_definition -> common lineage.  The crosswalk records
# T_NK / B_Plasma / Myeloid as pheno_resolution == "subclass": recovered from the
# per-marker gates, not from the broad label, which is exactly this column.
PANC_PHENO_IMMDEF_TO_COMMON = {
    "Imm-T-cell": "T_NK",
    "Imm-NK-cell": "T_NK",
    "Imm-B-cell": "B_Plasma",
    "Imm-plasma-cell": "B_Plasma",
    "Imm-macrophage": "Myeloid",
    "Imm-dendritic-cell": "Myeloid",
    "Imm-neutrophil": "Granulocyte",
}


def log(msg: str) -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# crosswalk
# ---------------------------------------------------------------------------
def load_crosswalk() -> pd.DataFrame:
    cw = pd.read_csv(CROSSWALK, dtype=str).fillna("")
    return cw[cw["kind"] == "lineage"].copy()


def comparable_lineages(cw: pd.DataFrame, *, immune_only: bool) -> tuple[list[str], dict[str, str]]:
    """The lineages a protein-vs-RNA comparison may use, and why each other one is out.

    A lineage is comparable only when the crosswalk gives it a name on BOTH sides
    and neither side's resolution disqualifies it.  Granulocyte is excluded by
    ``xenium_resolution == fine_only`` -- it has no broad class on the RNA side, so
    a Granulocyte share is not the same quantity in the two modalities.
    """
    keep: list[str] = []
    dropped: dict[str, str] = {}
    for _, r in cw.iterrows():
        name = r["common"]
        if immune_only and r["is_immune"] != "yes":
            dropped[name] = f"not an immune lineage (is_immune={r['is_immune']!r})"
            continue
        if not r["phenocycler"]:
            dropped[name] = f"no PhenoCycler label (pheno_resolution={r['pheno_resolution']!r})"
            continue
        if not r["xenium"]:
            dropped[name] = f"no Xenium label (xenium_resolution={r['xenium_resolution']!r})"
            continue
        if r["xenium_resolution"] == "fine_only":
            dropped[name] = "xenium_resolution=fine_only; no broad RNA class, so not the same quantity"
            continue
        if r["pheno_resolution"] in {"no", "never"}:
            dropped[name] = f"pheno_resolution={r['pheno_resolution']!r}"
            continue
        keep.append(name)
    return keep, dropped


def asymmetry_notes(cw: pd.DataFrame, lineages: list[str]) -> dict[str, str]:
    notes = {}
    for _, r in cw.iterrows():
        if r["common"] not in lineages:
            continue
        bits = []
        if r["pheno_resolution"] != "yes":
            bits.append(f"protein side is {r['pheno_resolution']} ({r['pheno_note']})")
        if r["xenium_resolution"] != "yes":
            bits.append(f"RNA side is {r['xenium_resolution']} ({r['xenium_note']})")
        if r["collapsed"]:
            bits.append(f"folds {r['collapsed']}")
        if bits:
            notes[r["common"]] = "; ".join(bits)
    return notes


# ---------------------------------------------------------------------------
# readers (retry: a concurrent workflow owns these trees and may be mid-write)
# ---------------------------------------------------------------------------
def read_parquet_retry(path: Path, columns=None, attempts: int = 5):
    import time
    last = None
    for i in range(attempts):
        try:
            return pd.read_parquet(path, columns=columns)
        except Exception as exc:  # partial write -> retry, never record a partial count
            last = exc
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"could not read {path} after {attempts} attempts: {last}")


def pln_protein_counts(donors_excluded: set[str]) -> tuple[pd.DataFrame, dict]:
    status = pd.read_csv(PLN_PHENO / "per_donor_status.csv", dtype={"donor_id": str})
    rows, states = [], {}
    for _, s in status.iterrows():
        d = str(s["donor_id"])
        if d in donors_excluded:
            states[d] = "donor centrally excluded"
            continue
        if s["state"] != "typed":
            states[d] = f"{s['state']}: {s['reason']}"
            continue
        p = PLN_PHENO / f"donor_id={d}" / "roi=pln_1" / "data.parquet"
        if not p.exists():
            states[d] = "typed in status table but no partition present"
            continue
        df = read_parquet_retry(p, columns=["object_id", "is_immune", "immune_class"])
        ic = df["immune_class"].astype(str)
        common = ic.map(PLN_PHENO_CLASS_TO_COMMON)
        vc = common.value_counts(dropna=True)
        rec = {"donor_id": d, "n_cells": len(df), "n_immune": int(df["is_immune"].sum())}
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rec["n_immune_unnamed"] = int(df["is_immune"].sum() - vc.sum())
        rows.append(rec)
        states[d] = "typed"
        log(f"  protein pLN {d}: {len(df):,} cells")
    return pd.DataFrame(rows), states


def release_columns(release, *roles: str) -> list[str]:
    """The columns that fill these roles in THIS release, skipping the absent ones.

    A release that does not carry an optional column is not a release with a zero
    in it, so the column is simply not read and the quantity it fed is reported as
    unavailable further down.
    """
    out = []
    for role in roles:
        name = release.column(role)
        if name is not None and name not in out:
            out.append(name)
    return out


def pln_rna_counts(donors_excluded: set[str]) -> tuple[pd.DataFrame, dict]:
    release = xenium_source("pln")
    if not release.available:
        # Named absence, not an empty comparison that reads as agreement.
        log(f"  RNA pLN: {release}")
        return pd.DataFrame(), {}
    lineage_col = release.require_column("common_lineage")
    wanted = release_columns(release, "cell_key", "common_lineage", "label",
                             "nucleus_count", "nucleus_measured", "depth")
    if "y_probe_detected" in release.calls_schema().names:
        wanted.append("y_probe_detected")
    rows, states = [], {}
    for d, files in sorted(release.donor_partitions().items(),
                           key=lambda kv: (len(kv[0]), kv[0])):
        if d in donors_excluded:
            states[d] = "donor centrally excluded"
            continue
        if not files:
            states[d] = f"no partition under {release.release_id}"
            continue
        df = pd.concat([read_parquet_retry(f, columns=wanted) for f in files],
                       ignore_index=True)
        if lineage_col != "common_lineage":
            df = df.rename(columns={lineage_col: "common_lineage"})
        n_raw = len(df)
        # A cell with a MEASURED absent nucleus is excluded from every calculation.
        # A cell whose nucleus was never measured is a missing measurement -- not a
        # measured zero, and not evidence that the cell has no nucleus.  Some sections
        # carry no nucleus table at all (nucleus_count is the -1 sentinel for every
        # cell); dropping those would discard whole sections on the strength of an
        # absent measurement.  Those cells are kept and the section is flagged, and
        # the flag is carried into every downstream table.
        if {"nucleus_measured", "nucleus_count"} <= set(df.columns):
            measured = df["nucleus_measured"].fillna(False)
            no_nuc = int((measured & (df["nucleus_count"] <= 0)).sum())
            unmeasured = int((~measured).sum())
            df = df[~(measured & (df["nucleus_count"] <= 0))]
        else:
            # The nucleus rule is UNMEASURABLE here; recorded as such, never as passed.
            measured = pd.Series(np.nan, index=df.index, dtype=float)
            no_nuc = unmeasured = 0
        vc = df["common_lineage"].replace("", np.nan).value_counts(dropna=True)
        rec = {
            "donor_id": d, "n_cells": n_raw, "n_retained": len(df),
            "n_no_nucleus_excluded": no_nuc, "n_nucleus_unmeasured": unmeasured,
            "nucleus_measured_frac": float(measured.mean()) if n_raw else np.nan,
            "median_depth": (float(df["depth"].median())
                             if len(df) and "depth" in df.columns else np.nan),
            "y_probe_rate": (float(df["y_probe_detected"].mean())
                             if len(df) and "y_probe_detected" in df.columns else np.nan),
        }
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rows.append(rec)
        states[d] = "called"
        frac = float(measured.mean()) if n_raw else float("nan")
        log(f"  RNA pLN {d}: {n_raw:,} cells, {len(df):,} retained, "
            + (f"nucleus measured on {100*frac:.0f}%" if frac == frac
               else "nucleus rule UNMEASURABLE in this release"))
    return pd.DataFrame(rows), states


def panc_protein_counts(donors: set[str], donors_excluded: set[str]) -> pd.DataFrame:
    rows = []
    for d in sorted(donors):
        if d in donors_excluded:
            continue
        p = PANC_PHENO / f"donor_id={d}" / "data.parquet"
        if not p.exists():
            continue
        df = read_parquet_retry(
            p, columns=["object_id", "cell_class", "immune_definition", "qc_analysis_eligible"]
        )
        df = df[df["qc_analysis_eligible"].fillna(False)]
        cc = df["cell_class"].astype(str)
        common = cc.map(PANC_PHENO_CLASS_TO_COMMON)
        imm = cc.eq("Immune")
        common = common.where(~imm, df["immune_definition"].astype(str).map(PANC_PHENO_IMMDEF_TO_COMMON))
        vc = common.value_counts(dropna=True)
        rec = {"donor_id": d, "n_cells": len(df),
               "n_immune_unnamed": int((imm & common.isna()).sum())}
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)

        # Every class that gets no common lineage, counted BY NAME and on every donor -- so a
        # class with no cells prints its zero instead of the row disappearing. The columns come
        # from the table above rather than from this donor's data, because a count that only
        # exists where it is non-zero cannot tell an empty class from an unasked question.
        for label in PANC_PHENO_CLASS_NOT_COMPARED:
            rec[not_compared_column(label)] = not_compared_count(cc, label)

        # Anything in neither table. Said out loud rather than absorbed: the reconciliation
        # below would otherwise be the only sign, and a number that does not add up is a worse
        # report than a named class nobody has decided about yet.
        named = (set(PANC_PHENO_CLASS_TO_COMMON) | set(PANC_PHENO_CLASS_NOT_COMPARED)
                 | set(PANC_PHENO_CLASS_COUNTED_ELSEWHERE))
        unnamed = sorted(set(cc.unique()) - named)
        rec["n_class_unnamed"] = int(cc.isin(unnamed).sum())
        rec["classes_unnamed"] = "; ".join(unnamed)

        # The parts have to sum to the total they are reported against. Every eligible cell is
        # in exactly one of: a mapped lineage, Immune with no definition, a named not-compared
        # class, an unnamed class. If that stops holding, a class is leaving the table again.
        # An unavailable bucket contributes nothing to the sum, and cannot: the class is not in
        # the column, so none of these cells is in it.  That is a different statement from "its
        # count is 0", which is why the column carries null and this line carries the 0.
        parts = (int(vc.sum()) + rec["n_immune_unnamed"] + rec["n_class_unnamed"]
                 + sum(rec[not_compared_column(l)] or 0
                       for l in PANC_PHENO_CLASS_NOT_COMPARED))
        if parts != len(df):
            raise RuntimeError(
                f"{d}: the per-class parts sum to {parts:,} and n_cells is {len(df):,}. Every "
                f"cell has to be counted exactly once, in a lineage or in a named bucket; a "
                f"difference means a class is being counted twice or dropped without a word. "
                f"classes seen: {sorted(set(cc.unique()))}"
            )
        rows.append(rec)
        not_compared = ", ".join(
            f"{l.split(' (')[0]} "
            + ("unavailable" if rec[not_compared_column(l)] is None
               else f"{rec[not_compared_column(l)]:,}")
            for l in PANC_PHENO_CLASS_NOT_COMPARED)
        log(f"  protein panc {d}: {len(df):,} eligible cells; "
            f"not compared: {not_compared}"
            + (f"; UNNAMED classes {unnamed}" if unnamed else ""))
    out = pd.DataFrame(rows)
    # Nullable integers, so an unavailable bucket stays a MISSING COUNT rather than becoming
    # a float NaN or an object column.  It matters the day a producer lands and some donors
    # carry the class while others do not: the column has to hold both without either being
    # read as the other.
    for _label in PANC_PHENO_CLASS_NOT_COMPARED:
        _col = not_compared_column(_label)
        out[_col] = out[_col].astype("Int64")
    return out


def panc_rna_counts(donors_excluded: set[str]) -> tuple[pd.DataFrame, Path | None, list[str]]:
    release = xenium_source("panc")
    if not release.available:
        return pd.DataFrame(), None, [
            str(release),
            f"roots searched: {[str(r) for r in release.roots]}",
            f"candidates rejected: {list(release.rejected)}",
        ]
    notes = [
        f"pancreas RNA release_id = {release.release_id}, lane = {release.lane}, "
        f"declared {release.generated_at_field} = {release.generated_at.isoformat()}, "
        f"root = {display_path(release.root)}",
        "resolved by scripts/senior/xenium_release.py from the release's own declaration; "
        "never from a directory mtime, and a tree with no calls/ is rejected by name rather "
        "than silently losing the comparison to an older one",
    ]
    lineage_col = release.require_column("common_lineage")
    wanted = release_columns(release, "common_lineage", "nucleus_count", "nucleus_measured")
    rows = []
    for d, files in sorted(release.donor_partitions().items(),
                           key=lambda kv: (len(kv[0]), kv[0])):
        if d in donors_excluded:
            continue
        if not files:
            continue
        df = pd.concat([read_parquet_retry(f, columns=wanted) for f in files],
                       ignore_index=True)
        if lineage_col != "common_lineage":
            df = df.rename(columns={lineage_col: "common_lineage"})
        cols = set(df.columns)
        if "nucleus_count" in cols and "nucleus_measured" in cols:
            measured = df["nucleus_measured"].fillna(False)
            n_unmeasured = int((~measured).sum())
            # Two different absences. A MEASURED absent nucleus excludes the cell.
            # A never-measured nucleus is a missing measurement: the rule never
            # applied to that cell, so dropping it would discard whole donors on
            # the strength of a measurement nobody made. Measured 2026-09-17 in
            # this release: every cell of 6347, 6579 and 6580.
            df = df[~(measured & (df["nucleus_count"] <= 0))]
            nuc = (f"nucleus filter applied; {n_unmeasured} cells had no nucleus "
                   "measurement at all and are retained and flagged")
        else:
            nuc = "NO nucleus column in this release; nucleated-only filter could NOT be applied"
        vc = df["common_lineage"].replace("", np.nan).value_counts(dropna=True)
        rec = {"donor_id": d, "n_cells": len(df), "nucleus_handling": nuc}
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rows.append(rec)
        log(f"  RNA panc {d}: {len(df):,} cells ({nuc})")
    return pd.DataFrame(rows), release.root, notes


# ---------------------------------------------------------------------------
# composition
# ---------------------------------------------------------------------------
def shares(counts: pd.DataFrame, lineages: list[str]) -> pd.DataFrame:
    if counts.empty:
        return pd.DataFrame()
    cols = [f"n_{l}" for l in lineages]
    present = [c for c in cols if c in counts.columns]
    sub = counts.reindex(columns=["donor_id"] + present).copy()
    for c in cols:
        if c not in sub.columns:
            sub[c] = np.nan  # the label does not exist in this modality: unavailable, not zero
    sub["n_comparable"] = sub[present].sum(axis=1)
    for l in lineages:
        sub[f"share_{l}"] = 100.0 * sub[f"n_{l}"] / sub["n_comparable"]
    return sub


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    ok = ~np.isnan(values)
    v, g = values[ok], groups[ok]
    if len(v) < 3 or len(np.unique(g)) < 2:
        return np.nan
    gm = v.mean()
    ss_tot = ((v - gm) ** 2).sum()
    if ss_tot <= 0:
        return np.nan
    ss_b = sum(len(v[g == k]) * (v[g == k].mean() - gm) ** 2 for k in np.unique(g))
    return float(ss_b / ss_tot)


def eta_with_permutation(values: np.ndarray, groups: np.ndarray, rng) -> tuple[float, float, float]:
    """eta-squared, its permutation-null median, and p.

    The null is the same values with the batch labels shuffled over donors, which
    is the only null that answers "is this more batch structure than a random
    partition of these donors into these group sizes would show".
    """
    obs = eta_squared(values, groups)
    if np.isnan(obs):
        return np.nan, np.nan, np.nan
    ok = ~np.isnan(values)
    v, g = values[ok], groups[ok].copy()
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        rng.shuffle(g)
        null[i] = eta_squared(v, g)
    p = (1 + int((null >= obs).sum())) / (1 + N_PERM)
    return obs, float(np.nanmedian(null)), float(p)


MIN_DONORS_FOR_VERDICT = 4


def verdict(med: float, iqr: float, rank_disagree: float, n_donors: int) -> tuple[str, str]:
    if np.isnan(med):
        return "UNAVAILABLE", "one modality does not produce this label"
    if n_donors < MIN_DONORS_FOR_VERDICT:
        return "INSUFFICIENT DONORS FOR A VERDICT", (
            f"{n_donors} donor(s); a median and an inter-quartile range over fewer than "
            f"{MIN_DONORS_FOR_VERDICT} donors is not a measurement of a cohort, and the "
            "numbers below are reported as observations only"
        )
    reasons = []
    if abs(med) > THRESHOLD_STRUCTURAL_PP:
        reasons.append(f"median offset {med:+.1f} pp exceeds the {THRESHOLD_STRUCTURAL_PP:.0f} pp structural threshold")
    if iqr > THRESHOLD_DISPERSION_IQR_PP:
        reasons.append(f"IQR {iqr:.1f} pp exceeds the {THRESHOLD_DISPERSION_IQR_PP:.0f} pp dispersion threshold; the offset is donor-dependent")
    if rank_disagree > THRESHOLD_RANK_DISAGREE_FRAC:
        reasons.append(f"the modalities disagree on composition rank order in {100*rank_disagree:.0f}% of donors")
    if reasons:
        return "MAY NOT CARRY A CROSS-MODALITY CLAIM", "; ".join(reasons)
    if abs(med) <= THRESHOLD_CONCORDANT_PP:
        return "MAY CARRY A CROSS-MODALITY CLAIM", f"median offset {med:+.1f} pp, IQR {iqr:.1f} pp, both inside the thresholds fixed in advance"
    return "MAY CARRY A CROSS-MODALITY CLAIM ONLY WITH THE OFFSET STATED", f"median offset {med:+.1f} pp is in the grey band ({THRESHOLD_CONCORDANT_PP:.0f}-{THRESHOLD_STRUCTURAL_PP:.0f} pp)"


def compare(tissue: str, prot: pd.DataFrame, rna: pd.DataFrame, lineages: list[str],
            batch: pd.DataFrame, rng) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if prot.empty or rna.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    sp = shares(prot, lineages).set_index("donor_id")
    sr = shares(rna, lineages).set_index("donor_id")
    donors = sorted(set(sp.index) & set(sr.index), key=lambda x: (len(x), x))
    long = []
    for d in donors:
        for l in lineages:
            a, b = sp.loc[d, f"share_{l}"], sr.loc[d, f"share_{l}"]
            long.append({
                "tissue": tissue, "donor_id": d, "lineage": l,
                "n_protein": sp.loc[d, f"n_{l}"], "n_rna": sr.loc[d, f"n_{l}"],
                "protein_comparable_total": sp.loc[d, "n_comparable"],
                "rna_comparable_total": sr.loc[d, "n_comparable"],
                "share_protein_pct": a, "share_rna_pct": b,
                "diff_pp": a - b,
            })
    long = pd.DataFrame(long)

    # Rank-order disagreement.  Two readings, both reported:
    #   whole_rank_frac  -- fraction of donors where the FULL ordering of the
    #                       comparable set differs between modalities.  This is the
    #                       criterion as first written; it is one number for the
    #                       whole composition, so it cannot separate one lineage from
    #                       another and it is reported as context, not as a verdict.
    #   per-lineage      -- fraction of donors in which THIS lineage sits at a
    #                       different position in its own modality's ordering.  This
    #                       is the lineage-specific reading of the same criterion and
    #                       is what the per-lineage verdict uses.  The 50% threshold
    #                       is unchanged.
    whole, per_lineage = {}, {l: [] for l in lineages}
    for d in donors:
        pa = np.array([sp.loc[d, f"share_{l}"] for l in lineages], dtype=float)
        rb = np.array([sr.loc[d, f"share_{l}"] for l in lineages], dtype=float)
        if np.isnan(pa).any() or np.isnan(rb).any():
            continue
        op, orr = list(np.argsort(pa)), list(np.argsort(rb))
        whole[d] = int(op != orr)
        rp, rr = np.argsort(np.argsort(pa)), np.argsort(np.argsort(rb))
        for i, l in enumerate(lineages):
            per_lineage[l].append(int(rp[i] != rr[i]))
    whole_rank_frac = float(np.mean(list(whole.values()))) if whole else np.nan

    bat = batch.set_index("donor_id")
    summary, eta_rows = [], []
    for l in lineages:
        sub = long[long.lineage == l]
        dif = sub["diff_pp"].to_numpy(dtype=float)
        med = float(np.nanmedian(dif)) if len(dif) else np.nan
        iqr = float(np.nanpercentile(dif, 75) - np.nanpercentile(dif, 25)) if np.isfinite(dif).sum() > 1 else np.nan
        rank_frac = float(np.mean(per_lineage[l])) if per_lineage[l] else np.nan
        v, why = verdict(med, iqr, rank_frac, int(np.isfinite(dif).sum()))
        summary.append({
            "tissue": tissue, "lineage": l, "n_donors": int(np.isfinite(dif).sum()),
            "median_share_protein_pct": float(np.nanmedian(sub["share_protein_pct"])),
            "median_share_rna_pct": float(np.nanmedian(sub["share_rna_pct"])),
            "median_diff_pp": med, "iqr_diff_pp": iqr,
            "min_diff_pp": float(np.nanmin(dif)) if np.isfinite(dif).any() else np.nan,
            "max_diff_pp": float(np.nanmax(dif)) if np.isfinite(dif).any() else np.nan,
            "rank_disagreement_frac": rank_frac,
            "whole_composition_rank_disagreement_frac": whole_rank_frac,
            "verdict": v, "verdict_reason": why,
        })
        # batch eta-squared with a permutation null, for every cross-donor quantity
        for qname, series in (("share_protein_pct", sub.set_index("donor_id")["share_protein_pct"]),
                              ("share_rna_pct", sub.set_index("donor_id")["share_rna_pct"]),
                              ("diff_pp", sub.set_index("donor_id")["diff_pp"])):
            for bvar in ("acq_run", "stain_group"):
                idx = [d for d in series.index if d in bat.index]
                if len(idx) < 3:
                    continue
                vals = series.loc[idx].to_numpy(dtype=float)
                grp = bat.loc[idx, bvar].astype(str).to_numpy()
                obs, nullmed, p = eta_with_permutation(vals, grp, rng)
                eta_rows.append({
                    "tissue": tissue, "lineage": l, "quantity": qname, "batch_var": bvar,
                    "n_donors": len(idx), "n_groups": int(len(set(grp))),
                    "eta_squared": obs, "perm_null_median": nullmed, "p_perm": p,
                    "n_perm": N_PERM,
                })
    return long, pd.DataFrame(summary), pd.DataFrame(eta_rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    excluded = set(DONOR_EXCLUSIONS)
    log(f"centrally excluded donors (from cohort.DONOR_EXCLUSIONS): {sorted(excluded) or 'none'}")

    check_pancreas_class_vocabulary()
    log("pancreas protein classes not compared (counted by name, never dropped): "
        + ", ".join(sorted(PANC_PHENO_CLASS_NOT_COMPARED)))

    cw = load_crosswalk()
    pln_lin, pln_dropped = comparable_lineages(cw, immune_only=True)
    panc_lin, panc_dropped = comparable_lineages(cw, immune_only=False)
    log(f"comparable pLN lineages: {pln_lin}")
    log(f"comparable pancreas lineages: {panc_lin}")

    batch = pd.read_csv(ACQ_RUNS, dtype={"donor_id": str})[["donor_id", "acq_run", "stain_group"]]
    batch = batch.drop_duplicates("donor_id")
    meta = pd.read_csv(DONOR_META, dtype={"case_id": str}).rename(columns={"case_id": "donor_id"})
    meta = meta[["donor_id", "donor_status"]].drop_duplicates("donor_id")

    log("reading pLN protein ...")
    pln_p, pln_p_state = pln_protein_counts(excluded)
    log("reading pLN RNA ...")
    pln_r, pln_r_state = pln_rna_counts(excluded)

    pln_long, pln_sum, pln_eta = compare("pln", pln_p, pln_r, pln_lin, batch, rng)

    # Controls, reported beside the comparison rather than used to gate it.
    # The Y-probe rate is a control only on donors with no Y template, so it is
    # reported against sex.  Transcript depth is the flag that matters for a 5K
    # panel call; a low-depth section is flagged on depth, not on a control rate.
    controls = pln_r.merge(meta, on="donor_id", how="left")
    sex = pd.read_csv(DONOR_META, dtype={"case_id": str}).rename(columns={"case_id": "donor_id"})
    controls = controls.merge(sex[["donor_id", "sex"]].drop_duplicates("donor_id"), on="donor_id", how="left")
    controls = controls[["donor_id", "sex", "donor_status", "n_cells", "n_retained",
                         "n_no_nucleus_excluded", "nucleus_measured_frac",
                         "median_depth", "y_probe_rate"]]
    controls.to_csv(out / "pln_rna_controls.csv", index=False)
    log("\n=== pLN RNA controls (Y-probe rate is a control only where there is no Y template) ===")
    log(controls.sort_values(["sex", "y_probe_rate"]).to_string(index=False))

    log("reading pancreas RNA ...")
    panc_r, panc_root, panc_notes = panc_rna_counts(excluded)
    panc_donors = set(panc_r["donor_id"]) if not panc_r.empty else set()
    log(f"pancreas RNA donors available: {sorted(panc_donors) or 'none'}")
    log("reading pancreas protein ...")
    panc_p = panc_protein_counts(panc_donors, excluded)

    # a lineage with no count column on one side is unavailable there, not zero
    panc_lin_eff = [l for l in panc_lin
                    if (f"n_{l}" in panc_p.columns) and (f"n_{l}" in panc_r.columns)]
    panc_unavail = {l: "no count produced by one modality on this tissue"
                    for l in panc_lin if l not in panc_lin_eff}
    panc_long, panc_sum, panc_eta = compare("pancreas", panc_p, panc_r, panc_lin_eff, batch, rng)

    for name, df in (("pln_per_donor_counts_protein", pln_p),
                     ("pln_per_donor_counts_rna", pln_r),
                     ("pln_per_donor_per_lineage", pln_long),
                     ("pln_lineage_summary", pln_sum),
                     ("pln_batch_eta_squared", pln_eta),
                     ("pancreas_per_donor_counts_protein", panc_p),
                     ("pancreas_per_donor_counts_rna", panc_r),
                     ("pancreas_per_donor_per_lineage", panc_long),
                     ("pancreas_lineage_summary", panc_sum),
                     ("pancreas_batch_eta_squared", panc_eta)):
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.to_csv(out / f"{name}.csv", index=False)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "script": "scripts/senior/modality_concordance.py",
        "latch": "development; no production typing is declared and nothing here is a release",
        "join": "donor_id ONLY. No cell identifier crosses modalities: PhenoCycler object_id is a "
                "QuPath UUID, Xenium cell_id is 10x segmentation unique only within a bundle.",
        "claim_scope": "COMPOSITIONAL ONLY. Different sections of different blocks. No distance, "
                       "adjacency, neighbourhood or shared scatter is computed or licensed here.",
        "thresholds_fixed_before_measurement": {
            "concordant_pp": THRESHOLD_CONCORDANT_PP,
            "structural_pp": THRESHOLD_STRUCTURAL_PP,
            "dispersion_iqr_pp": THRESHOLD_DISPERSION_IQR_PP,
            "rank_disagreement_frac": THRESHOLD_RANK_DISAGREE_FRAC,
        },
        "donor_exclusions": {d: DONOR_EXCLUSIONS[d] for d in sorted(excluded)},
        "crosswalk": str(CROSSWALK.relative_to(REPO)),
        "crosswalk_is_generated": "run_vocab in phenocycler/integration/vocab.py writes it; never hand-edited here",
        "pln": {
            "comparable_lineages": pln_lin,
            "excluded_lineages": pln_dropped,
            "asymmetries": asymmetry_notes(cw, pln_lin),
            "protein_source": display_path(PLN_PHENO),
            "rna_source": display_path(xenium_source("pln").root)
                          if xenium_source("pln").available else None,
            "rna_release": xenium_source("pln").describe(),
            "protein_donor_state": pln_p_state,
            "rna_donor_state": pln_r_state,
            "donors_compared": sorted(set(pln_long.donor_id)) if not pln_long.empty else [],
        },
        "pancreas": {
            "comparable_lineages": panc_lin_eff,
            "excluded_lineages": {**panc_dropped, **panc_unavail},
            "asymmetries": asymmetry_notes(cw, panc_lin_eff),
            "protein_source": display_path(PANC_PHENO),
            "rna_source": display_path(panc_root) if panc_root else None,
            "rna_release": xenium_source("panc").describe(),
            "rna_notes": panc_notes,
            "protein_classes_not_compared": PANC_PHENO_CLASS_NOT_COMPARED,
            "protein_classes_not_compared_note":
                "these cell_class values get no common lineage and are counted by name in "
                "pancreas_per_donor_counts_protein.csv (n_notcompared_*), on every donor, so "
                "an empty class prints a zero rather than vanishing. n_cells minus the "
                "comparable lineages, n_immune_unnamed, these and n_class_unnamed is zero by "
                "construction; the reader refuses the donor if it is not.",
            "donors_compared": sorted(set(panc_long.donor_id)) if not panc_long.empty else [],
        },
        "detection": "Xenium detection is lognorm >= log1p(1) as applied upstream by the calling "
                     "lane; no raw-count floor is introduced here.",
        "nucleus_rule": "A cell whose nucleus was MEASURED absent is excluded from every "
                        "calculation.  A cell whose nucleus was never measured is a missing "
                        "measurement, not a measured zero; it is retained and its section is "
                        "flagged by nucleus_measured_frac, because discarding it would record an "
                        "unavailable measurement as an absent nucleus.",
        "min_donors_for_verdict": MIN_DONORS_FOR_VERDICT,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")

    for title, df in (("pLN lineage summary", pln_sum), ("pancreas lineage summary", panc_sum)):
        log(f"\n=== {title} ===")
        if df.empty:
            log("  (nothing to compare)")
        else:
            log(df.to_string(index=False))
    for title, df in (("pLN batch eta-squared", pln_eta), ("pancreas batch eta-squared", panc_eta)):
        log(f"\n=== {title} ===")
        if df.empty:
            log("  (none)")
        else:
            log(df.to_string(index=False))
    log(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
