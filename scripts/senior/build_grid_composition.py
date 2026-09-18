#!/usr/bin/env python
"""The four-cell grid: per-donor composition for both tissues in both modalities.

WHAT THIS EXISTS FOR.  The app reads ONE of the four cells of the grid
(pancreas/protein) because every app composition artifact is built off
``data/development/typing/immune/``, which covers ten donors.  Three other cells
have calls and no reader at all.  This builder does not widen the pancreas
immune lane; it reads each of the four call trees directly, harmonises their
labels through the GENERATED vocabulary crosswalk, and writes ONE donor-level
layer that a card can read.

WHERE IT WRITES.  ``$ISLET_LAYERS_DIR`` defaulting to ``/data/islet_layers``,
because that is the root ``ieCore::ie_layers_dir()`` resolves and every card
with it.  Two builders have already shipped artifacts to a repo path no card
reads; a layer nothing can read is the whole failure, so there is no second
default and no search path here either.

EXHAUSTING THE DATA BEFORE DECLARING AN ABSENCE.  The pancreas protein immune
lineage split does NOT require ``typing/immune/``.  ``typing/cells`` carries
``immune_definition`` on all 20 donors -- the per-marker definition gate the
crosswalk describes as "recovered from per-marker gates, not the broad label".
``typing/immune`` adds a FINE subtype (Mac-M2-like, T-CD4-helper, ...) on ten
donors and is read here only for that finer level.  The broad grid is therefore
20 donors on the protein pancreas, not ten.

WHAT IS NEVER DONE HERE.
* A cell with no nucleus is excluded from every calculation.  It is counted in
  the presence table and excluded from every denominator -- never zero.
* A lineage a lane did not resolve gets a row with ``n_cells`` NULL and a named
  reason, never a zero.  A measured zero and a hole must not print the same.
* No distance, adjacency, neighbourhood or KDE crosses the pancreas/pLN
  boundary.  They are different blocks.  This layer is COMPOSITIONAL and joins
  on ``donor_id`` and ``disease_status`` alone.
* Nothing is tuned toward anything.  The cross-modality rows report agreement
  and disagreement identically; neither side is a target for the other.
* No cohort fact is written as a literal.  Every donor id, cohort size, status
  label and lineage name in the output is read out of the data.

Outputs under ``$ISLET_LAYERS_DIR/grid_composition``:

``composition.parquet``   (tissue, modality, donor_id, lineage_common) -> counts
``fine.parquet``          the fine call inside a common lineage, where a lane has one
``presence.parquet``      (tissue, modality, donor_id) -> one of the five states + reason
``crossmodal_tnk.csv``    per-donor T_NK share of immune, both modalities, per tissue
``batch_eta_squared.csv`` eta-squared of every cross-donor quantity with a permutation null
``manifest.json``         sources, domain root, builder, commit, totals
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cohort import DONOR_EXCLUSIONS, REINSTATED  # noqa: E402

sys.path.insert(0, str(REPO))

from scripts.senior import xenium_release  # noqa: E402

# --------------------------------------------------------------------------
# Inputs.  Every one is a call tree that already exists; nothing is re-typed.
# --------------------------------------------------------------------------
CELLS = REPO / "data" / "development" / "typing" / "cells"
IMMUNE_FINE = REPO / "data" / "development" / "typing" / "immune"
PLN_PHENO = REPO / "data" / "development" / "typing" / "pln" / "pheno"


@lru_cache(maxsize=None)
def xenium_source(tissue: str):
    """The current Xenium release for a tissue -- resolved, never typed in here.

    Both releases used to be pinned as literals in this file, one of them an
    absolute path outside the repository.  Five builders each held their own pin
    and only this one held the current pair.  ``scripts/senior/xenium_release.py``
    now selects on the release's own declared ``release_id``, write time and
    tissue -- never a directory mtime -- and the pins are gone from here.  It
    returns a NAMED unavailable rather than None, so an absent release cannot
    reach a count as a zero.
    """
    return xenium_release.resolve(tissue)

CROSSWALK = REPO / "core" / "phenocycler" / "data" / "integration" / "vocab_crosswalk.csv"
DONOR_META = REPO / "data" / "development" / "typing" / "app" / "donor_metadata.csv"
ACQ_RUNS = REPO / "data" / "donor_acquisition_runs.csv"
PLN_NODES = REPO / "core" / "phenocycler" / "data" / "integration" / "pln_node_annotations.csv"

#: The domain this layer belongs to.  The legacy CODEX build and the run cohort
#: share five donors and no vocabulary, so a matching donor id is not a licence
#: to read: the card compares this against domain_data_root() and refuses a
#: mismatch by name rather than rendering the wrong cells for a right-looking donor.
DOMAIN_ROOT = REPO / "data" / "development" / "typing" / "app"

OUT_DEFAULT = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")) / "grid_composition"

#: The crosswalk's scope token for the lymph node.  Read out of the file, never guessed.
PLN_SCOPE_TOKEN = "pancreatic_lymph_node"

#: bucket_kind values.  Only `lineage` rows are a lineage; the rest are the
#: book-keeping that makes the parts sum to the whole.
BK_LINEAGE = "lineage"
BK_PARTIAL = "partial"       # a positive call that stops short of a lineage
BK_NEGATIVE = "negative"     # a measured "not that" -- the pLN lanes' non-immune bucket
BK_UNASSIGNED = "unassigned"  # looked at, no call
BK_WITHDRAWN = "withdrawn"   # called, then the call was taken back

STATE_MEASURED = "measured"
STATE_BLANK = "blank"
STATE_UNAVAILABLE = "unavailable"
STATE_WITHHELD = "withheld"
STATE_FAILED = "failed"
STATES = (STATE_MEASURED, STATE_BLANK, STATE_UNAVAILABLE, STATE_WITHHELD, STATE_FAILED)


class GridError(RuntimeError):
    """A contract this builder will not paper over."""


def refuse_cross_tissue(quantity: str, tissues) -> None:
    """A lymph node has no islet, so a distance to one is undefined, not zero.

    Nothing in this module computes such a quantity; the guard is here so that a
    later edit that tries to has to delete a refusal rather than forget a rule.
    """
    if len(set(tissues)) > 1:
        raise GridError(
            f"{quantity!r} spans {sorted(set(tissues))}. The pancreas and the pLN are "
            "different blocks: no distance, adjacency, neighbourhood or KDE crosses that "
            "boundary. Join on donor_id and disease_status instead."
        )


# --------------------------------------------------------------------------
# The crosswalk.  GENERATED by run_vocab in phenocycler/integration/vocab.py --
# read, never hand-edited, and never filtered down to the lineages that happen
# to agree.  The pheno_resolution / xenium_resolution pair says which comparison
# is right for a lineage and how to read a disagreement, so both travel with
# every row instead of being consumed here.
# --------------------------------------------------------------------------
def load_crosswalk() -> pd.DataFrame:
    cw = pd.read_csv(CROSSWALK, dtype=str).fillna("")
    lin = cw[cw["kind"].eq("lineage")].copy()
    if lin.empty:
        raise GridError(f"{CROSSWALK} carries no lineage rows")
    return lin


#: Fine PhenoCycler calls -> the crosswalk `common` lineage they belong to.
#: These are FINE labels; the crosswalk names lineages, so it cannot name them,
#: and the map below is validated against it: every target must be a `common`
#: the crosswalk declares, or the build refuses.
PANC_DEFINITION_TO_COMMON = {
    "Imm-T-cell": "T_NK",
    "Imm-NK-cell": "T_NK",
    "Imm-B-cell": "B_Plasma",
    "Imm-plasma-cell": "B_Plasma",
    "Imm-macrophage": "Myeloid",
    "Imm-dendritic-cell": "Myeloid",
    "Imm-neutrophil": "Granulocyte",
}
PANC_CLASS_TO_COMMON = {
    "Beta": "Endocrine",
    "Alpha": "Endocrine",
    "Delta": "Endocrine",
    "Exocrine": "Exocrine",
    "Blood vessel": "Vascular",
    "Lymphatic": "Vascular",
    "Nerve": "Neural",
}
PLN_CLASS_TO_COMMON = {
    "T_cell": "T_NK",
    "NK_ILC": "T_NK",
    "B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}


def validate_maps(crosswalk: pd.DataFrame) -> None:
    known = set(crosswalk["common"])
    for name, mapping in (
        ("PANC_DEFINITION_TO_COMMON", PANC_DEFINITION_TO_COMMON),
        ("PANC_CLASS_TO_COMMON", PANC_CLASS_TO_COMMON),
        ("PLN_CLASS_TO_COMMON", PLN_CLASS_TO_COMMON),
    ):
        unknown = sorted(set(mapping.values()) - known)
        if unknown:
            raise GridError(
                f"{name} targets {unknown}, which the crosswalk does not declare as a "
                f"common lineage. The crosswalk is generated; fix the producer, not this map."
            )


# --------------------------------------------------------------------------
# Donor roster.  The universe is the union of the four call trees, minus the
# centrally excluded donors.  No lane restates the exclusion list.
# --------------------------------------------------------------------------
def partition_donors(root: Path, subdir: str = "") -> dict[str, list[Path]]:
    base = root / subdir if subdir else root
    out: dict[str, list[Path]] = {}
    for part in sorted(base.glob("donor_id=*")):
        donor = part.name.split("=", 1)[1]
        files = sorted(part.rglob("*.parquet"))
        out[donor] = files
    return out


def eligible(donors) -> list[str]:
    return sorted({d for d in donors if d not in DONOR_EXCLUSIONS}, key=lambda d: (len(d), d))


def load_status() -> pd.DataFrame:
    """One disease status per donor, from every roster that records one.

    Three files record it under two spellings ("Aab+" and "AAb+").  They are
    reconciled by case-folding the token, and the DISPLAY spelling is taken from
    the donor metadata of the domain this layer belongs to -- not typed in here.
    A genuine disagreement between sources is an error, not a vote.
    """
    frames = []
    meta = pd.read_csv(DONOR_META, dtype=str)
    frames.append(
        meta.rename(columns={"case_id": "donor_id", "donor_status": "disease_status"})[
            ["donor_id", "disease_status"]
        ].assign(source=str(DONOR_META), is_domain_roster=True)
    )
    for tissue, name in (("panc", "per_donor_summary.csv"), ("pln", "controls.csv")):
        release = xenium_source(tissue)
        if not release.available:
            continue
        path = release.sidecar(name)
        if path is None:
            continue
        tab = pd.read_csv(path, dtype=str)
        if "disease_status" not in tab.columns:
            continue
        frames.append(
            tab[["donor_id", "disease_status"]].assign(
                source=str(path), is_domain_roster=False
            )
        )
    roster = pd.concat(frames, ignore_index=True).dropna(subset=["donor_id", "disease_status"])
    roster["canon"] = roster["disease_status"].str.strip().str.casefold()

    # Display spelling: whatever the domain's own roster uses for that canonical token.
    domain = roster[roster["is_domain_roster"]]
    display = (
        domain.drop_duplicates(["canon"]).set_index("canon")["disease_status"].to_dict()
    )
    conflicts = []
    rows = []
    for donor, grp in roster.groupby("donor_id"):
        canons = sorted(set(grp["canon"]))
        if len(canons) > 1:
            conflicts.append((donor, canons))
        canon = canons[0]
        rows.append(
            {
                "donor_id": donor,
                "disease_status": display.get(canon, grp["disease_status"].iloc[0]),
                "disease_status_canon": canon,
                "status_sources": "|".join(sorted(set(grp["source"]))),
            }
        )
    if conflicts:
        raise GridError(f"disease status disagrees between sources: {conflicts}")
    return pd.DataFrame(rows)


def load_batch() -> pd.DataFrame:
    """The FINEST batch variable each source records, per donor.

    Acquisition batch is aliased with disease status at the coarse level.  The
    finer variables are what a correction could ever use, so those are what is
    carried and what eta-squared is reported on.  The Fusion version is never
    one of them.  Where a source records no batch variable at all, that is said
    out loud rather than filled from another source's.
    """
    acq = pd.read_csv(ACQ_RUNS, dtype=str)
    acq = acq[acq["excluded"].astype(str).str.lower() != "true"]
    panc_protein = (
        acq.groupby("donor_id")
        .agg(stain_group=("stain_group", lambda s: "|".join(sorted(set(s)))),
             acq_run=("acq_run", lambda s: "|".join(sorted(set(s)))))
        .reset_index()
        .assign(tissue="pancreas", modality="protein")
    )
    frames = [panc_protein]
    panc = xenium_source("panc")
    summary = panc.sidecar("per_donor_summary.csv") if panc.available else None
    if summary is not None:
        tab = pd.read_csv(summary, dtype=str)
        if "batch" in tab.columns:
            frames.append(
                tab[["donor_id", "batch"]]
                .rename(columns={"batch": "acq_run"})
                .assign(tissue="pancreas", modality="RNA", stain_group=pd.NA)
            )
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------
# Readers.  One per call tree.  Each returns a Counter over
# (bucket_kind, lineage_common, reason) plus the two book-keeping totals, and
# a Counter for the fine level.  No reader opens a source it does not own.
# --------------------------------------------------------------------------
def read_pheno_pancreas(files: list[Path]):
    """typing/cells -- the subpanel classification, 20 donors.

    ``immune_definition`` is what makes the immune lineage split reachable on
    every donor: it is the per-marker definition gate, not the broad label, and
    it is on the same rows as the class.  ``typing/immune`` is not required for
    the broad grid and is not opened here.
    """
    broad, fine = Counter(), Counter()
    acq = Counter()
    for path in files:
        tab = pq.read_table(
            path,
            columns=["cell_class", "immune_definition", "excluded_reason",
                     "acq_run", "stain_group"],
        ).to_pandas()
        pairs = (tab["stain_group"].astype(str) + "|" + tab["acq_run"].astype(str)).value_counts()
        for pair, n in pairs.items():
            sg, _, ar = str(pair).partition("|")
            acq[(sg, ar)] += int(n)
        cls = tab["cell_class"].astype(str)
        dfn = tab["immune_definition"].astype(str)
        for value, n in cls.value_counts().items():
            if value in PANC_CLASS_TO_COMMON:
                broad[(BK_LINEAGE, PANC_CLASS_TO_COMMON[value], "")] += int(n)
        immune = dfn[cls.eq("Immune")]
        for value, n in immune.value_counts().items():
            if value in PANC_DEFINITION_TO_COMMON:
                broad[(BK_LINEAGE, PANC_DEFINITION_TO_COMMON[value], "")] += int(n)
                fine[(PANC_DEFINITION_TO_COMMON[value], value, "definition")] += int(n)
            else:
                broad[(BK_PARTIAL, "Immune_unresolved",
                       "called Immune; no per-marker definition resolved a lineage")] += int(n)
        excl = tab.loc[cls.eq("Excluded"), "excluded_reason"].astype(str)
        for value, n in excl.value_counts().items():
            broad[(BK_WITHDRAWN, "Withdrawn", value or "no reason recorded")] += int(n)
    return broad, fine, 0, dict(acq)


def pheno_pancreas_ledger() -> pd.DataFrame:
    """``typing/cells`` drops the nucleus-less cells upstream and records how many.

    The count is carried into the presence table so that the pancreas protein cell
    of the grid reports the rule as applied rather than as never needed.
    """
    path = CELLS / "build_ledger.csv"
    if not path.exists():
        return pd.DataFrame(columns=["donor_id", "dropped_no_nucleus"])
    return pd.read_csv(path, dtype={"donor_id": str})[["donor_id", "dropped_no_nucleus"]]


def read_xenium(files: list[Path]):
    """A Xenium release: ``common_lineage`` is the producer's own crosswalk read.

    TWO DIFFERENT ABSENCES, and collapsing them is the specific error this reader
    exists to avoid.  ``nucleus_measured`` False means the nucleus was NEVER
    MEASURED on that cell -- measured 2026-09-17: every cell of donors 6347, 6579
    and 6580 in the pancreas release, 2,989,849 cells.  That is a missing
    measurement, which is unavailable; it is NOT a measured absence of a nucleus,
    and treating it as one would silently exclude three whole donors under a rule
    that never applied to them.  ``nucleus_count == 0`` on a cell whose nucleus WAS
    measured is the real thing, and that cell is excluded from every calculation.
    Both leave the denominator; they are counted apart and reported apart.
    """
    broad, fine = Counter(), Counter()
    no_nucleus = 0
    nucleus_not_measured = 0
    for path in files:
        tab = pq.read_table(
            path, columns=["lineage", "common_lineage", "nucleus_count", "nucleus_measured"]
        ).to_pandas()
        measured_nucleus = tab["nucleus_measured"].fillna(False).astype(bool)
        nucleus_not_measured += int((~measured_nucleus).sum())
        tab = tab[measured_nucleus]
        nuc = pd.to_numeric(tab["nucleus_count"], errors="coerce")
        dropped = tab[(nuc <= 0) | nuc.isna()
                      | tab["lineage"].astype(str).eq("Unavailable_no_nucleus")]
        no_nucleus += int(len(dropped))
        keep = tab.drop(dropped.index)
        common = keep["common_lineage"].astype(str)
        lineage = keep["lineage"].astype(str)
        for value, n in common[common.ne("")].value_counts().items():
            broad[(BK_LINEAGE, value, "")] += int(n)
        rest = lineage[common.eq("")]
        for value, n in rest.value_counts().items():
            if value == "Non_immune":
                broad[(BK_NEGATIVE, "Non_immune",
                       "measured as not immune; this lane resolves immune identity only")] += int(n)
            elif value == "Immune_uncalled":
                broad[(BK_PARTIAL, "Immune_unresolved",
                       "called immune; no lineage resolved")] += int(n)
            else:
                broad[(BK_UNASSIGNED, "Unassigned", value)] += int(n)
        named = pd.DataFrame({"common": common, "lineage": lineage})
        named = named[named["common"].ne("")]
        for (common_label, fine_label), n in named.value_counts().items():
            fine[(common_label, fine_label, "lineage")] += int(n)
    return broad, fine, no_nucleus, {"nucleus_not_measured": nucleus_not_measured}


def read_pheno_pln(files: list[Path]):
    """typing/pln/pheno -- the pLN protein lane.

    ``t_nk_class`` sits INSIDE ``immune_class == 'T_cell'`` (measured: the two
    totals are equal), so it is the fine level of T_NK and is never added to it.
    """
    broad, fine = Counter(), Counter()
    for path in files:
        tab = pq.read_table(path, columns=["is_immune", "immune_class", "t_nk_class"]).to_pandas()
        cls = tab["immune_class"].astype(str)
        for value, n in cls.value_counts().items():
            if value in PLN_CLASS_TO_COMMON:
                broad[(BK_LINEAGE, PLN_CLASS_TO_COMMON[value], "")] += int(n)
            elif value == "Uncalled":
                broad[(BK_PARTIAL, "Immune_unresolved",
                       "clustered and named immune; no class resolved")] += int(n)
        non_immune = tab.loc[cls.eq(""), "is_immune"]
        n_false = int((non_immune == False).sum())  # noqa: E712 -- object dtype, `is False` fails
        n_null = int(non_immune.isna().sum())
        if n_false:
            broad[(BK_NEGATIVE, "Non_immune",
                   "measured as not immune; this lane resolves immune identity only")] += n_false
        if n_null:
            broad[(BK_UNASSIGNED, "Unassigned",
                   "not clustered; is_immune was never answered")] += n_null
        tnk = tab["t_nk_class"].astype(str)
        for value, n in tnk[tnk.ne("")].value_counts().items():
            fine[("T_NK", value, "t_nk_class")] += int(n)
    return broad, fine, 0, {}


def read_immune_fine(donor: str):
    """typing/immune -- the FINE pancreas protein subtype, where it exists.

    This is the ten-donor tree.  It is read for the subtype level only; the
    broad lineage above it does not depend on it.
    """
    part = IMMUNE_FINE / f"donor_id={donor}"
    files = sorted(part.rglob("*.parquet")) if part.exists() else []
    fine = Counter()
    for path in files:
        tab = pq.read_table(path, columns=["immune_class", "immune_subtype"]).to_pandas()
        sub = tab["immune_subtype"].astype("object")
        cls = tab["immune_class"].astype("object")
        for (c, s), n in Counter(zip(cls, sub)).items():
            if s is None or (isinstance(s, float) and np.isnan(s)) or s == "":
                continue
            common = PANC_DEFINITION_TO_COMMON.get(c, "Immune_unresolved")
            fine[(common, str(s), "subtype")] += int(n)
    return fine, len(files)


# --------------------------------------------------------------------------
# Presence: one row per (tissue, modality, donor) over the whole universe.
# --------------------------------------------------------------------------
def pln_pheno_producer_status() -> pd.DataFrame:
    path = PLN_PHENO / "per_donor_status.csv"
    if not path.exists():
        return pd.DataFrame(columns=["donor_id", "state", "reason"])
    tab = pd.read_csv(path, dtype=str).fillna("")
    mapped = {"typed": STATE_MEASURED, "unavailable": STATE_UNAVAILABLE,
              "withheld": STATE_WITHHELD, "failed": STATE_FAILED, "blank": STATE_BLANK}
    tab["state"] = tab["state"].map(mapped)
    if tab["state"].isna().any():
        bad = sorted(set(pd.read_csv(path, dtype=str)["state"]) - set(mapped))
        raise GridError(f"{path} carries states this layer cannot map: {bad}")
    return tab[["donor_id", "state", "reason"]]


def pln_node_rows() -> pd.DataFrame:
    tab = pd.read_csv(PLN_NODES, dtype=str, comment="#").fillna("")
    return tab


# --------------------------------------------------------------------------
# Eta-squared with a permutation null.
# --------------------------------------------------------------------------
def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    ok = np.isfinite(values)
    values, groups = values[ok], groups[ok]
    if values.size < 3 or len(set(groups.tolist())) < 2:
        return float("nan")
    grand = values.mean()
    ss_total = float(((values - grand) ** 2).sum())
    if ss_total <= 0:
        return float("nan")
    ss_between = 0.0
    for level in set(groups.tolist()):
        sel = groups == level
        ss_between += sel.sum() * (values[sel].mean() - grand) ** 2
    return float(ss_between / ss_total)


def eta_with_null(values, groups, n_perm: int, seed: int = 0) -> dict:
    """With 14-21 donors over many acquisition levels a large eta-squared is the
    expected value of noise, not a finding, so the observed statistic is reported
    beside its own permutation null and never alone."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=object)
    obs = eta_squared(values, groups)
    if not np.isfinite(obs) or n_perm <= 0:
        return {"eta_squared": obs, "null_median": float("nan"),
                "p_permutation": float("nan"), "n_perm": 0, "n_donors": int(np.isfinite(values).sum())}
    ok = np.isfinite(values)
    v, g = values[ok], groups[ok]
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = eta_squared(rng.permutation(v), g)
    return {
        "eta_squared": obs,
        "null_median": float(np.nanmedian(null)),
        "p_permutation": float((np.nansum(null >= obs) + 1) / (n_perm + 1)),
        "n_perm": n_perm,
        "n_donors": int(ok.sum()),
    }


# --------------------------------------------------------------------------
# Build.
# --------------------------------------------------------------------------
def _xenium_entry(tissue: str, grid_tissue: str, method: str) -> dict:
    """One RNA grid cell, carrying the resolved release rather than a pinned path.

    An unresolved release is carried as the resolver's own named absence: the
    entry still exists, its donors are empty, and every presence row for it says
    which roots were searched and why each candidate lost.  It never becomes a
    path that silently globs to nothing.
    """
    release = xenium_source(tissue)
    return {
        "tissue": grid_tissue, "modality": "RNA", "release": release,
        "root": release.root if release.available else None,
        "source": (str(release.root) if release.available
                   else f"no current Xenium {tissue} release: {release.reason}"),
        "calls": "calls", "reader": read_xenium, "method": method,
    }


def sources_table() -> list[dict]:
    return [
        {"tissue": "pancreas", "modality": "protein", "root": CELLS, "source": str(CELLS),
         "release": None, "calls": "",
         "reader": read_pheno_pancreas, "method": "subpanel classification"},
        _xenium_entry("panc", "pancreas", "within-donor cluster naming"),
        {"tissue": "pln", "modality": "protein", "root": PLN_PHENO, "source": str(PLN_PHENO),
         "release": None, "calls": "",
         "reader": read_pheno_pln, "method": "cluster-named immune classes"},
        _xenium_entry("pln", "pln", "immune partition with a T/NK split"),
    ]


def build(n_perm: int, out: Path) -> dict:
    crosswalk = load_crosswalk()
    validate_maps(crosswalk)
    status = load_status()
    batch = load_batch()
    node_rows = pln_node_rows()
    pln_status = pln_pheno_producer_status()
    ledger = pheno_pancreas_ledger()

    sources = sources_table()
    donors_by_source = {}
    for src in sources:
        release = src["release"]
        if release is None:
            parts = partition_donors(src["root"], src["calls"])
        elif release.available:
            parts = release.donor_partitions()
        else:
            parts = {}
        donors_by_source[(src["tissue"], src["modality"])] = parts
    universe = eligible(
        {d for parts in donors_by_source.values() for d in parts}
        | set(status["donor_id"])
    )
    # Every donor in the universe must be one the cohort module admits.  The
    # exclusion list lives once; this asserts against it rather than restating it.
    blocked = sorted(set(universe) & set(DONOR_EXCLUSIONS))
    if blocked:
        raise GridError(f"excluded donors reached the roster: {blocked}")

    comp_rows, fine_rows, presence_rows = [], [], []
    acq_observed = {}

    for src in sources:
        tissue, modality = src["tissue"], src["modality"]
        parts = donors_by_source[(tissue, modality)]
        for donor in universe:
            files = parts.get(donor, [])
            key = dict(tissue=tissue, modality=modality, donor_id=donor)
            if not files:
                release = src["release"]
                if release is not None and not release.available:
                    # The absence is the RELEASE, not the donor. Saying "the lane
                    # was not run on this donor" when no release resolved at all
                    # would name the wrong absence.
                    state, reason = STATE_UNAVAILABLE, (
                        f"no current Xenium {release.tissue} release: {release.reason}. "
                        f"Roots searched: {[str(r) for r in release.roots]}"
                    )
                else:
                    state, reason = absence_reason(
                        tissue, modality, donor,
                        set(donors_by_source[("pancreas", "protein")]),
                        pln_status, node_rows, status,
                    )
                presence_rows.append({
                    **key, "state": state, "reason": reason, "n_cells_assessed": pd.NA,
                    "n_no_nucleus_excluded": pd.NA, "n_nucleus_not_measured": pd.NA,
                    "n_files": 0, "source": src["source"], "method": src["method"],
                })
                continue
            broad, fine, no_nucleus, extras = src["reader"](files)
            extras = dict(extras)
            not_measured = int(extras.pop("nucleus_not_measured", 0))
            if extras:
                acq_observed[(tissue, modality, donor)] = extras
            if modality == "protein" and tissue == "pancreas":
                extra_fine, _ = read_immune_fine(donor)
                fine.update(extra_fine)
                ledger_row = ledger[ledger["donor_id"].eq(donor)]
                if len(ledger_row) == 1:
                    no_nucleus = int(ledger_row["dropped_no_nucleus"].iloc[0])
            assessed = sum(broad.values())
            book = {"n_no_nucleus_excluded": no_nucleus,
                    "n_nucleus_not_measured": not_measured,
                    "n_files": len(files), "source": src["source"],
                    "method": src["method"]}
            if assessed <= 0:
                if not_measured > 0:
                    reason = (
                        f"nucleus_measured is False on all {not_measured} cells of this "
                        "partition: the nucleus was never measured here, so the nucleus rule "
                        "cannot be applied and no count from it may enter a calculation. This "
                        "is a missing measurement, not a measured absence of nuclei."
                    )
                else:
                    reason = "partition present but holds no assessable cell"
                presence_rows.append({**key, "state": STATE_UNAVAILABLE, "reason": reason,
                                      "n_cells_assessed": 0, **book})
                continue
            presence_rows.append({**key, "state": STATE_MEASURED, "reason": "",
                                  "n_cells_assessed": assessed, **book})
            seen = set()
            for (bucket, label, reason), n in sorted(broad.items()):
                seen.add(label)
                comp_rows.append({
                    **key, "bucket_kind": bucket, "lineage_common": label,
                    "n_cells": int(n), "denom_assessed": assessed,
                    "frac_of_assessed": n / assessed,
                    "state": STATE_MEASURED, "reason": reason,
                    "source": src["source"],
                })
            for _, cw in crosswalk.iterrows():
                label = cw["common"]
                if label in seen or label == "Unassigned":
                    continue
                state, reason = lineage_absence(cw, tissue, modality)
                comp_rows.append({
                    **key, "bucket_kind": BK_LINEAGE, "lineage_common": label,
                    "n_cells": pd.NA, "denom_assessed": assessed,
                    "frac_of_assessed": np.nan, "state": state, "reason": reason,
                    "source": src["source"],
                })
            for (common, label, level), n in sorted(fine.items()):
                fine_rows.append({
                    **key, "lineage_common": common, "lineage_fine": label,
                    "fine_level": level, "n_cells": int(n),
                    "source": src["source"],
                })

    comp = pd.DataFrame(comp_rows)
    fine = pd.DataFrame(fine_rows)
    presence = pd.DataFrame(presence_rows)

    cw_meta = crosswalk.set_index("common")
    for col in ("pheno_resolution", "xenium_resolution", "is_immune", "pheno_note", "xenium_note"):
        comp[col] = comp["lineage_common"].map(cw_meta[col]).fillna("")
    comp["in_crosswalk"] = comp["lineage_common"].isin(cw_meta.index)
    comp = comp.merge(status[["donor_id", "disease_status"]], on="donor_id", how="left")

    # The measured disagreement between the crosswalk and what a lane actually
    # emitted.  Named, not resolved: the crosswalk is generated and this builder
    # does not edit it, and a lineage is not dropped for disagreeing.
    comp["resolution_conflict"] = ""
    # Only a LINEAGE row can disagree with the crosswalk about a lineage. The
    # book-keeping buckets (Unassigned, Non_immune, Withdrawn, Immune_unresolved)
    # are not lineages and are not held to a lineage's resolution column.
    lineage_rows = comp["bucket_kind"].eq(BK_LINEAGE) & comp["in_crosswalk"]
    prot = lineage_rows & comp["modality"].eq("protein") & comp["state"].eq(STATE_MEASURED)
    comp.loc[prot & comp["pheno_resolution"].isin(["no", "never"]), "resolution_conflict"] = (
        "crosswalk pheno_resolution says this lineage is not resolvable on the protein "
        "side, and this lane emitted the call anyway"
    )
    # Restricted to the pancreas: on the pLN side the crosswalk's own
    # immune_tissue_scope already says which lineages are in scope, and a
    # non-immune lineage being absent from an immune lane is not a disagreement.
    prot_missing = (lineage_rows & comp["tissue"].eq("pancreas")
                    & comp["modality"].eq("protein") & comp["state"].ne(STATE_MEASURED))
    comp.loc[prot_missing & comp["pheno_resolution"].eq("yes"), "resolution_conflict"] = (
        "crosswalk pheno_resolution says this lineage is resolvable on the protein side, "
        "and this lane emitted no such class"
    )
    xen = lineage_rows & comp["modality"].eq("RNA") & comp["state"].eq(STATE_MEASURED)
    comp.loc[xen & comp["xenium_resolution"].isin(["no", "never"]), "resolution_conflict"] = (
        "crosswalk xenium_resolution says this lineage is not resolvable in RNA, and this "
        "lane emitted the call anyway"
    )

    fine = fine.merge(status[["donor_id", "disease_status"]], on="donor_id", how="left")
    presence = presence.merge(status[["donor_id", "disease_status"]], on="donor_id", how="left")

    presence = flag_stale_citations(presence, node_rows, pln_status)

    tnk = crossmodal_tnk(comp)
    eta = batch_table(comp, batch, status, n_perm)

    coherence = coherence_checks(comp, presence)

    out.mkdir(parents=True, exist_ok=True)
    comp.to_parquet(out / "composition.parquet", index=False)
    comp.to_csv(out / "composition.csv", index=False)
    fine.to_parquet(out / "fine.parquet", index=False)
    presence.to_parquet(out / "presence.parquet", index=False)
    presence.to_csv(out / "presence.csv", index=False)
    tnk.to_csv(out / "crossmodal_tnk.csv", index=False)
    eta.to_csv(out / "batch_eta_squared.csv", index=False)

    manifest = {
        "layer": "grid_composition",
        "written_utc": datetime.now(timezone.utc).isoformat(),
        "builder": "scripts/senior/build_grid_composition.py",
        "commit": git_commit(),
        "source_domain_root": str(DOMAIN_ROOT.resolve()),
        "domain_note": (
            "This layer belongs to one domain. The legacy CODEX build and the run cohort "
            "share five donors and no vocabulary, so a card must refuse this artifact when "
            "domain_data_root() does not match source_domain_root."
        ),
        "crosswalk": str(CROSSWALK),
        "crosswalk_note": "generated by run_vocab in phenocycler/integration/vocab.py; read, never hand-edited",
        "cohort_exclusions": dict(DONOR_EXCLUSIONS),
        "cohort_reinstated": dict(REINSTATED),
        "xenium_releases": {
            tissue: xenium_source(tissue).describe()
            for tissue in xenium_release.TISSUES
        },
        "xenium_release_rule": (
            "resolved by scripts/senior/xenium_release.py from each release's own declared "
            "release_id, write time and tissue, honouring a declared supersedes; never from a "
            "directory mtime, which a touch reorders"
        ),
        "n_donors_in_universe": int(len(universe)),
        "grid": [
            {
                "tissue": s["tissue"], "modality": s["modality"],
                "root": s["source"], "method": s["method"],
                "n_donors_measured": int(
                    ((presence["tissue"] == s["tissue"])
                     & (presence["modality"] == s["modality"])
                     & (presence["state"] == STATE_MEASURED)).sum()
                ),
                "n_cells_assessed": int(
                    presence.loc[(presence["tissue"] == s["tissue"])
                                 & (presence["modality"] == s["modality"]),
                                 "n_cells_assessed"].fillna(0).sum()
                ),
                "n_no_nucleus_excluded": int(
                    presence.loc[(presence["tissue"] == s["tissue"])
                                 & (presence["modality"] == s["modality"]),
                                 "n_no_nucleus_excluded"].fillna(0).sum()
                ),
                "n_nucleus_not_measured": int(
                    presence.loc[(presence["tissue"] == s["tissue"])
                                 & (presence["modality"] == s["modality"]),
                                 "n_nucleus_not_measured"].fillna(0).sum()
                ),
                "donors_measured": sorted(
                    presence.loc[(presence["tissue"] == s["tissue"])
                                 & (presence["modality"] == s["modality"])
                                 & (presence["state"] == STATE_MEASURED), "donor_id"]
                ),
            }
            for s in sources
        ],
        "nucleus_rule": "a cell with no nucleus is excluded from every calculation",
        "missing_rule": "a lineage a lane did not resolve carries n_cells NULL and a reason, never 0",
        "cross_tissue_rule": (
            "compositional only, joined on donor_id and disease_status; no distance, "
            "adjacency, neighbourhood or KDE crosses the pancreas/pLN boundary"
        ),
        "coherence": coherence,
        "acquisition_levels_observed": {
            f"{t}/{m}/{d}": {f"{sg}|{ar}": int(n) for (sg, ar), n in a.items()}
            for (t, m, d), a in acq_observed.items()
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    return manifest


def lineage_absence(cw_row, tissue: str, modality: str) -> tuple[str, str]:
    """Why this lane emitted no cell of this lineage -- never a zero.

    Reads the crosswalk's own columns. `immune_tissue_scope` is what says whether
    a lineage is even in scope for the lymph node, so an out-of-scope lineage is
    reported as out of scope rather than as an empty measurement.
    """
    label = cw_row["common"]
    scope = cw_row["immune_tissue_scope"]
    is_immune = cw_row["is_immune"]
    if tissue == "pln":
        if is_immune != "yes":
            return (STATE_BLANK,
                    f"{label} is not an immune lineage and this pLN lane resolves immune "
                    f"identity only; its non-immune cells are one measured bucket (Non_immune)")
        if scope and PLN_SCOPE_TOKEN not in scope.split("|"):
            return (STATE_BLANK,
                    f"crosswalk immune_tissue_scope ({scope}) does not extend {label} to the "
                    f"pancreatic lymph node")
    resolution = cw_row["pheno_resolution"] if modality == "protein" else cw_row["xenium_resolution"]
    note = cw_row["pheno_note"] if modality == "protein" else cw_row["xenium_note"]
    if resolution in ("no", "never"):
        return (STATE_BLANK,
                f"{label}: crosswalk "
                f"{('pheno' if modality == 'protein' else 'xenium')}_resolution={resolution}"
                f" -- {note or 'not resolvable in this modality'}")
    return (STATE_BLANK,
            f"this lane emitted no {label} call although the crosswalk records "
            f"{('pheno' if modality == 'protein' else 'xenium')}_resolution={resolution}"
            f"{': ' + note if note else ''}")


def absence_reason(tissue, modality, donor, parts_pancreas_protein, pln_status, node_rows, status):
    """Why a whole grid cell is absent for a donor -- with the producer's own words
    where the producer recorded them."""
    if tissue == "pln" and modality == "protein":
        row = pln_status[pln_status["donor_id"].eq(donor)]
        if len(row) == 1:
            return row["state"].iloc[0], row["reason"].iloc[0] or "recorded by the producer without a reason"
        if donor not in set(node_rows["donor_id"]):
            if donor not in parts_pancreas_protein:
                return (STATE_BLANK,
                        "no PhenoCycler acquisition of either tissue for this donor; it "
                        "reaches the cohort through the Xenium lane only")
            return (STATE_BLANK,
                    "the donor has a PhenoCycler pancreas section but no pLN image in "
                    "core/phenocycler/data/integration/pln_node_annotations.csv, so there "
                    "is no lymph node section to type")
        return (STATE_BLANK, "a pLN image exists but this lane recorded no status for the donor")
    if tissue == "pln" and modality == "RNA":
        return (STATE_BLANK,
                "no partition for this donor in the pLN Xenium release; the lane was not run "
                "on a lymph node section of this donor")
    if modality == "protein":
        if donor in REINSTATED:
            return (STATE_BLANK,
                    "no PhenoCycler pancreas partition for this donor. cohort.py records it "
                    f"as reinstated ({REINSTATED[donor]}), and the call tree predates that "
                    "notice, so the absence may be the lapsed exclusion rather than a missing "
                    "acquisition. Nothing here re-types it.")
        return (STATE_BLANK,
                "no PhenoCycler pancreas partition for this donor; this donor reaches the "
                "cohort through the Xenium lane only")
    return (STATE_BLANK,
            "no partition for this donor in the pancreas Xenium release; the lane was not run "
            "on a pancreas section of this donor")


def flag_stale_citations(presence: pd.DataFrame, node_rows: pd.DataFrame,
                         pln_status: pd.DataFrame) -> pd.DataFrame:
    """A producer's recorded reason can outlive the fact it cites.

    Measured here rather than corrected: this builder does not own the pLN protein
    lane and does not re-type a donor.  It says, on the row, that the reason no
    longer matches the file it names.
    """
    presence = presence.copy()
    presence["citation_conflict"] = ""
    nodes = node_rows.set_index("donor_id")
    sel = presence["tissue"].eq("pln") & presence["modality"].eq("protein")
    for idx in presence.index[sel]:
        reason = str(presence.at[idx, "reason"])
        donor = presence.at[idx, "donor_id"]
        if "pln_node_annotations" not in reason or donor not in nodes.index:
            continue
        excluded = str(nodes.at[donor, "excluded"]).strip().lower() == "true"
        cited_excluded = "excluded=true" in reason.lower()
        if cited_excluded and not excluded:
            presence.at[idx, "citation_conflict"] = (
                f"the reason cites pln_node_annotations.csv excluded=True; that file now "
                f"records excluded=False for this donor and {nodes.at[donor, 'n_cells']} cells. "
                f"The state stands as the producer recorded it; the citation does not."
            )
    return presence


def crossmodal_tnk(comp: pd.DataFrame) -> pd.DataFrame:
    """T_NK as a share of the immune cells a lane resolved, per donor, both modalities.

    The two modalities are different SECTIONS of different BLOCKS.  Agreement is
    informative and so is disagreement; neither side is tuned toward the other and
    neither is corrected by the other.  The offset is reported as it falls.
    """
    immune = comp[comp["is_immune"].eq("yes") & comp["state"].eq(STATE_MEASURED)
                  & comp["bucket_kind"].eq(BK_LINEAGE)]
    grouped = immune.groupby(["tissue", "modality", "donor_id"])["n_cells"].sum().rename("n_immune")
    tnk = (immune[immune["lineage_common"].eq("T_NK")]
           .groupby(["tissue", "modality", "donor_id"])["n_cells"].sum().rename("n_T_NK"))
    tab = pd.concat([grouped, tnk], axis=1).reset_index()
    tab["n_T_NK"] = tab["n_T_NK"].fillna(0).astype(int)
    tab["frac_T_NK_of_resolved_immune"] = tab["n_T_NK"] / tab["n_immune"]
    counts = tab.pivot_table(index=["tissue", "donor_id"], columns="modality",
                             values="n_immune", aggfunc="sum")
    wide = tab.pivot_table(index=["tissue", "donor_id"], columns="modality",
                           values="frac_T_NK_of_resolved_immune").reset_index()
    for modality in ("protein", "RNA"):
        if modality not in wide.columns:
            wide[modality] = np.nan
        wide[f"n_resolved_immune_{modality}"] = (
            counts[modality].reindex(
                pd.MultiIndex.from_frame(wide[["tissue", "donor_id"]])
            ).to_numpy() if modality in counts.columns else np.nan
        )
    wide["offset_pp_RNA_minus_protein"] = 100.0 * (wide["RNA"] - wide["protein"])
    wide["both_modalities"] = wide[["protein", "RNA"]].notna().all(axis=1)
    # Ranks ONLY among the donors a tissue measures in both modalities. A rank
    # taken over a different donor set in each modality would compare two
    # different cohorts and call the result a disagreement.
    wide["rank_protein"] = np.nan
    wide["rank_RNA"] = np.nan
    wide["rank_disagrees"] = pd.NA
    for tissue, grp in wide.groupby("tissue"):
        shared = grp[grp["both_modalities"]]
        if len(shared) < 2:
            continue
        rp = shared["protein"].rank(method="average")
        rr = shared["RNA"].rank(method="average")
        wide.loc[shared.index, "rank_protein"] = rp
        wide.loc[shared.index, "rank_RNA"] = rr
        wide.loc[shared.index, "rank_disagrees"] = (rp != rr)
    return wide.rename(columns={"protein": "frac_T_NK_protein", "RNA": "frac_T_NK_RNA"})


def batch_table(comp: pd.DataFrame, batch: pd.DataFrame, status: pd.DataFrame,
                n_perm: int) -> pd.DataFrame:
    """Eta-squared of every per-donor lineage fraction on disease status and on the
    finest batch variable each source records, each against its own permutation null.

    Coarse acquisition batch is aliased with disease status, so a correction is only
    ever defensible on the finer variables.  The Fusion version is not one of them
    and is not read here.
    """
    measured = comp[comp["state"].eq(STATE_MEASURED) & comp["bucket_kind"].eq(BK_LINEAGE)]
    rows = []
    keyed = batch.set_index(["tissue", "modality", "donor_id"])
    for (tissue, modality, lineage), grp in measured.groupby(
        ["tissue", "modality", "lineage_common"]
    ):
        grp = grp.merge(status[["donor_id", "disease_status"]], on="donor_id",
                        how="left", suffixes=("", "_roster"))
        values = grp["frac_of_assessed"].to_numpy(dtype=float)
        for variable in ("disease_status", "stain_group", "acq_run"):
            if variable == "disease_status":
                groups = grp["disease_status"].astype(str).to_numpy(dtype=object)
                available = not grp["disease_status"].isna().all()
                reason = ""
            else:
                lookup = []
                for donor in grp["donor_id"]:
                    try:
                        value = keyed.loc[(tissue, modality, donor), variable]
                    except KeyError:
                        value = None
                    if isinstance(value, pd.Series):
                        value = value.iloc[0]
                    lookup.append(None if pd.isna(value) else str(value))
                groups = np.asarray(lookup, dtype=object)
                available = any(v is not None for v in lookup)
                reason = ("" if available else
                          f"no {variable} is recorded for {tissue}/{modality} in this checkout")
            if not available:
                rows.append({"tissue": tissue, "modality": modality, "lineage_common": lineage,
                             "variable": variable, "state": STATE_BLANK, "reason": reason,
                             "eta_squared": np.nan, "null_median": np.nan,
                             "p_permutation": np.nan, "n_perm": 0,
                             "n_donors": int(len(grp)), "n_levels": 0})
                continue
            ok = np.array([g is not None for g in groups])
            stat = eta_with_null(values[ok], groups[ok], n_perm=n_perm)
            rows.append({"tissue": tissue, "modality": modality, "lineage_common": lineage,
                         "variable": variable, "state": STATE_MEASURED, "reason": "",
                         **stat, "n_levels": int(len(set(groups[ok].tolist())))})
    return pd.DataFrame(rows)


def coherence_checks(comp: pd.DataFrame, presence: pd.DataFrame) -> dict:
    """Internal coherence, which is legitimate: the subsets sum to their parent.

    No previous biological result is a target here.  What is checked is that the
    buckets of a grid cell add up to the number of cells that cell assessed.
    """
    measured = comp[comp["state"].eq(STATE_MEASURED)]
    summed = measured.groupby(["tissue", "modality", "donor_id"])["n_cells"].sum()
    declared = (presence[presence["state"].eq(STATE_MEASURED)]
                .set_index(["tissue", "modality", "donor_id"])["n_cells_assessed"])
    joined = pd.concat([summed.rename("summed"), declared.rename("declared")], axis=1).dropna()
    mismatched = joined[joined["summed"] != joined["declared"]]
    return {
        "rule": "every bucket of a grid cell sums to that cell's assessed cells",
        "n_grid_cells_checked": int(len(joined)),
        "n_mismatched": int(len(mismatched)),
        "mismatches": [
            {"key": list(k), "summed": int(v["summed"]), "declared": int(v["declared"])}
            for k, v in mismatched.iterrows()
        ],
    }


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        return "unavailable"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT,
                    help="layer directory; defaults to $ISLET_LAYERS_DIR/grid_composition")
    ap.add_argument("--n-perm", type=int, default=2000)
    args = ap.parse_args(argv)
    manifest = build(args.n_perm, args.out)
    print(json.dumps({k: manifest[k] for k in ("grid", "n_donors_in_universe", "coherence")},
                     indent=2, default=str))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
