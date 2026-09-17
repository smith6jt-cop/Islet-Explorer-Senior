#!/usr/bin/env python3
"""Aggregate current typing against UUID-keyed islet membership, never QuPath Parent.

The output is ``$ISLET_LAYERS_DIR/distance_profiles/data.parquet``. In addition
to the columns in ISLET_LAYERS_SPEC 1.3, ``state``/``reason`` distinguish absent
measurements, ``compartment`` preserves outline membership for the three
compartment card, and ``roi``/``tissue``/``modality`` name the section a row was
measured on. Core membership follows ``in_islet``; the peri-islet expansion is
outside the core with 0 < signed distance <= 20 um. Other pancreas cells remain
in the gradient but do not enter either compartment.

``n_cells`` is the eligible, non-excluded denominator, repeated for each immune
subtype/provenance pair. Do not sum denominators across those pairs. ``Immune``
with ``resolved_by=all`` is the total; other rows preserve the existing calls,
including supporting-tier cells whose subtype was not resolved. No typing or
marker threshold is fitted here. Missing UUID coverage withholds the whole
donor's profile instead of silently changing its denominator.

Cell files have no donor column: only their explicit donor_id= partition is
read. Membership files are read with partitioning=None and joined on object_id
alone, including when their optional donor column is absent. Root-level audit
files such as membership_diff.parquet are not cell partitions.

Two arms, and they are not the same kind of measurement
-------------------------------------------------------
* The **pancreas arm** (``roi=panc``, ``modality=phenocycler``) is spatial. It
  reads ``data/development/typing/cells`` and ``typing/islet_membership`` and
  carries ``frame=signed_distance_um``.
* The **pLN arm** (``roi=pln_*``, ``modality=xenium``) is **compositional only**
  and carries ``frame=none`` with a null ``distance_bin_um``. The pancreas and
  pLN sections mounted on one slide come from **different blocks**: they share a
  pixel grid and nothing anatomically. A signed distance to an islet is
  meaningless in a lymph node because there are no islets in it. So this arm
  joins to the pancreas arm on ``donor_id`` and ``donor_status`` and on nothing
  else. ``refuse_cross_tissue_distance`` makes that structural rather than
  conventional: pointing the distance path at a non-pancreas roi raises instead
  of emitting a number.

  Its immune fraction needs one thing the label alone cannot supply: a recorded
  statement of which lineages are immune. That is read from
  ``core/phenocycler/data/integration/vocab_crosswalk.csv`` (``is_immune``,
  ``immune_tissue_scope``), written by ``phenocycler.integration.vocab``.
  Maintainer ruling, 2026-09-17: "All immune cell calls from pancreas can be
  applied to pLN" -- the pancreas immune vocabulary is the pLN immune
  vocabulary. That licenses the *identity*, not the geometry, and it is not a
  licence to copy label rows: the panc and pLN bundles are segmented separately,
  their cell_id namespaces are disjoint, and the typer has to be run on the pLN
  cells. A lineage the crosswalk does not classify is withheld by name; it is
  never defaulted to non-immune.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from scripts.senior.build_sections_layer import filter_eligible_donors, norm_donor

TYPING = REPO / "data" / "development" / "typing"
INTEGRATION = REPO / "core" / "phenocycler" / "data" / "integration"
XEN_CELLS = INTEGRATION / "cells_xen"
MANIFEST = INTEGRATION / "manifest.csv"

# Where the two pLN modalities are typed. Both are read if present and neither
# is waited on: a missing one yields a named unavailable, never a zero and never
# a silently shorter artifact.
#
# cells_xen at roi=pln_1 is the INGEST, and its lineage columns are empty on
# every row -- the pLN cells were never typed there. Reading it for composition
# produces nothing but "not called", so the Xenium pLN arm reads a typed release
# instead. The panc and pLN bundles are segmented separately, so labels could
# never have been copied across in any case.
XEN_TYPING = REPO / "data" / "development" / "xenium_typing"
PLN_PHENO_CALLS = TYPING / "pln" / "pheno"
# The PhenoCycler pLN ingest. Used only to learn WHICH donors have a pLN section
# when the calls are absent, so the refusal names real donors rather than a
# hardcoded cohort list.
PLN_INGEST = INTEGRATION / "cells_pln"

# Hive partition values are read as strings, always. Every pLN donor id is all
# digits, so inferred partitioning turns '0112' into the integer 112 and the id
# stops joining against the string donor_id everywhere else.
PLN_PARTITIONING = ds.partitioning(
    pa.schema([("donor_id", pa.string()), ("roi", pa.string())]), flavor="hive"
)

PLN_TISSUE = "pancreatic_lymph_node"
# The release's own tissue tag, which is not the crosswalk's vocabulary.
RELEASE_TISSUE_ALIAS = {"pln": PLN_TISSUE, "pancreatic_lymph_node": PLN_TISSUE}

# A Xenium pLN release resolves immune membership in stage 1 and identity in
# stage 2. These are the stage-1 answers that resolve; anything else (Ambiguous,
# or an empty string) is not an answer and is never counted as either.
STAGE1_IMMUNE, STAGE1_NON_IMMUNE = "Immune", "Non_immune"

NO_PLN_PHENO = "phenocycler_pln_calls_absent"
NO_XEN_PLN_RELEASE = "xenium_pln_release_absent"

PANCREAS_ROI = "panc"
PANCREAS_TISSUE = "pancreas"
# The Xenium lineage column carried by cells_xen. An empty string is not a call.
XEN_LINEAGE = "lineage_common"

SUPERSEDED = (
    "Refusing superseded typing/subpanel/subpanel_calls.parquet: "
    "Beta 191,485 / Alpha 246,055 / Delta 62,797; "
    "current typing/cells has Beta 226,183 / Alpha 279,686 / Delta 90,658. "
    "Use typing/cells and typing/islet_membership."
)
CROSS_TISSUE = (
    "Refusing a signed distance to an islet for roi={roi!r}{tissue}: the pancreas "
    "and pLN sections mounted on one slide are from different blocks. They share a "
    "pixel grid and share nothing anatomically, so no distance, adjacency, "
    "neighbourhood, KDE or shared scatter may cross the tissue boundary, and there "
    "are no islets in a lymph node to measure a distance from. The pLN arm is a "
    "compositional comparison joined on donor_id and donor_status only."
)
# Which lineages are immune is read, never inferred from a label's spelling. The
# crosswalk (core/phenocycler/data/integration/vocab_crosswalk.csv) now records
# it per lineage in `is_immune`, with `immune_tissue_scope` naming the tissues
# that answer holds in -- written by phenocycler.integration.vocab, which is the
# only thing that may write that file. Every branch below still fires for
# anything the file does not record: this builder withholds a named number, it
# never invents a mapping and it never defaults an unrecorded lineage to
# non-immune.
VOCAB_CROSSWALK = INTEGRATION / "vocab_crosswalk.csv"
NO_IMMUNE_VOCAB = "xenium_immune_vocabulary_not_recorded"
# A lineage the data carries but the crosswalk does not classify. Named, and
# never silently treated as non-immune.
UNRECORDED_LINEAGE = "xenium_immune_lineage_not_recorded"
# The crosswalk records is-immune, but not for the tissue this section is.
OUT_OF_SCOPE_TISSUE = "xenium_immune_vocabulary_not_recorded_for_tissue"
IMMUNE_YES, IMMUNE_NO, IMMUNE_UNRESOLVABLE = "yes", "no", "unresolvable"


def load_immune_vocabulary(path: str | Path = VOCAB_CROSSWALK) -> dict | None:
    """Read is-immune per common lineage from the crosswalk, or None if unrecorded.

    Returns ``{"is_immune": {lineage: yes|no|unresolvable}, "scope": {tissue}}``.
    Marker rows are skipped; they are not lineages. A file without the column, or
    with the column present but no lineage classified, returns None so the caller
    withholds by name instead of guessing.
    """
    path = Path(path)
    if not path.exists():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "is_immune" not in rows[0]:
        return None
    is_immune: dict[str, str] = {}
    scope: set[str] = set()
    for row in rows:
        if (row.get("kind") or "").strip() != "lineage":
            continue
        lineage = (row.get("common") or "").strip()
        value = (row.get("is_immune") or "").strip()
        if not lineage or value not in (IMMUNE_YES, IMMUNE_NO, IMMUNE_UNRESOLVABLE):
            continue
        previous = is_immune.setdefault(lineage, value)
        if previous != value:
            raise ValueError(
                f"{path}: lineage {lineage!r} is recorded both {previous!r} and {value!r}"
            )
        if value == IMMUNE_YES:
            scope.update(t.strip() for t in
                         (row.get("immune_tissue_scope") or "").split("|") if t.strip())
    if not is_immune:
        return None
    return {"is_immune": is_immune, "scope": scope}

SCHEMA = pa.schema([
    ("donor_id", pa.string()), ("donor_status", pa.string()),
    ("roi", pa.string()), ("tissue", pa.string()), ("modality", pa.string()),
    ("distance_bin_um", pa.float64()), ("n_cells", pa.int64()),
    ("n_immune", pa.int64()), ("immune_fraction", pa.float64()),
    ("immune_subtype", pa.string()), ("resolved_by", pa.string()),
    ("frame", pa.string()), ("compartment", pa.string()),
    ("state", pa.string()), ("reason", pa.string()),
])

# donor_metadata.csv says "Aab+"; the hand-recorded integration manifest says
# "AAB". One vocabulary reaches the app, and it is the app's.
STATUS_ALIAS = {"ND": "ND", "AAB": "Aab+", "AAB+": "Aab+", "AAB_PLUS": "Aab+",
                "AAB+VE": "Aab+", "T1D": "T1D"}


def refuse_superseded(path: str | Path) -> None:
    """Reject the old input before any existence, schema, or output checks."""
    path = Path(path)
    if any(p.name == "subpanel_calls.parquet" for p in (path, path.resolve())):
        raise ValueError(SUPERSEDED)


def refuse_cross_tissue_distance(roi: str, tissue: str | None = None) -> None:
    """Refuse, by name, to put a non-pancreas section on the distance frame.

    Called on every roi the distance path is about to read. A pLN roi never
    reaches a signed distance, a bin, or a compartment: it raises here instead.
    """
    if norm_roi(roi) == PANCREAS_ROI:
        return
    suffix = f" (tissue {tissue!r})" if tissue else ""
    raise ValueError(CROSS_TISSUE.format(roi=roi, tissue=suffix))


def norm_roi(value: str) -> str:
    return str(value).strip().lower()


def is_pancreas_roi(value: str) -> bool:
    return norm_roi(value) == PANCREAS_ROI


def _roi_keys(root: Path) -> set[str]:
    """roi= partition keys present anywhere under a source root."""
    if not root.is_dir():
        return set()
    return {p.name.removeprefix("roi="): None for p in root.glob("**/roi=*") if p.is_dir()}.keys()


def _partitions(root: Path) -> list[Path]:
    refuse_superseded(root)
    if root.is_file():
        return [root]
    # The distance path is pancreas-only. A roi= partitioned tree (cells_xen)
    # would otherwise glob to nothing and be reported as an absent export.
    for roi in _roi_keys(root):
        refuse_cross_tissue_distance(roi)
    return sorted(root.glob("donor_id=*/*.parquet"))


def _dataset(files: list[Path], required: set[str]) -> ds.Dataset:
    for path in files:
        refuse_superseded(path)
        missing = required - set(pq.ParquetFile(path).schema_arrow.names)
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
    # In particular, never infer a hive donor key from the membership tree.
    return ds.dataset([str(p) for p in files], format="parquet", partitioning=None)


def _refuse_cross_tissue_columns(con: duckdb.DuckDBPyConnection, relation: str) -> None:
    """Refuse a non-pancreas section that arrives as a column, not a directory."""
    names = {row[0] for row in con.execute(f"DESCRIBE {relation}").fetchall()}
    for column in ("roi", "tissue"):
        if column not in names:
            continue
        for (value,) in con.execute(
            f"SELECT DISTINCT {column} FROM {relation} WHERE {column} IS NOT NULL"
        ).fetchall():
            if not str(value).strip():
                continue
            if column == "roi":
                refuse_cross_tissue_distance(value)
            elif str(value).strip().lower() != PANCREAS_TISSUE:
                refuse_cross_tissue_distance(PANCREAS_ROI + "?", value)


def _status_by_donor(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        key = "donor_id" if "donor_id" in fields else "case_id"
        if key not in fields or "donor_status" not in fields:
            raise ValueError(f"{path}: requires case_id (or donor_id) and donor_status")
        statuses = {}
        for row in reader:
            if not row[key] or not row["donor_status"]:
                continue
            donor, status = norm_donor(row[key]), row["donor_status"].strip()
            if donor in statuses and statuses[donor] != status:
                raise ValueError(f"{path}: conflicting donor_status for {donor}")
            statuses[donor] = status
        return statuses


def _manifest_statuses(path: Path, known: dict[str, str]) -> dict[str, str]:
    """Disease status for donors the app's donor_metadata.csv does not carry.

    The four pLN-only donors exist in no other tracked status table. The
    manifest is hand-recorded, so every donor it shares with donor_metadata.csv
    is checked for agreement rather than trusted; a disagreement raises.
    """
    if path is None or not Path(path).is_file():
        return {}
    statuses: dict[str, str] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            donor = (row.get("donor_id") or "").strip()
            raw = (row.get("disease_status") or "").strip()
            if not donor or not raw:
                continue
            donor = norm_donor(donor)
            status = STATUS_ALIAS.get(raw.upper())
            if status is None:
                raise ValueError(f"{path}: unrecognised disease_status {raw!r} for {donor}")
            if statuses.get(donor, status) != status:
                raise ValueError(f"{path}: conflicting disease_status for {donor}")
            statuses[donor] = status
    for donor, status in statuses.items():
        if donor in known and known[donor] != status:
            raise ValueError(
                f"{path}: disease_status {status!r} for {donor} disagrees with "
                f"donor_metadata.csv {known[donor]!r}"
            )
    return statuses


def _row(donor: str, status: str | None, **overrides) -> dict:
    row = {"donor_id": donor, "donor_status": status, "roi": PANCREAS_ROI,
           "tissue": PANCREAS_TISSUE, "modality": "phenocycler",
           "distance_bin_um": None, "n_cells": None, "n_immune": None,
           "immune_fraction": None, "immune_subtype": "Immune", "resolved_by": "all",
           "frame": "signed_distance_um", "compartment": None,
           "state": "unavailable", "reason": ""}
    row.update(overrides)
    return row


def _absence(donor: str, status: str | None, reason: str) -> dict:
    return _row(donor, status, state="unavailable", reason=reason)


def find_pln_release(root: Path) -> Path | None:
    """The newest Xenium release under ``root`` that declares a pLN tissue.

    Releases are siblings: the pancreas release declares roi=panc, and the pLN
    one declares a pLN tissue and roi=pln_*. Picking by the release's own
    declaration beats picking by directory name, and a tree with no pLN release
    returns None so the caller names the absence.
    """
    if not root.is_dir():
        return None
    candidates: list[tuple[float, Path]] = []
    for manifest in sorted(root.glob("*/release.json")):
        try:
            declared = json.loads(manifest.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        tissue = str(declared.get("tissue", "")).strip().lower()
        if RELEASE_TISSUE_ALIAS.get(tissue) != PLN_TISSUE:
            continue
        calls = manifest.parent / str(declared.get("calls_dir") or "calls")
        if calls.is_dir():
            candidates.append((manifest.stat().st_mtime, manifest.parent))
    if not candidates:
        return None
    return max(candidates)[1]


def _release_lineage_column(names: set[str]) -> str | None:
    """The release's common-lineage column, under either spelling.

    The pLN release writes ``common_lineage``; cells_xen writes
    ``lineage_common``. Guessing wrong would read every cell as uncalled, so the
    column is found rather than assumed and a release carrying neither is
    refused by name.
    """
    for candidate in ("common_lineage", XEN_LINEAGE):
        if candidate in names:
            return candidate
    return None


def _pln_composition(
    con: "duckdb.DuckDBPyConnection",
    paths: list[str],
    names: set[str],
    donor: str,
    roi: str,
    status: str | None,
    *,
    tissue: str | None,
    modality: str,
    vocabulary: dict | None,
) -> tuple[list[dict], dict]:
    """Compositional rows for ONE pLN section. Nothing spatial, in any modality.

    Both pLN lanes answer the same two questions in the same order -- immune vs
    non-immune from canonicals, then identity within the immune compartment from
    the subpanels -- so they are counted the same way here. What stays separate
    is the OUTPUT: every row carries its ``modality`` and the two are never
    pooled into one fraction.

    * A cell with no nucleus is excluded from every calculation; we cannot be
      sure it is a cell.
    * The DENOMINATOR is the cells stage 1 actually answered for (``Immune`` or
      ``Non_immune``). ``Ambiguous`` is not an answer: in the denominator it
      would quietly deflate every fraction, in the numerator inflate it. It is
      counted and named on its own row instead, so it stays visible.
    * The NUMERATOR is stage 1's ``Immune``, not the sum of the stage-2
      lineages. An immune cell stage 2 could not resolve is still immune, and
      dropping it would turn the immune fraction into a statement about stage
      2's success rate. The lineage rows therefore do not sum to the total, by
      construction.
    """
    lineage_column = _release_lineage_column(names)
    nucleated = ("nucleus_count IS NOT NULL AND nucleus_count > 0"
                 if "nucleus_count" in names else "TRUE")
    # read_parquet, not a registered Arrow dataset: these files carry
    # large_string columns and the Arrow scan path fails on a partition still
    # being written by the producing lane.
    total, with_nucleus = con.execute(
        f"SELECT count(*), count(*) FILTER (WHERE {nucleated}) FROM read_parquet(?)",
        [paths]).fetchone()
    by_stage1: dict[str, int] = {}
    if "stage1" in names:
        by_stage1 = {
            str(label).strip(): int(count)
            for label, count in con.execute(
                f"SELECT trim(coalesce(stage1, '')), count(*) FROM read_parquet(?) "
                f"WHERE {nucleated} GROUP BY 1", [paths]).fetchall()
        }
    by_lineage: dict[str, int] = {}
    if lineage_column is not None:
        by_lineage = {
            str(label).strip(): int(count)
            for label, count in con.execute(
                f"SELECT trim({lineage_column}), count(*) FROM read_parquet(?) "
                f"WHERE {nucleated} AND nullif(trim({lineage_column}), '') IS NOT NULL "
                f"GROUP BY 1", [paths]).fetchall()
        }

    immune = by_stage1.get(STAGE1_IMMUNE, 0)
    resolved = immune + by_stage1.get(STAGE1_NON_IMMUNE, 0)
    unresolved = with_nucleus - resolved
    audit = {"cells": total, "with_nucleus": with_nucleus, "stage1_resolved": resolved,
             "stage1_immune": immune, "stage1_unresolved": unresolved, "tissue": tissue,
             "donor_status": status, "modality": modality,
             "lineage_counts": dict(sorted(by_lineage.items())),
             "stage1_counts": dict(sorted(by_stage1.items()))}

    base = dict(roi=roi, tissue=tissue or "unrecorded", modality=modality, frame="none")
    prefix = f"{modality}_pln"
    state = reason = None
    if not status:
        state, reason = "unavailable", "donor_status_not_recorded"
    elif not total:
        state, reason = "unavailable", f"{prefix}_export_empty"
    elif not with_nucleus:
        state, reason = "blank", f"{prefix}_no_nucleus:{total}_of_{total}"
    elif lineage_column is None or "stage1" not in names:
        state, reason = "blank", f"{prefix}_missing_call_columns"
    elif not resolved:
        state, reason = "blank", f"{prefix}_stage1_unresolved:{with_nucleus}_of_{with_nucleus}"
    elif vocabulary is None:
        state, reason = "withheld", NO_IMMUNE_VOCAB
    else:
        unrecorded = sorted(lin for lin in by_lineage if lin not in vocabulary["is_immune"])
        if unrecorded:
            state, reason = "withheld", f"{UNRECORDED_LINEAGE}:{'|'.join(unrecorded)}"
        elif not tissue:
            state, reason = "withheld", f"{OUT_OF_SCOPE_TISSUE}:unrecorded"
        elif tissue not in vocabulary["scope"]:
            state, reason = "withheld", f"{OUT_OF_SCOPE_TISSUE}:{tissue}"
    if state is not None:
        return [_row(donor, status, state=state, reason=reason, **base)], audit

    rows = [_row(donor, status, n_cells=resolved, n_immune=immune,
                 immune_fraction=immune / resolved, immune_subtype="Immune",
                 resolved_by="all", state="measured", reason="", **base)]
    for lineage, count in sorted(by_lineage.items()):
        if vocabulary["is_immune"][lineage] != IMMUNE_YES:
            continue
        rows.append(_row(donor, status, n_cells=resolved, n_immune=count,
                         immune_fraction=count / resolved, immune_subtype=lineage,
                         resolved_by="identity_core", state="measured", reason="", **base))
    if unresolved:
        # Named, not absorbed: the cells stage 1 could not answer for.
        rows.append(_row(donor, status, state="blank", immune_subtype="Immune",
                         reason=f"{prefix}_stage1_unresolved:{unresolved}_of_{with_nucleus}",
                         resolved_by="stage1_unresolved", **base))
    return rows, audit


def _pln_xenium_release_arm(
    release_root: Path,
    statuses: dict[str, str],
    donors: set[str] | None,
    vocabulary: dict | None = None,
) -> tuple[list[dict], dict]:
    """Compositional rows from a typed Xenium pLN release. Nothing spatial.

    The release declares its own tissue and calls directory; both are read from
    release.json rather than assumed from the directory name.
    """
    rows: list[dict] = []
    audit: dict[str, dict] = {}
    declared = json.loads((release_root / "release.json").read_text())
    calls = release_root / str(declared.get("calls_dir") or "calls")
    tissue = RELEASE_TISSUE_ALIAS.get(str(declared.get("tissue", "")).strip().lower())

    partitions: dict[tuple[str, str], list[Path]] = {}
    for path in sorted(calls.glob("donor_id=*/roi=*/*.parquet")):
        roi = norm_roi(path.parent.name.removeprefix("roi="))
        if is_pancreas_roi(roi):
            continue
        donor = norm_donor(path.parents[1].name.removeprefix("donor_id="))
        partitions.setdefault((donor, roi), []).append(path)

    # A release may also use the declared A.4 DONOR-LEVEL layout, with shards directly below
    # donor_id=<d>/ and no roi= directory. `Release.calls_for` supports it (phenospace.py:128-130)
    # and donor-level releases exist in this codebase, so a roi-only scan would resolve such a
    # release, record it in `source_pln_by_modality`, and then emit ZERO rows from it -- the card
    # would read "not typed for this donor" while a release full of calls sat on disk. Take the
    # release's own declared roi for those shards, and only when it is not a pancreas one.
    declared_roi = norm_roi(str(declared.get("roi", "") or ""))
    if declared_roi and not is_pancreas_roi(declared_roi):
        for path in sorted(calls.glob("donor_id=*/*.parquet")):
            donor = norm_donor(path.parent.name.removeprefix("donor_id="))
            partitions.setdefault((donor, declared_roi), []).append(path)

    eligible = set(filter_eligible_donors({donor for donor, _ in partitions}))
    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        for (donor, roi), files in sorted(partitions.items()):
            if donor not in eligible:
                continue
            if donors is not None and donor not in donors:
                continue
            section, section_audit = _pln_composition(
                con, [str(f) for f in files], set(pq.read_schema(files[0]).names),
                donor, roi, statuses.get(donor),
                tissue=tissue, modality="xenium", vocabulary=vocabulary)
            rows.extend(section)
            section_audit["release"] = release_root.name
            audit[f"{donor}/{roi}"] = section_audit
    return rows, audit


def _pln_phenocycler_arm(
    calls_root: Path,
    ingest_root: Path,
    statuses: dict[str, str],
    donors: set[str] | None,
    vocabulary: dict | None = None,
) -> tuple[list[dict], dict]:
    """The PhenoCycler pLN arm: the donors that are called, and the ones that are not.

    This modality is typed separately from Xenium and the two are NOT
    interchangeable, so its rows stay distinct by ``modality`` and are never
    pooled into one fraction. A donor typed in one modality and not the other is
    exactly that, in both directions.

    The donor list comes from the pLN INGEST, not from a written-down cohort, so
    a donor whose calls have not landed is named rather than silently missing,
    and the set stays right as the producing lane fills in. Only
    ``donor_id=*/roi=*/*.parquet`` is a call partition -- the lane also writes
    ``_reference/`` tables beside them, and those are not cells.
    """
    audit: dict[str, dict] = {}
    by_donor: dict[tuple[str, str], list[Path]] = {}
    if calls_root.is_dir():
        for path in sorted(calls_root.glob("donor_id=*/roi=*/*.parquet")):
            roi = norm_roi(path.parent.name.removeprefix("roi="))
            if is_pancreas_roi(roi):
                continue
            donor = norm_donor(path.parents[1].name.removeprefix("donor_id="))
            by_donor.setdefault((donor, roi), []).append(path)
    audit["calls"] = {"root": str(calls_root),
                      "state": "present" if by_donor else
                               ("empty" if calls_root.is_dir() else "absent"),
                      "called_donors": sorted({d for d, _ in by_donor})}

    known: list[str] = []
    if ingest_root.is_dir():
        known = sorted({norm_donor(path.name.removeprefix("donor_id="))
                        for path in ingest_root.glob("donor_id=*") if path.is_dir()})
        audit["ingest"] = {"root": str(ingest_root), "donors": known}
    else:
        audit["ingest"] = {"root": str(ingest_root), "state": "absent"}

    rows: list[dict] = []
    called = {donor for donor, _ in by_donor}
    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        for (donor, roi), files in sorted(by_donor.items()):
            if donors is not None and donor not in donors:
                continue
            if donor not in set(filter_eligible_donors([donor])):
                continue
            paths = [str(f) for f in files]
            names = set(pq.read_schema(files[0]).names)
            section, section_audit = _pln_composition(
                con, paths, names, donor, roi, statuses.get(donor),
                tissue=PLN_TISSUE, modality="phenocycler", vocabulary=vocabulary)
            rows.extend(section)
            audit[f"{donor}/{roi}"] = section_audit

    # Sections that exist in the ingest and have no calls yet.
    reason = NO_PLN_PHENO if not called else "phenocycler_pln_calls_absent_for_donor"
    for donor in filter_eligible_donors([d for d in known if d not in called]):
        if donors is not None and donor not in donors:
            continue
        rows.append(_row(donor, statuses.get(donor), roi="pln_1", tissue=PLN_TISSUE,
                         modality="phenocycler", frame="none",
                         state="unavailable", reason=reason))
    return rows, audit


def _pln_arm(
    xen_cells_root: Path,
    statuses: dict[str, str],
    donors: set[str] | None,
    vocabulary: dict | None = None,
) -> tuple[list[dict], dict]:
    """Compositional rows for every non-pancreas section in cells_xen.

    No distance, bin or compartment is produced here, by construction: these
    rows carry ``frame=none`` and a null ``distance_bin_um``. The only quantity
    a lymph node section could contribute to this card is an immune fraction,
    and that needs a lineage call on every cell plus a recorded rule for which
    lineages are immune. Whatever is missing is named, never zeroed.

    The rule is now recorded. Maintainer ruling, 2026-09-17: "All immune cell
    calls from pancreas can be applied to pLN" -- the pancreas immune lineage
    vocabulary is the pLN immune lineage vocabulary, and the crosswalk carries
    it as ``is_immune`` with ``immune_tissue_scope``. That is a ruling about
    identity only: it changes what a lineage label means here, and nothing about
    anatomy or coordinates. These rows still cross no spatial relation, and a
    label row is still never copied from a pancreas cell onto a pLN cell -- the
    two ROIs are segmented separately, so the calls have to be produced on the
    pLN cells themselves.
    """
    rows: list[dict] = []
    audit: dict[str, dict] = {}
    parts: dict[tuple[str, str], list[Path]] = {}
    for path in sorted(xen_cells_root.glob("donor_id=*/roi=*/*.parquet")):
        roi = norm_roi(path.parent.name.removeprefix("roi="))
        if is_pancreas_roi(roi):
            continue  # the pancreas arm is PhenoCycler typing, not Xenium
        donor = norm_donor(path.parents[1].name.removeprefix("donor_id="))
        parts.setdefault((donor, roi), []).append(path)

    eligible = set(filter_eligible_donors({donor for donor, _ in parts}))
    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        for (donor, roi), files in sorted(parts.items()):
            if donor not in eligible:
                continue
            if donors is not None and donor not in donors:
                continue
            con.register("xen", ds.dataset([str(p) for p in files], format="parquet",
                                           partitioning=None))
            names = {row[0] for row in con.execute("DESCRIBE xen").fetchall()}
            tissue = None
            if "tissue" in names:
                values = sorted({
                    str(value).strip()
                    for (value,) in con.execute(
                        "SELECT DISTINCT tissue FROM xen WHERE tissue IS NOT NULL"
                    ).fetchall() if str(value).strip()
                })
                if len(values) > 1:
                    raise ValueError(f"{donor}/{roi}: mixed tissue values {values}")
                tissue = values[0] if values else None
            by_lineage: dict[str, int] = {}
            if XEN_LINEAGE in names:
                n, called = con.execute(
                    f"SELECT count(*), count(*) FILTER "
                    f"(WHERE nullif(trim({XEN_LINEAGE}), '') IS NOT NULL) FROM xen"
                ).fetchone()
                by_lineage = {
                    str(label).strip(): int(count)
                    for label, count in con.execute(
                        f"SELECT trim({XEN_LINEAGE}), count(*) FROM xen "
                        f"WHERE nullif(trim({XEN_LINEAGE}), '') IS NOT NULL GROUP BY 1"
                    ).fetchall()
                }
            else:
                n, called = con.execute("SELECT count(*), 0 FROM xen").fetchone()
            con.unregister("xen")

            status = statuses.get(donor)
            audit[f"{donor}/{roi}"] = {"cells": n, "lineage_called": called,
                                       "tissue": tissue, "donor_status": status,
                                       "lineage_counts": dict(sorted(by_lineage.items()))}
            # The state is the honest one. These cells were not dropped at QC:
            # they are present and were never typed, which ie_state_label calls
            # "Not typed" (blank). Reserve "unavailable" for an absent export and
            # "withheld" for a number this builder declines to compute.
            state = reason = None
            if not status:
                state, reason = "unavailable", "donor_status_not_recorded"
            elif not n:
                state, reason = "unavailable", "xenium_export_empty"
            elif called == 0:
                state, reason = "blank", f"xenium_lineage_not_called:{n}_of_{n}"
            elif called < n:
                state, reason = "blank", f"xenium_lineage_partially_called:{n - called}_of_{n}"
            elif vocabulary is None:
                # The branch that used to fire on every fully-called section. It
                # stays: an unrecorded vocabulary is still withheld, not guessed.
                state, reason = "withheld", NO_IMMUNE_VOCAB
            else:
                classes = vocabulary["is_immune"]
                unrecorded = sorted(lin for lin in by_lineage if lin not in classes)
                section_tissue = (tissue or "").strip()
                unresolvable = sum(count for lin, count in by_lineage.items()
                                   if classes.get(lin) == IMMUNE_UNRESOLVABLE)
                if unrecorded:
                    state = "withheld"
                    reason = f"{UNRECORDED_LINEAGE}:{'|'.join(unrecorded)}"
                elif not section_tissue:
                    state, reason = "withheld", f"{OUT_OF_SCOPE_TISSUE}:unrecorded"
                elif section_tissue not in vocabulary["scope"]:
                    state = "withheld"
                    reason = f"{OUT_OF_SCOPE_TISSUE}:{section_tissue}"
                elif unresolvable:
                    # Neither immune nor non-immune. Counting these either way
                    # would be inventing the numerator or the denominator.
                    state = "blank"
                    reason = f"xenium_lineage_unresolvable:{unresolvable}_of_{n}"
            if state is not None:
                rows.append(_row(
                    donor, status, roi=roi, tissue=tissue or "unrecorded", modality="xenium",
                    frame="none", state=state, reason=reason,
                ))
                continue

            immune_by_lineage = {lin: count for lin, count in sorted(by_lineage.items())
                                 if vocabulary["is_immune"][lin] == IMMUNE_YES}
            n_immune = sum(immune_by_lineage.values())
            base = dict(roi=roi, tissue=tissue, modality="xenium", frame="none",
                        n_cells=n, state="measured", reason="")
            rows.append(_row(donor, status, n_immune=n_immune,
                             immune_fraction=n_immune / n, immune_subtype="Immune",
                             resolved_by="all", **base))
            for lineage, count in immune_by_lineage.items():
                rows.append(_row(donor, status, n_immune=count,
                                 immune_fraction=count / n, immune_subtype=lineage,
                                 resolved_by=XEN_LINEAGE, **base))
    return rows, audit


def build_distance_profiles(
    cells_root: str | Path = TYPING / "cells",
    membership_root: str | Path = TYPING / "islet_membership",
    donor_metadata: str | Path = TYPING / "app" / "donor_metadata.csv",
    *,
    bin_width_um: float = 10.0,
    donors: list[str] | None = None,
    xen_cells_root: str | Path | None = None,
    xen_typing_root: str | Path | None = None,
    pln_pheno_root: str | Path | None = None,
    pln_ingest_root: str | Path | None = None,
    integration_manifest: str | Path | None = None,
    vocab_crosswalk: str | Path | None = VOCAB_CROSSWALK,
) -> pa.Table:
    """Build a small artifact in memory; no source or output file is modified."""
    for path in (cells_root, membership_root, donor_metadata):
        refuse_superseded(path)
    if not math.isfinite(bin_width_um) or bin_width_um <= 0:
        raise ValueError("bin_width_um must be finite and positive")
    cells_root, membership_root, donor_metadata = map(
        Path, (cells_root, membership_root, donor_metadata)
    )
    files_by_donor: dict[str, list[Path]] = {}
    for path in _partitions(cells_root):
        if not path.parent.name.startswith("donor_id="):
            raise ValueError(f"{path}: typing/cells requires an explicit donor_id= partition")
        donor = norm_donor(path.parent.name.removeprefix("donor_id="))
        files_by_donor.setdefault(donor, []).append(path)
    statuses = _status_by_donor(donor_metadata)
    requested = {norm_donor(d) for d in donors} if donors is not None else set(files_by_donor)
    if not requested:
        requested = set(statuses)
    eligible = sorted(filter_eligible_donors(requested))
    if not eligible:
        raise ValueError("unavailable(no_eligible_donors_in_cells_or_metadata)")

    member_files = _partitions(membership_root)
    rows, audit = [], {}
    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        if member_files:
            con.register("members", _dataset(member_files, {
                "object_id", "in_islet", "signed_dist_islet_um",
            }))
            _refuse_cross_tissue_columns(con, "members")
            n, unique, missing = con.execute("""
                SELECT count(*), count(DISTINCT object_id),
                       count(*) FILTER (WHERE object_id IS NULL OR trim(object_id) = '')
                FROM members
            """).fetchone()
            if missing or n != unique:
                raise ValueError("islet_membership has missing or duplicate object_id values")

        for donor in eligible:
            status = statuses.get(donor)
            if not status:
                rows.append(_absence(donor, None, "donor_status_not_recorded"))
                continue
            files = files_by_donor.get(donor, [])
            if not files or not member_files:
                reason = "typing_cells_absent" if not files else "islet_membership_export_absent"
                rows.append(_absence(donor, status, reason))
                continue
            con.register("cells", _dataset(files, {
                "object_id", "cell_class", "excluded_reason", "immune_definition",
                "resolved_by", "qc_analysis_eligible",
            }))
            _refuse_cross_tissue_columns(con, "cells")
            n, unique, missing = con.execute("""
                SELECT count(*), count(DISTINCT object_id),
                       count(*) FILTER (WHERE object_id IS NULL OR trim(object_id) = '')
                FROM cells
            """).fetchone()
            if missing or n != unique:
                raise ValueError(f"{donor}: typing/cells has missing or duplicate object_id values")
            con.execute("""
                CREATE OR REPLACE TEMP TABLE joined AS
                SELECT c.object_id, c.cell_class, c.immune_definition, c.resolved_by,
                       m.object_id AS member_id, m.in_islet, m.signed_dist_islet_um
                FROM cells c LEFT JOIN members m ON c.object_id = m.object_id
                WHERE c.qc_analysis_eligible IS TRUE
                  AND c.cell_class IS DISTINCT FROM 'Excluded'
                  AND coalesce(c.excluded_reason, '') = ''
            """)
            total, unmatched, invalid, untyped, disagree = con.execute("""
                SELECT count(*), count(*) FILTER (WHERE member_id IS NULL),
                       count(*) FILTER (WHERE member_id IS NOT NULL AND
                           (in_islet IS NULL OR signed_dist_islet_um IS NULL
                            OR NOT isfinite(signed_dist_islet_um))),
                       count(*) FILTER (WHERE cell_class IS NULL OR trim(cell_class) = ''
                           OR resolved_by IS NULL OR trim(resolved_by) = ''),
                       count(*) FILTER (WHERE in_islet != (signed_dist_islet_um <= 0))
                FROM joined
            """).fetchone()
            audit[donor] = {"source_cells": n, "eligible_cells": total,
                            "omitted_qc_or_excluded": n - total,
                            "missing_membership": unmatched, "invalid_membership": invalid,
                            "missing_typing": untyped, "membership_distance_disagreement": disagree}
            reason = None
            if not total:
                reason = "no_analysis_eligible_cells"
            elif unmatched:
                reason = f"islet_membership_missing_uuid:{unmatched}_of_{total}"
            elif invalid:
                reason = f"islet_membership_invalid_measurement:{invalid}_of_{total}"
            elif untyped:
                reason = f"typing_or_provenance_missing:{untyped}_of_{total}"
            if reason:
                rows.append(_absence(donor, status, reason))
                continue

            # Membership, rather than the sign of a possibly disagreeing distance,
            # decides the core. Exact zeros retain their own visible bin.
            con.execute("""
                CREATE OR REPLACE TEMP TABLE binned AS
                SELECT *, sign(signed_dist_islet_um) *
                          ceil(abs(signed_dist_islet_um) / ?) * ? AS distance_bin_um,
                       CASE WHEN in_islet THEN 'islet_core'
                            WHEN signed_dist_islet_um > 0 AND signed_dist_islet_um <= 20
                                THEN 'peri_islet'
                            ELSE 'other_pancreas' END AS compartment,
                       CASE WHEN coalesce(trim(immune_definition), '') = ''
                            THEN 'Immune (subtype unresolved)'
                            ELSE immune_definition END AS immune_subtype
                FROM joined
            """, [bin_width_um, bin_width_um])
            table = con.execute("""
                WITH denominators AS (
                    SELECT distance_bin_um, compartment, count(*) AS n_cells,
                           count(*) FILTER (WHERE cell_class = 'Immune') AS n_immune
                    FROM binned GROUP BY ALL
                ), pairs AS (
                    SELECT DISTINCT immune_subtype, resolved_by FROM binned
                    WHERE cell_class = 'Immune'
                ), numerators AS (
                    SELECT distance_bin_um, compartment, immune_subtype, resolved_by,
                           count(*) AS n_immune
                    FROM binned WHERE cell_class = 'Immune' GROUP BY ALL
                ), counts AS (
                    SELECT distance_bin_um, compartment, n_cells, n_immune,
                           'Immune' AS immune_subtype, 'all' AS resolved_by FROM denominators
                    UNION ALL
                    SELECT d.distance_bin_um, d.compartment, d.n_cells,
                           coalesce(n.n_immune, 0), p.immune_subtype, p.resolved_by
                    FROM denominators d CROSS JOIN pairs p LEFT JOIN numerators n
                      ON d.distance_bin_um = n.distance_bin_um AND d.compartment = n.compartment
                     AND p.immune_subtype = n.immune_subtype AND p.resolved_by = n.resolved_by
                )
                SELECT ? AS donor_id, ? AS donor_status, ? AS roi, ? AS tissue,
                       ? AS modality, distance_bin_um, n_cells, n_immune,
                       n_immune::DOUBLE / n_cells AS immune_fraction, immune_subtype,
                       resolved_by, 'signed_distance_um' AS frame, compartment,
                       'measured' AS state, '' AS reason
                FROM counts ORDER BY immune_subtype, resolved_by, compartment, distance_bin_um
            """, [donor, status, PANCREAS_ROI, PANCREAS_TISSUE, "phenocycler"]).to_arrow_table()
            rows.extend(table.to_pylist())

    pln_audit: dict[str, dict] = {}
    pln_sources: dict[str, str] = {}
    build_pln = any(root is not None for root in
                    (xen_cells_root, xen_typing_root, pln_pheno_root))
    if build_pln:
        pln_statuses = dict(statuses)
        pln_statuses.update({
            donor: status
            for donor, status in _manifest_statuses(integration_manifest, statuses).items()
            if donor not in pln_statuses
        })
        vocabulary = (load_immune_vocabulary(vocab_crosswalk)
                      if vocab_crosswalk is not None else None)
        wanted = {norm_donor(d) for d in donors} if donors is not None else None

        # Xenium. A typed release supersedes the untyped ingest: cells_xen at
        # roi=pln_1 carries no call on any row, so reading it can only ever say
        # "not called". The ingest is still read when no release exists, so the
        # arm degrades to a named absence instead of vanishing.
        release = find_pln_release(Path(xen_typing_root)) if xen_typing_root else None
        if release is not None:
            pln_sources["xenium"] = str(release.resolve())
            release_rows, release_audit = _pln_xenium_release_arm(
                release, pln_statuses, wanted, vocabulary)
            rows.extend(release_rows)
            pln_audit.update(release_audit)
        elif xen_cells_root is not None:
            xen_cells_root = Path(xen_cells_root)
            refuse_superseded(xen_cells_root)
            pln_sources["xenium"] = str(xen_cells_root.resolve())
            ingest_rows, ingest_audit = _pln_arm(
                xen_cells_root, pln_statuses, wanted, vocabulary)
            rows.extend(ingest_rows)
            pln_audit.update(ingest_audit)
        elif xen_typing_root is not None:
            pln_sources["xenium"] = NO_XEN_PLN_RELEASE

        # PhenoCycler, kept distinct. The two modalities are not interchangeable
        # and are never pooled into one fraction.
        if pln_pheno_root is not None:
            pheno_rows, pheno_audit = _pln_phenocycler_arm(
                Path(pln_pheno_root),
                Path(pln_ingest_root) if pln_ingest_root else PLN_INGEST,
                pln_statuses, wanted, vocabulary)
            rows.extend(pheno_rows)
            pln_audit["phenocycler"] = pheno_audit
            pln_sources["phenocycler"] = str(Path(pln_pheno_root).resolve())

    metadata = {
        "source_cells": str(cells_root.resolve()),
        "source_membership": str(membership_root.resolve()),
        "source_domain_root": str(donor_metadata.parent.resolve()),
        "source_xen_cells": str(xen_cells_root.resolve()) if xen_cells_root else "",
        "source_pln_by_modality": json.dumps(pln_sources, sort_keys=True),
        "bin_width_um": str(bin_width_um),
        "binning": "ceil(abs(distance)/width)*width with original sign; zero is its own bin",
        "denominator": "QC-eligible, non-excluded cells; repeated across subtype/provenance pairs",
        "source_immune_vocabulary": (str(Path(vocab_crosswalk).resolve())
                                     if vocab_crosswalk is not None else ""),
        "immune_vocabulary": (
            "is_immune per common lineage, read from the crosswalk; never inferred from a "
            "label. Maintainer ruling 2026-09-17: all immune cell calls from pancreas can "
            "be applied to pLN. Identity only -- no spatial relation crosses the tissues."
        ),
        "pln_semantics": (
            "roi != panc is compositional only, joined on donor_id and donor_status. "
            "frame=none, distance_bin_um is null, and no compartment is assigned. "
            "Different block from the pancreas section: no spatial relation crosses it."
        ),
        "audit": json.dumps(audit, sort_keys=True),
        "pln_audit": json.dumps(pln_audit, sort_keys=True),
    }
    table = pa.Table.from_pylist(rows, schema=SCHEMA)
    _assert_no_cross_tissue_distance(table)
    return table.replace_schema_metadata(metadata)


def _assert_no_cross_tissue_distance(table: pa.Table) -> None:
    """Last gate before the artifact exists: no non-pancreas row carries a distance.

    A non-pancreas row may now be ``measured`` -- the 2026-09-17 ruling makes the
    immune lineage vocabulary valid in the pLN, so a compositional immune
    fraction is a number this builder is licensed to produce. What the gate
    refuses is unchanged and is the spatial part: no bin, no compartment, and
    never the signed-distance frame.
    """
    for row in table.to_pylist():
        if is_pancreas_roi(row["roi"]):
            continue
        if (row["distance_bin_um"] is not None or row["compartment"] is not None
                or row["frame"] != "none"):
            raise ValueError(CROSS_TISSUE.format(roi=row["roi"],
                                                 tissue=f" (tissue {row['tissue']!r})"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells-root", "--cells", type=Path, default=TYPING / "cells")
    parser.add_argument("--membership-root", type=Path, default=TYPING / "islet_membership")
    parser.add_argument("--donor-metadata", type=Path, default=TYPING / "app" / "donor_metadata.csv")
    parser.add_argument("--xen-cells-root", type=Path, default=XEN_CELLS,
                        help="cells_xen INGEST, read only when no typed pLN release exists")
    parser.add_argument("--xen-typing-root", type=Path, default=XEN_TYPING,
                        help="tree of Xenium releases; the pLN one supersedes cells_xen")
    parser.add_argument("--pln-pheno-root", type=Path, default=PLN_PHENO_CALLS,
                        help="PhenoCycler pLN calls; absent yields a named unavailable")
    parser.add_argument("--pln-ingest-root", type=Path, default=PLN_INGEST,
                        help="PhenoCycler pLN ingest, used only to name the donors that have a section")
    parser.add_argument("--no-pln", action="store_true",
                        help="Build the pancreas arm alone")
    parser.add_argument("--integration-manifest", type=Path, default=MANIFEST)
    parser.add_argument("--vocab-crosswalk", type=Path, default=VOCAB_CROSSWALK,
                        help="crosswalk recording is_immune per common lineage")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))
    parser.add_argument("--bin-width-um", type=float, default=10.0)
    parser.add_argument("--donors", nargs="+")
    parser.add_argument("--dry-run", action="store_true", help="Measure without writing an artifact")
    args = parser.parse_args(argv)
    table = build_distance_profiles(
        args.cells_root, args.membership_root, args.donor_metadata,
        bin_width_um=args.bin_width_um, donors=args.donors,
        xen_cells_root=None if args.no_pln else args.xen_cells_root,
        xen_typing_root=None if args.no_pln else args.xen_typing_root,
        pln_pheno_root=None if args.no_pln else args.pln_pheno_root,
        pln_ingest_root=None if args.no_pln else args.pln_ingest_root,
        integration_manifest=args.integration_manifest,
        vocab_crosswalk=args.vocab_crosswalk,
    )
    measured = table.filter(pa.compute.equal(table["state"], "measured"))
    totals = measured.filter(pa.compute.equal(measured["resolved_by"], "all"))
    print(f"{len(set(measured['donor_id'].to_pylist()))} measured donors; "
          f"{sum(totals['n_cells'].to_pylist()):,} cells; "
          f"{sum(totals['n_immune'].to_pylist()):,} immune; {table.num_rows:,} artifact rows")
    for roi in sorted(set(table["roi"].to_pylist())):
        subset = table.filter(pa.compute.equal(table["roi"], roi))
        print(f"  roi={roi}: {subset.num_rows:,} rows; "
              f"{len(set(subset['donor_id'].to_pylist()))} donors; "
              f"frames {sorted(set(subset['frame'].to_pylist()))}")
    for row in table.filter(pa.compute.not_equal(table["state"], "measured")).to_pylist():
        print(f"{row['donor_id']}/{row['roi']}: {row['state']}({row['reason']})")
    if not args.dry_run:
        target = args.layers_dir / "distance_profiles" / "data.parquet"
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(".parquet.tmp")
        pq.write_table(table, temporary, compression="zstd")
        temporary.replace(target)
        print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
