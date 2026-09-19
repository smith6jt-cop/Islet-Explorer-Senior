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

The PhenoCycler pancreas subtype is ``immune_definition`` as it stands in
``typing/cells`` -- the per-marker immune definition, on all 20 donors. This
builder does NOT read ``data/development/typing/immune``: that second pass is a
finer subtype (rewritten 2026-09-17, 20 donors, 2,934,898 Immune cells), it
never modifies ``cell_class``, and no app artifact reads it. An Immune cell whose definition is empty is emitted as
``Immune (subtype unresolved)``, which is a named absence of a subtype and not a
zero.

``typing/cells`` here is the 2026-09-15 build: eight classes plus ``Excluded``.
It carries no CD99 question and no ``Unresolved`` class, because no donor has
been re-run with the 2026-09-18 islet chain. ``refined_class`` is not read --
``core/phenocycler/scripts/typing/refine.py`` is written and not wired, and
its 5th-percentile demotion is an open question rather than a decision.

Cell files have no donor column: only their explicit donor_id= partition is
read. Membership files are read with partitioning=None and joined on object_id
alone, including when their optional donor column is absent. Root-level audit
files such as membership_diff.parquet are not cell partitions.

Three arms, and they are not the same kind of measurement
---------------------------------------------------------
* The **pancreas PhenoCycler arm** (``roi=panc``, ``modality=phenocycler``) is
  spatial. It reads ``data/development/typing/cells`` and
  ``typing/islet_membership`` and carries ``frame=signed_distance_um``.
* The **pancreas Xenium arm** (``roi=panc``, ``modality=xenium``) is the same
  tissue in the other modality, and it is spatial ONLY where the four links
  below all hold. See "How a Xenium pancreas cell acquires a signed distance".
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

How a Xenium pancreas cell acquires a signed distance to an islet border
------------------------------------------------------------------------
It does not get one by being in the pancreas. Four links are required, each one
READ off disk rather than assumed, and the arm names the first one that is
missing instead of emitting a zero or a shorter artifact:

1. **A call for the cell.** From the current pancreas release, resolved by
   ``scripts/senior/xenium_release.py`` -- the release's own ``release_id``,
   declared write time and tissue, never a directory mtime. The release is keyed
   by ``cell_id`` and carries the call plus the nucleus columns. Measured on this
   checkout 2026-09-18, its calls schema is ``cell_id, section_key, cluster_id,
   subcluster_id, lineage, common_lineage, reason, ... nucleus_count,
   nucleus_measured, depth ...`` -- **no coordinate column at all**. A release
   alone can never reach a distance axis.
2. **Coordinates for that ``cell_id``.** The Xenium ingest
   ``core/phenocycler/data/integration/cells_xen/donor_id=<d>/roi=panc`` carries
   ``cell_id``, ``x_um`` and ``y_um`` in the **Xenium frame**. The join is on
   ``cell_id`` and on nothing else -- never row order, never a centroid match,
   never a filename. Measured on 6374: 1,012,696 of 1,012,696 release rows found
   their coordinates.
3. **An islet border expressed in the frame those coordinates live in.** The
   islet borders this cohort has are the QuPath ``Islet_N`` annotations, drawn on
   the **PhenoCycler** image. There is no islet outline in the Xenium frame
   anywhere in this checkout: ``integration/structures/xenium/donor_id=*/
   roi=panc/islets.parquet`` exists for two donors and carries ``x_um, y_um,
   n_cells, area_um2, eccentricity`` -- a centroid and an area, no geometry. A
   centroid plus an area is not a border, and turning one into a radius would be
   drawing the border rather than measuring against it.
4. **An acquisition design that licenses a per-cell claim between 2 and 3.**
   Reaching the PhenoCycler border from a Xenium cell means the two images are of
   ONE physical section. That is an acquisition fact, declared with evidence in
   ``core/phenocycler/data/integration/section_modes.csv`` and read here through
   ``phenocycler.integration.modes`` -- the one authority for it. A pair absent
   from that file is ``sequential``: two different physical sections from one
   block. Mapping a cell off one section onto a border drawn on another and
   calling the result that cell's distance to that border is the fabrication
   ``modes.py`` exists to prevent, so this arm refuses it by name.

Measured 2026-09-18 across the 21-donor pancreas release: link 3 fails for every
donor, and link 4 holds for exactly one declared pair (6533/panc/0080870), which
has no registration transform on disk here and whose lane is being developed
elsewhere. So every donor is a NAMED unavailable today and none is a zero. When a
same-slide pair does arrive with its border transported into the Xenium frame,
``panc_xenium_border_loader`` is where it enters and the arm measures without
another change here.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path
from typing import Callable, Sequence

import duckdb
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from scripts.senior import xenium_release
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
#
# WHICH release is no longer answered here. Until 2026-09-18 this file scanned
# ``data/development/xenium_typing`` and took the newest ``release.json`` by
# DIRECTORY MTIME -- a `touch` could reorder it, and the pancreas release could
# not be seen at all because it is written OUTSIDE the repository. Both roots and
# the selection rule now live once, in ``scripts/senior/xenium_release.py``,
# which selects on the release's own declared release_id, write time and tissue.
PLN_PHENO_CALLS = TYPING / "pln" / "pheno"

# Where the pancreas Xenium arm reads the three things the release does not carry.
REGISTRATION = INTEGRATION / "registration"
SECTION_MODES = INTEGRATION / "section_modes.csv"
XEN_STRUCTURES = INTEGRATION / "structures" / "xenium"
# The PhenoCycler pLN ingest. Used only to learn WHICH donors have a pLN section
# when the calls are absent, so the refusal names real donors rather than a
# hardcoded cohort list.
PLN_INGEST = INTEGRATION / "cells_pln"

# Hive partition values are read as strings, always -- every donor id in this
# cohort is all digits, so an inferring reader turns '0112' into the integer 112
# and it stops joining against the string donor_id everywhere else. That rule is
# now stated once, in xenium_release.Release.partitioning(), and the copy that
# used to sit here (defined, never called) is gone rather than left to drift.

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

# --------------------------------------------------------------------------
# The pancreas Xenium arm's named absences. Each one is the FIRST missing link
# of the four in the module docstring, so the reason says what to supply next.
# None of them is ever a zero, and none is ever left off the artifact.
# --------------------------------------------------------------------------
NO_XEN_PANC_RELEASE = "xenium_panc_release_absent"
#: A section exists in the Xenium ingest and the release carries no call for it.
#: Measured 2026-09-18: 6450, 6505, 6534 and 6623 -- four cohort donors whose
#: pancreas was imaged and is not in the 21-donor release. Named here rather than
#: left off, because a donor with no row reads as a donor with no section.
XEN_PANC_NOT_IN_RELEASE = "xenium_panc_release_has_no_calls_for_donor"
#: The release declares the donor but its partition holds no row.
XEN_PANC_EMPTY = "xenium_panc_release_partition_empty"
#: The nucleus was never measured on these cells. A missing measurement is not an
#: absent nucleus: the rule "a cell with no nucleus is excluded" never applied to
#: them, and counting them either way would invent the answer. Measured
#: 2026-09-18: every cell of 6347, 6579 and 6580 in the current release.
XEN_PANC_NUCLEUS_NOT_MEASURED = "xenium_panc_nucleus_not_measured:{missing}_of_{total}"
XEN_PANC_NOT_CALLED = "xenium_panc_lineage_not_called:{uncalled}_of_{total}"
#: Link 2. The release carries no coordinate column, so without the ingest
#: partition there is no position for the cell at all.
XEN_PANC_NO_COORDINATES = "xenium_panc_cell_coordinates_absent"
XEN_PANC_COORDINATES_PARTIAL = "xenium_panc_cell_coordinates_missing_cell_id:{missing}_of_{total}"
#: Link 3/4. The only islet borders this cohort has are drawn on the PhenoCycler
#: section, and this pair is two different physical sections.
XEN_PANC_BORDER_ON_ANOTHER_SECTION = "xenium_panc_islet_border_is_on_another_physical_section:{mode}"
#: Link 4, same-slide branch: the pair is one section imaged twice, but nothing
#: on disk maps between the two frames.
XEN_PANC_NO_REGISTRATION = "xenium_panc_registration_to_the_border_frame_absent"
#: Link 4, same-slide branch: the transport of the PhenoCycler outline into the
#: Xenium frame belongs to the same-slide lane, which owns that geometry. This
#: builder takes it as ``panc_xenium_border_loader`` and never invents one.
XEN_PANC_NO_BORDER_TRANSPORT = "xenium_panc_same_slide_border_transport_not_supplied"
#: The structures lane wrote an islets table for this section but it holds
#: centroids and areas, not outlines. Named rather than treated as no file.
XEN_PANC_STRUCTURES_WITHOUT_OUTLINE = "xenium_panc_islet_structures_carry_no_outline:{columns}"
#: Every spelling the Xenium structures lane might use for an outline. A
#: centroid column is deliberately NOT here: a point is not a border.
XEN_BORDER_COLUMNS = ("wkb", "geometry", "outline_wkb", "border_wkb")

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


def resolve_release(tissue: str, roots: Sequence[Path] | None):
    """The current Xenium release for a tissue, or a named ``Unavailable``.

    One line, because the whole selection belongs to
    ``scripts/senior/xenium_release.py`` and not to five builders: it reads both
    roots -- the repository one and the one outside it that the pancreas release
    is written to -- and selects on the release's OWN declared ``release_id``,
    write time and tissue. Nothing here reads a directory mtime, which is what
    this function replaced: a ``touch`` could reorder the previous scan, and a
    repo-relative scan could not see the pancreas release at all.

    It never returns ``None``. An absence comes back as ``Unavailable`` carrying
    its reason and every candidate it rejected, because a ``None`` that reaches a
    count reads as a zero and a missing measurement is unavailable, not zero.
    """
    return xenium_release.resolve(tissue, roots=roots)


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


def _release_tissue(release) -> str | None:
    """The pLN tissue a release declares, in the crosswalk's vocabulary.

    Read from the release's own ``tissue`` key, falling back to the canonical
    tissue the resolver derived from its declared ``roi``. Never the directory
    name: a release named ``anything_at_all`` that declares a pLN tissue is the
    pLN release, and one named ``release_pln_immune`` that declares pancreas is not.
    """
    declared = RELEASE_TISSUE_ALIAS.get(str(release.metadata.get("tissue", "")).strip().lower())
    if declared is not None:
        return declared
    return PLN_TISSUE if release.tissue == xenium_release.PLN else None


def _pln_xenium_release_arm(
    release,
    statuses: dict[str, str],
    donors: set[str] | None,
    vocabulary: dict | None = None,
) -> tuple[list[dict], dict]:
    """Compositional rows from a typed Xenium pLN release. Nothing spatial.

    ``release`` is a ``xenium_release.Release``: it has already been selected on
    its own declared release_id, write time and tissue, and it carries its calls
    directory and the hive keys read as STRINGS. A donor id that becomes an
    integer stops matching the donor it names, which is why nothing here infers a
    partition type.
    """
    rows: list[dict] = []
    audit: dict[str, dict] = {}
    tissue = _release_tissue(release)

    # `partition_files` handles BOTH declared layouts in one place: the
    # donor_id=/roi= tree, and the A.4 DONOR-LEVEL one with shards directly below
    # donor_id= and no roi directory. A roi-only scan would resolve a donor-level
    # release, record it as the source, and then emit ZERO rows from it -- the
    # card would read "not typed for this donor" while a release full of calls
    # sat on disk. A shard with no roi segment takes the release's own declared
    # roi, and only when that is not a pancreas one.
    declared_roi = next((norm_roi(value) for value in release.declared_roi
                         if value and not is_pancreas_roi(value)), "")
    partitions: dict[tuple[str, str], list[Path]] = {}
    for donor, roi, path in release.partition_files():
        roi = norm_roi(roi) if roi else declared_roi
        if not roi or is_pancreas_roi(roi):
            continue
        partitions.setdefault((norm_donor(donor), roi), []).append(path)

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
            section_audit["release"] = release.release_id
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


# --------------------------------------------------------------------------
# The pancreas Xenium arm
# --------------------------------------------------------------------------
def _xenium_serial(con: "duckdb.DuckDBPyConnection", paths: list[str],
                   names: set[str]) -> str:
    """The Xenium serial this section was acquired on, read off the ingest.

    ``section_modes.csv`` is keyed on donor + roi + SERIAL, because one donor+roi
    can hold both designs -- 6533/panc is a serial pair on 0041393 and a
    same-slide pair on 0080870. Asking without a serial gets the wildcard, which
    is the fail-safe (``sequential``) answer and would silently deny a same-slide
    pair its own row. So the serial is read from the ingest's own
    ``xenium_sample_key`` (measured: ``0041032__panc``), never parsed from a
    directory name and never guessed.
    """
    if "xenium_sample_key" not in names:
        return ""
    values = sorted({
        str(value).strip()
        for (value,) in con.execute(
            "SELECT DISTINCT xenium_sample_key FROM read_parquet(?) "
            "WHERE xenium_sample_key IS NOT NULL", [paths]).fetchall()
        if str(value).strip()
    })
    if len(values) != 1:
        return ""
    return values[0].split("__", 1)[0]


def _donor_partition(root: Path, donor: str, roi: str = PANCREAS_ROI) -> Path | None:
    """The ``donor_id=<d>/roi=<roi>`` directory for ``donor``, matched on the id.

    The hive segment is compared through ``norm_donor``, not as text, because
    ``0112`` and ``112`` are one donor and a tree that spells it the other way
    would otherwise read as an absent partition -- a named absence standing in
    for a donor whose data is on disk.
    """
    root = Path(root)
    if not root.is_dir():
        return None
    for path in sorted(root.glob("donor_id=*")):
        if norm_donor(path.name.removeprefix("donor_id=")) != donor:
            continue
        section = path / f"roi={roi}"
        if section.is_dir():
            return section
    return None


def _section_mode(donor: str, serial: str, section_modes_csv: Path) -> str:
    """The declared acquisition design for this pair, from the one authority.

    ``phenocycler.integration.modes`` owns this question; restating its rules
    here would be a second copy of a manifest that exists precisely so a flag
    cannot invent an acquisition fact. The global fallback is pinned to
    ``sequential`` -- the fail-safe direction -- because this builder has no
    config.ini and must never be the thing that promotes a pair to same_slide.
    """
    from types import SimpleNamespace

    from phenocycler.integration import modes  # imported here: only this arm needs it

    cfg = SimpleNamespace(section_modes_csv=Path(section_modes_csv),
                          integration_mode="sequential")
    return modes.mode_for(cfg, donor, PANCREAS_ROI, serial or None)


def xenium_frame_islet_border(
    donor: str,
    *,
    serial: str,
    structures_root: Path,
    registration_root: Path,
    section_modes_csv: Path,
    border_loader: Callable[..., object] | None = None,
) -> tuple[list | None, str, dict]:
    """(border geometries in the XENIUM micron frame, reason, audit).

    Link 3 of the four in the module docstring, and link 4 where it is needed.
    Two routes, and no third:

    * **The Xenium section's own outlines.** ``structures/xenium/donor_id=<d>/
      roi=panc/islets.parquet`` is where a structure segmented on the Xenium
      section lands. Measured 2026-09-18 it exists for two donors and carries
      ``x_um, y_um, n_cells, area_um2, eccentricity`` -- a centroid and an area.
      That is named, not read as an absent file, and it is not turned into a
      circle: an equivalent radius would be a border this builder drew itself.
    * **A same-slide pair's PhenoCycler outline, transported.** One physical
      section imaged twice means the QuPath ``Islet_N`` border IS a border on the
      Xenium tissue, and the registration transform expresses it in the Xenium
      frame. Both prerequisites are checked here -- the ``same_slide``
      declaration and the transform on disk -- and the transport itself is the
      same-slide lane's geometry, taken as ``border_loader`` rather than
      reimplemented. A sequential pair never reaches this route.

    Returns ``(None, reason, audit)`` when no border exists, so the caller emits
    a named unavailable. It never returns an empty border list as if it were one.
    """
    audit: dict = {}
    section = _donor_partition(structures_root, donor)
    structures = (section / "islets.parquet") if section is not None else None
    audit["xenium_structures"] = str(structures) if structures else ""
    if structures is not None and structures.is_file():
        columns = list(pq.read_schema(structures).names)
        audit["xenium_structures_columns"] = columns
        outline = next((c for c in XEN_BORDER_COLUMNS if c in columns), None)
        if outline is not None:
            from shapely import wkb as shapely_wkb  # only when a border is really there

            table = pq.read_table(structures, columns=[outline])
            border = [shapely_wkb.loads(bytes(value)) for value in table[outline].to_pylist()
                      if value is not None]
            audit["border_source"] = f"{structures}:{outline}"
            audit["n_border_outlines"] = len(border)
            if border:
                return border, "", audit
        else:
            audit["xenium_structures_state"] = "no outline column"

    mode = _section_mode(donor, serial, section_modes_csv)
    audit["section_mode"] = mode
    audit["xenium_serial"] = serial
    if mode != "same_slide":
        reason = XEN_PANC_BORDER_ON_ANOTHER_SECTION.format(mode=mode)
        if "xenium_structures_state" in audit:
            reason = XEN_PANC_STRUCTURES_WITHOUT_OUTLINE.format(
                columns="|".join(audit["xenium_structures_columns"]))
        return None, reason, audit

    registered = _donor_partition(registration_root, donor)
    transform = (registered / "transform.json") if registered is not None else None
    audit["transform"] = str(transform) if transform else ""
    if transform is None or not transform.is_file():
        return None, XEN_PANC_NO_REGISTRATION, audit
    if border_loader is None:
        return None, XEN_PANC_NO_BORDER_TRANSPORT, audit
    border = list(border_loader(donor=donor, transform=transform, serial=serial) or [])
    audit["border_source"] = "same_slide_transport"
    audit["n_border_outlines"] = len(border)
    if not border:
        return None, XEN_PANC_NO_BORDER_TRANSPORT, audit
    return border, "", audit


def signed_distance_to_border_um(border: list, x, y):
    """Signed distance to the nearest islet border: negative inside, positive out.

    Containment decides the sign, the distance to the BOUNDARY gives the
    magnitude, and a point exactly on the border keeps its own zero. This is the
    same decision order ``typing/islet_membership.py`` uses on the PhenoCycler
    side -- membership from the outline, distance as the magnitude -- so the two
    modalities mean the same thing by "inside".
    """
    import numpy as np
    from shapely import STRtree, points

    coordinates = np.column_stack([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    pts = points(coordinates)
    boundaries = [geometry.boundary for geometry in border]
    nearest = STRtree(boundaries).nearest(pts)
    magnitude = np.array([boundaries[int(j)].distance(p) for p, j in zip(pts, nearest)],
                         dtype=float)
    inside = np.zeros(len(pts), dtype=bool)
    point_index, _ = STRtree(border).query(pts, predicate="within")
    inside[point_index] = True
    return np.where(inside, -magnitude, magnitude), inside


def _panc_xenium_measure(
    con: "duckdb.DuckDBPyConnection",
    call_paths: list[str],
    coordinate_paths: list[str],
    lineage_column: str,
    nucleated: str,
) -> tuple[list[str], list[str], list[float], list[float], int, int]:
    """(cell_id, lineage, x, y, release rows kept, rows with no coordinate).

    The join is ``cell_id`` to ``cell_id`` and nothing else. Row order, file
    order and a centroid match are all excluded by construction: an inner join on
    the key either finds the cell or does not, and the ones it does not find are
    counted and named rather than dropped.
    """
    kept, missing = con.execute(
        f"SELECT count(*), count(*) FILTER (WHERE x.cell_id IS NULL) "
        f"FROM read_parquet(?) r LEFT JOIN read_parquet(?) x ON r.cell_id = x.cell_id "
        f"WHERE {nucleated}", [call_paths, coordinate_paths]).fetchone()
    rows = con.execute(
        f"SELECT r.cell_id, trim(coalesce(r.{lineage_column}, '')) AS lineage, "
        f"       x.x_um, x.y_um "
        f"FROM read_parquet(?) r JOIN read_parquet(?) x ON r.cell_id = x.cell_id "
        f"WHERE {nucleated}", [call_paths, coordinate_paths]).fetchall()
    return ([r[0] for r in rows], [r[1] for r in rows],
            [r[2] for r in rows], [r[3] for r in rows], int(kept), int(missing))


def _panc_xenium_arm(
    release,
    xen_cells_root: Path,
    statuses: dict[str, str],
    donors: set[str] | None,
    vocabulary: dict | None,
    *,
    bin_width_um: float,
    structures_root: Path,
    registration_root: Path,
    section_modes_csv: Path,
    border_loader: Callable[..., object] | None = None,
) -> tuple[list[dict], dict]:
    """The pancreas in the OTHER modality: spatial where the four links hold.

    Every donor the release carries gets a row, whatever state it is in. A donor
    the arm cannot place on the distance axis is a named unavailable carrying the
    first missing link; it is never a zero, never omitted, and never quietly
    given the PhenoCycler donor's number.

    These rows are ``roi=panc`` because they ARE the pancreas. Keeping them there
    is what leaves ``refuse_cross_tissue_distance`` meaning what it says -- a
    tissue boundary, not a modality one. They carry ``modality=xenium`` and are
    never pooled with the PhenoCycler arm: the two sections are segmented
    separately, typed by different lanes against different panels, and their
    denominators are not interchangeable.
    """
    rows: list[dict] = []
    audit: dict[str, dict] = {}
    tallies = xenium_release.donor_tallies(release)
    partitions = release.donor_partitions()
    eligible = set(filter_eligible_donors({norm_donor(d) for d in partitions}))
    label_column = release.require_column("common_lineage")

    def absent(donor, status, reason, **extra):
        return _row(donor, status, roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE,
                    modality="xenium", frame="none", state="unavailable",
                    reason=reason, **extra)

    with duckdb.connect(config={"threads": 2, "memory_limit": "2GB"}) as con:
        for raw_donor, files in sorted(partitions.items()):
            donor = norm_donor(raw_donor)
            if donor not in eligible:
                continue
            if donors is not None and donor not in donors:
                continue
            tally = tallies[raw_donor]
            status = statuses.get(donor)
            entry: dict = {"release": release.release_id, "donor_status": status,
                           "cells": tally.n_rows, "called": tally.n_called,
                           "uncalled": tally.n_uncalled,
                           "label_unavailable": tally.n_label_unavailable,
                           "no_nucleus_excluded": tally.n_no_nucleus_excluded,
                           "nucleus_not_measured": tally.n_nucleus_not_measured}
            audit[f"{donor}/{PANCREAS_ROI}"] = entry

            if not status:
                rows.append(absent(donor, None, "donor_status_not_recorded"))
                continue
            if not tally.n_rows:
                rows.append(absent(donor, status, XEN_PANC_EMPTY))
                continue
            # The nucleus rule and a missing nucleus measurement are two different
            # things and are never collapsed. A section where the nucleus was
            # never measured cannot have the rule applied to it at all.
            if tally.nucleus_rule_applied and tally.n_nucleus_not_measured >= tally.n_rows:
                rows.append(_row(donor, status, roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE,
                                 modality="xenium", frame="none", state="blank",
                                 reason=XEN_PANC_NUCLEUS_NOT_MEASURED.format(
                                     missing=tally.n_nucleus_not_measured, total=tally.n_rows)))
                continue
            if not tally.n_called:
                rows.append(_row(donor, status, roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE,
                                 modality="xenium", frame="none", state="blank",
                                 reason=XEN_PANC_NOT_CALLED.format(
                                     uncalled=tally.n_rows, total=tally.n_rows)))
                continue

            section = _donor_partition(xen_cells_root, donor)
            coordinates = sorted(section.glob("*.parquet")) if section is not None else []
            entry["coordinates"] = [str(p) for p in coordinates]
            if not coordinates:
                rows.append(absent(donor, status, XEN_PANC_NO_COORDINATES))
                continue
            coordinate_names = set(pq.read_schema(coordinates[0]).names)
            if not {"cell_id", "x_um", "y_um"} <= coordinate_names:
                rows.append(absent(donor, status, XEN_PANC_NO_COORDINATES))
                continue
            coordinate_paths = [str(p) for p in coordinates]
            serial = _xenium_serial(con, coordinate_paths, coordinate_names)

            border, reason, border_audit = xenium_frame_islet_border(
                donor, serial=serial, structures_root=structures_root,
                registration_root=registration_root,
                section_modes_csv=section_modes_csv, border_loader=border_loader)
            entry.update(border_audit)
            if border is None:
                rows.append(absent(donor, status, reason))
                continue
            if vocabulary is None:
                rows.append(_row(donor, status, roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE,
                                 modality="xenium", frame="none", state="withheld",
                                 reason=NO_IMMUNE_VOCAB))
                continue
            if PANCREAS_TISSUE not in vocabulary["scope"]:
                rows.append(_row(donor, status, roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE,
                                 modality="xenium", frame="none", state="withheld",
                                 reason=f"{OUT_OF_SCOPE_TISSUE}:{PANCREAS_TISSUE}"))
                continue

            call_paths = [str(p) for p in files]
            nucleus_count = release.column("nucleus_count")
            nucleus_measured = release.column("nucleus_measured")
            if nucleus_count is None:
                nucleated = "TRUE"
            elif nucleus_measured is None:
                nucleated = f"r.{nucleus_count} IS NOT NULL AND r.{nucleus_count} > 0"
            else:
                nucleated = (f"coalesce(r.{nucleus_measured}, FALSE) "
                             f"AND r.{nucleus_count} IS NOT NULL AND r.{nucleus_count} > 0")
            _ids, lineages, xs, ys, kept, missing = _panc_xenium_measure(
                con, call_paths, coordinate_paths, label_column, nucleated)
            entry["nucleated_release_rows"] = kept
            entry["coordinates_missing"] = missing
            if missing:
                rows.append(absent(donor, status, XEN_PANC_COORDINATES_PARTIAL.format(
                    missing=missing, total=kept)))
                continue
            if not xs:
                rows.append(absent(donor, status, XEN_PANC_EMPTY))
                continue
            unrecorded = sorted({lin for lin in lineages
                                 if lin and lin not in vocabulary["is_immune"]})
            if unrecorded:
                rows.append(_row(donor, status, roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE,
                                 modality="xenium", frame="none", state="withheld",
                                 reason=f"{UNRECORDED_LINEAGE}:{'|'.join(unrecorded)}"))
                continue
            rows.extend(_panc_xenium_rows(donor, status, lineages, xs, ys, border,
                                          vocabulary, bin_width_um, entry))

    # Sections the ingest has and the release does not call. The donor list comes
    # from the data on both sides, never from a written-down cohort, so a section
    # the release has not reached yet is named instead of silently absent -- and
    # the set stays right as the typing lane fills in.
    called = {norm_donor(d) for d in partitions}
    ingested = sorted({norm_donor(path.name.removeprefix("donor_id="))
                       for path in Path(xen_cells_root).glob("donor_id=*")
                       if (path / f"roi={PANCREAS_ROI}").is_dir()})
    for donor in filter_eligible_donors([d for d in ingested if d not in called]):
        if donors is not None and donor not in donors:
            continue
        audit[f"{donor}/{PANCREAS_ROI}"] = {"release": release.release_id,
                                            "state": XEN_PANC_NOT_IN_RELEASE}
        rows.append(_row(donor, statuses.get(donor), roi=PANCREAS_ROI,
                         tissue=PANCREAS_TISSUE, modality="xenium", frame="none",
                         state="unavailable", reason=XEN_PANC_NOT_IN_RELEASE))
    return rows, audit


def _panc_xenium_rows(donor, status, lineages, xs, ys, border, vocabulary,
                      bin_width_um: float, entry: dict) -> list[dict]:
    """One donor's measured spatial rows, on the same axis and cut as the protein arm."""
    import math as _math

    distance, inside = signed_distance_to_border_um(border, xs, ys)
    entry["measured_cells"] = len(lineages)
    entry["frame"] = "signed_distance_um"

    denominators: dict[tuple[float, str], int] = {}
    numerators: dict[tuple[float, str, str], int] = {}
    totals: dict[tuple[float, str], int] = {}
    for lineage, signed, is_inside in zip(lineages, distance.tolist(), inside.tolist()):
        # The protein arm's SQL: sign(d) * ceil(|d|/w) * w. Exact zero keeps its
        # own bin and does not acquire a sign from which side it was reached on.
        edge = (0.0 if signed == 0 else
                _math.copysign(_math.ceil(abs(signed) / bin_width_um) * bin_width_um, signed))
        if is_inside:
            compartment = "islet_core"
        elif 0 < signed <= 20:
            compartment = "peri_islet"
        else:
            compartment = "other_pancreas"
        key = (edge, compartment)
        denominators[key] = denominators.get(key, 0) + 1
        if lineage and vocabulary["is_immune"].get(lineage) == IMMUNE_YES:
            totals[key] = totals.get(key, 0) + 1
            numerators[(edge, compartment, lineage)] = \
                numerators.get((edge, compartment, lineage), 0) + 1

    subtypes = sorted({lineage for (_edge, _compartment, lineage) in numerators})
    rows: list[dict] = []
    for (edge, compartment), n_cells in sorted(denominators.items()):
        base = dict(roi=PANCREAS_ROI, tissue=PANCREAS_TISSUE, modality="xenium",
                    frame="signed_distance_um", distance_bin_um=edge,
                    compartment=compartment, n_cells=n_cells, state="measured", reason="")
        n_immune = totals.get((edge, compartment), 0)
        rows.append(_row(donor, status, n_immune=n_immune,
                         immune_fraction=n_immune / n_cells,
                         immune_subtype="Immune", resolved_by="all", **base))
        for lineage in subtypes:
            count = numerators.get((edge, compartment, lineage), 0)
            rows.append(_row(donor, status, n_immune=count, immune_fraction=count / n_cells,
                             immune_subtype=lineage, resolved_by="common_lineage", **base))
    return rows


def build_distance_profiles(
    cells_root: str | Path = TYPING / "cells",
    membership_root: str | Path = TYPING / "islet_membership",
    donor_metadata: str | Path = TYPING / "app" / "donor_metadata.csv",
    *,
    bin_width_um: float = 10.0,
    donors: list[str] | None = None,
    xen_cells_root: str | Path | None = None,
    xenium_roots: Sequence[str | Path] | None = None,
    pln_pheno_root: str | Path | None = None,
    pln_ingest_root: str | Path | None = None,
    integration_manifest: str | Path | None = None,
    vocab_crosswalk: str | Path | None = VOCAB_CROSSWALK,
    pln: bool = True,
    panc_xenium: bool = True,
    panc_xenium_border_loader: Callable[..., object] | None = None,
    registration_root: str | Path = REGISTRATION,
    section_modes_csv: str | Path = SECTION_MODES,
    xenium_structures_root: str | Path = XEN_STRUCTURES,
) -> pa.Table:
    """Build a small artifact in memory; no source or output file is modified.

    ``xenium_roots`` is the list of roots the Xenium releases are searched under,
    and passing ``None`` builds no Xenium arm at all. That is deliberate: the
    resolver's own default reaches ``/data/xenium_typing``, so a caller that says
    nothing would read real releases. The CLI names the default roots explicitly;
    a caller handing in its own tree gets that tree and nothing else.
    """
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
    panc_xenium_audit: dict[str, dict] = {}
    panc_xenium_source = ""
    if xen_cells_root is not None:
        xen_cells_root = Path(xen_cells_root)
        refuse_superseded(xen_cells_root)
    xenium_roots = ([Path(root) for root in xenium_roots]
                    if xenium_roots is not None else None)
    build_pln = pln and any(root is not None for root in
                            (xen_cells_root, xenium_roots, pln_pheno_root))
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
        release = resolve_release(xenium_release.PLN, xenium_roots) if xenium_roots else None
        if release is not None and release.available:
            pln_sources["xenium"] = str(release.root.resolve())
            release_rows, release_audit = _pln_xenium_release_arm(
                release, pln_statuses, wanted, vocabulary)
            rows.extend(release_rows)
            pln_audit.update(release_audit)
            pln_audit["xenium_release"] = release.describe()
        elif xen_cells_root is not None:
            pln_sources["xenium"] = str(xen_cells_root.resolve())
            ingest_rows, ingest_audit = _pln_arm(
                xen_cells_root, pln_statuses, wanted, vocabulary)
            rows.extend(ingest_rows)
            pln_audit.update(ingest_audit)
            if release is not None:
                pln_audit["xenium_release"] = release.describe()
        elif release is not None:
            pln_sources["xenium"] = NO_XEN_PLN_RELEASE
            pln_audit["xenium_release"] = release.describe()

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

    # The pancreas in the other modality. It is its own arm because it is the
    # same tissue and a different measurement, and it runs whether or not the pLN
    # arm produced anything: a donor with no pLN section still has a pancreas.
    if panc_xenium and xenium_roots is not None:
        panc_release = resolve_release(xenium_release.PANC, xenium_roots)
        panc_xenium_audit["release"] = panc_release.describe()
        if not panc_release.available:
            panc_xenium_source = NO_XEN_PANC_RELEASE
        elif xen_cells_root is None:
            # Link 2 for the whole arm, not one donor: without the ingest there
            # is no coordinate anywhere, and the release carries none.
            panc_xenium_source = XEN_PANC_NO_COORDINATES
        else:
            panc_xenium_source = str(panc_release.root.resolve())
            panc_statuses = dict(statuses)
            panc_statuses.update({
                donor: status
                for donor, status in _manifest_statuses(integration_manifest, statuses).items()
                if donor not in panc_statuses
            })
            panc_rows, donor_audit = _panc_xenium_arm(
                panc_release, xen_cells_root, panc_statuses,
                {norm_donor(d) for d in donors} if donors is not None else None,
                (load_immune_vocabulary(vocab_crosswalk)
                 if vocab_crosswalk is not None else None),
                bin_width_um=bin_width_um,
                structures_root=Path(xenium_structures_root),
                registration_root=Path(registration_root),
                section_modes_csv=Path(section_modes_csv),
                border_loader=panc_xenium_border_loader)
            rows.extend(panc_rows)
            panc_xenium_audit.update(donor_audit)

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
        "source_xenium_panc": panc_xenium_source,
        "panc_xenium_semantics": (
            "roi=panc, modality=xenium is the pancreas in the RNA modality. It is "
            "spatial only where the section's own islet border exists in the frame "
            "its coordinates live in; a serial pair's PhenoCycler border is on "
            "another physical section and is refused by name, never mapped across. "
            "Its counts are never pooled with the PhenoCycler pancreas arm: the two "
            "sections are segmented separately and their denominators differ."
        ),
        "audit": json.dumps(audit, sort_keys=True),
        "pln_audit": json.dumps(pln_audit, sort_keys=True),
        "panc_xenium_audit": json.dumps(panc_xenium_audit, sort_keys=True, default=str),
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
                        help="cells_xen INGEST: the pLN fallback when no typed pLN release "
                             "exists, and the only source of a pancreas cell's coordinates")
    parser.add_argument("--xenium-root", type=Path, action="append", dest="xenium_roots",
                        help="repeatable root to search for Xenium releases; the default is "
                             "both roots from scripts/senior/xenium_release.py, which selects "
                             "on each release's own declared id, write time and tissue")
    parser.add_argument("--no-panc-xenium", action="store_true",
                        help="Drop the pancreas Xenium arm. Its rows say per donor why the "
                             "release cannot reach the distance axis, so dropping it removes "
                             "named absences rather than adding anything")
    parser.add_argument("--pln-pheno-root", type=Path, default=PLN_PHENO_CALLS,
                        help="PhenoCycler pLN calls; absent yields a named unavailable")
    parser.add_argument("--pln-ingest-root", type=Path, default=PLN_INGEST,
                        help="PhenoCycler pLN ingest, used only to name the donors that have a section")
    parser.add_argument("--no-pln", action="store_true",
                        help="Drop both pLN arms. The pancreas arms, in both modalities, "
                             "are unaffected")
    parser.add_argument("--integration-manifest", type=Path, default=MANIFEST)
    parser.add_argument("--registration-root", type=Path, default=REGISTRATION,
                        help="per-pair registration transforms, read only for a same-slide pair")
    parser.add_argument("--section-modes", type=Path, default=SECTION_MODES,
                        help="the acquisition-design manifest; a pair absent from it is serial")
    parser.add_argument("--xenium-structures-root", type=Path, default=XEN_STRUCTURES,
                        help="structures segmented on the Xenium section, where an islet "
                             "outline in the Xenium frame would live")
    parser.add_argument("--vocab-crosswalk", type=Path, default=VOCAB_CROSSWALK,
                        help="crosswalk recording is_immune per common lineage")
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))
    parser.add_argument("--bin-width-um", type=float, default=10.0)
    parser.add_argument("--donors", nargs="+")
    parser.add_argument("--dry-run", action="store_true", help="Measure without writing an artifact")
    args = parser.parse_args(argv)
    # Named explicitly rather than left to the resolver's own default, so the
    # roots this build read are in the artifact's provenance and a caller that
    # hands in its own tree gets that tree and nothing else.
    roots = args.xenium_roots or list(xenium_release.default_roots(REPO))
    table = build_distance_profiles(
        args.cells_root, args.membership_root, args.donor_metadata,
        bin_width_um=args.bin_width_um, donors=args.donors,
        xen_cells_root=args.xen_cells_root,
        xenium_roots=roots,
        pln_pheno_root=None if args.no_pln else args.pln_pheno_root,
        pln_ingest_root=None if args.no_pln else args.pln_ingest_root,
        integration_manifest=args.integration_manifest,
        vocab_crosswalk=args.vocab_crosswalk,
        pln=not args.no_pln,
        panc_xenium=not args.no_panc_xenium,
        registration_root=args.registration_root,
        section_modes_csv=args.section_modes,
        xenium_structures_root=args.xenium_structures_root,
    )
    measured = table.filter(pa.compute.equal(table["state"], "measured"))
    totals = measured.filter(pa.compute.equal(measured["resolved_by"], "all"))
    print(f"{len(set(measured['donor_id'].to_pylist()))} measured donors; "
          f"{sum(totals['n_cells'].to_pylist()):,} cells; "
          f"{sum(totals['n_immune'].to_pylist()):,} immune; {table.num_rows:,} artifact rows")
    # Broken out by roi AND modality. The pancreas now carries two, and they are
    # two measurements of one tissue, not one measurement: a summary that added
    # them would be the first place the two got pooled.
    keys = sorted(set(zip(table["roi"].to_pylist(), table["modality"].to_pylist())))
    for roi, modality in keys:
        subset = table.filter(pa.compute.and_(pa.compute.equal(table["roi"], roi),
                                              pa.compute.equal(table["modality"], modality)))
        print(f"  roi={roi} modality={modality}: {subset.num_rows:,} rows; "
              f"{len(set(subset['donor_id'].to_pylist()))} donors; "
              f"frames {sorted(set(subset['frame'].to_pylist()))}")
    for row in table.filter(pa.compute.not_equal(table["state"], "measured")).to_pylist():
        print(f"{row['donor_id']}/{row['roi']}/{row['modality']}: "
              f"{row['state']}({row['reason']})")
    # The reader contract, said out loud, because the pancreas arm now carries
    # two modalities and a reader that filters it by donor alone will read one
    # modality's absence as the other modality's.
    panc = table.filter(pa.compute.equal(table["roi"], PANCREAS_ROI))
    modalities = sorted(set(panc["modality"].to_pylist()))
    if len(modalities) > 1:
        print(f"NOTE: roi=panc carries {modalities}. A reader of this arm must select a "
              f"modality before it selects a donor; the two are separate measurements of "
              f"one tissue and are never pooled or substituted for one another.")
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
