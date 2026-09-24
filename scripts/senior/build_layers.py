#!/usr/bin/env python
"""One driver for the layer artifacts, and the staleness gate that was missing.

Nine layer builders -- eight under ``scripts/senior/`` and ``pln``, which is
``core/phenocycler/scripts/pln/build_pln_layer.py`` -- and the ``typing/app`` artifact builders
under ``core/phenocycler/scripts/typing/`` were each invoked by hand, in an order that existed
only in prose.  Every artifact therefore sat at whatever donor set and whatever input generation it
happened to be built against, and nothing compared the two.  This file declares the graph once --
``(name, command, inputs, outputs)`` per node -- so the comparison can be made by a command rather
than by remembering.

Three modes, and only the third writes anything::

    .venv/bin/python scripts/senior/build_layers.py --check
    .venv/bin/python scripts/senior/build_layers.py --plan
    .venv/bin/python scripts/senior/build_layers.py --run [--only NODE ...]

``--check`` compares each node's OLDEST output against its NEWEST input and exits non-zero if any
node is stale, has a missing output or has a missing input.  A node whose outputs were left by a
run that did not finish is stale too (AN INTERRUPTED RUN, below).  That is the gate.  It reads
mtimes and nothing else; it opens no parquet and writes no file.

``--plan`` prints the topologically ordered commands and runs nothing.

``--run`` executes in dependency order, streaming each builder's output, and REFUSES to start a
node whose own input is the output of a node that is stale and is not being rebuilt first.  A
half-fresh artifact built on a stale input is the failure this driver exists to stop, so it is
refused rather than reported afterwards.

It also SKIPS any node that is already fresh.  Until 2026-09-21 it did not: it rebuilt every
scheduled node, including the ones ``--check`` had reported fresh in the same breath -- the
staleness was computed, printed, and then thrown away.  On this cohort that is 2 nodes of 26,
but the cost is not the point; a driver that cannot do less work when less has changed is a
driver that makes every change cost a full rebuild.  ``--only`` still forces a named node, which
is the escape hatch for one whose inputs are not fully declared.

Fresh is judged WHEN A NODE IS REACHED, not when the run starts.  Until 2026-09-24 the skip read the
state measured before the first node ran, so a node that was fresh then was skipped even after an
upstream node in the same run had rewritten its input -- and was left stale.  ``islet_leiden`` was
that case: ``build_workbook`` rewrites ``islet_celltype_proportions.csv``, and the Leiden table
would have stayed the one built on the previous proportions.

WHAT COUNTS AS AN INPUT
-----------------------
The builder's own source file is declared as an input alongside its data.  A builder edited after
its artifact was written means the artifact no longer describes that builder, which is the same
kind of staleness as a newer data file and is reported the same way.  The table names which input
won, so code staleness and data staleness are told apart at a glance, and the lines under it list
EVERY input newer than the node's oldest output: a node stale for its code and for its data says
both, where the table alone would show whichever is newer.

A directory input is scanned recursively and the newest file in it wins; a ``*`` in a declared path
is a glob.  ``typing/geometry`` holds two unrelated families written by two different builders
(``endocrine_union_*.parquet`` and ``nonendocrine_*.geojson``), so those are declared as globs and
not as the directory -- declaring the directory would make each builder its own input.

AN INTERRUPTED RUN
------------------
``--run`` runs one node at a time, so a kill interrupts at most one builder: every node before it
finished, and none after it started.  What that builder leaves depends on how it writes.

* ``sections``, ``availability`` and ``distance_profiles`` write a ``.parquet.tmp`` and rename it
  (build_sections_layer.py:266-276, build_distance_profiles.py:1544-1548), so a kill leaves their
  previous artifact whole.  Every other builder writes in place, and a kill during a write can
  leave a truncated file at the path the app reads.  So a DATA run copies the files it will
  replace before it starts (Step 7 of the 2026-09-21 plan), and a rerun that cannot finish is
  undone by copying them back.
* A builder that writes per donor -- ``islet_geometry``'s geojsons, the four ``app_*`` nodes, the
  donor partitions of ``pln``, ``restain`` and ``pln_node_membership`` -- is left with some donors
  rewritten and the rest not.  ``--check`` reports that as stale: the donors not yet rewritten
  predate the input that made the node stale.  The one case that rule cannot see is a run forced
  with ``--only`` on a node that was fresh, where every output, old or new, postdates every input.
  A node that declares ``written_last`` -- the one file its builder writes after all the others --
  is stale while any output is newer than it; ``islet_geometry`` declares its parquet, which
  build_islet_geometry.py writes after the donor loop that writes the geojsons (:257, :270).
* The remedy is the same every time: rerun the node through ``--run``.  It is red, so it is not
  skipped, and every builder rewrites its whole cohort; nothing resumes where it stopped.  Do not
  finish a killed run by hand with ``--donors``.  Twelve builders take that flag.  Eleven write a
  cohort file from the donors given alone (``build_parquet_views.py`` refuses a write that would
  drop islets it already holds), and ``build_app_data.py`` writes per-donor files only.  The
  nodes whose numbers are computed across donors -- ``islet_leiden``'s pooled ordering,
  ``grid_composition``, ``compartments``' batch eta-squared, ``trajectory`` -- are wrong on a
  subset even where nothing is dropped.

ONE WRITER OF islets.parquet, SINCE 2026-09-22
-----------------------------------------------
``build_parquet_views.py`` is the only node that writes ``app/parquet/islets.parquet``.  The two
islet-level rollups -- ``summarise_crossed_fraction.py`` (14 ``crossed_*`` columns) and
``summarise_singlecell_trajectory.py`` (7 ``sc_*`` columns) -- write only their own summary
files, which ``build_parquet_views.py`` merges through its ``ROLLUPS`` table.  Until 2026-09-22
each rollup also re-merged itself into ``islets.parquet`` afterwards, so three nodes wrote one
path and the graph needed a declared ``after`` chain to keep the last writer from dropping the
other two's columns.  The dependency is now the two summary files themselves, declared as inputs
of ``parquet_views``, and ``build_graph`` still REFUSES any graph in which two nodes write one
path without an ordering between them.

ONE DECLARED CYCLE BREAK
------------------------
``build_nonendocrine_outline.py:97`` reads ``app/islet_geometry.parquet`` into a local named
``geometry`` and never uses it; the cross-check below it runs off ``cells["cell_area"]``.  The real
flow is one-way -- ``build_nonendocrine_outline`` writes ``endocrine_union_*.parquet`` and
``build_islet_geometry`` reads them -- so that dead read is declared as a ``soft_input``: recorded,
printed by ``--plan``, and deliberately not an edge.  Reported as a finding rather than edited here;
that file belongs to another lane.

WHAT IS DELIBERATELY NOT A NODE
-------------------------------
``nonendocrine_outline`` was a node until 2026-09-21 and is no longer one.  It streams
``data/redsea_scratch/geojson`` -- 61 GB across 22 single-line FeatureCollections, up to 5.5 GB
each, parsed one feature at a time through ``ijson`` -- to union each islet's endocrine cell
outlines.  It cost 987 s of a 26-node rebuild, and that cost is JSON parsing, not geometry.

Running it on every rebuild buys nothing, because the numbers it feeds are already measured per
cell:

* the non-endocrine AREA never came from the polygon.  ``build_islet_geometry.py``'s own
  docstring says so -- "its area is the summed segmented area of the non-endocrine cells in it
  -- a measurement, not the region area minus an estimate" -- and every row it writes already
  carries ``cells_um2 = subset["cell_area"].sum()``;
* the STRUCTURE areas are the same evidence twice.  ``data/development/typing/crossed/
  crossed_<donor>.parquet`` already holds ``object_id, cell_um2, crossed_frac_vessel,
  crossed_frac_lymphatic, crossed_frac_nerve`` -- the per-cell fraction each classified object
  covers -- so the area over any cell set is a join and a sum.

The one thing only the polygon gives is the drawn outline.  The script is kept and still runs on
demand for that; it is simply not part of the routine build, so a typing change no longer costs
sixteen minutes of JSON parsing to re-derive what the matrix already holds.

The upstream typing lane (``typing/cells``, ``typing/immune``, ``typing/crossed``,
``typing/islet_membership``, ``celltype/uniform``) and the two export ingests
(``rebuild_donor_from_revised_export.py``, ``ingest_pln_cells.py``) are SOURCES here, not nodes.
They have no single declared command -- the rebuild takes a donor list, the ingest reads a 43.7 GB
TSV -- so this driver names them as the newest input when they are the reason something is stale,
and does not offer to run them.  ``SOURCE_PRODUCERS`` records what does produce each, so the check
table can say so out loud.

A source can itself be stale against what it is built from, and a node reading it would never
show it.  ``typing/immune`` is the case measured on 2026-09-21: ``immune_subtype.py`` subtyped
2,934,898 Immune cells on 2026-09-17, ``typing/cells`` was rewritten on 2026-09-19 with 2,875,751,
and ``grid_composition`` reads both.  ``SOURCE_CHECKS`` declares such a source with its own inputs,
and ``--check`` reports its state under the table with the nodes that read it.  It is reported,
never run, and never part of the exit code: rerunning it is the typing lane's step, on its own go.

The refinement pass is part of that lane, not a node.  Since 2026-09-23 (maintainer's ruling)
``typing/cells`` is written by ``refine.py --out data/development/typing/refine`` followed by
``build_cells.py --refine data/development/typing/refine``, which joins the pass's verdict beside
``cell_class``.  ``refine.py`` READS ``typing/cells`` and the join REWRITES it, so a ``refine``
node here would go stale every time its own result was joined back, and once rerun it would
report fresh while the served table still carried an older join.  The pass stays with the lane
that consumes it.

``--include-marginal`` is spelled into the app-artifact commands because the cohort those artifacts
were rebuilt on is ``cohort.CLEAN`` (17) + ``cohort.MARGINAL`` (3) = 20, each marginal donor
carrying its own measured defect (docs/APP_COVERAGE_2026-09-17.md, "The --include-marginal
decision").  It is a recorded choice, not a default.
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
#: The repo venv, never a bare `python`: /usr/bin/python3 lacks pyarrow and duckdb here.
VENV_PYTHON = REPO / ".venv" / "bin" / "python"

#: The historical content-addressed run tree, spelled ONCE and taken from the one table that
#: knows it is historical -- core/phenocycler/scripts/typing/superseded.py -- rather than typed
#: out again here. `build_islet_geometry.py` really does read it (the 20-donor cohort output
#: lives nowhere else), so it is a declared input of `islet_geometry` and has to stay one; what
#: it must NOT be is a second, unattributed spelling of a path this repository has already
#: recorded as superseded, drifting from it.
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
from superseded import REGISTRY as _SUPERSEDED_REGISTRY  # noqa: E402

_RUNS_ROOT = next(r for r in _SUPERSEDED_REGISTRY if r.key == "content_addressed_runs").fragments[0]
#: The single run directory on disk, and the only stage of it any node here declares.
HISTORICAL_RUN_CELLS = f"{_RUNS_ROOT}/30f6689e4d27d295f890/cells"

FRESH = "fresh"
STALE = "stale"
MISSING_OUTPUT = "missing_output"
MISSING_INPUT = "missing_input"
#: Every state but `fresh` makes --check exit non-zero. A gate that passed on a missing
#: artifact would be a gate that passes before the build has happened.
RED_STATES = (STALE, MISSING_OUTPUT, MISSING_INPUT)

#: Directory entries never counted when scanning a tree for its newest file.
SKIP_DIR_NAMES = {"__pycache__", ".ipynb_checkpoints"}
SKIP_SUFFIXES = (".tmp", ".parquet.tmp")

#: Sources with no node here, and what does write them. Printed beside the check table whenever
#: one of them is the newest input, so "stale against typing/cells" names its producer rather
#: than leaving the reader to grep for it.
SOURCE_PRODUCERS = {
    "data/development/typing/geometry":
        "core/phenocycler/scripts/typing/build_nonendocrine_outline.py -- run ON DEMAND, not as "
        "a node: it streams 61 GB of per-cell GeoJSON for numbers the cell matrix already holds "
        "(see WHAT IS DELIBERATELY NOT A NODE). Only the drawn outline needs it",
    "data/development/typing/cells":
        "core/phenocycler/scripts/typing/build_cells.py (the subpanel classification lane; "
        "see core/phenocycler/scripts/typing/RECONSTRUCTION.md), with the refinement joined: "
        "refine.py --out data/development/typing/refine, then build_cells.py --include-marginal "
        "--refine data/development/typing/refine (a plain build_cells.py run first after new "
        "subpanel.py calls, because refine.py reads typing/cells; see REFINEMENT.md)",
    "data/development/typing/islet_membership":
        "core/phenocycler/scripts/typing/islet_membership.py",
    "data/development/typing/immune":
        "core/phenocycler/scripts/typing/immune_subtype.py",
    "data/development/typing/crossed":
        "core/phenocycler/scripts/typing/build_crossed_fraction.py",
    "data/development/celltype/uniform":
        "the UniFORM alignment lane",
    "data/development/typing/trajectory/singlecell_pseudotime.parquet":
        "core/phenocycler/scripts/typing/trajectory_singlecell.R (R 4.5.1)",
    "data/development/typing/trajectory/singlecell_features.parquet":
        "core/phenocycler/scripts/typing/trajectory_singlecell_features.py",
    "core/phenocycler/data/integration/cells_xen":
        "the Xenium ingest",
    "core/phenocycler/data/integration/cells_pln":
        "core/phenocycler/scripts/pln/ingest_pln_cells.py (reads /data/IO60pLN/results/"
        "IO60pLN_cells.tsv; no node here -- it takes a 43.7 GB TSV and a duckdb temp dir)",
    "data/development/typing/revised_cells":
        "core/phenocycler/scripts/typing/rebuild_donor_from_revised_export.py (takes a donor "
        "list; no node here)",
    "data/development/typing/pln/pheno":
        "core/phenocycler/scripts/typing/build_pln_immune_pheno.py",
    "data/development/restain":
        "the post-Xenium IF registration lane: core/phenocycler/scripts/restain/ and "
        "phenocycler.integration.postxen",
    "data/development/xenium_typing":
        "the Xenium typing releases (pln_immune_typing_v3.py and its predecessors); selected by "
        "their declared release_id, never by directory mtime",
    "/data/xenium_typing":
        "the Xenium typing releases written outside the repo (panc_xenium_cluster_name.py)",
    "data/development/broadtype/scratch/cosine_scale.parquet":
        "the withdrawn broad-typing producer comparison's scratch; still read by "
        "subpanel_io.load_scale()",
    "data/redsea_scratch/geojson":
        "the REDSEA scratch per-cell outline export -- the only cell outlines in this checkout",
    HISTORICAL_RUN_CELLS:
        "the historical content-addressed run tree; nothing has written to it since the run "
        "framework was removed 2026-09-16, and reading it announces itself through "
        "core/phenocycler/scripts/typing/superseded.py",
    "results/Geojson/Islets":
        "QuPath, as an export made by hand",
    "results/AnnotationsMeasurements.csv":
        "QuPath, as an export made by hand",
    "/data/IO60pLN/results/node_geojson":
        "core/phenocycler/scripts/groovy/export_pln_node_annotations.groovy",
    "/data/islet_layers/pln/provenance/qupath_annotation_audit.log":
        "the captured stdout of export_pln_node_annotations.groovy; it lives under the pln "
        "artifact root but it is an input, not something build_pln_layer.py writes",
}

#: Sources whose OWN staleness --check reports, measured exactly as a node's is (oldest file of
#: the source against the newest file of each thing it is built from), and which --run never
#: starts. Repo-relative. A stale source does not change --check's exit code: the gate covers
#: what this driver can rebuild, and rerunning a source is the typing lane's own step, taken on
#: its own go. What it does is stop a node that reads the source from looking sound while it
#: pairs two generations of the typing -- grid_composition reads typing/immune beside
#: typing/cells (build_grid_composition.py:102, :585).
SOURCE_CHECKS = {
    "data/development/typing/immune": (
        "core/phenocycler/scripts/typing/immune_subtype.py",
        "data/development/typing/cells",
        "data/development/celltype/uniform",
        "data/development/broadtype/scratch/cosine_scale.parquet",
    ),
}


@dataclass(frozen=True)
class Node:
    """One builder, and the files it reads and writes.

    ``inputs`` and ``outputs`` are repo-relative or absolute paths; a ``*`` makes one a glob and a
    directory is scanned recursively.  ``after`` names nodes this one must follow for a reason the
    files do not express.  ``soft_inputs`` are files the builder opens that are deliberately NOT
    edges, each with its reason in ``why``.  ``written_last`` names the one output file a complete
    run writes after every other output; any output newer than it was written by a run that
    stopped partway (see AN INTERRUPTED RUN in the module docstring).
    """

    name: str
    command: tuple[str, ...]
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    after: tuple[str, ...] = ()
    soft_inputs: tuple[str, ...] = ()
    why: str = ""
    written_last: str = ""

    def __post_init__(self) -> None:
        if self.written_last and ("*" in self.written_last
                                  or self.written_last not in self.outputs):
            raise ValueError(
                f"{self.name}: written_last must be one of its declared outputs and one file, "
                f"not a glob; got {self.written_last!r}")


@dataclass
class NodeStatus:
    """What ``--check`` measured for one node."""

    name: str
    state: str
    newest_input: str | None = None
    newest_input_mtime: float | None = None
    oldest_output: str | None = None
    oldest_output_mtime: float | None = None
    missing: tuple[str, ...] = field(default_factory=tuple)
    #: Every declared input newer than the oldest output, as (newest file in it, mtime), newest
    #: first. ``newest_input`` is only the first of these; a node stale for its code AND its data
    #: lists both here.
    stale_inputs: tuple[tuple[str, float], ...] = field(default_factory=tuple)
    #: Output files newer than the node's ``written_last``: left by a run that did not finish.
    unfinished: tuple[str, ...] = field(default_factory=tuple)

    @property
    def is_red(self) -> bool:
        return self.state in RED_STATES


# --------------------------------------------------------------------------------------------
# The declared graph
# --------------------------------------------------------------------------------------------

def _tokens(repo: Path, layers_dir: Path) -> dict[str, str]:
    typing = repo / "data/development/typing"
    return {
        "repo": str(repo),
        "layers": str(layers_dir),
        "typing": str(typing),
        "app": str(typing / "app"),
        "integration": str(repo / "core/phenocycler/data/integration"),
        "senior": str(repo / "scripts/senior"),
        "typingscripts": str(repo / "core/phenocycler/scripts/typing"),
        "plnscripts": str(repo / "core/phenocycler/scripts/pln"),
        # The QuPath export tree sits beside the checkout, exactly as subpanel_io.py resolves it.
        "qupath_results": str(repo.parent / "results"),
        "py": str(VENV_PYTHON),
    }


def _expand(values, tokens: dict[str, str]) -> tuple[str, ...]:
    """Command arguments: substituted, otherwise left exactly as declared."""
    return tuple(str(v).format(**tokens) for v in values)


def _expand_paths(values, tokens: dict[str, str]) -> tuple[str, ...]:
    """Paths: substituted and normalised, so a declared `..` prints as the path it names.

    normpath, never resolve(): resolve() would follow the symlinks some data roots are mounted
    through and report a path no declaration contains.
    """
    return tuple(os.path.normpath(str(v).format(**tokens)) for v in values)


def _declaration() -> list[dict]:
    """The graph as written, before token expansion. Read this to know what runs."""
    return [
        # ``nonendocrine_outline`` IS NOT A NODE, as of 2026-09-21 -- see the section of that
        # name in the module docstring. It streams 61 GB of per-cell GeoJSON to derive numbers
        # the cell matrix already carries.
        dict(
            name="islet_geometry",
            command=("{py}", "{typingscripts}/build_islet_geometry.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/build_islet_geometry.py",
                "{typing}/cells",
                "{typing}/geometry/endocrine_union_*.parquet",
                "{repo}/core/phenocycler/data/qupath_manifest.json",
                "{repo}/" + HISTORICAL_RUN_CELLS,
                "{qupath_results}/Geojson/Islets",
            ),
            outputs=(
                "{app}/islet_geometry.parquet",
                "{typing}/geometry/nonendocrine_*.geojson",
            ),
            written_last="{app}/islet_geometry.parquet",
            why="The islet floor every per-islet artifact meets: 18,404 core outlines, 20 donors. "
                "Each donor's nonendocrine_<donor>.geojson is written inside the donor loop "
                "(build_islet_geometry.py:257) and the parquet once, after it (:270), so a "
                "geojson newer than the parquet is a run that stopped partway.",
        ),

        # -- the typing/app artifacts, in the order docs/APP_COVERAGE_2026-09-17.md gives -------
        dict(
            name="aggregate_islets",
            command=("{py}", "{typingscripts}/aggregate_islets.py", "--include-marginal"),
            inputs=("{typingscripts}/aggregate_islets.py", "{typing}/cells"),
            outputs=(
                "{app}/islet_composition_long.parquet",
                "{app}/islet_hormone_composition.parquet",
            ),
            why="Drops the exocrine-inside-an-islet-border set from the composition denominator.",
        ),
        dict(
            name="build_workbook",
            command=("{py}", "{typingscripts}/build_workbook.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/build_workbook.py",
                "{app}/donor_metadata.csv",
                "{app}/islet_geometry.parquet",
                "{app}/islet_composition_long.parquet",
                "{app}/islet_hormone_composition.parquet",
            ),
            outputs=(
                "{app}/master_results_run10.xlsx",
                "{app}/islet_celltype_proportions.csv",
                "{app}/islet_spatial_lookup.csv",
            ),
        ),
        dict(
            name="neighborhood",
            command=("{py}", "{typingscripts}/neighborhood.py", "--include-marginal"),
            inputs=("{typingscripts}/neighborhood.py", "{typing}/cells"),
            outputs=("{app}/neighborhood_metrics.csv",),
        ),
        dict(
            name="islet_marker_means",
            command=("{py}", "{typingscripts}/build_islet_marker_means.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/build_islet_marker_means.py",
                "{typing}/cells",
                "{repo}/data/development/celltype/uniform",
            ),
            outputs=("{app}/islet_marker_means.parquet",),
            why="Takes islet membership from typing/cells joined on object_id, never from the "
                "aligned table's own islet_num -- the wrong-islet join found 2026-09-17.",
        ),
        dict(
            name="islet_leiden",
            command=("{py}", "{typingscripts}/trajectory_islet.py", "--include-marginal",
                     "--leiden-only"),
            inputs=(
                "{typingscripts}/trajectory_islet.py",
                "{app}/islet_celltype_proportions.csv",
                "{app}/donor_metadata.csv",
            ),
            outputs=("{app}/islet_leiden.parquet",),
            why="--leiden-only (2026-09-24): the Leiden table from the graph the pooled ordering "
                "is built on, without the ordering -- no root, no pseudotime, no method check, "
                "no H5AD. The full run writes the identical table (measured cell for cell on "
                "2026-09-24) and then exits 1 when its ordering fails pseudotime_is_not_mostly_"
                "donor, which stopped every node after this one. The pooled ordering and its H5AD "
                "are run by hand, not by this graph; adata_ins_root.h5ad has its own producer "
                "(scripts/rebuild_adata_ins_root.py).",
        ),
        dict(
            name="galectin",
            command=("{py}", "{typingscripts}/build_galectin_layer.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/build_galectin_layer.py",
                "{qupath_results}/AnnotationsMeasurements.csv",
                "{typing}/cells",
            ),
            outputs=("{app}/islet_galectin.parquet",),
        ),
        dict(
            name="parquet_views",
            command=("{py}", "{typingscripts}/build_parquet_views.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/build_parquet_views.py",
                "{app}/islet_celltype_proportions.csv",
                "{app}/islet_geometry.parquet",
                "{app}/donor_metadata.csv",
                "{app}/neighborhood_metrics.csv",
                "{app}/islet_leiden.parquet",
                "{app}/islet_galectin.parquet",
                "{app}/master_results_run10.xlsx",
                "{typing}/crossed/crossed_islet_summary.parquet",
                "{typing}/trajectory/singlecell_islet_summary.parquet",
            ),
            outputs=(
                "{app}/parquet/islets.parquet",
                "{app}/parquet/uns_targets.parquet",
                "{app}/parquet/uns_composition.parquet",
            ),
            why="The one writer of islets.parquet. The two rollup summaries above are merged in "
                "here (ROLLUPS in the builder); each is recorded present or absent by name.",
        ),
        dict(
            name="crossed_fraction",
            command=("{py}", "{typingscripts}/summarise_crossed_fraction.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/summarise_crossed_fraction.py",
                "{typing}/crossed",
                "{typing}/cells",
            ),
            outputs=("{typing}/crossed/crossed_islet_summary.parquet",),
            why="The 14 crossed_* columns per islet, as a summary file only; parquet_views "
                "reads it.",
        ),
        dict(
            name="singlecell_trajectory",
            command=("{py}", "{typingscripts}/summarise_singlecell_trajectory.py"),
            inputs=(
                "{typingscripts}/summarise_singlecell_trajectory.py",
                "{typing}/trajectory/singlecell_pseudotime.parquet",
            ),
            outputs=("{typing}/trajectory/singlecell_islet_summary.parquet",),
            why="The 7 sc_* columns per islet, as a summary file only; parquet_views reads it. "
                "Takes no --include-marginal; its cohort comes from singlecell_pseudotime.parquet.",
        ),

        # -- the per-cell app artifacts the drill-down and the Spatial tab read ------------------
        dict(
            name="app_tissue",
            command=("{py}", "{typingscripts}/build_app_data.py", "--include-marginal",
                     "--only", "tissue"),
            inputs=("{typingscripts}/build_app_data.py", "{typing}/cells"),
            outputs=("{app}/parquet/tissue",),
        ),
        dict(
            name="app_drilldown_cells",
            command=("{py}", "{typingscripts}/build_app_data.py", "--include-marginal",
                     "--only", "cells"),
            inputs=(
                "{typingscripts}/build_app_data.py",
                "{typing}/cells",
                "{repo}/data/development/celltype/uniform",
                "{repo}/data/development/broadtype/scratch/cosine_scale.parquet",
            ),
            outputs=("{app}/cells",),
            why="One CSV per islet, identity by filename. Its file count should match "
                "islet_geometry's islet count.",
        ),
        dict(
            name="app_outlines",
            command=("{py}", "{typingscripts}/build_app_data.py", "--include-marginal",
                     "--only", "outlines"),
            inputs=(
                "{typingscripts}/build_app_data.py",
                "{repo}/core/phenocycler/data/qupath_manifest.json",
                "{qupath_results}/Geojson/Islets",
                "{typing}/geometry/nonendocrine_*.geojson",
            ),
            outputs=("{app}/json", "{app}/gson"),
        ),
        dict(
            name="app_cell_distance",
            command=("{py}", "{typingscripts}/build_app_data.py", "--include-marginal",
                     "--only", "distance"),
            inputs=(
                "{typingscripts}/build_app_data.py",
                "{typing}/cells",
                "{app}/donor_metadata.csv",
            ),
            outputs=("{app}/parquet/cell_distance",),
            why="The writer is spelled `distance`, not `cell_distance`; the OUTPUT root is "
                "parquet/cell_distance.",
        ),

        # -- the pLN node membership the pln layer is built from ---------------------------------
        dict(
            name="pln_node_membership",
            command=("{py}", "{senior}/derive_pln_node_membership.py",
                     "--centroids", "/data/tmp/duckdb_pln/pln_centroids.tsv"),
            inputs=(
                "{senior}/derive_pln_node_membership.py",
                "{integration}/pln_node_annotations.csv",
                "/data/IO60pLN/results/node_geojson",
                "/data/tmp/duckdb_pln/pln_centroids.tsv",
            ),
            outputs=(
                "{repo}/data/development/pln/node_membership",
                "{repo}/data/development/pln/node_geometry_audit.csv",
                "{repo}/data/development/pln/node_membership_summary.csv",
            ),
            why="The centroid cut is not kept in the checkout; the duckdb COPY that remakes it is "
                "in the script's docstring. Membership is point-in-polygon, never the TSV's Parent.",
        ),

        # -- the nine layer artifacts -------------------------------------------------------------
        dict(
            name="sections",
            command=("{py}", "{senior}/build_sections_layer.py"),
            inputs=(
                "{senior}/build_sections_layer.py",
                "{integration}/restain_sections.csv",
                "{integration}/xenium_paths.csv",
            ),
            outputs=("{layers}/sections/data.parquet",),
        ),
        dict(
            name="availability",
            command=("{py}", "{senior}/build_availability.py"),
            inputs=(
                "{senior}/build_availability.py",
                "{senior}/build_sections_layer.py",
                "{integration}/restain_sections.csv",
                "{integration}/xenium_paths.csv",
                "{integration}/cells_xen",
                "{typing}/cells",
                "{typing}/geometry/endocrine_union_*.parquet",
                "{typing}/trajectory/singlecell_islet_summary.parquet",
                "{repo}/data/development/restain",
            ),
            outputs=("{layers}/availability/data.parquet",),
            why="Four capabilities, each read by a card: phenocycler_panc and "
                "islet_outline (card_x1_gradient.R:284), restain_if (card_x2_overlay.R:47), "
                "trajectory (card_s2_trajectory.R:56). Until 2026-09-23 it also wrote "
                "phenocycler_pln, xenium_panc, xenium_pln and pln_nodes, which no card read; with "
                "them went its reads of both Xenium releases and pln_node_annotations.csv, so "
                "none is declared here. cells_xen is read for the IF pairing only.",
        ),
        dict(
            name="distance_profiles",
            command=("{py}", "{senior}/build_distance_profiles.py"),
            inputs=(
                "{senior}/build_distance_profiles.py",
                "{senior}/build_sections_layer.py",
                "{typing}/cells",
                "{typing}/islet_membership",
                "{app}/donor_metadata.csv",
                "{integration}/cells_xen",
                "{integration}/cells_pln",
                "{integration}/manifest.csv",
                "{integration}/vocab_crosswalk.csv",
                "{typing}/pln/pheno",
                "/data/xenium_typing/panc_2026-09-17/release",
                "{repo}/data/development/xenium_typing",
            ),
            outputs=("{layers}/distance_profiles/data.parquet",),
            why="Carries the study question: the immune fraction against signed distance to the "
                "islet edge. Refuses subpanel_calls.parquet by path. Every Xenium call comes from "
                "the two declared releases; cells_xen is declared for a pancreas cell's "
                "coordinates only.",
        ),
        dict(
            name="restain",
            command=("{py}", "{senior}/build_restain_layer.py"),
            inputs=(
                "{senior}/build_restain_layer.py",
                "{repo}/data/development/restain",
                "{integration}/restain_sections.csv",
                "{integration}/xenium_paths.csv",
                "{integration}/cells_xen",
            ),
            outputs=(
                "{layers}/restain/availability.parquet",
                "{layers}/restain/cells",
            ),
        ),
        dict(
            name="islet_summary",
            command=("{py}", "{senior}/build_islet_summary.py"),
            inputs=(
                "{senior}/build_islet_summary.py",
                "{app}/islet_geometry.parquet",
                "{app}/islet_composition_long.parquet",
                "{typing}/geometry/endocrine_union_*.parquet",
            ),
            outputs=(
                "{layers}/islet_summary/islets.parquet",
                "{layers}/islet_summary/composition.parquet",
                "{layers}/islet_summary/donors.parquet",
                "{layers}/islet_summary/availability.parquet",
            ),
        ),
        dict(
            name="trajectory",
            command=("{py}", "{senior}/build_trajectory_layer.py"),
            inputs=(
                "{senior}/build_trajectory_layer.py",
                "{typing}/trajectory/singlecell_islet_summary.parquet",
                "{typing}/trajectory/singlecell_pseudotime.parquet",
                "{typing}/trajectory/singlecell_trajectory.json",
            ),
            outputs=("{layers}/trajectory/data.parquet",),
        ),
        dict(
            name="pln",
            command=("{py}", "{plnscripts}/build_pln_layer.py"),
            inputs=(
                "{plnscripts}/build_pln_layer.py",
                "{integration}/pln_node_annotations.csv",
                "{integration}/cells_pln_ingest_report.json",
                "{integration}/cells_pln",
                "{repo}/data/development/pln/node_membership",
                "{repo}/data/development/pln/node_geometry_audit.csv",
                "/data/IO60pLN/results/node_geojson",
                "{layers}/pln/provenance/qupath_annotation_audit.log",
            ),
            outputs=(
                "{layers}/pln/MANIFEST.json",
                "{layers}/pln/nodes.parquet",
                "{layers}/pln/donor_id=*/pln_node_membership.parquet",
            ),
            why="Its provenance log lives under the artifact root and is an INPUT, so the outputs "
                "are declared file by file rather than as the directory.",
        ),
        dict(
            name="compartments",
            command=("{py}", "{senior}/build_compartment_gradient.py"),
            inputs=(
                "{senior}/build_compartment_gradient.py",
                "{typing}/cells",
                "{typing}/islet_membership",
                "{typing}/pln/pheno",
                "{repo}/data/development/xenium_typing/release_pln_immune_v3",
                "{app}/donor_metadata.csv",
                "{repo}/data/donor_acquisition_runs.csv",
                "{integration}/vocab_crosswalk.csv",
            ),
            outputs=(
                "{layers}/compartments/compartment_composition.csv",
                "{layers}/compartments/islet_size_strata.csv",
                "{layers}/compartments/islet_size_summary.csv",
                "{layers}/compartments/batch_eta_squared.csv",
                "{layers}/compartments/manifest.json",
            ),
            why="No longer pinned: since 2026-09-18 it resolves the pLN release through "
                "scripts/senior/xenium_release.py, which honours _v3's declaration that it "
                "supersedes _v2, and pLN is the only release it resolves "
                "(build_compartment_gradient.py:870). The declared input is the v3 root it "
                "reads, so a new v3 write makes it stale; until 2026-09-23 it declared the "
                "superseded _v2 and the staleness check watched the wrong tree. Not in "
                "docs/ISLET_LAYERS_SPEC.md; it exists anyway.",
        ),
        dict(
            name="grid_composition",
            command=("{py}", "{senior}/build_grid_composition.py"),
            inputs=(
                "{senior}/build_grid_composition.py",
                "{typing}/cells",
                "{typing}/immune",
                "{typing}/pln/pheno",
                "/data/xenium_typing/panc_2026-09-17/release",
                "{repo}/data/development/xenium_typing/release_pln_immune_v3",
                "{integration}/vocab_crosswalk.csv",
                "{integration}/pln_node_annotations.csv",
                "{app}/donor_metadata.csv",
                "{repo}/data/donor_acquisition_runs.csv",
            ),
            outputs=(
                "{layers}/grid_composition/composition.csv",
                "{layers}/grid_composition/fine.csv",
                "{layers}/grid_composition/presence.csv",
                "{layers}/grid_composition/crossmodal_tnk.csv",
                "{layers}/grid_composition/batch_eta_squared.csv",
                "{layers}/grid_composition/manifest.json",
            ),
            why="The one builder reading both current Xenium releases, and the only place all four "
                "tissue x modality cells are populated. Not in ISLET_LAYERS_SPEC.md; it exists "
                "anyway. One of its two pins is absolute and outside the repo.",
        ),
    ]


#: The nine layer artifacts, named so --check can summarise them separately from the app lane.
LAYER_NODES = (
    "availability", "sections", "distance_profiles", "restain", "islet_summary",
    "trajectory", "pln", "compartments", "grid_composition",
)


def build_graph(repo: Path = REPO, layers_dir: Path | None = None) -> dict[str, Node]:
    """Expand the declaration and refuse a graph that cannot be ordered safely."""
    layers_dir = Path(layers_dir or os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers"))
    tokens = _tokens(Path(repo), layers_dir)
    nodes: dict[str, Node] = {}
    for spec in _declaration():
        node = Node(
            name=spec["name"],
            command=_expand(spec["command"], tokens),
            inputs=_expand_paths(spec["inputs"], tokens),
            outputs=_expand_paths(spec["outputs"], tokens),
            after=tuple(spec.get("after", ())),
            soft_inputs=_expand_paths(spec.get("soft_inputs", ()), tokens),
            why=spec.get("why", ""),
            written_last=(_expand_paths((spec["written_last"],), tokens)[0]
                          if spec.get("written_last") else ""),
        )
        if node.name in nodes:
            raise ValueError(f"Duplicate node name: {node.name}")
        nodes[node.name] = node
    for node in nodes.values():
        for earlier in node.after:
            if earlier not in nodes:
                raise ValueError(f"{node.name}: `after` names an unknown node {earlier!r}")
    _assert_shared_outputs_are_ordered(nodes)
    topological_order(nodes)          # raises on a cycle
    return nodes


def producers(nodes: dict[str, Node]) -> dict[str, list[str]]:
    """output path -> the node names that write it."""
    out: dict[str, list[str]] = {}
    for node in nodes.values():
        for path in node.outputs:
            out.setdefault(path, []).append(node.name)
    return out


def edges(nodes: dict[str, Node]) -> dict[str, set[str]]:
    """node -> the nodes it must follow. A file edge, or an explicit `after`.

    A path a node both reads and writes is a rewrite in place, not a dependency on the other
    nodes that rewrite it: islets.parquet is read and rewritten by both merge steps, so inferring
    an edge from every writer of it would make the three-writer chain circular.  The order among
    rewriters is declared with `after`, which is exactly where that decision should be visible.
    """
    by_output = producers(nodes)
    result = {name: set(node.after) for name, node in nodes.items()}
    for name, node in nodes.items():
        for path in node.inputs:
            if path in node.outputs:
                continue
            for upstream in by_output.get(path, ()):
                if upstream != name:
                    result[name].add(upstream)
    return result


def topological_order(nodes: dict[str, Node]) -> list[str]:
    """Kahn, taking the first READY node in declaration order each step.

    Taking one at a time rather than a whole ready set keeps the printed order as close to the
    order this file declares as the dependencies allow, so --plan reads the way the declaration
    reads.  It is deterministic, which matters because the test pins it.
    """
    dependencies = edges(nodes)
    declared = list(nodes)
    remaining = set(declared)
    order: list[str] = []
    done: set[str] = set()
    while remaining:
        nxt = next((n for n in declared
                    if n in remaining and not (dependencies[n] - done)), None)
        if nxt is None:
            raise ValueError(
                "The declared graph has a cycle among: " + ", ".join(sorted(remaining)))
        order.append(nxt)
        done.add(nxt)
        remaining.discard(nxt)
    return order


def _assert_shared_outputs_are_ordered(nodes: dict[str, Node]) -> None:
    """Two nodes writing one path must be totally ordered, or the last writer is luck.

    This is the islets.parquet hazard: build_parquet_views.py rewrites it from scratch and drops
    the 21 columns the other two merge in.  The graph refuses to build unless the three are in a
    declared order, so build_parquet_views cannot be last.
    """
    dependencies = edges(nodes)

    def reaches(later: str, earlier: str, seen: set[str] | None = None) -> bool:
        seen = seen if seen is not None else set()
        if later in seen:
            return False
        seen.add(later)
        for upstream in dependencies[later]:
            if upstream == earlier or reaches(upstream, earlier, seen):
                return True
        return False

    for path, writers in producers(nodes).items():
        if len(writers) < 2:
            continue
        for i, a in enumerate(writers):
            for b in writers[i + 1:]:
                if not (reaches(a, b) or reaches(b, a)):
                    raise ValueError(
                        f"{a} and {b} both write {path} and nothing orders them. A path with two "
                        f"writers and no order means whichever ran last decides what is in it; "
                        f"declare the order with `after`.")


# --------------------------------------------------------------------------------------------
# Measuring
# --------------------------------------------------------------------------------------------

def _files_under(root: Path):
    """Every file under a directory, skipping caches and half-written temporaries."""
    stack = [root]
    while stack:
        current = stack.pop()
        try:
            entries = list(os.scandir(current))
        except (NotADirectoryError, PermissionError, FileNotFoundError):
            continue
        for entry in entries:
            if entry.name in SKIP_DIR_NAMES or entry.name.startswith("."):
                continue
            if entry.is_dir(follow_symlinks=False):
                stack.append(Path(entry.path))
            elif entry.is_file(follow_symlinks=False):
                if entry.name.endswith(SKIP_SUFFIXES):
                    continue
                yield entry


def _matches(path: str) -> list[tuple[float, str]]:
    """(mtime, file) for every file a declared path names: one file, a directory tree or a glob."""
    target = Path(path)
    if "*" in path:
        parent = Path(path).parent
        pattern = Path(path).name
        if "*" in str(parent):
            # a glob in the directory part, e.g. donor_id=*/file.parquet
            base = Path(path)
            anchor = base
            while "*" in str(anchor):
                anchor = anchor.parent
            candidates = list(anchor.glob(str(base.relative_to(anchor))))
        else:
            candidates = list(parent.glob(pattern))
        return [(c.stat().st_mtime, str(c)) for c in candidates if c.is_file()]
    if target.is_dir():
        return [(e.stat().st_mtime, e.path) for e in _files_under(target)]
    if target.is_file():
        return [(target.stat().st_mtime, str(target))]
    return []


def _extremum(path: str, newest: bool) -> tuple[float, str] | None:
    """(mtime, the file carrying it) over a file, a directory tree or a glob; None if absent."""
    stats = _matches(path)
    if not stats:
        return None
    return max(stats) if newest else min(stats)


def _unfinished(node: Node) -> tuple[str, ...]:
    """Output files newer than the one a complete run of this node writes last.

    The oldest-output rule cannot see this.  A run forced with ``--only`` on a node that was
    fresh, and killed partway, leaves every output newer than every input -- some from the killed
    run, the rest from the run before it -- and reads as fresh.  The builder's last write is the
    one file that says whether a run finished, so anything newer than it did not.
    """
    if not node.written_last:
        return ()
    last = _extremum(node.written_last, newest=True)
    if last is None:
        return ()
    later = {path for output in node.outputs if output != node.written_last
             for mtime, path in _matches(output) if mtime > last[0] and path != last[1]}
    return tuple(sorted(later))


def check_node(node: Node) -> NodeStatus:
    """Compare the node's oldest output against every input it reads. Reads mtimes only."""
    missing_inputs = tuple(p for p in node.inputs if _extremum(p, newest=True) is None)
    if missing_inputs:
        return NodeStatus(node.name, MISSING_INPUT, missing=missing_inputs)
    missing_outputs = tuple(p for p in node.outputs if _extremum(p, newest=False) is None)
    if missing_outputs:
        newest = max((_extremum(p, newest=True) for p in node.inputs), default=None)
        return NodeStatus(node.name, MISSING_OUTPUT, missing=missing_outputs,
                          newest_input=newest[1] if newest else None,
                          newest_input_mtime=newest[0] if newest else None)
    newest_per_input = [_extremum(p, newest=True) for p in node.inputs]
    newest = max(newest_per_input)
    oldest = min(_extremum(p, newest=False) for p in node.outputs)
    stale_inputs = tuple((path, mtime) for mtime, path in
                         sorted((x for x in newest_per_input if x[0] > oldest[0]), reverse=True))
    unfinished = _unfinished(node)
    state = STALE if stale_inputs or unfinished else FRESH
    return NodeStatus(node.name, state, newest_input=newest[1], newest_input_mtime=newest[0],
                      oldest_output=oldest[1], oldest_output_mtime=oldest[0],
                      stale_inputs=stale_inputs, unfinished=unfinished)


def check(nodes: dict[str, Node], only: list[str] | None = None) -> list[NodeStatus]:
    """One NodeStatus per node, in topological order. The public entry point for Lane V."""
    order = topological_order(nodes)
    if only:
        order = [n for n in order if n in set(only)]
    return [check_node(nodes[name]) for name in order]


def check_sources(repo: Path = REPO, checks: dict | None = None) -> list[NodeStatus]:
    """One NodeStatus per reported source (SOURCE_CHECKS), named by its repo-relative path."""
    checks = SOURCE_CHECKS if checks is None else checks
    return [check_node(Node(name=source, command=(),
                            inputs=tuple(str(Path(repo) / p) for p in inputs),
                            outputs=(str(Path(repo) / source),)))
            for source, inputs in checks.items()]


def readers_of(path: str, nodes: dict[str, Node]) -> list[str]:
    """The nodes that declare ``path``, or a file under it, as an input."""
    return [name for name, node in nodes.items()
            if any(i == path or i.startswith(path.rstrip("/") + "/") for i in node.inputs)]


# --------------------------------------------------------------------------------------------
# Printing
# --------------------------------------------------------------------------------------------

def _stamp(mtime: float | None) -> str:
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(mtime)) if mtime else "-"


def _short(path: str | None, repo: Path) -> str:
    if not path:
        return "-"
    try:
        return str(Path(path).relative_to(repo))
    except ValueError:
        return path


def _row(status: NodeStatus, width: int, repo: Path) -> str:
    mark = "  " if status.state == FRESH else "* "
    by = _short(status.newest_input, repo)
    if status.state == MISSING_INPUT:
        by = "MISSING: " + ", ".join(_short(p, repo) for p in status.missing)
    elif status.state == MISSING_OUTPUT:
        by = f"absent: {', '.join(_short(p, repo) for p in status.missing)}"
    elif status.state == STALE and not status.stale_inputs:
        by = f"a run that did not finish: {len(status.unfinished)} outputs newer than its last write"
    return (f"{mark + status.name:<{width}}  {status.state:<14} "
            f"{_stamp(status.oldest_output_mtime):<16} "
            f"{_stamp(status.newest_input_mtime):<16} {by}")


def _why(status: NodeStatus, node: Node | None, repo: Path) -> list[str]:
    """Every reason one stale node is stale: each newer input, then an unfinished run."""
    lines = [f"    {_stamp(mtime)}  {_short(path, repo)}" for path, mtime in status.stale_inputs]
    if status.unfinished:
        last = _short(node.written_last, repo) if node else "its last write"
        named = ", ".join(Path(p).name for p in status.unfinished[:3])
        more = f" and {len(status.unfinished) - 3} more" if len(status.unfinished) > 3 else ""
        lines.append(f"    a run that did not finish: {len(status.unfinished)} outputs are newer "
                     f"than {last}, which a complete run writes last ({named}{more})")
    return lines


def print_check(statuses: list[NodeStatus], repo: Path, sources: list[NodeStatus] = (),
                nodes: dict[str, Node] | None = None) -> None:
    width = max(len(s.name) for s in statuses) + 2
    print(f"{'  node':<{width}}  {'state':<14} {'output':<16} {'newest input':<16} by")
    print("-" * (width + 90))
    for status in statuses:
        print(_row(status, width, repo))
    red = [s for s in statuses if s.is_red]
    layers_red = [s.name for s in red if s.name in LAYER_NODES]
    print()
    print(f"{len(statuses) - len(red)} fresh, {len(red)} not fresh "
          f"({len(layers_red)} of the nine layers: {', '.join(layers_red) or 'none'})")
    stale = [s for s in red if s.state == STALE]
    if stale:
        # The table names the newest input only. A node stale for its code AND its data would
        # show just one of the two, so every reason is listed here.
        print("\nWhy, node by node: every declared input newer than the node's oldest output, "
              "newest first.")
        for status in stale:
            print(f"  {status.name}")
            for line in _why(status, (nodes or {}).get(status.name), repo):
                print(line)
    named = {}
    for status in red:
        paths = [p for p, _ in status.stale_inputs] or [status.newest_input]
        for by in filter(None, paths):
            for source, producer in SOURCE_PRODUCERS.items():
                if source in by and source not in named:
                    named[source] = producer
    if named:
        print("\nThe sources named above have no node here; what writes each:")
        for source, producer in named.items():
            print(f"  {source} is written by {producer}")
    if sources:
        swidth = max(len(s.name) for s in sources) + 2
        print("\nSources this driver reports and never runs (a stale one does not change the "
              "exit code):")
        for status in sources:
            print(_row(status, swidth, repo))
            readers = readers_of(str(repo / status.name), nodes or {})
            producer = SOURCE_PRODUCERS.get(status.name, "no producer recorded")
            print(f"    written by {producer}; read by {', '.join(readers) or 'no node'}")
            for line in _why(status, None, repo) if status.state == STALE else ():
                print(line)
    if red:
        print("\nStale means the artifact was built before an input it reads, or that a run of it "
              "did not finish. It is not a\nclaim about the numbers in it; it is a claim that they "
              "were computed from a different generation of the\ninputs than the ones on disk now.")


def print_plan(nodes: dict[str, Node], repo: Path, only: list[str] | None = None) -> None:
    order = topological_order(nodes)
    if only:
        order = [n for n in order if n in set(only)]
    print(f"# {len(order)} nodes, in dependency order. Run from {repo}.")
    for position, name in enumerate(order, 1):
        node = nodes[name]
        print(f"\n# [{position}/{len(order)}] {name}")
        if node.why:
            for line in _wrapped(node.why):
                print(f"#   {line}")
        for path in node.soft_inputs:
            print(f"#   soft input (read, deliberately not an edge): {_short(path, repo)}")
        print(" ".join(shlex.quote(part) for part in node.command))


def _wrapped(text: str, width: int = 92):
    words, line = text.split(), ""
    for word in words:
        if len(line) + len(word) + 1 > width:
            yield line
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        yield line


# --------------------------------------------------------------------------------------------
# Running
# --------------------------------------------------------------------------------------------

def run(nodes: dict[str, Node], repo: Path, only: list[str] | None = None) -> int:
    """Execute in dependency order, refusing any node whose own input is stale."""
    order = topological_order(nodes)
    scheduled = [n for n in order if not only or n in set(only)]
    if not scheduled:
        print("Nothing selected.", file=sys.stderr)
        return 2
    before = {s.name: s for s in check(nodes)}
    dependencies = edges(nodes)
    rebuilt: set[str] = set()
    skipped: list[str] = []
    durations: list[tuple[str, float, int]] = []
    refused: list[tuple[str, str]] = []
    for name in scheduled:
        # Iteration is topological, so any scheduled upstream has already been attempted and is
        # in `rebuilt` only if it actually succeeded. An upstream that was refused, that failed,
        # or that was never selected is therefore still red here, which is the whole point.
        stale_upstream = sorted(
            up for up in dependencies[name] if before[up].is_red and up not in rebuilt
        )
        if stale_upstream:
            reason = ("it reads the output of a node that is not fresh and is not being rebuilt "
                      f"first: {', '.join(stale_upstream)}")
            print(f"\n== {name}: REFUSED -- {reason}")
            refused.append((name, reason))
            continue
        # Measured now, not taken from `before`: an upstream node rebuilt earlier in this run may
        # have rewritten this node's input since the run began.
        now = check_node(nodes[name])
        missing = now.missing if now.state == MISSING_INPUT else ()
        if missing:
            reason = "a declared input does not exist: " + ", ".join(missing)
            print(f"\n== {name}: REFUSED -- {reason}")
            refused.append((name, reason))
            continue

        # A node whose output already postdates every input it reads has nothing to do. Until
        # 2026-09-21 `--run` rebuilt every scheduled node, including the ones `--check` had just
        # reported fresh; until 2026-09-24 it judged fresh from `before`, the state measured
        # before any node ran. `--only` is still honoured, so naming a fresh node explicitly
        # forces it -- the escape hatch for a node whose inputs are not fully declared.
        if not now.is_red and not (only and name in set(only)):
            print(f"\n== {name}: SKIPPED -- fresh (output postdates every declared input)")
            skipped.append(name)
            rebuilt.add(name)      # downstream may proceed: this node is not stale
            continue

        command = list(nodes[name].command)
        print(f"\n== {name}\n$ {' '.join(shlex.quote(p) for p in command)}", flush=True)
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=str(repo), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
        code = process.wait()
        elapsed = time.monotonic() - started
        durations.append((name, elapsed, code))
        print(f"-- {name}: exit {code} in {elapsed:,.1f}s", flush=True)
        if code != 0:
            print(f"\n{name} failed; stopping before anything downstream of it runs.",
                  file=sys.stderr)
            _print_durations(durations, refused, skipped)
            return code
        rebuilt.add(name)
    _print_durations(durations, refused, skipped)
    return 1 if refused else 0


def _print_durations(durations, refused, skipped=()) -> None:
    if durations:
        print("\nper-node duration")
        width = max(len(n) for n, _, _ in durations)
        for name, elapsed, code in durations:
            print(f"  {name:<{width}}  {elapsed:>9,.1f}s  exit {code}")
    if skipped:
        print(f"\nskipped as fresh ({len(skipped)}): {', '.join(skipped)}")
    for name, reason in refused:
        print(f"  REFUSED {name}: {reason}")


# --------------------------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true",
                      help="compare every output's mtime against its inputs' newest; "
                           "exit non-zero if anything is stale")
    mode.add_argument("--plan", action="store_true",
                      help="print the topologically ordered commands and run nothing")
    mode.add_argument("--run", action="store_true",
                      help="execute in dependency order. SKIPS a node that is already fresh, "
                           "and refuses one whose input is stale. --only forces a named node "
                           "even when it is fresh.")
    parser.add_argument("--only", nargs="+", metavar="NODE",
                        help="restrict to these nodes (dependency order is kept)")
    parser.add_argument("--repo-root", type=Path, default=REPO)
    parser.add_argument("--layers-dir", type=Path, default=None,
                        help="default: $ISLET_LAYERS_DIR, else /data/islet_layers")
    parser.add_argument("--list", action="store_true", help="print the node names and exit")
    args = parser.parse_args(argv)

    nodes = build_graph(args.repo_root, args.layers_dir)
    if args.only:
        unknown = sorted(set(args.only) - set(nodes))
        if unknown:
            parser.error(f"unknown node(s): {', '.join(unknown)}. "
                         f"Known: {', '.join(topological_order(nodes))}")
    if args.list:
        for name in topological_order(nodes):
            print(name)
        return 0
    if not (args.check or args.plan or args.run):
        parser.error("choose one of --check, --plan, --run (or --list)")
    if args.plan:
        print_plan(nodes, args.repo_root, args.only)
        return 0
    if args.check:
        statuses = check(nodes, args.only)
        print_check(statuses, args.repo_root, sources=check_sources(args.repo_root), nodes=nodes)
        # The nodes alone decide the exit code; the reported sources never do (SOURCE_CHECKS).
        return 1 if any(s.is_red for s in statuses) else 0
    return run(nodes, args.repo_root, args.only)


if __name__ == "__main__":
    raise SystemExit(main())
