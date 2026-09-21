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
node is stale, has a missing output or has a missing input.  That is the gate.  It reads mtimes and
nothing else; it opens no parquet and writes no file.

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

WHAT COUNTS AS AN INPUT
-----------------------
The builder's own source file is declared as an input alongside its data.  A builder edited after
its artifact was written means the artifact no longer describes that builder, which is the same
kind of staleness as a newer data file and is reported the same way.  The table names which input
won, so code staleness and data staleness are told apart at a glance.

A directory input is scanned recursively and the newest file in it wins; a ``*`` in a declared path
is a glob.  ``typing/geometry`` holds two unrelated families written by two different builders
(``endocrine_union_*.parquet`` and ``nonendocrine_*.geojson``), so those are declared as globs and
not as the directory -- declaring the directory would make each builder its own input.

THE ORDERING HAZARD THIS ENCODES
--------------------------------
``build_parquet_views.py`` rewrites ``app/parquet/islets.parquet`` FROM SCRATCH and silently drops
the 21 ``crossed_*`` and ``sc_*`` columns that ``summarise_crossed_fraction.py`` and
``summarise_singlecell_trajectory.py`` merge in afterwards (docs/APP_COVERAGE_2026-09-17.md:200-212).
All three declare that file as an output, and ``build_graph`` REFUSES a graph in which two nodes
write one path without an ordering between them.  ``build_parquet_views`` therefore cannot be last:
the graph will not build unless ``crossed_fraction`` runs after it and ``singlecell_trajectory``
after that.

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
        "see core/phenocycler/scripts/typing/RECONSTRUCTION.md)",
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
        "exported from QuPath",
    "results/AnnotationsMeasurements.csv":
        "exported from QuPath",
    "/data/IO60pLN/results/node_geojson":
        "core/phenocycler/scripts/groovy/export_pln_node_annotations.groovy",
    "/data/islet_layers/pln/provenance/qupath_annotation_audit.log":
        "the captured stdout of export_pln_node_annotations.groovy; it lives under the pln "
        "artifact root but it is an input, not something build_pln_layer.py writes",
}


@dataclass(frozen=True)
class Node:
    """One builder, and the files it reads and writes.

    ``inputs`` and ``outputs`` are repo-relative or absolute paths; a ``*`` makes one a glob and a
    directory is scanned recursively.  ``after`` names nodes this one must follow for a reason the
    files do not express.  ``soft_inputs`` are files the builder opens that are deliberately NOT
    edges, each with its reason in ``why``.
    """

    name: str
    command: tuple[str, ...]
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    after: tuple[str, ...] = ()
    soft_inputs: tuple[str, ...] = ()
    why: str = ""


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
            why="The islet floor every per-islet artifact meets: 18,404 core outlines, 20 donors.",
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
                "{app}/annotations.tsv",
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
            command=("{py}", "{typingscripts}/trajectory_islet.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/trajectory_islet.py",
                "{app}/islet_celltype_proportions.csv",
                "{app}/donor_metadata.csv",
            ),
            outputs=("{app}/islet_leiden.parquet",),
            why="Leiden is written at :481, BEFORE the pseudotime gate at :598, so a rejected "
                "pooled ordering does not withhold the clustering. adata_ins_root.h5ad is NOT "
                "declared here: it has a second producer (scripts/rebuild_adata_ins_root.py) and "
                "which one owns it is another lane's decision.",
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
            ),
            outputs=(
                "{app}/parquet/islets.parquet",
                "{app}/parquet/uns_targets.parquet",
                "{app}/parquet/uns_composition.parquet",
            ),
            why="Rewrites islets.parquet FROM SCRATCH. The two nodes below merge 21 columns back "
                "in, so this one must never run last -- build_graph refuses a graph where it can.",
        ),
        dict(
            name="crossed_fraction",
            command=("{py}", "{typingscripts}/summarise_crossed_fraction.py", "--include-marginal"),
            inputs=(
                "{typingscripts}/summarise_crossed_fraction.py",
                "{typing}/crossed",
                "{typing}/cells",
                "{app}/parquet/islets.parquet",
            ),
            outputs=(
                "{typing}/crossed/crossed_islet_summary.parquet",
                "{app}/parquet/islets.parquet",
            ),
            after=("parquet_views",),
            why="Merges the 14 crossed_* columns back into islets.parquet after "
                "build_parquet_views.py has dropped them.",
        ),
        dict(
            name="singlecell_trajectory",
            command=("{py}", "{typingscripts}/summarise_singlecell_trajectory.py"),
            inputs=(
                "{typingscripts}/summarise_singlecell_trajectory.py",
                "{typing}/trajectory/singlecell_pseudotime.parquet",
                "{app}/parquet/islets.parquet",
            ),
            outputs=(
                "{typing}/trajectory/singlecell_islet_summary.parquet",
                "{app}/parquet/islets.parquet",
            ),
            after=("crossed_fraction",),
            why="Merges the 7 sc_* columns back in. Takes no --include-marginal; its cohort comes "
                "from singlecell_pseudotime.parquet.",
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
            name="app_donors",
            command=("{py}", "{typingscripts}/build_app_data.py", "--include-marginal",
                     "--only", "donors"),
            inputs=("{typingscripts}/build_app_data.py", "{typing}/cells"),
            outputs=("{app}/donors",),
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
                "{integration}/pln_node_annotations.csv",
                "{integration}/cells_xen",
                "{typing}/cells",
                "{typing}/geometry/endocrine_union_*.parquet",
                "{typing}/trajectory/singlecell_islet_summary.parquet",
                "{repo}/data/development/restain",
            ),
            outputs=("{layers}/availability/data.parquet",),
            why="The table every other card consults. Since 2026-09-18 it answers the two Xenium "
                "capabilities from the typed releases, resolved by scripts/senior/xenium_release.py; "
                "cells_xen, the UNTYPED ingest, is read only for the IF pairing and as the fallback "
                "when no release resolves. The declared inputs below still name only cells_xen and "
                "do not name the release roots, so this node's staleness check does not see a new "
                "release. Reported, not edited here.",
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
                "{repo}/data/development/xenium_typing",
            ),
            outputs=("{layers}/distance_profiles/data.parquet",),
            why="Carries the study question: the immune fraction against signed distance to the "
                "islet edge. Refuses subpanel_calls.parquet by path.",
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
                "{repo}/data/development/xenium_typing/release_pln_immune_v2",
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
                "supersedes _v2. The declared input below is still the superseded _v2 path, so "
                "the staleness check watches the wrong tree. Reported, not edited here. Not in "
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
                "{layers}/grid_composition/composition.parquet",
                "{layers}/grid_composition/composition.csv",
                "{layers}/grid_composition/fine.parquet",
                "{layers}/grid_composition/presence.parquet",
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


def _extremum(path: str, newest: bool) -> tuple[float, str] | None:
    """(mtime, the file carrying it) over a file, a directory tree or a glob; None if absent."""
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
        stats = [(c.stat().st_mtime, str(c)) for c in candidates if c.is_file()]
    elif target.is_dir():
        stats = [(e.stat().st_mtime, e.path) for e in _files_under(target)]
    elif target.is_file():
        stats = [(target.stat().st_mtime, str(target))]
    else:
        return None
    if not stats:
        return None
    return max(stats) if newest else min(stats)


def check_node(node: Node) -> NodeStatus:
    """Compare the node's oldest output against its newest input. Reads mtimes only."""
    missing_inputs = tuple(p for p in node.inputs if _extremum(p, newest=True) is None)
    if missing_inputs:
        return NodeStatus(node.name, MISSING_INPUT, missing=missing_inputs)
    missing_outputs = tuple(p for p in node.outputs if _extremum(p, newest=False) is None)
    if missing_outputs:
        newest = max((_extremum(p, newest=True) for p in node.inputs), default=None)
        return NodeStatus(node.name, MISSING_OUTPUT, missing=missing_outputs,
                          newest_input=newest[1] if newest else None,
                          newest_input_mtime=newest[0] if newest else None)
    newest = max(_extremum(p, newest=True) for p in node.inputs)
    oldest = min(_extremum(p, newest=False) for p in node.outputs)
    state = STALE if newest[0] > oldest[0] else FRESH
    return NodeStatus(node.name, state, newest_input=newest[1], newest_input_mtime=newest[0],
                      oldest_output=oldest[1], oldest_output_mtime=oldest[0])


def check(nodes: dict[str, Node], only: list[str] | None = None) -> list[NodeStatus]:
    """One NodeStatus per node, in topological order. The public entry point for Lane V."""
    order = topological_order(nodes)
    if only:
        order = [n for n in order if n in set(only)]
    return [check_node(nodes[name]) for name in order]


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


def print_check(statuses: list[NodeStatus], repo: Path) -> None:
    width = max(len(s.name) for s in statuses) + 2
    print(f"{'  node':<{width}}  {'state':<14} {'output':<16} {'newest input':<16} by")
    print("-" * (width + 90))
    for status in statuses:
        mark = "  " if status.state == FRESH else "* "
        by = _short(status.newest_input, repo)
        if status.state == MISSING_INPUT:
            by = "MISSING: " + ", ".join(_short(p, repo) for p in status.missing)
        elif status.state == MISSING_OUTPUT:
            by = f"absent: {', '.join(_short(p, repo) for p in status.missing)}"
        print(f"{mark + status.name:<{width}}  {status.state:<14} "
              f"{_stamp(status.oldest_output_mtime):<16} "
              f"{_stamp(status.newest_input_mtime):<16} {by}")
    red = [s for s in statuses if s.is_red]
    layers_red = [s.name for s in red if s.name in LAYER_NODES]
    print()
    print(f"{len(statuses) - len(red)} fresh, {len(red)} not fresh "
          f"({len(layers_red)} of the nine layers: {', '.join(layers_red) or 'none'})")
    named = set()
    for status in red:
        by = status.newest_input
        if not by:
            continue
        for source, producer in SOURCE_PRODUCERS.items():
            if source in by and source not in named:
                named.add(source)
                print(f"  {source} is written by {producer}")
    if red:
        print("\nStale means the artifact was built before an input it reads. It is not a claim "
              "about the numbers in it;\nit is a claim that they were computed from a different "
              "generation of the inputs than the ones on disk now.")


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
        missing = before[name].missing if before[name].state == MISSING_INPUT else ()
        if missing:
            reason = "a declared input does not exist: " + ", ".join(missing)
            print(f"\n== {name}: REFUSED -- {reason}")
            refused.append((name, reason))
            continue

        # A node whose output already postdates every input it reads has nothing to do.
        # `check()` computed exactly that a few lines above and the answer was thrown away:
        # before this, `--run` rebuilt every scheduled node, including the ones `--check`
        # had just reported fresh. `--only` is still honoured, so naming a fresh node
        # explicitly forces it -- that is the escape hatch for a node whose inputs are not
        # fully declared.
        if not before[name].is_red and not (only and name in set(only)):
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
        print_check(statuses, args.repo_root)
        return 1 if any(s.is_red for s in statuses) else 0
    return run(nodes, args.repo_root, args.only)


if __name__ == "__main__":
    raise SystemExit(main())
