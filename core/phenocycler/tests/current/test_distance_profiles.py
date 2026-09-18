"""UUID, denominator, availability, and app checks for the distance layer.

Run from the repository root with .venv/bin/python -m pytest <this file>.
Fixtures are synthetic measurements, never expected biological proportions.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import time
import subprocess
import sys
from pathlib import Path
from uuid import UUID

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
BUILDER = REPO / "scripts/senior/build_distance_profiles.py"
spec = importlib.util.spec_from_file_location("build_distance_profiles", BUILDER)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


@pytest.fixture
def sources(tmp_path):
    cells = tmp_path / "cells" / "donor_id=0112" / "data.parquet"
    members = tmp_path / "islet_membership" / "donor_id=unrelated_directory" / "data.parquet"
    metadata = tmp_path / "app" / "donor_metadata.csv"
    metadata.parent.mkdir()
    metadata.write_text("case_id,donor_status\n112,ND\n")
    classes = ["Beta", "Immune", "Immune", "Exocrine", "Immune", "Exocrine",
               "Immune", "Immune", "Excluded", "Immune"]
    definitions = ["", "Imm-T-cell", "", "", "Imm-B-cell", "", "Imm-T-cell",
                   "Imm-T-cell", "", "Imm-T-cell"]
    provenance = ["hormone", "definition", "supporting", "fallthrough", "definition",
                  "fallthrough", "definition", "definition", "fallthrough", "definition"]
    distances = [-12., -3., 2., 1., 20., 21., 0., -8., -9., -10.]
    cell_rows, member_rows = [], []
    for i, cell_class in enumerate(classes):
        object_id = str(UUID(int=i + 1))
        cell_rows.append(dict(
            object_id=object_id, cell_class=cell_class, immune_definition=definitions[i],
            resolved_by=provenance[i], qc_analysis_eligible=False if i == 7 else None if i == 9 else True,
            excluded_reason="exocrine_inside_islet_border" if i == 8 else "",
            # Deliberately wrong cached geometry and Parent must have no effect.
            Parent="Annotation (Tissue)", signed_dist_islet_um=9000.,
        ))
        member_rows.append(dict(object_id=object_id, in_islet=i < 3 or i >= 7,
                                signed_dist_islet_um=distances[i]))
    write_rows(cells, cell_rows)
    write_rows(members, list(reversed(member_rows)))
    # This real-tree hazard must not enter cell-dataset schema inference.
    write_rows(members.parents[1] / "membership_diff.parquet", [{"donor_id": "112", "gained": 1}])
    return dict(cells=cells, members=members, metadata=metadata,
                cell_rows=cell_rows, member_rows=member_rows)


def build(sources, **kwargs):
    return builder.build_distance_profiles(sources["cells"].parents[1],
                                            sources["members"].parents[1],
                                            sources["metadata"], **kwargs)


@pytest.mark.parametrize("root", [
    Path("data/development/typing/subpanel/subpanel_calls.parquet"),
    Path("/does/not/exist/typing/subpanel/subpanel_calls.parquet"),
])
def test_builder_refuses_superseded_path_with_measured_counts(root):
    """The acceptance check: rejection precedes even checking that input exists."""
    with pytest.raises(ValueError) as exc:
        builder.build_distance_profiles(cells_root=root)
    message = str(exc.value)
    for counts in ("Beta 191,485", "Alpha 246,055", "Delta 62,797",
                   "Beta 226,183", "Alpha 279,686", "Delta 90,658"):
        assert counts in message


def test_superseded_symlink_and_cli_are_also_refused(tmp_path):
    target = tmp_path / "subpanel_calls.parquet"
    target.touch()
    alias = tmp_path / "cells.parquet"
    alias.symlink_to(target)
    with pytest.raises(ValueError, match="Beta 191,485"):
        builder.build_distance_profiles(alias)
    result = subprocess.run([sys.executable, "-B", str(BUILDER), "--cells-root", str(target),
                             "--layers-dir", str(tmp_path / "output")],
                            capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert "Beta 191,485 / Alpha 246,055 / Delta 62,797" in result.stderr
    assert not (tmp_path / "output").exists()


def test_uuid_join_qc_provenance_and_exact_compartment_boundaries(sources):
    table = build(sources)
    data = table.to_pandas()
    assert set(data.donor_id) == {"112"}
    assert set(data.donor_status) == {"ND"}
    assert set(data.frame) == {"signed_distance_um"}
    assert set(data.state) == {"measured"}
    assert set(data.reason) == {""}
    total = data[data.resolved_by == "all"]
    assert total.n_cells.sum() == 7
    assert total.n_immune.sum() == 4
    assert set(total.distance_bin_um) == {-20., -10., 0., 10., 20., 30.}
    counts = total.groupby("compartment")[["n_cells", "n_immune"]].sum().to_dict("index")
    assert counts == {"islet_core": {"n_cells": 3, "n_immune": 2},
                      "peri_islet": {"n_cells": 2, "n_immune": 1},
                      "other_pancreas": {"n_cells": 2, "n_immune": 1}}
    assert set(data[data.resolved_by != "all"].resolved_by) == {"definition", "supporting"}
    assert data.loc[data.resolved_by == "supporting", "n_immune"].sum() == 1
    assert "Immune (subtype unresolved)" in set(data.immune_subtype)
    assert all(data.immune_fraction == data.n_immune / data.n_cells)
    # Zero immune cells in an observed bin is a measurement, not an absence.
    assert ((total.n_immune == 0) & (total.immune_fraction == 0)).any()
    audit = json.loads(table.schema.metadata[b"audit"])["112"]
    assert audit["omitted_qc_or_excluded"] == 3
    assert audit["membership_distance_disagreement"] == 2


def test_exclusions_use_central_cohort_after_normalizing(sources):
    from phenocycler.cohort import DONOR_EXCLUSIONS

    excluded = next(iter(DONOR_EXCLUSIONS))
    excluded_file = sources["cells"].parents[1] / f"donor_id=0{excluded}" / "data.parquet"
    # Intentionally malformed: an excluded donor must never be read for analysis.
    write_rows(excluded_file, [{"unrelated": 1}])
    assert set(build(sources)["donor_id"].to_pylist()) == {"112"}


@pytest.mark.parametrize("defect,reason", [
    ("missing_uuid", "islet_membership_missing_uuid:1_of_7"),
    ("missing_distance", "islet_membership_invalid_measurement:1_of_7"),
    ("infinite_distance", "islet_membership_invalid_measurement:1_of_7"),
    ("missing_typing", "typing_or_provenance_missing:1_of_7"),
    ("missing_status", "donor_status_not_recorded"),
    ("missing_export", "islet_membership_export_absent"),
])
def test_missing_measurement_is_named_and_never_zero(sources, defect, reason):
    members = sources["member_rows"]
    if defect == "missing_uuid":
        members = members[1:]
    elif defect == "missing_distance":
        members[0]["signed_dist_islet_um"] = None
    elif defect == "infinite_distance":
        members[0]["signed_dist_islet_um"] = float("inf")
    elif defect == "missing_typing":
        sources["cell_rows"][0]["resolved_by"] = None
        write_rows(sources["cells"], sources["cell_rows"])
    elif defect == "missing_status":
        sources["metadata"].unlink()
    write_rows(sources["members"], members)
    if defect == "missing_export":
        sources["members"].unlink()
    rows = build(sources).to_pylist()
    assert len(rows) == 1
    assert rows[0]["state"] == "unavailable"
    assert rows[0]["reason"] == reason
    for name in ("n_cells", "n_immune", "immune_fraction", "distance_bin_um"):
        assert rows[0][name] is None


@pytest.mark.parametrize("source", ["cells", "members"])
def test_duplicate_uuid_is_refused(sources, source):
    rows = sources["cell_rows" if source == "cells" else "member_rows"]
    write_rows(sources[source], rows + rows[:1])
    with pytest.raises(ValueError, match="duplicate object_id"):
        build(sources)


def test_cli_writes_only_one_layer_root(sources, tmp_path):
    layers = tmp_path / "layers"
    assert builder.main([
        "--cells-root", str(sources["cells"].parents[1]),
        "--membership-root", str(sources["members"].parents[1]),
        "--donor-metadata", str(sources["metadata"]), "--layers-dir", str(layers),
    ]) == 0
    assert list(layers.rglob("*.parquet")) == [layers / "distance_profiles" / "data.parquet"]
    result = pq.ParquetFile(layers / "distance_profiles" / "data.parquet").read()
    assert result.schema.metadata[b"source_domain_root"].decode() == str(sources["metadata"].parent)
    assert not list(layers.rglob("*.tmp"))


def test_cards_render_measured_and_absent_states_and_refuse_mixed_frames(sources, tmp_path):
    rscript = Path("/opt/R/4.5.1/bin/Rscript")
    if not rscript.is_file():
        discovered = shutil.which("Rscript")
        if discovered is None:
            pytest.skip("Rscript is not installed")
        rscript = Path(discovered)
    dependencies = subprocess.run([
        str(rscript), "--vanilla", "-e",
        'quit(status = if (all(vapply(c("shiny", "arrow", "ggplot2", "scales"), '
        'requireNamespace, logical(1), quietly = TRUE))) 0 else 1)',
    ], capture_output=True, text=True, check=False)
    if dependencies.returncode:
        pytest.skip("R card dependencies are not installed")
    layers = tmp_path / "layers"
    target = layers / "distance_profiles" / "data.parquet"
    target.parent.mkdir(parents=True)
    pq.write_table(build(sources), target)
    write_rows(layers / "availability" / "data.parquet", [
        dict(donor_id="0112", capability=c, state="measured", reason="", n=7, source="fixture")
        for c in ("phenocycler_panc", "islet_outline")
    ])
    env = dict(os.environ, ISLET_LAYERS_DIR=str(layers),
               DISTANCE_TEST_DOMAIN=str(sources["metadata"].parent))
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
source("apps/R/ieCore/R/ie_explain.R")
source("apps/R/ieCore/R/ie_frames.R")
source("apps/islet/R/card_x1_gradient.R")
source("apps/islet/R/card_s3_compartments.R")
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
table_text <- function(html) {
  trimws(gsub("[[:space:]]+", " ", gsub("<[^>]+>", " ", html)))
}
layer <- .ie_distance_read()
stopifnot(ie_is_measured(layer$value))
stopifnot(ie_is_measured(.ie_distance_select(layer, "112")$value))
stopifnot(identical(.ie_distance_donors(layer, list(donor_ids = function() "00112")), "112"))
stopifnot(length(.ie_distance_donors(layer, list(donor_ids = function() character()))) == 0)
for (state in setdiff(IE_STATES, "measured")) {
  absent <- layer
  absent$ledger$state[[1]] <- state
  absent$ledger$reason[[1]] <- "measured_fixture_reason"
  value <- .ie_distance_select(absent, "112")$value
  stopifnot(identical(value$state, state))
  stopifnot(grepl("guard_test", as.character(.ie_distance_note(value, "guard_test"))))
}
domain_data_root <- function() "/different/domain"
stopifnot(identical(.ie_distance_read()$value$reason, "distance_profile_domain_mismatch"))
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
shiny::testServer(function(input, output, session) {
  result <- card_x1_gradient_server(input, output, session, active_tab = shiny::reactive("CrossModal"))
}, {
  session$setInputs(x1_gradient_donor = "112", x1_gradient_pair = "Immune | all")
  stopifnot(sum(result$data()$n_cells) == 7, sum(result$data()$n_immune) == 4)
  plot <- result$plot()
  stopifnot(inherits(ggplot2::ggplot_build(plot), "ggplot_built"))
  stopifnot(0 %in% plot$data$distance_bin_um)
  # LINEAR in micrometres, not pseudo_log. The card no longer draws past X1_MAX_DIST_UM, and
  # pseudo_log spent most of the panel on distances it now excludes while compressing the
  # peri-islet band -- the part being read -- into the left third. Zero stays visible on a linear
  # axis on its own, which is the only thing pseudo_log was buying here.
  stopifnot(identical(plot$scales$get_scales("x")$trans$name, "identity"))
  # The cap is a hard ceiling, and the slider may only narrow it.
  stopifnot(X1_MAX_DIST_UM == 500)
  stopifnot(all(abs(plot$data$distance_bin_um) <= X1_MAX_DIST_UM))
  # No counts table under the plot: it restated the plotted data row for row.
  stopifnot(!grepl("x1_gradient_counts",
                   paste(as.character(output$x1_gradient_body$html), collapse = " "),
                   fixed = TRUE))
})
shiny::testServer(function(input, output, session) {
  result <- card_s3_compartments_server(input, output, session, active_tab = shiny::reactive("LymphNode"))
}, {
  session$setInputs(s3_compartments_donor = "112", s3_compartments_pair = "Immune | all")
  data <- result$data()
  stopifnot(sum(data$n_cells) == 7, sum(data$n_immune) == 4)
  stopifnot(data$immune_fraction[data$compartment == "Islet core"] == 2/3)
  stopifnot(data$immune_fraction[data$compartment == "Peri-islet expansion"] == 1/2)
  stopifnot(data$immune_fraction[data$compartment == "Other pancreas"] == 1/2)
  built <- ggplot2::ggplot_build(result$plot())
  stopifnot(nrow(built$data[[1]]) == 3)
  stopifnot(identical(sort(built$data[[2]]$label), c("1 / 2", "1 / 2", "2 / 3")))
  counts <- table_text(output$s3_compartments_counts)
  for (text in c("Eligible cells", "Immune cells", "Islet core 3 2",
                 "Peri-islet expansion 2 1", "Other pancreas 2 1")) {
    stopifnot(grepl(text, counts, fixed = TRUE))
  }
  stopifnot(grepl("unavailable(guard_pln_export_absent)",
                 output$s3_compartments_pln$html, fixed = TRUE))
  # A subtype changes only the numerator, including a measured zero in other pancreas.
  session$setInputs(s3_compartments_pair = "Imm-B-cell | definition")
  data <- result$data()
  stopifnot(sum(data$n_cells) == 7, sum(data$n_immune) == 1)
  stopifnot(data$immune_fraction[data$compartment == "Other pancreas"] == 0)
  stopifnot(grepl("Other pancreas 2 0", table_text(output$s3_compartments_counts), fixed = TRUE))
})
actual_reader <- .ie_distance_read
.ie_distance_read <- function() {
  value <- actual_reader()
  value$data <- value$data[value$data$compartment == "other_pancreas", , drop = FALSE]
  value
}
shiny::testServer(function(input, output, session) {
  result <- card_s3_compartments_server(input, output, session, active_tab = shiny::reactive("LymphNode"))
}, {
  session$setInputs(s3_compartments_donor = "112", s3_compartments_pair = "Immune | all")
  data <- result$data()
  stopifnot(is.data.frame(data), nrow(data) == 1, as.character(data$compartment) == "Other pancreas")
  stopifnot(data$n_cells == 2, data$n_immune == 1)
  stopifnot(inherits(ggplot2::ggplot_build(result$plot()), "ggplot_built"))
  stopifnot(grepl("Other pancreas 2 1", table_text(output$s3_compartments_counts), fixed = TRUE))
  for (name in c("Islet core", "Peri-islet expansion")) {
    stopifnot(grepl(paste0("no_eligible_cells_in_", name), output$s3_compartments_body$html, fixed = TRUE))
  }
  stopifnot(grepl("guard_s3_compartment_empty", output$s3_compartments_body$html, fixed = TRUE))
})
.ie_distance_read <- function() {
  value <- actual_reader()
  value$data <- value$data[value$data$compartment != "other_pancreas", , drop = FALSE]
  value
}
shiny::testServer(function(input, output, session) {
  result <- card_s3_compartments_server(input, output, session, active_tab = shiny::reactive("LymphNode"))
}, {
  session$setInputs(s3_compartments_donor = "112", s3_compartments_pair = "Immune | all")
  stopifnot(sum(result$data()$n_cells) == 5, sum(result$data()$n_immune) == 3)
  stopifnot(grepl("no_eligible_cells_in_Other pancreas", output$s3_compartments_body$html, fixed = TRUE))
  stopifnot(grepl("guard_s3_compartment_empty", output$s3_compartments_body$html, fixed = TRUE))
  stopifnot(!grepl("Other pancreas", table_text(output$s3_compartments_counts), fixed = TRUE))
})
.ie_distance_read <- function() {
  bad <- actual_reader()
  bad$data$frame[which(bad$data$resolved_by == "all")[[1]]] <- "phenocycler_um:112:panc"
  bad
}
shiny::testServer(function(input, output, session) {
  result <- card_x1_gradient_server(input, output, session, active_tab = shiny::reactive("CrossModal"))
}, {
  session$setInputs(x1_gradient_donor = "112", x1_gradient_pair = "Immune | all")
  error <- tryCatch({result$plot(); NULL}, error = identity)
  stopifnot(inherits(error, "error"), grepl("frame", conditionMessage(error)))
})
.ie_distance_read <- actual_reader
Sys.setenv(ISLET_LAYERS_DIR = file.path(Sys.getenv("ISLET_LAYERS_DIR"), "absent"))
shiny::testServer(function(input, output, session) {
  card_x1_gradient_server(input, output, session, active_tab = shiny::reactive("CrossModal"))
}, {
  stopifnot(grepl("distance_profiles_export_absent", output$x1_gradient_body$html))
})
cat("R distance cards: measured, unavailable, domain and frame checks passed\n")
'''
    result = subprocess.run([str(rscript), "--vanilla", "-"], input=code, cwd=REPO,
                            env=env, capture_output=True, text=True, check=False, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr


# --------------------------------------------------------------------------- #
# The pLN arm: a recorded immune vocabulary, and nothing spatial
# --------------------------------------------------------------------------- #

IMMUNE_HEADER = ("tissue,kind,common,phenocycler,xenium,pheno_resolution,xenium_resolution,"
                 "pheno_note,xenium_note,collapsed,identity_core_genes,is_immune,"
                 "immune_tissue_scope,immune_note\n")
IMMUNE_SCOPE = "pancreas|pancreatic_lymph_node"


def write_crosswalk(path, lineages, scope=IMMUNE_SCOPE, header=IMMUNE_HEADER):
    """A crosswalk fixture: {lineage: is_immune}. Never a copy of the real file."""
    lines = [header]
    for lineage, value in lineages.items():
        lines.append(f"pancreas,lineage,{lineage},,,yes,yes,,,,0,{value},"
                     f"{scope if value == 'yes' else ''},ruling 2026-09-17\n")
    lines.append("pancreas,marker,,CD3e,CD3E,TYPE,on_panel,,,,,not_applicable,,\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines))
    return path


RECORDED = {"Endocrine": "no", "Exocrine": "no", "Stromal": "no", "Vascular": "no",
            "Neural": "no", "T_NK": "yes", "B_Plasma": "yes", "Myeloid": "yes",
            "Granulocyte": "yes", "Unassigned": "unresolvable"}


@pytest.fixture
def pln(tmp_path):
    """One pLN section for donor 112, fully called. Counts are synthetic."""
    calls = (["T_NK"] * 4 + ["B_Plasma"] * 3 + ["Myeloid"] * 2 + ["Stromal"] * 1)
    path = tmp_path / "cells_xen" / "donor_id=0112" / "roi=pln_1" / "data.parquet"
    write_rows(path, [dict(cell_id=f"c{i}", lineage_common=call,
                           tissue="pancreatic_lymph_node")
                      for i, call in enumerate(calls)])
    return dict(root=path.parents[2], calls=calls)


def build_pln(sources, pln, crosswalk, **kwargs):
    return build(sources, xen_cells_root=pln["root"], vocab_crosswalk=crosswalk, **kwargs)


def pln_rows(table):
    return [r for r in table.to_pylist() if r["roi"] != "panc"]


def test_pln_immune_fraction_is_measured_once_the_vocabulary_is_recorded(sources, pln, tmp_path):
    table = build_pln(sources, pln, write_crosswalk(tmp_path / "vocab.csv", RECORDED))
    rows = pln_rows(table)
    total = [r for r in rows if r["resolved_by"] == "all"]
    assert len(total) == 1
    assert total[0]["state"] == "measured" and total[0]["reason"] == ""
    # 4 T_NK + 3 B_Plasma + 2 Myeloid immune of 10 cells; Stromal is recorded not-immune.
    assert (total[0]["n_cells"], total[0]["n_immune"]) == (10, 9)
    assert total[0]["immune_fraction"] == 0.9
    by_lineage = {r["immune_subtype"]: r["n_immune"] for r in rows if r["resolved_by"] != "all"}
    assert by_lineage == {"T_NK": 4, "B_Plasma": 3, "Myeloid": 2}
    assert all(r["resolved_by"] == "lineage_common" for r in rows if r["resolved_by"] != "all")
    assert set(json.loads(table.schema.metadata[b"pln_audit"].decode())) == {"112/pln_1"}


def test_a_measured_pln_row_still_carries_nothing_spatial(sources, pln, tmp_path):
    rows = pln_rows(build_pln(sources, pln, write_crosswalk(tmp_path / "vocab.csv", RECORDED)))
    assert rows
    for row in rows:
        assert row["frame"] == "none"
        assert row["distance_bin_um"] is None
        assert row["compartment"] is None
        assert row["tissue"] == "pancreatic_lymph_node"
        assert row["modality"] == "xenium"


def test_an_unrecorded_lineage_is_withheld_by_name_not_treated_as_non_immune(sources, pln, tmp_path):
    partial = {k: v for k, v in RECORDED.items() if k != "Stromal"}
    rows = pln_rows(build_pln(sources, pln, write_crosswalk(tmp_path / "vocab.csv", partial)))
    assert len(rows) == 1
    assert rows[0]["state"] == "withheld"
    assert rows[0]["reason"] == "xenium_immune_lineage_not_recorded:Stromal"
    # Withheld means withheld: no numerator, no denominator, no zero.
    for name in ("n_cells", "n_immune", "immune_fraction"):
        assert rows[0][name] is None


def test_an_unclassifiable_lineage_is_blank_and_never_counted_either_way(sources, pln, tmp_path):
    path = pln["root"] / "donor_id=0112" / "roi=pln_1" / "data.parquet"
    write_rows(path, [dict(cell_id="c0", lineage_common="Unassigned",
                           tissue="pancreatic_lymph_node")])
    rows = pln_rows(build_pln(sources, pln, write_crosswalk(tmp_path / "vocab.csv", RECORDED)))
    assert (rows[0]["state"], rows[0]["reason"]) == ("blank", "xenium_lineage_unresolvable:1_of_1")
    assert rows[0]["n_immune"] is None


def test_the_recorded_tissue_scope_is_what_licenses_the_lymph_node(sources, pln, tmp_path):
    """Scope is read, never assumed: pancreas-only scope withholds on a pLN section."""
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED, scope="pancreas")
    rows = pln_rows(build_pln(sources, pln, crosswalk))
    assert rows[0]["state"] == "withheld"
    assert rows[0]["reason"] == (
        "xenium_immune_vocabulary_not_recorded_for_tissue:pancreatic_lymph_node")
    # The same section with the ruling's scope recorded is measured.
    assert pln_rows(build_pln(sources, pln, write_crosswalk(
        tmp_path / "scoped.csv", RECORDED)))[0]["state"] == "measured"


@pytest.mark.parametrize("crosswalk_state", ["absent", "no_column"])
def test_the_unrecorded_vocabulary_branch_stays_in_place(sources, pln, tmp_path, crosswalk_state):
    if crosswalk_state == "absent":
        crosswalk = tmp_path / "missing.csv"
    else:
        # The crosswalk as it stood before the ruling was recorded: lineages
        # named, is-immune not among the columns.
        crosswalk = write_crosswalk(tmp_path / "old.csv", RECORDED)
        kept = [line.rsplit(",", 3)[0] for line in crosswalk.read_text().splitlines()]
        assert "is_immune" not in kept[0]
        crosswalk.write_text("\n".join(kept) + "\n")
    rows = pln_rows(build_pln(sources, pln, crosswalk))
    assert rows[0]["state"] == "withheld"
    assert rows[0]["reason"] == "xenium_immune_vocabulary_not_recorded"


def test_the_shipped_crosswalk_records_exactly_the_four_immune_lineages():
    """The file the builder reads by default, read the way the builder reads it."""
    vocabulary = builder.load_immune_vocabulary()
    assert vocabulary is not None, f"{builder.VOCAB_CROSSWALK} records no is_immune column"
    assert {lin for lin, v in vocabulary["is_immune"].items() if v == "yes"} == {
        "T_NK", "B_Plasma", "Myeloid", "Granulocyte"}
    assert vocabulary["is_immune"]["Unassigned"] == "unresolvable"
    assert vocabulary["scope"] == {"pancreas", "pancreatic_lymph_node"}


# --------------------------------------------------------------------------- #
# The last gate, exercised directly on hand-made rows
# --------------------------------------------------------------------------- #
#
# The builder cannot itself produce a pLN row carrying a distance, so going through
# build() can only ever show the gate not firing. These construct the rows by hand,
# which is how a reviewer demonstrated that the gate refused a pLN row carrying only
# a COUNT -- making the compositional branch, the one thing the pLN arm is allowed to
# emit, unreachable by construction.

def gate_row(**overrides):
    """One artifact row. Defaults are the compositional pLN shape: nothing spatial."""
    row = dict(donor_id="112", donor_status="ND", roi="pln_1",
               tissue="pancreatic_lymph_node", modality="xenium",
               distance_bin_um=None, n_cells=5, n_immune=4, immune_fraction=0.8,
               immune_subtype="", resolved_by="all", frame="none", compartment=None,
               state="measured", reason="")
    row.update(overrides)
    return pa.Table.from_pylist([row], schema=builder.SCHEMA)


def test_a_pln_row_carrying_only_a_count_passes_the_gate():
    """The reviewer's row, exactly: measured, five cells, no spatial quantity."""
    builder._assert_no_cross_tissue_distance(
        gate_row(distance_bin_um=None, compartment=None, frame="none",
                 state="measured", n_cells=5))


@pytest.mark.parametrize("state", ["measured", "withheld", "blank", "absent"])
def test_the_gate_judges_the_spatial_quantity_not_the_state(state):
    """State is not what the gate is about; a pLN row in any state may carry counts."""
    builder._assert_no_cross_tissue_distance(gate_row(state=state))


@pytest.mark.parametrize("spatial", [
    dict(distance_bin_um=20.0),
    dict(distance_bin_um=-10.0),
    dict(distance_bin_um=0.0),          # zero is its own bin, not the absence of one
    dict(compartment="Islet_N"),
    dict(compartment="Islet_N_exp20um"),
    dict(frame="signed_distance"),
])
def test_a_pln_row_carrying_a_spatial_quantity_is_still_refused(spatial):
    with pytest.raises(ValueError) as excinfo:
        builder._assert_no_cross_tissue_distance(gate_row(**spatial))
    message = str(excinfo.value)
    assert "pln_1" in message and "pancreatic_lymph_node" in message


def test_the_pancreas_arm_still_carries_distances_through_the_gate():
    """The gate is about the tissue boundary, not about distances as such."""
    builder._assert_no_cross_tissue_distance(
        gate_row(roi="panc", tissue="pancreas", modality="phenocycler",
                 distance_bin_um=-20.0, compartment="Islet_N", frame="signed_distance"))


# --------------------------------------------------------------------------- #
# Every recorded source path, and a count that describes what the card draws
# --------------------------------------------------------------------------- #

# The real layout, as strings: the run cohort's domain root and the legacy CODEX
# root are not in one another's tree, which is what the containment check turns
# on. Nothing is created on disk -- the reader compares recorded paths, and the
# active root arrives through the environment.
DOMAIN_TREE = "/repo/data/development/typing"
DOMAIN_ROOT = DOMAIN_TREE + "/app"
LEGACY_ROOT = "/repo/data/app_data"


def patch_metadata(table, **overrides):
    """Re-record the artifact's source paths at the real layout, then override."""
    metadata = {k.decode(): v.decode() for k, v in (table.schema.metadata or {}).items()}
    metadata.update({"source_domain_root": DOMAIN_ROOT,
                     "source_cells": DOMAIN_TREE + "/cells",
                     "source_membership": DOMAIN_TREE + "/islet_membership"})
    for key, value in overrides.items():
        metadata.pop(key) if value is None else metadata.update({key: value})
    return table.replace_schema_metadata(metadata)


def layers_dir(tmp_path, name, table):
    root = tmp_path / "layers" / name
    target = root / "distance_profiles" / "data.parquet"
    target.parent.mkdir(parents=True)
    pq.write_table(table, target)
    write_rows(root / "availability" / "data.parquet", [
        dict(donor_id="0112", capability=c, state="measured", reason="", n=7, source="fixture")
        for c in ("phenocycler_panc", "islet_outline")
    ])
    return root


@pytest.fixture
def pln_two_donors(tmp_path):
    """Donor 112 (in the cohort) and 6347 (a pLN section with no pancreas here)."""
    root = tmp_path / "cells_xen"
    write_rows(root / "donor_id=0112" / "roi=pln_1" / "data.parquet",
               [dict(cell_id=f"a{i}", lineage_common=call, tissue="pancreatic_lymph_node")
                for i, call in enumerate(["T_NK"] * 4 + ["B_Plasma"] * 3 + ["Stromal"] * 3)])
    write_rows(root / "donor_id=6347" / "roi=pln_1" / "data.parquet",
               [dict(cell_id=f"b{i}", lineage_common="T_NK", tissue="pancreatic_lymph_node")
                for i in range(5)])
    return dict(root=root)


def r_or_skip():
    rscript = Path("/opt/R/4.5.1/bin/Rscript")
    if not rscript.is_file():
        discovered = shutil.which("Rscript")
        if discovered is None:
            pytest.skip("Rscript is not installed")
        rscript = Path(discovered)
    probe = subprocess.run([
        str(rscript), "--vanilla", "-e",
        'quit(status = if (all(vapply(c("shiny", "arrow", "ggplot2", "scales"), '
        'requireNamespace, logical(1), quietly = TRUE))) 0 else 1)',
    ], capture_output=True, text=True, check=False)
    if probe.returncode:
        pytest.skip("R card dependencies are not installed")
    return rscript


#: (metadata override, the reason the reader must give). Each names the source
#: that escaped, so the card's note says which one -- not just "mismatch".
SOURCE_DEFECTS = [
    ({"source_cells": LEGACY_ROOT + "/cells"},
     "distance_profile_source_cells_outside_domain"),
    ({"source_membership": LEGACY_ROOT + "/islet_membership"},
     "distance_profile_source_membership_outside_domain"),
    ({"source_cells": ""}, "distance_profile_source_cells_not_recorded"),
    ({"source_cells": None}, "distance_profile_source_cells_not_recorded"),
    ({"source_membership": None}, "distance_profile_source_membership_not_recorded"),
    ({"source_cells": "data/development/typing/cells"},
     "distance_profile_source_cells_not_absolute"),
    ({"source_xen_cells": ""}, "distance_profile_source_xen_cells_not_recorded"),
    ({"source_xen_cells": None}, "distance_profile_source_xen_cells_not_recorded"),
    ({"source_immune_vocabulary": ""},
     "distance_profile_source_immune_vocabulary_not_recorded"),
    ({"source_immune_vocabulary": None},
     "distance_profile_source_immune_vocabulary_not_recorded"),
]


def test_cards_validate_every_source_path_not_only_the_domain_root(
        sources, pln_two_donors, tmp_path):
    """Defect 3: a correct source_domain_root beside cells from the other cohort.

    The reader used to check only source_domain_root, so an artifact that named
    the right domain and resolved its cells in the legacy CODEX root passed --
    the wrong cells for a right-looking islet.
    """
    rscript = r_or_skip()
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED)
    table = build(sources, xen_cells_root=pln_two_donors["root"], vocab_crosswalk=crosswalk)
    good = layers_dir(tmp_path, "good", patch_metadata(table))
    cases = []
    for index, (override, reason) in enumerate(SOURCE_DEFECTS):
        root = layers_dir(tmp_path, f"bad{index}", patch_metadata(table, **override))
        cases.append((str(root), reason))
    env = dict(os.environ, ISLET_LAYERS_DIR=str(good), DISTANCE_TEST_DOMAIN=DOMAIN_ROOT,
               DISTANCE_TEST_CASES=json.dumps(cases))
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
source("apps/R/ieCore/R/ie_explain.R")
source("apps/R/ieCore/R/ie_frames.R")
source("apps/islet/R/card_x1_gradient.R")
source("apps/islet/R/card_s3_compartments.R")
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
cohort <- list(donor_ids = function() c("00112"))

# The artifact whose every recorded source is inside the active domain tree.
layer <- .ie_distance_read()
stopifnot(ie_is_measured(layer$value))

# Defect 4: n describes the rows the card draws, not the whole artifact.
stopifnot(identical(layer$value$n, nrow(layer$data)))
stopifnot(layer$n_pln_rows > 0, nrow(layer$pln) == layer$n_pln_rows)
stopifnot(layer$n_artifact_rows == nrow(layer$data) + layer$n_pln_rows)
stopifnot(layer$value$n < layer$n_artifact_rows)  # they really do differ here

# Defect 2: the pancreas arm narrowed to the cohort, and the pLN-only donor
# named as withheld rather than unioned in and then intersected away.
stopifnot(identical(.ie_distance_donors(layer, cohort), "112"))
stopifnot(identical(.ie_distance_pln_outside_cohort(layer, cohort), "6347"))
stopifnot(!"6347" %in% .ie_distance_donors(layer, cohort))
note <- as.character(.ie_pln_outside_cohort_note(layer, cohort))
stopifnot(grepl("pln_donor_outside_active_domain_cohort", note, fixed = TRUE))
stopifnot(grepl("6347", note, fixed = TRUE))
# With no cohort to narrow against, nothing is claimed to be outside one.
stopifnot(length(.ie_distance_pln_outside_cohort(layer, NULL)) == 0)
stopifnot(is.null(.ie_pln_outside_cohort_note(layer, NULL)))

# Defect 3: every recorded source path, each refused by its own name.
cases <- jsonlite::fromJSON(Sys.getenv("DISTANCE_TEST_CASES"), simplifyVector = FALSE)
for (case in cases) {
  Sys.setenv(ISLET_LAYERS_DIR = case[[1]])
  value <- .ie_distance_read()$value
  if (!identical(value$reason, case[[2]])) {
    stop(sprintf("expected %s, got %s(%s)", case[[2]], value$state, value$reason))
  }
  stopifnot(identical(value$state, "failed"))
}
Sys.setenv(ISLET_LAYERS_DIR = Sys.getenv("DISTANCE_TEST_GOOD"))
cat("R source-path, donor and count checks passed\n")
'''
    env["DISTANCE_TEST_GOOD"] = str(good)
    result = subprocess.run([str(rscript), "--vanilla", "-"], input=code, cwd=REPO,
                            env=env, capture_output=True, text=True, check=False, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_pln_column_names_the_donor_it_will_not_show(sources, pln_two_donors, tmp_path):
    """The card renders the withheld pLN donor rather than dropping it silently."""
    rscript = r_or_skip()
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED)
    table = build(sources, xen_cells_root=pln_two_donors["root"], vocab_crosswalk=crosswalk)
    root = layers_dir(tmp_path, "card", patch_metadata(table))
    env = dict(os.environ, ISLET_LAYERS_DIR=str(root), DISTANCE_TEST_DOMAIN=DOMAIN_ROOT)
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
source("apps/R/ieCore/R/ie_explain.R")
source("apps/R/ieCore/R/ie_frames.R")
source("apps/islet/R/card_x1_gradient.R")
source("apps/islet/R/card_s3_compartments.R")
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
shiny::testServer(function(input, output, session) {
  result <- card_s3_compartments_server(input, output, session,
                                        active_tab = shiny::reactive("LymphNode"),
                                        cohort = list(donor_ids = function() c("00112")))
}, {
  session$setInputs(s3_compartments_donor = "112", s3_compartments_pair = "Immune | all")
  html <- output$s3_compartments_pln$html
  stopifnot(grepl("pln_donor_outside_active_domain_cohort", html, fixed = TRUE))
  stopifnot(grepl("6347", html, fixed = TRUE))
  # Donor 112's own pLN section is compositional and measured: 7 immune of 10.
  stopifnot(grepl("Compositional only", html, fixed = TRUE))
  stopifnot(grepl("T_NK", html, fixed = TRUE))
  # And it never grows a distance axis.
  for (word in c("distance_bin", "signed_distance", "Islet core")) {
    stopifnot(!grepl(word, html, fixed = TRUE))
  }
})
cat("R pLN column checks passed\n")
'''
    result = subprocess.run([str(rscript), "--vanilla", "-"], input=code, cwd=REPO,
                            env=env, capture_output=True, text=True, check=False, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr


def test_a_shared_source_is_required_only_where_a_number_depended_on_it(
        sources, pln_two_donors, tmp_path):
    """The pancreas arm is not refused over a pLN-arm metadata gap.

    The artifact built on 2026-09-17 records no source_immune_vocabulary and its
    pLN rows are all blank -- the cells were never typed, so no number came from
    the vocabulary. Demanding it there would make the whole card unavailable,
    pancreas arm included. Demanding it once a pLN row is measured is the case
    where a number really did come from it.
    """
    rscript = r_or_skip()
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED)
    # Uncalled pLN cells: the rows exist and are blank, exactly as shipped.
    uncalled = tmp_path / "uncalled"
    write_rows(uncalled / "donor_id=0112" / "roi=pln_1" / "data.parquet",
               [dict(cell_id=f"u{i}", lineage_common="",
                     tissue="pancreatic_lymph_node") for i in range(6)])
    blank = build(sources, xen_cells_root=uncalled, vocab_crosswalk=crosswalk)
    assert {r["state"] for r in blank.to_pylist() if r["roi"] != "panc"} == {"blank"}
    measured = build(sources, xen_cells_root=pln_two_donors["root"],
                     vocab_crosswalk=crosswalk)
    assert "measured" in {r["state"] for r in measured.to_pylist() if r["roi"] != "panc"}

    no_vocab = dict(source_immune_vocabulary=None)
    roots = {
        "blank_no_vocab": (layers_dir(tmp_path, "b1", patch_metadata(blank, **no_vocab)),
                           "measured"),
        "measured_no_vocab": (layers_dir(tmp_path, "b2", patch_metadata(measured, **no_vocab)),
                              "distance_profile_source_immune_vocabulary_not_recorded"),
        # source_xen_cells is load-bearing whatever state the rows are in.
        "blank_no_xen": (layers_dir(tmp_path, "b3",
                                    patch_metadata(blank, source_xen_cells=None)),
                         "distance_profile_source_xen_cells_not_recorded"),
    }
    cases = [[str(root), expected] for root, expected in roots.values()]
    env = dict(os.environ, ISLET_LAYERS_DIR=str(roots["blank_no_vocab"][0]),
               DISTANCE_TEST_DOMAIN=DOMAIN_ROOT, DISTANCE_TEST_CASES=json.dumps(cases))
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
source("apps/R/ieCore/R/ie_explain.R")
source("apps/R/ieCore/R/ie_frames.R")
source("apps/islet/R/card_x1_gradient.R")
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
cases <- jsonlite::fromJSON(Sys.getenv("DISTANCE_TEST_CASES"), simplifyVector = FALSE)
for (case in cases) {
  Sys.setenv(ISLET_LAYERS_DIR = case[[1]])
  layer <- .ie_distance_read()
  if (identical(case[[2]], "measured")) {
    if (!ie_is_measured(layer$value)) {
      stop(sprintf("expected measured, got %s(%s)", layer$value$state, layer$value$reason))
    }
    # The pancreas arm is intact and counted on its own.
    stopifnot(nrow(layer$data) > 0, identical(layer$value$n, nrow(layer$data)))
  } else if (!identical(layer$value$reason, case[[2]])) {
    stop(sprintf("expected %s, got %s(%s)", case[[2]], layer$value$state, layer$value$reason))
  }
}
cat("R shared-source dependence checks passed\n")
'''
    result = subprocess.run([str(rscript), "--vanilla", "-"], input=code, cwd=REPO,
                            env=env, capture_output=True, text=True, check=False, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr


# --------------------------------------------------------------------------- #
# The pLN arm points at typed calls, and keeps the two modalities apart
# --------------------------------------------------------------------------- #

def write_release(root, donors, *, tissue="pln", lineage_column="common_lineage",
                  calls_dir="calls", release_id=None, generated_at="2026-09-17T20:52:22+00:00",
                  roi="pln_1", supersedes=None):
    """A Xenium pLN release fixture. {donor: [(stage1, lineage, nucleus_count)]}.

    It declares a ``release_id`` and a write time because that is what the shared
    resolver selects on. A release that declares neither is not a candidate, which
    is deliberate: a directory mtime is not a release's own word for when it was
    written, and a ``touch`` must not be able to reorder the cohort's sources.
    """
    root.mkdir(parents=True, exist_ok=True)
    declared = {"schema_version": 1, "modality": "xenium", "tissue": tissue,
                "calls_dir": calls_dir, "cell_key": "cell_id", "roi": [roi],
                "release_id": release_id or f"fixture-{root.name}",
                "generated_at_utc": generated_at}
    if supersedes is not None:
        declared["producer"] = {"supersedes": supersedes}
    (root / "release.json").write_text(json.dumps(declared))
    for donor, cells in donors.items():
        write_rows(root / calls_dir / f"donor_id={donor}" / f"roi={roi}" / "data.parquet",
                   [{"cell_id": f"{donor}-{i}", "stage1": stage1,
                     lineage_column: lineage, "nucleus_count": nucleus}
                    for i, (stage1, lineage, nucleus) in enumerate(cells)])
    return root


#: One section: 4 immune (3 resolved to a lineage, 1 not), 2 non-immune,
#: 2 ambiguous, and 1 cell with no nucleus.
SECTION = ([("Immune", "T_NK", 1)] * 2 + [("Immune", "B_Plasma", 1)]
           + [("Immune", "", 1)] + [("Non_immune", "", 1)] * 2
           + [("Ambiguous", "", 1)] * 2 + [("Immune", "T_NK", 0)])


def test_the_xenium_pln_arm_reads_a_typed_release_not_the_untyped_ingest(
        sources, tmp_path):
    """cells_xen at roi=pln_1 carries no call on any row; the release does."""
    typing_root = tmp_path / "xenium_typing"
    write_release(typing_root / "release_pln_immune", {"0112": SECTION})
    # A pancreas release beside it must not be mistaken for the pLN one.
    write_release(typing_root / "release", {"0112": SECTION}, tissue="pancreas")
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED)
    table = build(sources, xenium_roots=[typing_root], vocab_crosswalk=crosswalk)
    rows = [r for r in table.to_pylist() if r["roi"] != "panc"]
    assert rows and all(r["modality"] == "xenium" for r in rows)

    total = [r for r in rows if r["resolved_by"] == "all"]
    assert len(total) == 1 and total[0]["state"] == "measured"
    # Denominator: stage 1 answered for 6 nucleated cells (4 Immune + 2
    # Non_immune). The 2 Ambiguous are excluded, and so is the cell with no
    # nucleus -- including the one that carried a call.
    assert (total[0]["n_cells"], total[0]["n_immune"]) == (6, 4)
    assert total[0]["immune_fraction"] == 4 / 6

    # Stage 2's lineages do not sum to the immune total: one immune cell was
    # never resolved to a lineage, and that is not a reason to drop it.
    lineages = {r["immune_subtype"]: r["n_immune"]
                for r in rows if r["resolved_by"] == "identity_core"}
    assert lineages == {"T_NK": 2, "B_Plasma": 1}
    assert sum(lineages.values()) < total[0]["n_immune"]

    # The cells stage 1 could not answer for are named, never absorbed.
    unresolved = [r for r in rows if r["resolved_by"] == "stage1_unresolved"]
    assert len(unresolved) == 1
    assert unresolved[0]["state"] == "blank"
    assert unresolved[0]["reason"] == "xenium_pln_stage1_unresolved:2_of_8"
    for name in ("n_cells", "n_immune", "immune_fraction"):
        assert unresolved[0][name] is None

    # Still nothing spatial, whatever the source.
    for row in rows:
        assert row["frame"] == "none"
        assert row["distance_bin_um"] is None and row["compartment"] is None
        assert row["tissue"] == "pancreatic_lymph_node"


def test_a_release_supersedes_the_ingest_and_the_ingest_is_the_fallback(
        sources, pln, tmp_path):
    """Both roots given: the typed release wins, and its source is recorded."""
    typing_root = tmp_path / "xenium_typing"
    release = write_release(typing_root / "release_pln_immune", {"0112": SECTION})
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED)
    table = build(sources, xenium_roots=[typing_root], xen_cells_root=pln["root"],
                  vocab_crosswalk=crosswalk)
    rows = [r for r in table.to_pylist() if r["roi"] != "panc"]
    total = [r for r in rows if r["resolved_by"] == "all"][0]
    assert (total["n_cells"], total["n_immune"]) == (6, 4)   # the release, not the ingest
    sources_by_modality = json.loads(
        table.schema.metadata[b"source_pln_by_modality"].decode())
    assert sources_by_modality["xenium"] == str(release.resolve())

    # With no release present, the ingest is read and named as the source.
    fallback = build(sources, xenium_roots=[tmp_path / "empty"],
                     xen_cells_root=pln["root"], vocab_crosswalk=crosswalk)
    fallback_rows = [r for r in fallback.to_pylist() if r["roi"] != "panc"]
    assert [r for r in fallback_rows if r["resolved_by"] == "all"][0]["n_cells"] == 10
    assert json.loads(fallback.schema.metadata[b"source_pln_by_modality"].decode())[
        "xenium"] == str(pln["root"].resolve())


def test_the_two_pln_modalities_stay_distinct_and_are_never_pooled(sources, tmp_path):
    typing_root = tmp_path / "xenium_typing"
    write_release(typing_root / "release_pln_immune", {"0112": SECTION})
    pheno = tmp_path / "pln" / "pheno"
    write_rows(pheno / "donor_id=0112" / "roi=pln_1" / "data.parquet",
               [{"cell_id": f"p{i}", "stage1": stage1, "common_lineage": lineage,
                 "nucleus_count": nucleus}
                for i, (stage1, lineage, nucleus) in enumerate(
                    [("Immune", "Myeloid", 1)] * 3 + [("Non_immune", "", 1)])])
    # The lane writes reference tables beside the calls; they are not cells.
    write_rows(pheno / "_reference" / "pancreas_reference_levels.parquet",
               [{"marker": "CD3e", "level": 1.0}])
    ingest = tmp_path / "cells_pln"
    for donor in ("0112", "6539"):
        (ingest / f"donor_id={donor}").mkdir(parents=True)
    crosswalk = write_crosswalk(tmp_path / "vocab.csv", RECORDED)
    table = build(sources, xenium_roots=[typing_root], pln_pheno_root=pheno,
                  pln_ingest_root=ingest, vocab_crosswalk=crosswalk)
    rows = [r for r in table.to_pylist() if r["roi"] != "panc"]

    by_modality = {}
    for row in rows:
        if row["resolved_by"] == "all" and row["state"] == "measured":
            by_modality[row["modality"]] = (row["n_cells"], row["n_immune"])
    # Two sections for one donor, counted separately and never added together.
    assert by_modality == {"xenium": (6, 4), "phenocycler": (4, 3)}
    assert {r["modality"] for r in rows} == {"xenium", "phenocycler"}

    # The donor whose PhenoCycler calls have not landed is named, not missing.
    waiting = [r for r in rows if r["modality"] == "phenocycler" and r["donor_id"] == "6539"]
    assert len(waiting) == 1
    assert waiting[0]["state"] == "unavailable"
    assert waiting[0]["reason"] == "phenocycler_pln_calls_absent_for_donor"
    assert waiting[0]["n_cells"] is None


def test_an_absent_phenocycler_lane_names_every_donor_that_has_a_section(
        sources, tmp_path):
    """The refusal comes from the ingest, never from a written-down donor list."""
    ingest = tmp_path / "cells_pln"
    for donor in ("0112", "6539", "6380"):
        (ingest / f"donor_id={donor}").mkdir(parents=True)
    table = build(sources, pln_pheno_root=tmp_path / "never_written",
                  pln_ingest_root=ingest,
                  vocab_crosswalk=write_crosswalk(tmp_path / "vocab.csv", RECORDED))
    rows = [r for r in table.to_pylist() if r["modality"] == "phenocycler"
            and r["roi"] != "panc"]
    assert {r["donor_id"] for r in rows} == {"112", "6539", "6380"}
    assert {r["reason"] for r in rows} == {"phenocycler_pln_calls_absent"}
    assert all(r["state"] == "unavailable" and r["n_cells"] is None for r in rows)


def test_hive_partitions_keep_donor_ids_as_strings(tmp_path):
    """Every donor id here is all digits, so an inferring reader makes them ints.

    The rule lives in the shared resolver now; the builder's own copy of it was
    defined and never called, so it is gone rather than left to drift out of step
    with the one the releases are actually read through.
    """
    assert not hasattr(builder, "PLN_PARTITIONING")
    root = write_release(tmp_path / "xenium_typing" / "release_pln", {"0112": SECTION})
    release = builder.resolve_release("pln", [root.parent])
    schema = release.partitioning().schema
    for key in ("donor_id", "roi"):
        assert pa.types.is_string(schema.field(key).type) or \
               pa.types.is_large_string(schema.field(key).type)

    # And the donors come back as the text they were written as, not as integers.
    assert [donor for donor, _roi, _path in release.partition_files()] == ["0112"]


def test_a_release_is_chosen_by_its_declared_tissue_not_its_directory_name(tmp_path):
    """The selection moved to scripts/senior/xenium_release.py; the rule did not."""
    root = tmp_path / "xenium_typing"
    write_release(root / "release", {"0112": SECTION}, tissue="pancreas", roi="panc")
    missing = builder.resolve_release("pln", [root])
    assert missing.available is False
    assert "declares this tissue" in missing.reason

    pln_release = write_release(root / "anything_at_all", {"0112": SECTION}, tissue="pln")
    found = builder.resolve_release("pln", [root])
    assert found.available is True and found.root == pln_release

    absent = builder.resolve_release("pln", [tmp_path / "does_not_exist"])
    assert absent.available is False
    # Never None: a None that reaches a count reads as a zero.
    assert absent is not None and absent.reason


def test_the_release_is_chosen_by_its_declared_time_not_a_directory_mtime(tmp_path):
    """The defect this replaced: a `touch` could reorder which release was read."""
    root = tmp_path / "xenium_typing"
    older = write_release(root / "v1", {"0112": SECTION}, release_id="fixture-v1",
                          generated_at="2026-09-16T00:00:00+00:00")
    newer = write_release(root / "v2", {"0112": SECTION}, release_id="fixture-v2",
                          generated_at="2026-09-17T00:00:00+00:00")
    assert builder.resolve_release("pln", [root]).root == newer

    # Make the SUPERSEDED one the newest thing on the filesystem. The answer
    # must not move, because a file's mtime is not a release's own word.
    future = time.time() + 10_000
    os.utime(older / "release.json", (future, future))
    os.utime(older, (future, future))
    assert builder.resolve_release("pln", [root]).root == newer

    # And a release that declares it supersedes another retires it whatever the
    # timestamps say. v3 is the OLDEST declaration here, so it does not win by
    # being newest -- what it does is take v2 out, leaving v1 current. Superseded
    # means not current, not "current unless something newer exists".
    write_release(root / "v3", {"0112": SECTION}, release_id="fixture-v3",
                  generated_at="2026-09-15T00:00:00+00:00", supersedes="fixture-v2")
    current = builder.resolve_release("pln", [root])
    assert current.release_id != "fixture-v2"
    assert current.root == older


def test_the_pln_release_is_searched_under_every_root_it_is_given(tmp_path):
    """The pancreas release lives OUTSIDE the repository; one root can never see both."""
    inside, outside = tmp_path / "repo_tree", tmp_path / "external_tree"
    write_release(outside / "release_pln", {"0112": SECTION})
    assert builder.resolve_release("pln", [inside]).available is False
    assert builder.resolve_release("pln", [inside, outside]).available is True


def test_the_card_shows_both_modalities_without_conflating_them(sources, tmp_path):
    """Two sections of one node, typed by two lanes, must stay two measurements."""
    rscript = r_or_skip()
    typing_root = tmp_path / "xenium_typing"
    write_release(typing_root / "release_pln_immune", {"0112": SECTION})
    pheno = tmp_path / "pln" / "pheno"
    write_rows(pheno / "donor_id=0112" / "roi=pln_1" / "data.parquet",
               [{"cell_id": f"p{i}", "stage1": stage1, "common_lineage": lineage,
                 "nucleus_count": nucleus}
                for i, (stage1, lineage, nucleus) in enumerate(
                    [("Immune", "Myeloid", 1)] * 3 + [("Non_immune", "", 1)])])
    ingest = tmp_path / "cells_pln"
    (ingest / "donor_id=0112").mkdir(parents=True)
    table = build(sources, xenium_roots=[typing_root], pln_pheno_root=pheno,
                  pln_ingest_root=ingest,
                  vocab_crosswalk=write_crosswalk(tmp_path / "vocab.csv", RECORDED))
    root = layers_dir(tmp_path, "both", patch_metadata(
        table, source_xen_cells="/repo/core/phenocycler/data/integration/cells_xen",
        source_immune_vocabulary="/repo/core/phenocycler/data/integration/vocab_crosswalk.csv"))
    env = dict(os.environ, ISLET_LAYERS_DIR=str(root), DISTANCE_TEST_DOMAIN=DOMAIN_ROOT)
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
source("apps/R/ieCore/R/ie_explain.R")
source("apps/R/ieCore/R/ie_frames.R")
source("apps/islet/R/card_x1_gradient.R")
source("apps/islet/R/card_s3_compartments.R")
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
text_of <- function(html) trimws(gsub("[[:space:]]+", " ", gsub("<[^>]+>", " ", html)))
shiny::testServer(function(input, output, session) {
  result <- card_s3_compartments_server(input, output, session,
                                        active_tab = shiny::reactive("LymphNode"),
                                        cohort = list(donor_ids = function() c("00112")))
}, {
  session$setInputs(s3_compartments_donor = "112", s3_compartments_pair = "Immune | all")
  body <- text_of(output$s3_compartments_pln$html)
  # Both modalities are named, and their own denominators are shown.
  stopifnot(grepl("xenium", body, fixed = TRUE), grepl("phenocycler", body, fixed = TRUE))
  stopifnot(grepl("Modality", body, fixed = TRUE))
  # The Xenium section resolved 6 cells, the PhenoCycler one 4. Neither is 10.
  stopifnot(grepl("xenium Immune all 4 6", body, fixed = TRUE))
  stopifnot(grepl("phenocycler Immune all 3 4", body, fixed = TRUE))
  stopifnot(!grepl(" 10 ", body, fixed = TRUE))
  # The don't-conflate guarantee must be REACHABLE, not printed above every panel. It moved into
  # the "how this is measured" modal when the tab was measured at ~470 words of standing prose;
  # what the reader must still get is the one-line claim and a way to the full statement.
  stopifnot(grepl("one modality at a time", body, fixed = TRUE))
  stopifnot(grepl("how this is measured", body, fixed = TRUE))
  # And the statement itself still exists, in the modal that link opens.
  stopifnot(grepl("Compare within one, never across",
                  paste(readLines("apps/islet/R/card_s3_compartments.R"), collapse = " "),
                  fixed = TRUE))
  # Also moved into the modal, and checked there for the same reason.
  stopifnot(grepl("do not sum to the Immune total",
                  paste(readLines("apps/islet/R/card_s3_compartments.R"), collapse = " "),
                  fixed = TRUE))
  # And the node column still grows no distance axis.
  for (word in c("signed_distance", "distance_bin", "Islet core")) {
    stopifnot(!grepl(word, body, fixed = TRUE))
  }
})
cat("R two-modality pLN column checks passed\n")
'''
    result = subprocess.run([str(rscript), "--vanilla", "-"], input=code, cwd=REPO,
                            env=env, capture_output=True, text=True, check=False, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr


# --------------------------------------------------------------------------- #
# The pancreas Xenium arm: the same tissue, the other modality
#
# Four links have to hold before a Xenium pancreas cell has a signed distance to
# an islet border (see the builder's module docstring). These fixtures supply
# them one at a time so each refusal is the FIRST missing link and nothing is
# ever a zero standing in for an absence.
# --------------------------------------------------------------------------- #

#: (native lineage, common lineage, nucleus_count, nucleus_measured)
PANC_CELL = ("Immune_T", "T_NK", 1, True)


def write_panc_release(root, donors, *, release_id="fixture-panc",
                       generated_at="2026-09-17T20:36:54+00:00"):
    """A Xenium pancreas release fixture: {donor: [(native, common, count, measured)]}.

    It carries NO coordinate column, which is the fact the whole arm turns on:
    the release says what a cell is, never where it is. The measured release does
    the same -- schema read 2026-09-18.
    """
    root.mkdir(parents=True, exist_ok=True)
    (root / "release.json").write_text(json.dumps(
        {"schema_version": 1, "modality": "xenium", "roi": "panc", "calls_dir": "calls",
         "cell_key": "cell_id", "release_id": release_id, "written_utc": generated_at,
         "columns": {"label": "lineage", "common_lineage": "common_lineage"}}))
    for donor, cells in donors.items():
        write_rows(root / "calls" / f"donor_id={donor}" / "roi=panc" / "data.parquet",
                   [{"cell_id": f"{donor}-{i}", "lineage": native,
                     "common_lineage": common, "nucleus_count": count,
                     "nucleus_measured": measured}
                    for i, (native, common, count, measured) in enumerate(cells)])
    return root


def write_panc_coordinates(root, donor, positions, *, sample_key="0041032__panc",
                           first_index=0):
    """cells_xen for one pancreas section: the Xenium-frame coordinates, by cell_id."""
    write_rows(root / f"donor_id={donor}" / "roi=panc" / "data.parquet",
               [{"cell_id": f"{donor}-{i + first_index}", "x_um": float(x), "y_um": float(y),
                 "xenium_sample_key": sample_key, "tissue": "pancreas"}
                for i, (x, y) in enumerate(positions)])
    return root


def write_section_modes(path, rows=()):
    """The acquisition-design manifest. A pair absent from it is serial."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["donor_id,roi,xenium_serial,mode,evidence,note\n"]
    for donor, roi, serial, mode, evidence in rows:
        lines.append(f"{donor},{roi},{serial},{mode},{evidence},\n")
    path.write_text("".join(lines))
    return path


def write_square_border(root, donor, squares):
    """Islet outlines in the XENIUM frame, as the structures lane would write them."""
    shapely_wkb = pytest.importorskip("shapely.wkb")
    from shapely.geometry import box
    write_rows(root / f"donor_id={donor}" / "roi=panc" / "islets.parquet",
               [{"structure_id": f"islet_{i}", "kind": "islet", "x_um": (x0 + x1) / 2,
                 "y_um": (y0 + y1) / 2, "area_um2": (x1 - x0) * (y1 - y0),
                 "wkb": shapely_wkb.dumps(box(x0, y0, x1, y1))}
                for i, (x0, y0, x1, y1) in enumerate(squares)])
    return root


@pytest.fixture
def panc_xenium(tmp_path):
    """Every input the pancreas Xenium arm reads, each one switchable on its own."""
    return dict(
        releases=tmp_path / "xenium_releases",
        cells=tmp_path / "cells_xen",
        structures=tmp_path / "structures_xenium",
        registration=tmp_path / "registration",
        modes=write_section_modes(tmp_path / "section_modes.csv"),
        crosswalk=write_crosswalk(tmp_path / "vocab.csv", RECORDED),
    )


def build_panc(sources, panc_xenium, **kwargs):
    return build(sources, xenium_roots=[panc_xenium["releases"]],
                 xen_cells_root=panc_xenium["cells"],
                 xenium_structures_root=panc_xenium["structures"],
                 registration_root=panc_xenium["registration"],
                 section_modes_csv=panc_xenium["modes"],
                 vocab_crosswalk=panc_xenium["crosswalk"], **kwargs)


def panc_xenium_rows(table):
    return [r for r in table.to_pylist()
            if r["roi"] == "panc" and r["modality"] == "xenium"]


def test_no_xenium_arm_is_built_unless_a_root_is_named(sources):
    """A build that names no root must not reach a release on this machine."""
    table = build(sources)
    assert panc_xenium_rows(table) == []
    metadata = table.schema.metadata
    assert metadata[b"source_xenium_panc"].decode() == ""


def test_the_release_alone_can_never_reach_the_distance_axis(sources, panc_xenium):
    """Link 2: the release carries no coordinate column, so without the ingest
    there is no position for any cell -- and that is named, not counted as zero."""
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 4})
    assert "x_um" not in set(pq.read_schema(
        panc_xenium["releases"] / "panc_release" / "calls" / "donor_id=0112" /
        "roi=panc" / "data.parquet").names)

    table = build(sources, xenium_roots=[panc_xenium["releases"]],
                  vocab_crosswalk=panc_xenium["crosswalk"])
    assert panc_xenium_rows(table) == []
    assert (table.schema.metadata[b"source_xenium_panc"].decode()
            == builder.XEN_PANC_NO_COORDINATES)


def test_the_arm_names_the_first_missing_link_for_each_donor(sources, panc_xenium):
    write_panc_release(panc_xenium["releases"] / "panc_release", {
        "0112": [PANC_CELL] * 4,                          # has coordinates below
        "6347": [("Immune_T", "T_NK", None, False)] * 3,  # nucleus never measured
        "6380": [PANC_CELL] * 2,                          # no ingest partition
    })
    write_panc_coordinates(panc_xenium["cells"], "0112", [(1, 1), (2, 2), (3, 3), (4, 4)])
    # The Xenium cohort is not the PhenoCycler cohort: these two donors have an
    # RNA section and no pancreas partition here, and they still get a row.
    with sources["metadata"].open("a") as handle:
        handle.write("6347,ND\n6380,T1D\n")
    table = build_panc(sources, panc_xenium)
    by_donor = {r["donor_id"]: r for r in panc_xenium_rows(table)}
    assert set(by_donor) == {"112", "6347", "6380"}

    # A nucleus that was never measured is a MISSING MEASUREMENT: the rule "a
    # cell with no nucleus is excluded" never applied to these cells, and the
    # count is named rather than folded into an absent-nucleus tally.
    assert by_donor["6347"]["state"] == "blank"
    assert by_donor["6347"]["reason"] == "xenium_panc_nucleus_not_measured:3_of_3"
    # Link 2 for one donor.
    assert by_donor["6380"]["reason"] == builder.XEN_PANC_NO_COORDINATES
    # Link 3/4: coordinates exist, and the only islet border is on another section.
    assert by_donor["112"]["reason"] == (
        "xenium_panc_islet_border_is_on_another_physical_section:sequential")
    for row in by_donor.values():
        assert row["state"] != "measured"
        assert row["reason"]
        for column in ("n_cells", "n_immune", "immune_fraction", "distance_bin_um",
                       "compartment"):
            assert row[column] is None


def test_a_section_the_release_does_not_call_is_named_not_left_off(sources, panc_xenium):
    """A donor with no row reads as a donor with no section, and it has one.

    Measured 2026-09-18: 6450, 6505, 6534 and 6623 have a pancreas Xenium
    partition in the ingest and no call in the 21-donor release.
    """
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 2})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(1, 1), (2, 2)])
    write_panc_coordinates(panc_xenium["cells"], "6450", [(1, 1)])
    with sources["metadata"].open("a") as handle:
        handle.write("6450,T1D\n")
    by_donor = {r["donor_id"]: r for r in panc_xenium_rows(build_panc(sources, panc_xenium))}
    assert set(by_donor) == {"112", "6450"}
    assert by_donor["6450"]["reason"] == builder.XEN_PANC_NOT_IN_RELEASE
    assert by_donor["6450"]["state"] == "unavailable"
    assert by_donor["6450"]["donor_status"] == "T1D"
    assert by_donor["6450"]["n_cells"] is None


def test_a_serial_pair_is_refused_by_the_declared_design_not_by_a_guess(
        sources, panc_xenium):
    """The acquisition design is read from the manifest that owns it.

    Declaring the pair same_slide does NOT produce a number. It moves the refusal
    to the next missing link, which is the point: the declaration is what makes a
    per-cell claim admissible, and the geometry still has to exist.
    """
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 4})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(1, 1), (2, 2), (3, 3), (4, 4)],
                           sample_key="0080870__panc")
    serial = build_panc(sources, panc_xenium)
    assert panc_xenium_rows(serial)[0]["reason"].endswith(":sequential")

    write_section_modes(panc_xenium["modes"], [
        ("112", "panc", "0080870", "same_slide",
         "one physical section imaged twice; fixture")])
    declared = build_panc(sources, panc_xenium)
    row = panc_xenium_rows(declared)[0]
    assert row["reason"] == builder.XEN_PANC_NO_REGISTRATION
    assert row["state"] == "unavailable" and row["n_cells"] is None

    # With the transform on disk the transport is still not this builder's to
    # invent: it belongs to the same-slide lane and arrives as a loader.
    transform = (panc_xenium["registration"] / "donor_id=112" / "roi=panc" / "transform.json")
    transform.parent.mkdir(parents=True)
    transform.write_text(json.dumps({"affine": [[1, 0, 0], [0, 1, 0]]}))
    assert panc_xenium_rows(build_panc(sources, panc_xenium))[0]["reason"] == (
        builder.XEN_PANC_NO_BORDER_TRANSPORT)


def test_a_structures_table_without_an_outline_is_named_not_read_as_absent(
        sources, panc_xenium):
    """A centroid and an area are not a border, and nothing turns one into a radius."""
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 4})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(1, 1), (2, 2), (3, 3), (4, 4)])
    write_rows(panc_xenium["structures"] / "donor_id=112" / "roi=panc" / "islets.parquet",
               [{"structure_id": "islet_0", "x_um": 5.0, "y_um": 5.0, "area_um2": 100.0,
                 "eccentricity": 0.1}])
    row = panc_xenium_rows(build_panc(sources, panc_xenium))[0]
    assert row["state"] == "unavailable"
    assert row["reason"].startswith("xenium_panc_islet_structures_carry_no_outline:")
    assert "area_um2" in row["reason"] and "x_um" in row["reason"]


def test_a_border_in_the_xenium_frame_is_measured_on_the_protein_arm_axis(
        sources, panc_xenium):
    """The one route that measures today: an outline segmented on the Xenium section.

    Distances are checked against the geometry, not against any previous result.
    """
    pytest.importorskip("shapely")
    # An islet occupying (0,0)-(10,10). Cell positions are chosen so each signed
    # distance is exact: -5 at the centre, +2 just outside, +30 well away.
    cells = [("Immune_T", "T_NK", 1, True), ("Exocrine", "Exocrine", 1, True),
             ("Immune_T", "T_NK", 1, True), ("Immune_B", "B_Plasma", 1, True),
             ("Immune_T", "T_NK", 0, True)]
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": cells})
    write_panc_coordinates(panc_xenium["cells"], "0112",
                           [(5, 5), (5, 5), (5, 12), (5, 40), (5, 5)])
    write_square_border(panc_xenium["structures"], "112", [(0, 0, 10, 10)])

    rows = panc_xenium_rows(build_panc(sources, panc_xenium))
    assert rows and all(r["state"] == "measured" for r in rows)
    assert {r["frame"] for r in rows} == {"signed_distance_um"}

    totals = {(r["distance_bin_um"], r["compartment"]): (r["n_cells"], r["n_immune"])
              for r in rows if r["resolved_by"] == "all"}
    # The cell with no nucleus is excluded from every calculation, so the core
    # bin holds two cells, not three.
    assert totals == {(-10.0, "islet_core"): (2, 1),
                      (10.0, "peri_islet"): (1, 1),
                      (30.0, "other_pancreas"): (1, 1)}

    # The lineage rows read the crosswalk's common vocabulary, not the release's
    # native label: nothing in the crosswalk records "Immune_T".
    subtypes = {(r["distance_bin_um"], r["immune_subtype"]): r["n_immune"]
                for r in rows if r["resolved_by"] == "common_lineage"}
    assert subtypes[(-10.0, "T_NK")] == 1
    assert subtypes[(30.0, "B_Plasma")] == 1
    assert subtypes[(-10.0, "B_Plasma")] == 0     # a measured zero, in a real bin
    for row in rows:
        assert abs(row["immune_fraction"] - row["n_immune"] / row["n_cells"]) < 1e-12


def test_a_cell_on_the_border_is_cut_the_same_way_in_both_modalities(
        sources, panc_xenium):
    """Zero keeps its own bin, and it does not take a sign from either side.

    Containment decides the core in both arms, and a centroid exactly on the
    outline is not within it -- the same answer `typing/islet_membership.py`
    gives on the protein side. The two arms are not made to agree by a rule
    written here; they agree because both read containment.
    """
    pytest.importorskip("shapely")
    write_panc_release(panc_xenium["releases"] / "panc_release",
                       {"0112": [PANC_CELL, ("Exocrine", "Exocrine", 1, True)]})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(0, 5), (0, 5)])
    write_square_border(panc_xenium["structures"], "112", [(0, 0, 10, 10)])
    rows = [r for r in panc_xenium_rows(build_panc(sources, panc_xenium))
            if r["resolved_by"] == "all"]
    assert len(rows) == 1
    assert rows[0]["distance_bin_um"] == 0.0
    assert str(rows[0]["distance_bin_um"]) != "-0.0"
    assert rows[0]["compartment"] == "other_pancreas"
    assert (rows[0]["n_cells"], rows[0]["n_immune"]) == (2, 1)


def test_the_two_pancreas_modalities_are_never_pooled(sources, panc_xenium):
    """One tissue, two measurements. Their denominators are not interchangeable."""
    pytest.importorskip("shapely")
    write_panc_release(panc_xenium["releases"] / "panc_release",
                       {"0112": [PANC_CELL] * 3})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(5, 5), (5, 5), (5, 40)])
    write_square_border(panc_xenium["structures"], "112", [(0, 0, 10, 10)])
    table = build_panc(sources, panc_xenium)

    panc = [r for r in table.to_pylist() if r["roi"] == "panc"]
    assert {r["modality"] for r in panc} == {"phenocycler", "xenium"}
    # 3 Xenium cells; 7 QC-eligible, non-excluded PhenoCycler cells of the 10 in
    # the fixture. Neither denominator is the other's, and neither is 10.
    for modality, expected in (("xenium", 3), ("phenocycler", 7)):
        totals = [r for r in panc if r["modality"] == modality
                  and r["resolved_by"] == "all" and r["state"] == "measured"]
        assert sum(r["n_cells"] for r in totals) == expected
    # Nothing in the artifact carries both modalities' cells in one row.
    assert not any(r["n_cells"] == 10 for r in panc if r["n_cells"] is not None)


def test_the_cell_id_join_is_the_only_join_and_a_miss_is_named(sources, panc_xenium):
    """Row order, file order and a centroid match are excluded by construction."""
    pytest.importorskip("shapely")
    write_panc_release(panc_xenium["releases"] / "panc_release",
                       {"0112": [PANC_CELL] * 4})
    # Three of the four release cells have coordinates, and they are written in
    # the opposite order: only the key can match them.
    write_rows(panc_xenium["cells"] / "donor_id=0112" / "roi=panc" / "data.parquet",
               [{"cell_id": f"0112-{i}", "x_um": 5.0, "y_um": 5.0,
                 "xenium_sample_key": "0041032__panc"} for i in (2, 1, 0)])
    write_square_border(panc_xenium["structures"], "112", [(0, 0, 10, 10)])
    row = panc_xenium_rows(build_panc(sources, panc_xenium))[0]
    assert row["state"] == "unavailable"
    assert row["reason"] == "xenium_panc_cell_coordinates_missing_cell_id:1_of_4"
    assert row["n_cells"] is None

    # Supply the fourth and the same three still match by key, not by position.
    write_rows(panc_xenium["cells"] / "donor_id=0112" / "roi=panc" / "data.parquet",
               [{"cell_id": f"0112-{i}", "x_um": 5.0, "y_um": 5.0,
                 "xenium_sample_key": "0041032__panc"} for i in (3, 2, 1, 0)])
    measured = panc_xenium_rows(build_panc(sources, panc_xenium))
    assert all(r["state"] == "measured" for r in measured)
    assert sum(r["n_cells"] for r in measured if r["resolved_by"] == "all") == 4


def test_an_unrecorded_lineage_is_withheld_by_name_in_the_pancreas_arm(
        sources, panc_xenium):
    """Never defaulted to non-immune, exactly as in the pLN arm."""
    pytest.importorskip("shapely")
    write_panc_release(panc_xenium["releases"] / "panc_release",
                       {"0112": [("Immune_T", "T_NK", 1, True),
                                 ("Something_New", "Tuft", 1, True)]})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(5, 5), (5, 5)])
    write_square_border(panc_xenium["structures"], "112", [(0, 0, 10, 10)])
    row = panc_xenium_rows(build_panc(sources, panc_xenium))[0]
    assert row["state"] == "withheld"
    assert row["reason"] == "xenium_immune_lineage_not_recorded:Tuft"
    assert row["n_cells"] is None and row["n_immune"] is None


def test_the_pancreas_xenium_rows_reach_the_final_cross_tissue_gate_unchanged(
        sources, panc_xenium):
    """(c): the gate judges the TISSUE, and a pancreas row is a pancreas row.

    The arm's rows are roi=panc because they are the pancreas. Had they been
    given a different roi to dodge a reader, `refuse_cross_tissue_distance` would
    have started meaning "a different modality" instead of "a different tissue".
    """
    pytest.importorskip("shapely")
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 2})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(5, 5), (5, 40)])
    write_square_border(panc_xenium["structures"], "112", [(0, 0, 10, 10)])
    table = build_panc(sources, panc_xenium)
    builder._assert_no_cross_tissue_distance(table)      # does not raise
    assert {r["tissue"] for r in panc_xenium_rows(table)} == {"pancreas"}
    with pytest.raises(ValueError, match="different blocks"):
        builder.refuse_cross_tissue_distance("pln_1", "pancreatic_lymph_node")


def test_the_arm_records_what_it_read_and_what_each_donor_held(sources, panc_xenium):
    """The refusal carries the counts, so a reader is not left with only a verdict."""
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 4})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(1, 1), (2, 2), (3, 3), (4, 4)])
    table = build_panc(sources, panc_xenium)
    audit = json.loads(table.schema.metadata[b"panc_xenium_audit"].decode())
    assert audit["release"]["release_id"] == "fixture-panc"
    assert audit["release"]["selected_by"].endswith("never a directory mtime")
    entry = audit["112/panc"]
    assert (entry["cells"], entry["called"]) == (4, 4)
    assert entry["section_mode"] == "sequential"
    assert entry["xenium_serial"] == "0041032"


# --------------------------------------------------------------------------- #
# The reader contract the second pancreas modality creates
# --------------------------------------------------------------------------- #

@pytest.mark.xfail(
    strict=False,
    reason="card_x1_gradient.R:.ie_distance_select selects the pancreas arm by donor "
           "alone, so it reads the Xenium row's absence as the PhenoCycler gradient's. "
           "Measured 2026-09-18: donor 112 with 7 measured PhenoCycler rows and one "
           "unavailable Xenium row returns unavailable(xenium_panc_islet_border_is_on_"
           "another_physical_section:sequential). The card must filter on modality "
           "before donor. That file is outside this lane; this test is the request.")
def test_a_donors_protein_gradient_survives_the_other_modalitys_absence(
        sources, panc_xenium, tmp_path):
    """One tissue, two measurements: one modality's hole is not the other's.

    The pancreas arm now carries both modalities, and a reader that narrows it to
    a donor without first narrowing it to a modality will report the RNA
    section's refusal as though the protein section had none.
    """
    rscript = r_or_skip()
    write_panc_release(panc_xenium["releases"] / "panc_release", {"0112": [PANC_CELL] * 4})
    write_panc_coordinates(panc_xenium["cells"], "0112", [(1, 1), (2, 2), (3, 3), (4, 4)])
    table = build_panc(sources, panc_xenium)
    panc = [r for r in table.to_pylist() if r["roi"] == "panc"]
    assert any(r["modality"] == "phenocycler" and r["state"] == "measured" for r in panc)
    assert any(r["modality"] == "xenium" and r["state"] != "measured" for r in panc)

    root = layers_dir(tmp_path, "two_modalities", patch_metadata(table))
    env = dict(os.environ, ISLET_LAYERS_DIR=str(root), DISTANCE_TEST_DOMAIN=DOMAIN_ROOT)
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
source("apps/R/ieCore/R/ie_explain.R")
source("apps/R/ieCore/R/ie_frames.R")
source("apps/islet/R/card_x1_gradient.R")
domain_data_root <- function() Sys.getenv("DISTANCE_TEST_DOMAIN")
layer <- .ie_distance_read()
stopifnot(ie_is_measured(layer$value))
selected <- .ie_distance_select(layer, "112", "Immune | all")
if (!ie_is_measured(selected$value)) {
  stop(sprintf("the protein gradient was withheld over the RNA section: %s(%s)",
               selected$value$state, selected$value$reason))
}
stopifnot(all(selected$data$modality == "phenocycler"))
cat("R two-pancreas-modality selection check passed\n")
'''
    result = subprocess.run([str(rscript), "--vanilla", "-"], input=code, cwd=REPO,
                            env=env, capture_output=True, text=True, check=False, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
