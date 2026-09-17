"""UUID, denominator, availability, and app checks for the distance layer.

Run from the repository root with .venv/bin/python -m pytest <this file>.
Fixtures are synthetic measurements, never expected biological proportions.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
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
  stopifnot(identical(plot$scales$get_scales("x")$trans$name, "pseudo_log"))
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
                  calls_dir="calls"):
    """A Xenium pLN release fixture. {donor: [(stage1, lineage, nucleus_count)]}."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "release.json").write_text(json.dumps(
        {"schema_version": 1, "modality": "xenium", "tissue": tissue,
         "calls_dir": calls_dir, "cell_key": "cell_id", "roi": ["pln_1"]}))
    for donor, cells in donors.items():
        write_rows(root / calls_dir / f"donor_id={donor}" / "roi=pln_1" / "data.parquet",
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
    table = build(sources, xen_typing_root=typing_root, vocab_crosswalk=crosswalk)
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
    table = build(sources, xen_typing_root=typing_root, xen_cells_root=pln["root"],
                  vocab_crosswalk=crosswalk)
    rows = [r for r in table.to_pylist() if r["roi"] != "panc"]
    total = [r for r in rows if r["resolved_by"] == "all"][0]
    assert (total["n_cells"], total["n_immune"]) == (6, 4)   # the release, not the ingest
    sources_by_modality = json.loads(
        table.schema.metadata[b"source_pln_by_modality"].decode())
    assert sources_by_modality["xenium"] == str(release.resolve())

    # With no release present, the ingest is read and named as the source.
    fallback = build(sources, xen_typing_root=tmp_path / "empty",
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
    table = build(sources, xen_typing_root=typing_root, pln_pheno_root=pheno,
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


def test_hive_partitions_keep_donor_ids_as_strings():
    """Every pLN donor id is all digits, so inferred partitioning makes them ints."""
    schema = builder.PLN_PARTITIONING.schema
    assert schema.field("donor_id").type == pa.string()
    assert schema.field("roi").type == pa.string()


def test_a_release_is_chosen_by_its_declared_tissue_not_its_directory_name(tmp_path):
    root = tmp_path / "xenium_typing"
    write_release(root / "release", {"0112": SECTION}, tissue="pancreas")
    assert builder.find_pln_release(root) is None
    pln_release = write_release(root / "anything_at_all", {"0112": SECTION}, tissue="pln")
    assert builder.find_pln_release(root) == pln_release
    assert builder.find_pln_release(tmp_path / "does_not_exist") is None


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
    table = build(sources, xen_typing_root=typing_root, pln_pheno_root=pheno,
                  pln_ingest_root=ingest,
                  vocab_crosswalk=write_crosswalk(tmp_path / "vocab.csv", RECORDED))
    root = layers_dir(tmp_path, "both", patch_metadata(
        table, source_xen_cells="/repo/core/phenocycler/data/integration/cells_xen",
        source_immune_vocabulary="/repo/core/phenocycler/data/integration/vocab_crosswalk.csv"))
    env = dict(os.environ, ISLET_LAYERS_DIR=str(root), DISTANCE_TEST_DOMAIN=DOMAIN_ROOT)
    code = r'''
source("apps/R/ieCore/R/ie_availability.R")
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
  stopifnot(grepl("compare within one, never across", body, fixed = TRUE))
  stopifnot(grepl("do not sum to the Immune total", body, fixed = TRUE))
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
