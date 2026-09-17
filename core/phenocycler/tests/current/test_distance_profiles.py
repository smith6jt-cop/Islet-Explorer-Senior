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
