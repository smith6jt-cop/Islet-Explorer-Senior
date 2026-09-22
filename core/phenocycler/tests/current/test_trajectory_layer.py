"""Within-donor trajectory acceptance, absence semantics, and card smoke checks.

Run from the repository root with:
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
        core/phenocycler/tests/current/test_trajectory_layer.py

Parquet round trips use memory buffers; the source trees are never modified.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
SPEC = importlib.util.spec_from_file_location(
    "build_trajectory_layer", REPO / "scripts/senior/build_trajectory_layer.py"
)
layer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(layer)


@pytest.fixture
def inputs():
    summary = pd.DataFrame([
        {"case_id": donor, "islet_key": key, "sc_pseudotime_core": core,
         "sc_pseudotime_mean_core": core, "sc_n_core": 5,
         "sc_pseudotime_peri": peri, "sc_pseudotime_mean_peri": peri, "sc_n_peri": 6,
         layer.DELTA: peri - core}
        for donor, key, core, peri in [
            ("0112", "Islet_1", 0.9, 0.8), ("112", "Islet_2", 0.8, 0.9),
            ("6533", "Islet_1", 0.1, 0.2), ("6533", "Islet_2", 0.2, 0.2),
            ("6533", "Islet_3", 0.2, 0.3),
        ]
    ])
    cells = pd.DataFrame([
        {"donor_id": donor, "object_id": f"00000000-0000-4000-8000-{i:012d}",
         "donor_status": "ND", "stain_group": batch, "pseudotime": value}
        for i, (donor, batch, value) in enumerate([
            ("0112", "S0", 0.0), ("112", "S0", 2.0),
            ("6533", "S1", 2.0), ("6533", "S1", 4.0),
        ])
    ])
    return summary, cells


def build(inputs):
    summary, cells = inputs
    return layer.trajectory_table(
        summary, cells, {"n_ordered": len(cells), "n_population": 100},
        domain_root=str((REPO / "data/development/typing/app").resolve()),
    )


def test_acceptance_builder_groups_by_donor_and_emits_batch_eta_sq(inputs):
    table = build(inputs)
    assert "batch_eta_sq_within_status" in table.column_names
    assert table.schema.field("donor_id").type == pa.string()
    rows = table.to_pandas().set_index(["donor_id", "islet_key"])
    # Both donors restart at 1 despite disjoint core-pseudotime ranges. Identical
    # islet keys in different donors cannot merge, and padding cannot split a donor.
    assert rows.loc[("112", "Islet_2"), "within_donor_rank"] == 1
    assert rows.loc[("112", "Islet_1"), "within_donor_rank"] == 2
    assert rows.loc[("6533", "Islet_1"), "within_donor_rank"] == 1
    assert rows.loc[("6533", "Islet_2"), "within_donor_rank"] == 2.5
    assert rows.loc[("6533", "Islet_3"), "within_donor_rank"] == 2.5
    assert set(rows.n_islets_ordered) == {2, 3}
    assert rows.batch_eta_sq_within_status.eq(0.5).all()
    assert rows.frame.eq("none").all()


def test_other_donors_and_disease_labels_cannot_change_a_donors_order(inputs):
    baseline = build(inputs).to_pandas()
    summary, cells = (part.copy() for part in inputs)
    other = summary.case_id.eq("6533")
    for column in ("sc_pseudotime_core", "sc_pseudotime_mean_core",
                   "sc_pseudotime_peri", "sc_pseudotime_mean_peri"):
        summary.loc[other, column] += 100
    cells["donor_status"] = "changed audit label"
    altered = build((summary.sample(frac=1, random_state=5), cells)).to_pandas()
    columns = ["donor_id", "islet_key", "within_donor_rank"]
    pd.testing.assert_frame_equal(baseline[columns], altered[columns])


def test_batch_diagnostic_stratifies_status_instead_of_pooling(inputs):
    cells = inputs[1].copy()
    other = cells.copy()
    other["donor_status"] = "T1D"
    other["pseudotime"] += 1000
    both = pd.concat([cells, other], ignore_index=True)
    diagnostics = layer.batch_eta_sq_within_status(both)
    assert set(diagnostics) == {"ND", "T1D"}
    for result in diagnostics.values():
        assert result["batch_eta_sq_within_status"] == pytest.approx(0.5)
        assert result["batch_n_cells"] == 4
        assert result["batch_n_groups"] == 2
        assert result["batch_eta_sq_state"] == "measured"


@pytest.mark.parametrize("condition,reason", [
    ("single_batch", "fewer_than_two_stain_groups_within_status"),
    ("constant", "constant_pseudotime_within_status"),
    ("missing_batch", "sampled_cells_missing_stain_group"),
    ("nonfinite", "no_finite_sampled_pseudotime"),
])
def test_unestimable_batch_effect_is_unavailable_not_zero(inputs, condition, reason):
    cells = inputs[1].copy()
    if condition == "single_batch":
        cells["stain_group"] = "S0"
    elif condition == "constant":
        cells["pseudotime"] = 7.0
    elif condition == "missing_batch":
        cells.loc[0, "stain_group"] = None
    else:
        cells["pseudotime"] = np.nan
    diagnostic = layer.batch_eta_sq_within_status(cells)["ND"]
    assert diagnostic["batch_eta_sq_within_status"] is None
    assert diagnostic["batch_eta_sq_state"] == "unavailable"
    assert diagnostic["batch_eta_sq_reason"] == reason


def test_measured_zero_batch_effect_and_zero_paired_difference_are_preserved(inputs):
    summary, cells = inputs
    cells = cells.copy()
    cells["pseudotime"] = [0.0, 2.0, 0.0, 2.0]
    table = build((summary, cells)).to_pandas()
    assert table.batch_eta_sq_within_status.eq(0).all()
    assert table.batch_eta_sq_state.eq("measured").all()
    zero = table[table.donor_id.eq("6533") & table.islet_key.eq("Islet_2")].iloc[0]
    assert zero[layer.DELTA] == 0
    assert zero.delta_state == "measured"


def test_exclusions_are_imported_and_applied_before_ordering_and_diagnostic(inputs):
    from phenocycler.cohort import DONOR_EXCLUSIONS

    summary, cells = (part.copy() for part in inputs)
    for i, donor in enumerate(DONOR_EXCLUSIONS):
        row = summary.iloc[[0]].copy()
        row["case_id"] = "0" + donor
        summary = pd.concat([summary, row], ignore_index=True)
        cell = cells.iloc[[0]].copy()
        cell["donor_id"] = "0" + donor
        cell["object_id"] = f"excluded-{i}"
        cell["stain_group"], cell["pseudotime"] = "another_batch", 100000.0
        cells = pd.concat([cells, cell], ignore_index=True)
    result = build((summary, cells)).to_pandas()
    assert set(result.donor_id) == {"112", "6533"}
    assert result.batch_eta_sq_within_status.eq(0.5).all()
    assert result.batch_n_cells.eq(4).all()


def test_missing_core_peri_and_delta_remain_explicit_unavailable(inputs):
    summary, cells = (part.copy() for part in inputs)
    summary.loc[0, ["sc_pseudotime_core", "sc_pseudotime_mean_core", layer.DELTA]] = np.nan
    summary.loc[0, "sc_n_core"] = 2
    summary.loc[1, ["sc_pseudotime_peri", "sc_pseudotime_mean_peri", "sc_n_peri", layer.DELTA]] = np.nan
    result = build((summary, cells)).to_pylist()
    by_key = {(row["donor_id"], row["islet_key"]): row for row in result}
    core = by_key[("112", "Islet_1")]
    assert core["state"] == core["delta_state"] == "unavailable"
    assert core["reason"] == "core_pseudotime_not_reported(n_sampled=2)"
    assert core["within_donor_rank"] is core[layer.DELTA] is core["sc_pseudotime_core"] is None
    assert core["sc_n_core"] == 2
    peri = by_key[("112", "Islet_2")]
    assert peri["state"] == "measured"
    assert peri["within_donor_rank"] == 1
    assert peri["peri_reason"] == peri["delta_reason"] == "no_sampled_peri_cells"
    assert peri["sc_n_peri"] is peri[layer.DELTA] is None


def test_no_measured_islets_keeps_donor_rows_and_reasons(inputs):
    summary, cells = (part.copy() for part in inputs)
    columns = [column for column in summary if column.startswith("sc_pseudotime")]
    summary[columns] = np.nan
    result = build((summary, cells)).to_pylist()
    assert len(result) == len(summary)
    assert all(row["state"] == "unavailable" and row["reason"] for row in result)
    assert all(row["within_donor_rank"] is None and row["n_islets_ordered"] is None for row in result)


def test_missing_summary_keeps_sampled_donor_with_named_reason(inputs):
    summary, cells = inputs
    result = build((summary.iloc[:0], cells)).to_pylist()
    assert {row["donor_id"] for row in result} == {"112", "6533"}
    assert all(row["reason"] == "no_islet_summary_rows" for row in result)
    assert all(row["islet_key"] is None for row in result)


def test_absent_sources_emit_one_typed_unavailable_row():
    missing = REPO / "data/development/typing/nonexistent_trajectory_test_source"
    assert not missing.exists()
    table = layer.build_trajectory_layer(missing)
    assert table.num_rows == 1
    row = table.to_pylist()[0]
    assert row["state"] == "unavailable"
    assert row["reason"].startswith("source_not_found:")
    assert row["batch_eta_sq_within_status"] is None
    assert row["n_ordered_cells"] is None


@pytest.mark.parametrize("problem,message", [
    ("duplicate_islet", "duplicate .* after donor normalization"),
    ("duplicate_uuid", "unique, nonmissing object_id"),
    ("floating_donor", "integer strings"),
    ("bad_delta", "peri-minus-core differs"),
    ("fractional_count", "positive integer"),
    ("conflicting_status", "conflicting donor_status"),
])
def test_inconsistent_source_contracts_fail_loudly(inputs, problem, message):
    summary, cells = (part.copy() for part in inputs)
    if problem == "duplicate_islet":
        row = summary.iloc[[0]].copy()
        row["case_id"] = "112"
        summary = pd.concat([summary, row], ignore_index=True)
    elif problem == "duplicate_uuid":
        cells.loc[1, "object_id"] = cells.loc[0, "object_id"]
    elif problem == "floating_donor":
        summary.loc[0, "case_id"] = "112.0"
    elif problem == "bad_delta":
        summary.loc[0, layer.DELTA] = 100.0
    elif problem == "fractional_count":
        summary["sc_n_core"] = summary.sc_n_core.astype(float)
        summary.loc[0, "sc_n_core"] = 0.5
    elif problem == "conflicting_status":
        cells.loc[0, "donor_status"] = "T1D"
    with pytest.raises(ValueError, match=message):
        build((summary, cells))


def test_parquet_output_preserves_nullable_values_and_states(inputs):
    original = build(inputs)
    sink = pa.BufferOutputStream()
    pq.write_table(original, sink)
    restored = pq.ParquetFile(pa.BufferReader(sink.getvalue())).read()
    assert original.equals(restored)
    for column in restored.column_names:
        if column == "state" or column.endswith("_state"):
            assert set(restored[column].to_pylist()) <= layer.STATES


def test_sampling_comes_from_provenance_and_disagreement_is_rejected(inputs):
    summary, cells = inputs
    result = build(inputs).to_pandas()
    assert result.n_ordered_cells.eq(4).all()
    assert result.n_population_cells.eq(100).all()
    missing = layer.trajectory_table(summary, cells, None).to_pylist()
    assert all(row["sampling_state"] == "unavailable" for row in missing)
    assert all(row["n_population_cells"] is None for row in missing)
    with pytest.raises(ValueError, match="sampling metadata disagrees"):
        layer.trajectory_table(summary, cells, {"n_ordered": 5, "n_population": 100})


@pytest.mark.data
def test_real_data_builder_and_in_memory_acceptance_check(monkeypatch, capsys):
    if not (layer.DEFAULT_SOURCE / layer.SUMMARY_FILE).is_file():
        pytest.skip("trajectory source data is not installed")
    observed = []
    dataset = layer.ds.dataset

    def no_hive(*args, **kwargs):
        assert "partitioning" in kwargs and kwargs["partitioning"] is None
        observed.append(args[0])
        return dataset(*args, **kwargs)

    monkeypatch.setattr(layer.ds, "dataset", no_hive)
    assert layer.main(["--check"]) == 0
    assert len(observed) == 2
    report = capsys.readouterr().out
    assert "batch_eta_sq_within_status" in report
    assert "ordered islets" in report and "paired differences" in report


def test_card_single_donor_guards_caption_and_server(inputs):
    rscript = Path("/opt/R/4.5.1/bin/Rscript")
    if not rscript.is_file():
        pytest.skip("the supported R 4.5.1 runtime is not installed")
    code = r'''
    source("apps/R/ieCore/R/ie_availability.R")
    source("apps/R/ieCore/R/ie_frames.R")
    source("apps/islet/R/card_s2_trajectory.R")
    data <- jsonlite::fromJSON(file("stdin"))
    domain_data_root <- function() unique(data$domain_data_root)
    theme_islet <- function() ggplot2::theme_minimal()
    one <- data[data$donor_id == "112", , drop = FALSE]
    one$n_ordered_cells <- 49999L
    one$n_population_cells <- 800562L
    # .s2_trajectory_caption() HAS NOT EXISTED since the card split it in two, and this
    # test had called it ever since -- it is the known red the app lane inherited. The two
    # functions that do exist are checked instead: the PNG caption, which travels with a
    # figure that has no modal, and the one-line on-screen claim.
    caption <- .s2_trajectory_plot_caption(one)
    stopifnot(grepl("49,999 of 800,562", caption, fixed = TRUE),
              grepl("fit within donor", caption, fixed = TRUE),
              # No batch eta-squared, deliberately: .s2_trajectory_modal() states why --
              # every donor is rank-uniform on its own axis, so a pooled figure beside the
              # word "batch" would be arithmetic about the construction, and a zero there
              # would read as "batch was removed".
              !grepl("batch", caption, ignore.case = TRUE))
    claim <- paste(as.character(.s2_trajectory_claim(one, "s2_trajectory_method")), collapse = "")
    stopifnot(grepl("49,999 of 800,562", claim, fixed = TRUE),
              grepl("how this is measured", claim, fixed = TRUE))
    stopifnot(is.null(.s2_trajectory_domain_note(data)))
    stopifnot(inherits(.s2_trajectory_domain_note(data, "/different/domain"), "ie_guard_note"))
    stopifnot(inherits(.s2_trajectory_read("/absent/trajectory/test", "trajectory", character())$note,
                       "ie_guard_note"))
    framed <- .s2_trajectory_plot_data(one)
    stopifnot(identical(ie_frame_of(framed), "none"))
    stopifnot(inherits(try(.s2_trajectory_plot_data(data), silent = TRUE), "try-error"))
    wrong <- one; wrong$frame[[1L]] <- "phenocycler_um:112:panc"
    stopifnot(inherits(try(.s2_trajectory_plot_data(wrong), silent = TRUE), "try-error"))
    ledger <- list(data = data.frame(donor_id = c("112", "6533"), capability = "trajectory",
                                    state = "measured", reason = ""), note = NULL)
    for (state in IE_STATES) {
      trial <- ledger
      trial$data$state <- state
      trial$data$reason <- if (state == "measured") "" else "recorded_absence"
      note <- .s2_trajectory_availability_note(trial, "112")
      stopifnot(if (state == "measured") is.null(note) else inherits(note, "ie_guard_note"))
    }
    .s2_trajectory_read <- function(root, name, required) {
      if (name == "trajectory") list(data = data, note = NULL) else ledger
    }
    shiny::testServer(function(input, output, session) {
      card <- card_s2_trajectory_server(input, output, session)
    }, {
      session$setInputs(s2_trajectory_donor = "112")
      stopifnot(is.null(card$data()$note), all(card$data()$data$donor_id == "112"))
      p <- card$plot()
      stopifnot(inherits(p, "ggplot"), all(p$data$donor_id == "112"))
      session$setInputs(s2_trajectory_donor = "6533")
      stopifnot(all(card$plot()$data$donor_id == "6533"))
    })
    ledger$data$state <- "withheld"
    ledger$data$reason <- "recorded_withholding"
    shiny::testServer(function(input, output, session) {
      card <- card_s2_trajectory_server(input, output, session)
    }, {
      session$setInputs(s2_trajectory_donor = "112")
      stopifnot(inherits(card$data()$note, "ie_guard_note"),
                is.null(card$data()$data), is.null(card$plot()))
    })
    domain_data_root <- function() "/different/domain"
    shiny::testServer(function(input, output, session) {
      card <- card_s2_trajectory_server(input, output, session)
    }, {
      stopifnot(inherits(card$data()$note, "ie_guard_note"),
                is.null(card$data()$data), is.null(card$plot()))
    })
    cat("trajectory card smoke passed\n")
    '''
    result = subprocess.run(
        [str(rscript), "--vanilla", "-e", code], cwd=REPO,
        input=json.dumps(build(inputs).to_pylist()), text=True, capture_output=True,
        env={**os.environ, "R_LIBS_USER": str(Path.home() / "R/x86_64-pc-linux-gnu-library/4.5")},
        timeout=60, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "trajectory card smoke passed" in result.stdout
