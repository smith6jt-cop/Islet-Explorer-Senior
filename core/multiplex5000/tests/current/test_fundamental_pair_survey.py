from __future__ import annotations

import json
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from scripts.review import fundamental_pair_survey as survey
from scripts.review.qupath_fundamental_inputs import (
    METHOD_VERSION as PROJECTION_METHOD_VERSION,
)
from scripts.review.qupath_fundamental_inputs import (
    SOURCE_MEASUREMENTS,
)
from scripts.review.qupath_fundamental_inputs import (
    WORKSPACE_MARKER as PROJECTION_WORKSPACE_MARKER,
)


def _declare_projection(root, donor: str = "D", row_count: int = 8) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / PROJECTION_WORKSPACE_MARKER).write_text(
        PROJECTION_METHOD_VERSION + "\n", encoding="utf-8"
    )
    (root / "projection_summary.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "method_version": PROJECTION_METHOD_VERSION,
                "source_kind": "upstream QuPath full-cell measurement export",
                "source_measurements": SOURCE_MEASUREMENTS,
                "restore_outputs_read": False,
                "donors": {donor: row_count},
                "donor_images": {donor: "image"},
            }
        ),
        encoding="utf-8",
    )


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "donor_id": ["D"] * 8,
            "image": ["image"] * 8,
            "object_id": [f"cell-{index}" for index in range(8)],
            "CD3e": [900.0, 100.0, 120.0, 700.0, 80.0, 600.0, np.nan, 90.0],
            "Iba1_norm": [-0.4, -0.2, 0.8, 0.1, -0.1, 1.2, 0.0, np.nan],
            "Ker8_18": [100.0, 120.0, 110.0, 130.0, 900.0, 800.0, 105.0, 115.0],
        }
    )


def test_roster_is_two_unique_upstream_qupath_pairs_without_sparse_cd79a():
    pairs = survey.ROSTER.pairs
    assert len(pairs) == 2
    assert [pair.target_column for pair in pairs] == ["CD3e", "Iba1_norm"]
    assert "CD79a" not in survey.ROSTER.source_columns
    assert {pair.reference_column for pair in pairs} == {"Ker8_18"}
    assert {pair.target for pair in pairs}.isdisjoint({"CD68", "MPO", "EpCAM"})
    assert survey.ROSTER.method_status == "development_mutable"
    assert not hasattr(survey.ROSTER, "normalization_threshold")
    iba1 = pairs[-1]
    assert iba1.target_source_measurement == "Cell: Iba1: Mean (norm)"
    assert iba1.target_source_kind == "qupath_acquisition_normalized_cell_mean"
    assert iba1.reference_source_measurement == "Cell: Ker8-18: Mean"


def test_roster_accepts_reduction_but_rejects_source_aliasing():
    payload = json.loads(survey.ROSTER_PATH.read_text(encoding="utf-8"))
    reduced = deepcopy(payload)
    reduced["pairs"] = reduced["pairs"][:1]
    assert len(survey.roster_from_dict(reduced).pairs) == 1

    fallback = deepcopy(payload)
    fallback["pairs"][1]["target_source_measurement"] = "Cell: Iba1: Mean"
    with pytest.raises(ValueError, match=r"Cell: Iba1: Mean \(norm\)"):
        survey.roster_from_dict(fallback)

    restore = deepcopy(payload)
    restore["pairs"][1]["target_source_kind"] = "restore_normalized"
    with pytest.raises(ValueError, match="qupath_acquisition_normalized"):
        survey.roster_from_dict(restore)

    duplicate = deepcopy(payload)
    duplicate["pairs"].append(deepcopy(duplicate["pairs"][0]))
    with pytest.raises(ValueError, match="duplicate pairs"):
        survey.roster_from_dict(duplicate)


def test_pair_diagnostic_reports_coverage_without_creating_calls():
    row = survey.qupath_pair_diagnostics(
        _frame(),
        donor="D",
        pair=survey.ROSTER.pairs[-1],
    )

    assert row["target_column"] == "Iba1_norm"
    assert row["reference_column"] == "Ker8_18"
    assert row["measured_n"] == 7
    assert row["target_negative_n"] == 3
    assert row["reference_negative_n"] == 0
    assert row["target_zero_n"] == 1
    assert not any("positive" in column for column in row)
    assert "threshold" not in row


def test_input_summary_retains_signed_qupath_norm():
    summary = survey.summarize_input_features(_frame(), donor="D", roster=survey.ROSTER)
    iba1 = summary.set_index("column").loc["Iba1_norm"]

    assert iba1["source_measurement"] == "Cell: Iba1: Mean (norm)"
    assert iba1["finite_n"] == 7
    assert iba1["negative_n"] == 3
    assert iba1["minimum"] == pytest.approx(-0.4)
    assert iba1["maximum"] == pytest.approx(1.2)


def test_missing_qupath_iba1_norm_never_falls_back(tmp_path):
    _declare_projection(tmp_path / "input")
    donor = tmp_path / "input" / "donor_id=D"
    donor.mkdir(parents=True)
    frame = _frame().drop(columns="Iba1_norm")
    frame["Iba1"] = 99.0
    frame["Iba1__restore_normalized"] = 2.0
    frame.to_parquet(donor / "data.parquet", index=False)

    with pytest.raises(ValueError, match="Iba1_norm"):
        survey.load_qupath_donor(tmp_path / "input", "D", survey.ROSTER)


def test_lookalike_legacy_norm_table_without_projection_contract_is_refused(tmp_path):
    donor = tmp_path / "legacy" / "donor_id=D"
    donor.mkdir(parents=True)
    _frame().to_parquet(donor / "data.parquet", index=False)

    with pytest.raises(ValueError, match="not a declared upstream QuPath"):
        survey.load_qupath_donor(tmp_path / "legacy", "D", survey.ROSTER)


@pytest.mark.parametrize("dtype", [np.uint16, np.int32, np.float32, np.float64])
def test_donor_loader_canonicalizes_numeric_partner_storage_to_float64(tmp_path, dtype):
    frame = _frame().dropna().reset_index(drop=True)
    _declare_projection(tmp_path / "input", row_count=len(frame))
    donor = tmp_path / "input" / "donor_id=D"
    donor.mkdir(parents=True)
    frame["Ker8_18"] = frame["Ker8_18"].astype(dtype)
    frame.to_parquet(donor / "data.parquet", index=False)

    loaded, _ = survey.load_qupath_donor(tmp_path / "input", "D", survey.ROSTER)

    assert all(
        loaded[column].dtype == np.dtype("float64") for column in survey.ROSTER.source_columns
    )
    assert loaded["Iba1_norm"].tolist() == frame["Iba1_norm"].astype(float).tolist()
    assert loaded["Ker8_18"].tolist() == frame["Ker8_18"].astype(float).tolist()


@pytest.mark.parametrize("bad", [True, "12"])
def test_donor_loader_rejects_non_numeric_partner_storage(tmp_path, bad):
    _declare_projection(tmp_path / "input")
    donor = tmp_path / "input" / "donor_id=D"
    donor.mkdir(parents=True)
    frame = _frame()
    frame["Ker8_18"] = bad
    frame.to_parquet(donor / "data.parquet", index=False)

    with pytest.raises(ValueError, match="integer or floating numeric storage"):
        survey.load_qupath_donor(tmp_path / "input", "D", survey.ROSTER)


def test_donor_loader_rejects_infinite_partner(tmp_path):
    _declare_projection(tmp_path / "input")
    donor = tmp_path / "input" / "donor_id=D"
    donor.mkdir(parents=True)
    frame = _frame()
    frame.loc[0, "Ker8_18"] = np.inf
    frame.to_parquet(donor / "data.parquet", index=False)

    with pytest.raises(ValueError, match="contains infinite"):
        survey.load_qupath_donor(tmp_path / "input", "D", survey.ROSTER)


def test_mutable_workspace_replaces_same_label_without_touching_siblings(tmp_path):
    root = tmp_path / "development"
    output = survey.prepare_development_workspace(root, "diagnostic")
    (output / "old.txt").write_text("stale", encoding="utf-8")
    sibling = root / "keep"
    sibling.mkdir()
    (sibling / "user.txt").write_text("keep", encoding="utf-8")

    refreshed = survey.prepare_development_workspace(root, "diagnostic")

    assert refreshed == output
    assert not (refreshed / "old.txt").exists()
    assert (sibling / "user.txt").read_text(encoding="utf-8") == "keep"
    with pytest.raises(ValueError, match="one nonempty path component"):
        survey.prepare_development_workspace(root, "../escape")


def test_cli_writes_input_diagnostics_but_no_assignments(tmp_path):
    input_root = tmp_path / "input"
    _declare_projection(input_root)
    donor = input_root / "donor_id=D"
    donor.mkdir(parents=True)
    _frame().to_parquet(donor / "data.parquet", index=False)
    development = tmp_path / "development"
    argv = [
        "--donors",
        "D",
        "--input-root",
        str(input_root),
        "--development-root",
        str(development),
        "--label",
        "current",
    ]

    assert survey.main(argv) == 0
    output = development / "current"
    (output / "stale.txt").write_text("old", encoding="utf-8")
    assert survey.main(argv) == 0

    assert not (output / "stale.txt").exists()
    assert not (output / "assignments").exists()
    assert not (output / "fundamental_type_summary.parquet").exists()
    diagnostics = pd.read_parquet(output / "fundamental_pair_diagnostics.parquet")
    assert diagnostics["target_column"].tolist() == ["CD3e", "Iba1_norm"]
    assert not any("positive" in column for column in diagnostics)
    assert (output / "fundamental_input_summary.parquet").is_file()
    metadata = json.loads((output / "development_summary.json").read_text(encoding="utf-8"))
    assert metadata["restore_outputs_read"] is False
    assert metadata["creates_cell_calls"] is False
    assert "Cell: Iba1: Mean (norm)" in metadata["source_semantics"]
