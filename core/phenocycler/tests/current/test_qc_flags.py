"""Synthetic tests for the v3 reason-coded QC flags (sibling of geometry_qc)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.qc_flags import (
    QC_FLAG_REASON_ORDER,
    QCFlagConfig,
    evaluate_qc_flags,
    summarize_qc_flags,
)


def _table(n: int = 20, **overrides) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    values = {
        "object_id": [f"c{i}" for i in range(n)],
        "cell_area_um2": rng.normal(60.0, 5.0, n),
        "cell_solidity": rng.normal(0.9, 0.02, n),
        "nucleus_cell_area_ratio": rng.normal(0.4, 0.05, n),
        "dapi_mean": rng.normal(500.0, 50.0, n),
        "ecadherin_mean": rng.normal(200.0, 20.0, n),
        "channel_a_max": np.full(n, 1000.0),
        "channel_b_max": np.full(n, 2000.0),
        "channel_a_mean": rng.normal(300.0, 30.0, n),
        "channel_b_mean": rng.normal(150.0, 15.0, n),
        "channel_c_mean": rng.normal(80.0, 8.0, n),
        "signed_distance_to_tissue_um": rng.normal(200.0, 50.0, n),
        "cell_region": ["tissue"] * n,
    }
    values.update(overrides)
    return pd.DataFrame(values)


def _config(**overrides) -> QCFlagConfig:
    return QCFlagConfig(**overrides)


def _evaluate(table, **kwargs):
    defaults = dict(
        donor_id="D1",
        config=_config(),
        cell_area_column="cell_area_um2",
        solidity_column="cell_solidity",
        nucleus_cell_ratio_column="nucleus_cell_area_ratio",
        dapi_column="dapi_mean",
        ecadherin_column="ecadherin_mean",
        channel_max_columns=["channel_a_max", "channel_b_max"],
        channel_mean_columns=["channel_a_mean", "channel_b_mean", "channel_c_mean"],
        distance_to_tissue_column="signed_distance_to_tissue_um",
        region_column="cell_region",
    )
    defaults.update(kwargs)
    return evaluate_qc_flags(table, **defaults)


def test_every_reason_has_its_own_switchable_boolean_column():
    result = _evaluate(_table())
    frame = result.to_frame()
    for name in QC_FLAG_REASON_ORDER:
        assert name in frame.columns
        assert frame[name].dtype == "boolean"
    assert "qc_flag_reasons" in frame.columns
    assert (frame["donor_id"] == "D1").all()


def test_a_planted_area_outlier_is_caught_by_the_tukey_fence():
    table = _table(n=30)
    table.loc[0, "cell_area_um2"] = 5000.0  # far above any reasonable donor fence
    result = _evaluate(table)
    frame = result.to_frame()
    assert frame.loc[frame.object_id == "c0", "tukey_outlier_log_cell_area"].iloc[0]
    assert (
        "tukey_outlier_log_cell_area"
        in frame.loc[frame.object_id == "c0", "qc_flag_reasons"].iloc[0]
    )
    # a cell near the donor's own mean is not flagged
    assert not frame.loc[frame.object_id == "c1", "tukey_outlier_log_cell_area"].iloc[0]


def test_solidity_and_nucleus_ratio_outliers_are_caught_independently():
    table = _table(n=30)
    table.loc[1, "cell_solidity"] = 0.01
    table.loc[2, "nucleus_cell_area_ratio"] = 5.0
    frame = _evaluate(table).to_frame()
    assert frame.loc[frame.object_id == "c1", "tukey_outlier_solidity"].iloc[0]
    assert frame.loc[frame.object_id == "c2", "tukey_outlier_nucleus_cell_ratio"].iloc[0]
    # neither planted outlier trips the other's flag
    assert not frame.loc[frame.object_id == "c1", "tukey_outlier_nucleus_cell_ratio"].iloc[0]
    assert not frame.loc[frame.object_id == "c2", "tukey_outlier_solidity"].iloc[0]


def test_dapi_and_ecadherin_top_quantile_flags_are_independent():
    table = _table(n=100)
    table.loc[0, "dapi_mean"] = 1e6
    frame = _evaluate(table, config=_config(dapi_top_quantile=0.95)).to_frame()
    assert frame.loc[frame.object_id == "c0", "dapi_top_quantile"].iloc[0]
    assert not frame.loc[frame.object_id == "c0", "ecadherin_top_quantile"].iloc[0]


def test_saturated_cell_is_flagged_and_saturation_is_unavailable_without_max_columns():
    table = _table(n=10)
    table.loc[3, "channel_a_max"] = 65535
    result = _evaluate(table)
    frame = result.to_frame()
    assert frame.loc[frame.object_id == "c3", "channel_saturated"].iloc[0]
    assert frame.loc[frame.object_id != "c3", "channel_saturated"].eq(False).all()

    no_max = _evaluate(table, channel_max_columns=[]).to_frame()
    assert no_max["channel_saturated"].isna().all()


def test_a_custom_saturation_value_and_edge_sign_are_honoured():
    table = _table(n=5)
    table.loc[0, "channel_a_max"] = 4095
    flagged = _evaluate(table, config=_config(saturation_value=4095)).to_frame()
    assert flagged.loc[flagged.object_id == "c0", "channel_saturated"].iloc[0]


def test_edge_proximal_cells_are_flagged_within_the_configured_distance():
    table = _table(n=5, signed_distance_to_tissue_um=[2.0, 9.9, 10.1, 50.0, -1.0])
    frame = _evaluate(table).to_frame()
    assert frame.loc[0, "edge_proximal"]
    assert frame.loc[1, "edge_proximal"]
    assert not frame.loc[2, "edge_proximal"]
    assert not frame.loc[3, "edge_proximal"]
    # inside_sign=1 means negative distance is outside the tissue, not "inside" edge-proximal
    assert not frame.loc[4, "edge_proximal"]


def test_edge_distance_and_inside_sign_are_parameterizable():
    table = _table(n=3, signed_distance_to_tissue_um=[-2.0, -30.0, 100.0])
    frame = _evaluate(table, config=_config(edge_distance_um=5.0, inside_sign=-1)).to_frame()
    assert frame.loc[0, "edge_proximal"]
    assert not frame.loc[1, "edge_proximal"]
    assert not frame.loc[2, "edge_proximal"]


def test_artifact_region_membership_is_read_verbatim_from_parent():
    table = _table(n=4, cell_region=["tissue", "artifact", "pLN", "artifact"])
    frame = _evaluate(table).to_frame()
    assert not frame.loc[0, "artifact_region"]
    assert frame.loc[1, "artifact_region"]
    assert not frame.loc[2, "artifact_region"]
    assert frame.loc[3, "artifact_region"]


def test_a_null_region_is_unavailable_not_false():
    table = _table(n=3, cell_region=["tissue", None, "artifact"])
    frame = _evaluate(table).to_frame()
    assert not frame.loc[0, "artifact_region"]
    assert pd.isna(frame.loc[1, "artifact_region"])
    assert frame.loc[2, "artifact_region"]


def test_pan_channel_autofluorescence_flags_a_cell_bright_on_every_channel():
    table = _table(n=50)
    # Bright on every channel at once, not just one -- the autofluorescence signature.
    table.loc[0, "channel_a_mean"] = 5000.0
    table.loc[0, "channel_b_mean"] = 5000.0
    table.loc[0, "channel_c_mean"] = 5000.0
    frame = _evaluate(table, config=_config(autofluorescence_quantile=0.9)).to_frame()
    assert frame.loc[frame.object_id == "c0", "pan_channel_autofluorescence"].iloc[0]


def test_disabling_a_flag_makes_it_unavailable_everywhere():
    table = _table(n=10)
    table.loc[0, "cell_area_um2"] = 1e6
    frame = _evaluate(table, config=_config(enable_tukey_log_cell_area=False)).to_frame()
    assert frame["tukey_outlier_log_cell_area"].isna().all()


def test_a_missing_optional_column_is_unavailable_never_false():
    table = _table(n=10).drop(columns=["dapi_mean"])
    frame = _evaluate(table, dapi_column="dapi_mean").to_frame()
    assert frame["dapi_top_quantile"].isna().all()
    # every other flag remains assessable
    assert frame["tukey_outlier_log_cell_area"].notna().any()


def test_missing_object_id_column_is_a_contract_error():
    with pytest.raises(ContractError, match="object_id"):
        evaluate_qc_flags(_table().drop(columns=["object_id"]), donor_id="D1", config=_config())


def test_empty_table_is_a_contract_error():
    with pytest.raises(ContractError, match="empty"):
        evaluate_qc_flags(_table().iloc[0:0], donor_id="D1", config=_config())


def test_config_rejects_a_non_positive_tukey_k():
    with pytest.raises(ContractError, match="tukey_k"):
        QCFlagConfig(tukey_k=0.0)


def test_config_rejects_an_out_of_range_quantile():
    with pytest.raises(ContractError, match="dapi_top_quantile"):
        QCFlagConfig(dapi_top_quantile=1.0)


def test_summary_table_reports_counts_and_fences_per_donor():
    table = _table(n=30)
    table.loc[0, "cell_area_um2"] = 5000.0
    result = _evaluate(table)
    summary = summarize_qc_flags(result)
    assert summary.loc[0, "donor_id"] == "D1"
    assert summary.loc[0, "n_cells"] == 30
    assert summary.loc[0, "tukey_outlier_log_cell_area_flagged"] >= 1
    assert "fence_log_cell_area_high" in summary.columns

    stacked = summarize_qc_flags([result, result])
    assert len(stacked) == 2


def test_result_is_deterministic_across_repeated_evaluation():
    table = _table(n=200)
    table.loc[0, "cell_area_um2"] = 9000.0
    first = _evaluate(table).to_frame()
    second = _evaluate(table).to_frame()
    pd.testing.assert_frame_equal(first, second)


def test_row_order_does_not_change_the_per_donor_fence():
    table = _table(n=40)
    table.loc[5, "cell_area_um2"] = 4000.0
    forward = _evaluate(table).to_frame().sort_values("object_id").reset_index(drop=True)
    shuffled = _evaluate(table.sample(frac=1.0, random_state=1).reset_index(drop=True)).to_frame()
    shuffled = shuffled.sort_values("object_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(forward, shuffled)
