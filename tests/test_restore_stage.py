import numpy as np
import pandas as pd
from pathlib import Path

from phenocycler.marker_registry import load_registry
from phenocycler.simultaneous_typing import type_simultaneous_cells
from phenocycler.restore_stage import (
    CELL_NEGATIVE,
    CELL_POSITIVE,
    CELL_UNAVAILABLE,
    RestoreStageConfig,
    build_restore_input_frame,
    cell_calls_from_evaluation,
    production_pairs,
    restore_marker_evidence,
)
from phenocycler.restore_validation import (
    ExpressionReview,
    PairState,
    load_expression_reviews,
)


def test_production_restore_pairs_come_from_active_registry():
    pairs = production_pairs(load_registry())

    assert len(pairs) == 24
    assert len({target for target, _ in pairs}) == 24
    assert ("CD11b", "EpCAM") in pairs
    assert ("CD11b", "E_cadherin") not in pairs


def test_versioned_expression_review_preserves_accepted_6566_ins_absence():
    path = (
        Path(__file__).resolve().parents[1]
        / "phenocycler"
        / "restore_expression_reviews.csv"
    )
    reviews, audit = load_expression_reviews(path)

    assert reviews == {("6566", "INS"): ExpressionReview.CONFIRMED_ABSENT}
    assert audit.loc[0, "evidence"].strip()


def test_build_restore_input_frame_respects_geometry_and_raw_compartments():
    cells = pd.DataFrame(
        {
            "object_id": ["a", "b", "c"],
            "image": ["slide"] * 3,
            "X_centroid": [1.0, 2.0, 3.0],
            "Y_centroid": [4.0, 5.0, 6.0],
            "E_cadherin": [10.0, 20.0, 30.0],
            "CD31": [40.0, 50.0, 60.0],
        }
    )
    redsea = pd.DataFrame(
        {
            "object_id": ["c", "a", "b"],
            "Cell__E_cadherin": [3.0, 1.0, 2.0],
            "Cell__CD31": [6.0, 4.0, 5.0],
        }
    )
    geometry = pd.DataFrame(
        {
            "object_id": ["b", "c", "a"],
            "analysis_eligible": [True, False, True],
        }
    )
    sample = pd.DataFrame(
        {
            "object_id": ["a", "b"],
            "Membrane__E_cadherin": [7.0, 8.0],
            "Membrane__CD31": [9.0, 10.0],
        }
    )

    result = build_restore_input_frame(
        cells,
        redsea,
        geometry,
        donor_id="115",
        pairs=(("E_cadherin", "CD31"),),
        compartment_sample=sample,
    )

    assert result["object_id"].tolist() == ["a", "b", "c"]
    assert result["redsea__E_cadherin"].tolist()[:2] == [1.0, 2.0]
    assert result.loc[0, "Membrane__CD31"] == 9.0
    assert np.isnan(result.loc[2, "raw__CD31"])
    assert np.isnan(result.loc[2, "redsea__CD31"])
    assert np.isnan(result.loc[2, "Cell__CD31"])
    assert result.attrs["analysis_eligible_n"] == 2
    assert result.attrs["compartment_sample_covered_n"] == 2


def test_cell_calls_are_tri_state_and_never_make_failure_negative():
    valid = {
        "n_input": 3,
        "state": PairState.VALID.value,
        "valid_idx": np.array([0, 1]),
        "full_lineage_positive": np.array([False, True, False]),
        "full_normalized": np.array([0.5, 2.0, np.nan]),
    }
    positive, normalized, state = cell_calls_from_evaluation(
        valid, analysis_eligible=np.array([True, True, False])
    )
    assert positive.tolist() == [False, True, False]
    assert normalized[:2].tolist() == [0.5, 2.0]
    assert state.tolist() == [
        CELL_NEGATIVE,
        CELL_POSITIVE,
        CELL_UNAVAILABLE,
    ]

    failure = {"n_input": 2, "state": PairState.MODEL_UNSTABLE.value}
    _, _, failed_state = cell_calls_from_evaluation(failure)
    assert failed_state.tolist() == [CELL_UNAVAILABLE, CELL_UNAVAILABLE]


def test_confirmed_absence_never_calls_an_unmeasured_cell_negative():
    """A reviewed absence is a statement about the donor, not about a missing measurement.

    Cell 2 is outside ``valid_idx``: it had no finite measurement for this
    pair.  A confirmed absence of target expression says nothing about it, so
    it must stay ``unavailable`` rather than becoming a confident negative.
    """

    confirmed = {
        "n_input": 3,
        "state": PairState.NO_TARGET_EXPRESSION_CONFIRMED.value,
        "valid_idx": np.array([0, 1]),
        # normalize_for_assessment returns all-NaN / all-False for this state.
        "full_lineage_positive": np.array([False, False, False]),
        "full_normalized": np.array([np.nan, np.nan, np.nan]),
    }
    positive, normalized, state = cell_calls_from_evaluation(
        confirmed, analysis_eligible=np.array([True, True, True])
    )
    assert state.tolist() == [
        CELL_NEGATIVE,
        CELL_NEGATIVE,
        CELL_UNAVAILABLE,
    ]
    assert not positive.any()
    assert np.isnan(normalized).all()


def test_restore_evidence_drives_rules_typing_without_labels_or_model(monkeypatch):
    registry = load_registry()
    pairs = production_pairs(registry)
    donor_frame = pd.DataFrame(
        {
            "object_id": ["immune", "other"],
            "qc_analysis_eligible": [True, True],
            **{
                f"redsea__{target}": [2.0, 0.5]
                for target, _ in pairs
            },
        }
    )

    def evaluator(frame, donor, target, reference, **kwargs):
        is_cd3 = target == "CD3e"
        return {
            "n_input": 2,
            "state": PairState.VALID.value,
            "state_reason": "",
            "valid_idx": np.array([0, 1]),
            "full_lineage_positive": np.array([is_cd3, False]),
            "full_normalized": np.array(
                [2.0 if is_cd3 else 0.5, 0.5]
            ),
        }

    def compact_metrics(evaluation):
        return {
            "state": evaluation["state"],
            "canonical_divisor": 1.0,
            "reference_n": 20,
        }

    monkeypatch.setattr(
        "phenocycler.restore_stage.locked_pair_metrics", compact_metrics
    )
    result = restore_marker_evidence(
        donor_frame,
        donor_id="115",
        registry=registry,
        evaluator=evaluator,
    )

    assert result.evidence_wide["CD3e__state"].astype(str).tolist() == [
        CELL_POSITIVE,
        CELL_NEGATIVE,
    ]
    assert result.model_table.shape[0] == 24

    assignments = type_simultaneous_cells(
        result.evidence_wide,
        marker_registry=registry,
    ).set_index("object_id")
    assert assignments.loc["immune", "broad_label"] == "Immune"
    assert assignments.loc["immune", "broad_assignment_status"] == "anchor"
    assert assignments.loc["other", "broad_label"] == "Other"
    assert assignments.loc["other", "broad_assignment_status"] == "other"


def test_pair_exception_becomes_unavailable_when_configured():
    registry = load_registry()
    pairs = production_pairs(registry)
    donor_frame = pd.DataFrame(
        {
            "object_id": ["a"],
            "qc_analysis_eligible": [True],
            **{f"redsea__{target}": [1.0] for target, _ in pairs},
        }
    )

    def broken(*args, **kwargs):
        raise ValueError("synthetic fit failure")

    result = restore_marker_evidence(
        donor_frame,
        donor_id="115",
        registry=registry,
        config=RestoreStageConfig(on_technical_failure="unavailable"),
        evaluator=broken,
    )

    assert set(result.model_table["state"]) == {
        PairState.TECHNICAL_FAILURE.value
    }
    assert all(
        result.evidence_wide[f"{target}__state"].iloc[0]
        == CELL_UNAVAILABLE
        for target, _ in pairs
    )
