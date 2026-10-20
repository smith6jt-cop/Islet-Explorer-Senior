from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.restore_seeded import (
    FEATURE_TRANSFORMS,
    _fit_indices,
    _reference_stratified_fit_indices,
    evaluate_reference_seeded_pair,
    preprocess_nnmf_features,
)
from phenocycler.restore_validation import (
    PairReview,
    PairState,
    PairValidationConfig,
    negative_control_statistics,
)


def _pair_config(**overrides: object) -> PairValidationConfig:
    values: dict[str, object] = {
        "feature_compartments": ("Nucleus", "Cytoplasm", "Membrane"),
        "seeds": (0, 1, 2),
        "fit_cap": 5_000,
        "nnmf_fit_basis": "all",
        "control_membership": "component",
        "control_definition": "reference_only",
        "divisor_statistic": "empirical_alpha",
    }
    values.update(overrides)
    return PairValidationConfig(**values)


def _synthetic_pair(n_reference: int = 1_300, n_target: int = 900) -> pd.DataFrame:
    n = n_reference + n_target
    rng = np.random.default_rng(17)
    reference_group = np.arange(n) < n_reference
    target = np.where(
        reference_group,
        rng.gamma(shape=2.0, scale=0.4, size=n),
        rng.gamma(shape=8.0, scale=2.5, size=n),
    )
    reference = np.where(
        reference_group,
        rng.gamma(shape=8.0, scale=2.5, size=n),
        rng.gamma(shape=2.0, scale=0.4, size=n),
    )
    frame: dict[str, object] = {
        "object_id": [f"cell-{index:05d}" for index in range(n)],
        "qc_analysis_eligible": np.ones(n, dtype=bool),
        "qc_estimation_eligible": np.ones(n, dtype=bool),
        "selected__Target": target,
        "selected__Reference": reference,
        # Deliberately absurd legacy values: generic evaluation must not reach
        # behind its selected-input boundary.
        "legacy_raw__Target": rng.lognormal(mean=10.0, sigma=3.0, size=n),
        "legacy_raw__Reference": rng.lognormal(mean=-5.0, sigma=4.0, size=n),
    }
    for offset, compartment in enumerate(("Nucleus", "Cytoplasm", "Membrane")):
        jitter = rng.uniform(0.95, 1.05, size=n)
        frame[f"{compartment}__Target"] = (target + 0.2 + offset) * jitter
        frame[f"{compartment}__Reference"] = (
            reference + 0.2 + offset
        ) * jitter
    result = pd.DataFrame(frame)
    result.attrs["measurement_source_by_marker"] = {
        "Target": "qupath_norm",
        "Reference": "redsea",
    }
    result.attrs["input_transform_by_marker"] = {
        "Target": "qupath_norm_left_censor_zero",
        "Reference": "none",
    }
    return result


def _fixed_components(monkeypatch: pytest.MonkeyPatch, n_reference: int) -> None:
    def fake_fit(features, **kwargs):
        labels = np.zeros(len(features), dtype=int)
        labels[n_reference:] = 1
        return (
            labels,
            np.ones((2, features.shape[1]), dtype=float),
            {
                "converged": True,
                "n_iter": 10,
                "fit_n": len(kwargs["fit_indices"]),
                "reconstruction_error": 0.01,
            },
        )

    monkeypatch.setattr("phenocycler.restore_seeded.fit_nnmf_predict", fake_fit)


def test_feature_transforms_are_explicit_nonnegative_and_scale_bounded():
    features = np.array(
        [
            [1.0, 10.0],
            [2.0, 20.0],
            [4.0, 40.0],
            [8.0, 80.0],
        ]
    )
    for transform in FEATURE_TRANSFORMS:
        scaled, lower, upper = preprocess_nnmf_features(
            features,
            transform=transform,
            lower_quantile=0.0,
            upper_quantile=1.0,
        )
        assert scaled.shape == features.shape
        assert np.isfinite(scaled).all()
        assert ((scaled >= 0) & (scaled <= 1)).all()
        assert np.all(upper > lower)


def test_linear_robust_minmax_is_invariant_to_channelwise_units():
    rng = np.random.default_rng(4)
    features = rng.gamma(shape=3.0, scale=2.0, size=(500, 4))
    first, _, _ = preprocess_nnmf_features(
        features,
        transform="linear_robust_minmax",
        lower_quantile=0.005,
        upper_quantile=0.995,
    )
    second, _, _ = preprocess_nnmf_features(
        features * np.array([0.01, 2.0, 100.0, 8.0]),
        transform="linear_robust_minmax",
        lower_quantile=0.005,
        upper_quantile=0.995,
    )
    np.testing.assert_allclose(first, second, atol=1e-12)


def test_linear_upper_scale_does_not_shift_the_zero_anchor():
    features = np.array([[0.0, 2.0], [1.0, 4.0], [3.0, 8.0]])
    scaled, lower, upper = preprocess_nnmf_features(
        features,
        transform="linear_upper_scale",
        lower_quantile=0.25,
        upper_quantile=1.0,
    )
    np.testing.assert_array_equal(lower, np.zeros(2))
    np.testing.assert_allclose(upper, [3.0, 8.0])
    assert scaled[0, 0] == 0.0


def test_nested_fit_samples_change_only_sample_size():
    small = _fit_indices(10_000, 1_000, 7, pool_cap=4_000)
    medium = _fit_indices(10_000, 2_000, 7, pool_cap=4_000)
    large = _fit_indices(10_000, 4_000, 7, pool_cap=4_000)

    assert set(small).issubset(set(medium))
    assert set(medium).issubset(set(large))
    np.testing.assert_array_equal(small, np.sort(small))
    with pytest.raises(ValueError, match="at least the fit cap"):
        _fit_indices(10_000, 4_001, 7, pool_cap=4_000)


def test_reference_stratified_sample_uses_no_target_values():
    reference_seed = np.zeros(10_000, dtype=bool)
    reference_seed[:2_000] = True
    selected = _reference_stratified_fit_indices(reference_seed, 5_000, 7)

    assert len(selected) == 5_000
    assert reference_seed[selected].sum() == 2_000
    assert (~reference_seed[selected]).sum() == 3_000
    np.testing.assert_array_equal(selected, np.sort(selected))
    np.testing.assert_array_equal(
        selected,
        _reference_stratified_fit_indices(reference_seed, 5_000, 7),
    )


def test_empirical_tail_refuses_censored_zero_and_unresolved_alpha():
    with pytest.raises(ValueError, match="strictly positive"):
        negative_control_statistics(
            np.array([0.0, 1.0, 2.0]),
            "empirical_alpha",
            marker_alpha=0.1,
        )
    with pytest.raises(ValueError, match="cannot resolve marker alpha"):
        negative_control_statistics(
            np.linspace(1.0, 2.0, 200),
            "empirical_alpha",
            marker_alpha=0.001,
        )


def test_evaluator_rejects_negative_input_that_bypasses_restore_boundary():
    frame = _synthetic_pair()
    frame.loc[0, "selected__Target"] = -0.01

    with pytest.raises(ValueError, match="frame construction"):
        evaluate_reference_seeded_pair(
            frame,
            "D1",
            "Target",
            "Reference",
            config=_pair_config(),
            marker_alpha=0.001,
            reference_positive=np.ones(len(frame), dtype=bool),
        )


def test_reference_seeded_pair_uses_only_selected_inputs_and_explicit_seeds(
    monkeypatch,
):
    n_reference = 1_300
    frame = _synthetic_pair(n_reference=n_reference)
    _fixed_components(monkeypatch, n_reference)
    reference_seed = np.zeros(len(frame), dtype=bool)
    reference_seed[:1_100] = True
    # A clamped target value inside the seed/component intersection is censored
    # and excluded from both the control and the call universe.
    frame.loc[0, "selected__Target"] = 0.0

    first = evaluate_reference_seeded_pair(
        frame,
        "D1",
        "Target",
        "Reference",
        config=_pair_config(),
        marker_alpha=0.001,
        reference_positive=reference_seed,
        pair_review=PairReview.ACCEPTED,
    )
    changed_raw = frame.copy()
    changed_raw["legacy_raw__Target"] = np.linspace(1e-30, 1e30, len(frame))
    changed_raw["legacy_raw__Reference"] = changed_raw[
        "legacy_raw__Target"
    ][::-1].to_numpy()
    second = evaluate_reference_seeded_pair(
        changed_raw,
        "D1",
        "Target",
        "Reference",
        config=_pair_config(),
        marker_alpha=0.001,
        reference_positive=reference_seed,
        pair_review=PairReview.ACCEPTED,
    )

    assert first["state"] == PairState.VALID.value
    assert first["canonical_divisor"] == second["canonical_divisor"]
    np.testing.assert_array_equal(
        first["full_lineage_positive"], second["full_lineage_positive"]
    )
    control = np.asarray(first["full_reference_control"], dtype=bool)
    assert control.sum() == 1_099
    assert not control[0]
    assert not control[1_100:n_reference].any()
    assert bool(first["full_censored_below"][0])
    assert 0 not in set(np.asarray(first["call_idx"], dtype=int))
    model = first["model_row"]
    assert model["target_input_floor_method"] == "not_used"
    assert model["call_intensity_source"] == "selected_expression"
    assert model["target_measurement_source"] == "qupath_norm"
    assert model["reference_measurement_source"] == "redsea"
    assert model["target_input_transform"] == "qupath_norm_left_censor_zero"
    assert model["reference_input_transform"] == "none"
    assert model["nnmf_input_source"] == "source-routed selected compartments"
    assert model["control_source"] == "reference_controls_intersect_nnmf_component"
    assert model["divisor_statistic"] == "empirical_alpha"
    assert model["nnmf_n_components"] == 2
    assert model["nnmf_beta_loss"] == "frobenius"
    assert model["nnmf_solver"] == "cd"
    assert model["nnmf_max_iter"] == 2_000
    assert model["nnmf_fit_sample_strategy"] == "uniform"
    assert model["nnmf_fit_sample_seed_policy"] == "coupled_to_model_seed"


def test_review_fit_controls_are_forwarded_and_audited(monkeypatch):
    frame = _synthetic_pair()
    calls: list[dict[str, object]] = []

    def fake_fit(features, **kwargs):
        calls.append(dict(kwargs))
        labels = np.zeros(len(features), dtype=int)
        labels[1_300:] = 1
        return (
            labels,
            np.ones((2, features.shape[1]), dtype=float),
            {
                "converged": True,
                "n_iter": 12,
                "fit_n": len(kwargs["fit_indices"]),
                "reconstruction_error": 0.01,
            },
        )

    monkeypatch.setattr("phenocycler.restore_seeded.fit_nnmf_predict", fake_fit)
    reference_seed = np.zeros(len(frame), dtype=bool)
    reference_seed[:1_100] = True
    result = evaluate_reference_seeded_pair(
        frame,
        "D1",
        "Target",
        "Reference",
        config=_pair_config(fit_cap=1_000),
        marker_alpha=0.001,
        reference_positive=reference_seed,
        nnmf_max_iter=7_500,
        fit_sample_seeds=(9, 9, 9),
        fit_sample_pool_cap=2_000,
    )

    assert result["state"] == PairState.PAIR_REVIEW_PENDING.value
    assert [call["seed"] for call in calls] == [0, 1, 2]
    assert {call["max_iter"] for call in calls} == {7_500}
    assert all(len(call["fit_indices"]) == 1_000 for call in calls)
    for call in calls[1:]:
        np.testing.assert_array_equal(calls[0]["fit_indices"], call["fit_indices"])
    assert [row["sample_seed"] for row in result["stability_rows"]] == [9, 9, 9]
    assert {row["n_iter"] for row in result["stability_rows"]} == {12}
    model = result["model_row"]
    assert model["nnmf_max_iter"] == 7_500
    assert model["nnmf_fit_sample_seed_policy"] == "explicit"
    assert model["nnmf_fit_sample_seeds"] == "[9,9,9]"
    assert model["nnmf_fit_sample_pool_cap"] == 2_000


def test_reference_stratified_review_fit_is_audited(monkeypatch):
    frame = _synthetic_pair()
    calls: list[dict[str, object]] = []

    def fake_fit(features, **kwargs):
        calls.append(dict(kwargs))
        labels = np.zeros(len(features), dtype=int)
        labels[1_300:] = 1
        return (
            labels,
            np.ones((2, features.shape[1]), dtype=float),
            {
                "converged": True,
                "n_iter": 10,
                "fit_n": len(kwargs["fit_indices"]),
                "reconstruction_error": 0.01,
            },
        )

    monkeypatch.setattr("phenocycler.restore_seeded.fit_nnmf_predict", fake_fit)
    reference_seed = np.zeros(len(frame), dtype=bool)
    reference_seed[:1_100] = True
    result = evaluate_reference_seeded_pair(
        frame,
        "D1",
        "Target",
        "Reference",
        config=_pair_config(fit_cap=1_000),
        marker_alpha=0.001,
        reference_positive=reference_seed,
        fit_sample_strategy="reference_stratified",
    )

    for call in calls:
        fit_indices = np.asarray(call["fit_indices"], dtype=int)
        assert reference_seed[fit_indices].sum() == 500
    assert result["model_row"]["nnmf_fit_sample_strategy"] == "reference_stratified"


def test_excess_reference_control_contamination_fails_closed(monkeypatch):
    n_reference = 1_300
    frame = _synthetic_pair(n_reference=n_reference)
    _fixed_components(monkeypatch, n_reference)
    reference_seed = np.zeros(len(frame), dtype=bool)
    reference_seed[:1_100] = True
    frame.loc[:109, "selected__Target"] = 100.0

    result = evaluate_reference_seeded_pair(
        frame,
        "D1",
        "Target",
        "Reference",
        config=_pair_config(),
        marker_alpha=0.001,
        reference_positive=reference_seed,
    )

    assert result["state"] == PairState.TECHNICAL_FAILURE.value
    assert "contamination" in result["state_reason"]
    assert result["canonical_divisor"] is None
