"""Characterization of the input floors and of the method that replaces them.

The floors had NO test coverage at all before this file, which is why the
p99.5 clip defect survived. These tests pin current behaviour -- including the
behaviour that is wrong -- so that retiring the floors is a measurable change
rather than an assertion.

Two halves:

* ``TestOtsuFloorAsShipped`` pins ``raw_otsu_floor`` as it behaves today,
  including the defect. Tests marked DEFECT assert the WRONG answer on purpose;
  when the floors are retired these tests go with them.
* ``TestEmpiricalUpperTail`` pins the successor primitives in
  ``marker_calibration``: a finite-sample upper-tail order statistic at a
  per-marker alpha, which refuses when the control cannot resolve that alpha.
"""

from __future__ import annotations

import numpy as np
import pytest

from phenocycler.marker_calibration import (
    _empirical_order_index,
    empirical_tail_probability,
    empirical_upper_threshold,
)
from phenocycler.marker_registry import load_registry
from phenocycler.restore_validation import (
    LINEAR_OTSU_FLOOR,
    MARKER_INPUT_FLOOR_METHODS,
    MAX_LINEAR_OTSU_TRIANGLE_FLOOR,
    estimate_marker_input_floor,
    pair_input_qc,
    raw_otsu_floor,
    raw_triangle_floor,
)

N_CELLS = 100_000
BACKGROUND_SCALE = 50.0
POSITIVE_SCALE = 200.0


def _bimodal(fraction_positive: float, seed: int = 0):
    """A marker with an explicit, well-separated positive population."""
    rng = np.random.default_rng(seed)
    n_positive = int(N_CELLS * fraction_positive)
    background = rng.gamma(2.0, BACKGROUND_SCALE, N_CELLS - n_positive)
    positive = rng.gamma(20.0, POSITIVE_SCALE, n_positive)
    return np.concatenate([background, positive]), n_positive


class TestOtsuFloorAsShipped:
    def test_otsu_separates_an_abundant_population(self):
        """At 30% prevalence Otsu is correct, and that is why it survived."""
        values, n_positive = _bimodal(0.30)
        floor = raw_otsu_floor(values)
        called = int((values > floor).sum())
        assert 0.9 < called / n_positive < 1.1

    def test_otsu_separates_a_five_percent_population(self):
        values, n_positive = _bimodal(0.05)
        floor = raw_otsu_floor(values)
        called = int((values > floor).sum())
        assert 0.9 < called / n_positive < 1.2

    def test_DEFECT_otsu_collapses_a_population_rarer_than_the_clip(self):
        """The p99.5 clip erases a sparse positive mode, so Otsu splits background.

        ``raw_otsu_floor`` clips its input at ``upper_quantile=0.995`` and then
        histograms over ``range=(0, upper)``. A population rarer than 0.5% is
        therefore entirely inside the terminal bin, contributes no between-class
        structure, and the argmax lands inside the negative mode.

        This asserts the WRONG answer deliberately. Retire it with the floors.
        """
        values, n_positive = _bimodal(0.003)
        floor = raw_otsu_floor(values)
        clip_point = float(np.quantile(values, 0.995))
        positive_mode_start = float(np.quantile(values[-n_positive:], 0.01))

        # The entire positive population sits above the clip point.
        assert clip_point < positive_mode_start
        # And the floor lands in the background, admitting an order of magnitude
        # more cells than exist. Measured at ~89x; asserted loosely as >20x so
        # the test pins the failure mode, not a fragile constant.
        called = int((values > floor).sum())
        assert called > 20 * n_positive

    def test_triangle_is_not_clipped_and_survives_the_sparse_case(self):
        """Which is why ten markers were moved to max(otsu, triangle)."""
        values, n_positive = _bimodal(0.003)
        assert (values > raw_triangle_floor(values)).sum() < 10 * n_positive

    def test_floor_method_selection_matches_the_registry_policy(self):
        values, _ = _bimodal(0.05)
        otsu_marker = next(
            m for m, k in MARKER_INPUT_FLOOR_METHODS.items() if k == LINEAR_OTSU_FLOOR
        )
        both_marker = next(
            m
            for m, k in MARKER_INPUT_FLOOR_METHODS.items()
            if k == MAX_LINEAR_OTSU_TRIANGLE_FLOOR
        )
        otsu_only = estimate_marker_input_floor(values, otsu_marker)
        assert otsu_only.value == otsu_only.linear_otsu
        assert otsu_only.linear_triangle is None

        combined = estimate_marker_input_floor(values, both_marker)
        assert combined.value == max(combined.linear_otsu, combined.linear_triangle)

    def test_an_unconfigured_marker_fails_closed(self):
        values, _ = _bimodal(0.05)
        with pytest.raises(ValueError, match="no input-floor method"):
            estimate_marker_input_floor(values, "not_a_marker")

    def test_pair_input_qc_retains_the_union_of_both_arms(self):
        target, _ = _bimodal(0.05, seed=1)
        reference, _ = _bimodal(0.05, seed=2)
        retained, target_floor, reference_floor = pair_input_qc(
            target, reference, target_marker="CD3e", reference_marker="E_cadherin"
        )
        expected = (target > target_floor.value) | (reference > reference_floor.value)
        assert np.array_equal(retained, expected)


class TestEmpiricalUpperTail:
    """The successor: an order statistic on the pair's own control population."""

    def test_the_threshold_is_an_order_statistic_of_the_control(self):
        rng = np.random.default_rng(0)
        control = np.log1p(rng.gamma(2.0, BACKGROUND_SCALE, 30_000))
        threshold = empirical_upper_threshold(control, 0.001)
        assert threshold in set(control.tolist())

    def test_it_refuses_when_the_control_cannot_resolve_alpha(self):
        """No observation can reach p <= alpha when alpha < 1/(n+1)."""
        assert _empirical_order_index(200, 0.001) is None
        assert _empirical_order_index(999, 0.001) == 998

        rng = np.random.default_rng(0)
        control = np.log1p(rng.gamma(2.0, BACKGROUND_SCALE, 200))
        with pytest.raises(ValueError, match="cannot resolve marker alpha"):
            empirical_upper_threshold(control, 0.001)

    def test_it_holds_across_prevalences_that_break_otsu(self):
        """0.3%, 5% and 30% positives, one rule, no per-case tuning."""
        rng = np.random.default_rng(0)
        control = np.log1p(rng.gamma(2.0, BACKGROUND_SCALE, 30_000))
        for fraction in (0.003, 0.05, 0.30):
            values, n_positive = _bimodal(fraction)
            threshold = empirical_upper_threshold(control, 0.001)
            called = int((np.log1p(values) > threshold).sum())
            assert 0.8 < called / n_positive < 1.5, f"failed at {fraction:.1%}"

    def test_tail_probability_is_never_zero_and_is_monotone(self):
        rng = np.random.default_rng(0)
        control = np.sort(np.log1p(rng.gamma(2.0, BACKGROUND_SCALE, 5_000)))
        probabilities = empirical_tail_probability(
            np.array([control[0] - 1.0, control[len(control) // 2], 1e6]), control
        )
        assert (probabilities > 0).all()
        assert probabilities[0] > probabilities[1] > probabilities[2]

    def test_alpha_is_a_familywise_allocation_not_a_tuned_constant(self):
        registry = load_registry()
        for name, family in registry.families.items():
            assert family.member_alpha == family.familywise_alpha / len(family.members)
        # Each broad gate shares one 0.01 familywise budget.
        assert registry.marker_alpha("CD3e") == pytest.approx(0.01 / 10)
        assert registry.marker_alpha("CD31") == pytest.approx(0.01 / 3)

    def test_DRIFT_the_immune_family_still_carries_CD11b(self):
        """CD11b left the Immune gate in 3192a26 but not the alpha family.

        The nine markers that do gate Immune are therefore corrected as if there
        were ten, i.e. slightly more conservatively than the 0.01 familywise
        budget intends. Conservative, so not urgent -- but it is drift, and the
        allocation is a scientific decision, so it is pinned here rather than
        silently changed.
        """
        from phenocycler.simultaneous_typing import load_broad_gate_policy

        registry = load_registry()
        policy = load_broad_gate_policy(registry=registry)
        family = set(registry.families["broad_immune"].members)
        assert family - set(policy.gates["Immune"]) == {"CD11b"}
