"""Small invariant tests for the roster-independent empirical tail helpers."""

from __future__ import annotations

import numpy as np
import pytest

from phenocycler.empirical_tail import (
    empirical_tail_probability,
    empirical_upper_threshold,
)


def test_empirical_threshold_is_an_observed_order_statistic() -> None:
    values = np.arange(1.0, 101.0)
    threshold = empirical_upper_threshold(values[::-1], alpha=0.05)
    assert threshold == 96.0
    assert threshold in values


def test_empirical_threshold_fails_when_controls_cannot_resolve_alpha() -> None:
    with pytest.raises(ValueError, match="cannot resolve"):
        empirical_upper_threshold(np.arange(20.0), alpha=0.001)


def test_tail_probability_has_finite_sample_plus_one_correction() -> None:
    controls = np.array([1.0, 2.0, 3.0, 4.0])
    observed = empirical_tail_probability([0.0, 2.0, 5.0], controls)
    np.testing.assert_allclose(observed, [1.0, 0.8, 0.2])

