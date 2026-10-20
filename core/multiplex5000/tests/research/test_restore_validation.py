import numpy as np

from phenocycler.restore_validation import (
    PairValidationConfig,
    divisor_reproducibility_across_seeds,
)


def test_reproducibility_is_the_min_over_max_of_usable_seed_divisors():
    count, value = divisor_reproducibility_across_seeds([100.0, 95.0, 90.0])
    assert count == 3
    assert value == 0.9


def test_seeds_without_a_usable_divisor_are_excluded_from_the_ratio():
    """A collapsed seed contributes nothing; it must not widen or narrow the ratio."""

    count, value = divisor_reproducibility_across_seeds(
        [100.0, np.nan, 90.0, 0.0, -5.0]
    )
    assert count == 2
    assert value == 0.9


def test_fewer_than_two_usable_seeds_is_not_measured_and_fails_the_gate():
    """The gate gets NaN, not 1.0.

    Reproducibility across seeds is undefined with one seed.  Reporting 1.0
    there is perfect agreement asserted by fiat, and it let a pair whose fit
    collapsed on two of three seeds reach production on the third seed's
    divisor.  NaN is the honest value and fails the stability comparison.
    """

    threshold = PairValidationConfig().min_divisor_reproducibility

    for divisors in ([200.0, np.nan, np.nan], [np.nan, np.nan, np.nan], []):
        count, value = divisor_reproducibility_across_seeds(divisors)
        assert count < 2
        assert np.isnan(value)
        # This is the expression evaluate_locked_pair gates `stable` on.
        assert not (np.isfinite(value) and value >= threshold)


def test_perfectly_reproducible_seeds_still_pass():
    count, value = divisor_reproducibility_across_seeds([133.8, 133.8, 133.8])
    assert count == 3
    assert value == 1.0
    assert value >= PairValidationConfig().min_divisor_reproducibility
