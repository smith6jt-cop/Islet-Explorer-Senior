"""The claims features.py rests on, tested on data whose answer is known."""
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import features                                     # noqa: E402


def _two_donors(seed=0, donor_offset=0.0, n=4000, n_markers=6):
    """Two donors, same biology; ``donor_offset`` adds a per-donor shift in marker 0."""
    rng = np.random.default_rng(seed)
    donors = np.array(["A"] * n + ["B"] * n)
    values = rng.lognormal(mean=1.0, sigma=0.8, size=(2 * n, n_markers))
    values[:, 1] *= 500.0                            # one very bright marker
    values[donors == "B", 0] += donor_offset
    return values, donors


def test_standard_scaler_is_feature_weighting_not_donor_correction():
    """The operation KINTSUGI applies must leave the donor-median spread invariant."""
    from sklearn.preprocessing import StandardScaler

    values, donors = _two_donors(donor_offset=3.0)
    scaled = StandardScaler().fit_transform(values)
    report = features.scaler_is_feature_weighting(values, donors, scaled)
    assert report["passes"], report
    assert report["max_relative_change"] < 1e-6


def test_the_check_catches_actual_donor_correction():
    """A per-donor centring is what INTENSITY_SCALE.md forbids, and the check must fail it."""
    values, donors = _two_donors(donor_offset=3.0)
    corrected = values.copy()
    for d in np.unique(donors):
        m = donors == d
        corrected[m] = corrected[m] - corrected[m].mean(axis=0)
    report = features.scaler_is_feature_weighting(values, donors, corrected)
    assert not report["passes"], "per-donor centring must be refused, not passed"


def test_within_donor_rank_is_monotone_inside_each_donor():
    values, donors = _two_donors()
    ranked = features.within_donor_rank(values[:, 0], donors)
    for d in np.unique(donors):
        m = donors == d
        order_before = values[m, 0].argsort()
        order_after = ranked[m].argsort()
        assert np.array_equal(order_before, order_after)
        assert ranked[m].min() == 0.0 and ranked[m].max() == 1.0


def test_donor_purity_flags_a_donor_pure_cluster():
    donors = np.array(["A"] * 500 + ["B"] * 500)
    mixed = np.zeros(1000, dtype=int)                # one cluster, mirrors the cohort
    assert features.donor_purity(mixed, donors)["purity"].iloc[0] == pytest.approx(0.0)

    split = np.where(donors == "A", 0, 1)            # clusters ARE the donors
    assert features.donor_purity(split, donors)["purity"].min() == pytest.approx(1.0)


def test_feature_block_excludes_what_the_registry_and_shifts_declare():
    classes = features.shift_classes()
    for arm in features.ARMS:
        markers = features.feature_markers(arm)
        assert "IAPP" not in markers                 # registry: failed_marker
        assert "DAPI" not in markers                 # registry: structural
        assert "CD38" not in markers and "b_Catenin1" not in markers   # panel variants
        for refused in classes["refused"]:
            assert refused not in markers, f"{refused} is a reagent change, never shifted"
        # CD11b is deliberately KEPT: the rule table's Myeloid gate depends on it.
        if arm != "C":
            assert "CD11b" in markers
    assert set(features.feature_markers("A")) - set(features.feature_markers("C")) \
        == set(classes["over_bar"])
