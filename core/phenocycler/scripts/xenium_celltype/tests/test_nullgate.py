"""The gate is only worth having if it refuses. These check that it does."""
import sys, pathlib
import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from nullgate import (NullGateError, annotate_groups, detection_deciles,  # noqa: E402
                      MIN_GROUP_CELLS)

GENES = [f"g{i}" for i in range(200)]


def _counts(n_cells=400, seed=0, hot=None, hot_rows=None):
    """Sparse counts like Xenium: ~96% zeros, a few genes hot in a few cells."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.04, size=(n_cells, len(GENES))).astype(float)
    if hot is not None:
        X[np.ix_(hot_rows, hot)] += rng.poisson(6.0, size=(len(hot_rows), len(hot)))
    return X


def test_a_real_lineage_clears_its_null():
    hot = list(range(0, 10))
    rows = np.arange(0, 200)
    X = _counts(hot=hot, hot_rows=rows)
    groups = np.where(np.arange(400) < 200, "A", "B")
    cores = {"Real": [GENES[i] for i in hot], "Other": [GENES[i] for i in range(100, 110)]}
    out = annotate_groups(X, GENES, groups, cores, n_draws=40, seed=1).set_index("group")
    assert out.loc["A", "lineage"] == "Real"
    assert out.loc["A", "real"] > out.loc["A", "null_p95"]


def test_a_group_with_no_real_signal_is_refused_not_argmaxed():
    """The defect that motivated this module: nearest-core instead of no-call."""
    X = _counts(n_cells=600, seed=3)          # pure background, no lineage present
    groups = np.where(np.arange(600) < 300, "A", "B")   # split is arbitrary by construction
    cores = {"A": GENES[:10], "B": GENES[20:30], "C": GENES[40:48]}
    out = annotate_groups(X, GENES, groups, cores, n_draws=40, seed=2)
    assert (out.lineage == "").all()
    assert (out.reason == "no_core_beats_its_matched_null").all()


def test_ubiquitous_core_does_not_win_on_magnitude_alone():
    """A core bright in EVERY cell distinguishes no group and must clear no null.

    This is why the statistic is a contrast rather than a level: on a level, the brightest
    core wins every group. Here the real lineage must be awarded and the ubiquitous one
    refused, from the same data.
    """
    rng = np.random.default_rng(5)
    X = rng.poisson(0.04, size=(600, len(GENES))).astype(float)
    ubiq = list(range(0, 10))
    X[:, ubiq] += rng.poisson(5.0, size=(600, len(ubiq)))        # hot everywhere
    real = list(range(50, 60))
    X[np.ix_(np.arange(300), real)] += rng.poisson(6.0, size=(300, len(real)))  # hot in A only
    groups = np.where(np.arange(600) < 300, "A", "B")
    out = annotate_groups(
        X, GENES, groups,
        {"Ubiquitous": [GENES[i] for i in ubiq], "Real": [GENES[i] for i in real]},
        n_draws=40, seed=4).set_index("group")
    assert out.loc["A", "lineage"] == "Real"
    assert "Ubiquitous" not in out.loc["A", "runners_up"]


def test_one_group_cannot_be_a_contrast_and_says_so():
    """With a single group there is nothing to be specific against. Raise, never fall back
    to scoring the level — that is precisely the statistic this module rejects."""
    X = _counts()
    with pytest.raises(NullGateError, match="both sides of the group boundary"):
        annotate_groups(X, GENES, np.full(400, "only"), {"A": GENES[:10]}, n_draws=5)


def test_small_groups_are_reported_not_called():
    X = _counts(n_cells=300)
    groups = np.where(np.arange(300) < MIN_GROUP_CELLS - 1, "tiny", "big")
    out = annotate_groups(X, GENES, groups, {"A": GENES[:10]}, n_draws=5).set_index("group")
    assert out.loc["tiny", "reason"] == "too_few_cells"
    assert out.loc["tiny", "lineage"] == ""


def test_a_core_absent_from_the_panel_raises_rather_than_being_skipped():
    X = _counts()
    groups = np.where(np.arange(400) < 200, "A", "B")
    with pytest.raises(NullGateError, match="no identity-core gene on this panel"):
        annotate_groups(X, GENES, groups, {"Ghost": ["not_on_panel"]}, n_draws=5)


def test_detection_deciles_survive_massive_ties():
    X = np.zeros((100, 200)); X[0, :50] = 1     # 150 genes share detection rate 0.0
    d = detection_deciles(X)
    assert len(np.unique(d)) == 10
