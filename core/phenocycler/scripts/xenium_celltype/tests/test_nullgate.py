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


def test_a_group_reading_lower_than_the_tissue_is_never_awarded_that_core():
    """`real > p95` alone awarded cores on NEGATIVE contrasts.

    If a group reads lower on a core than the rest of the tissue, and the random cores read
    lower still, the old gate handed it the lineage anyway. Measured in the shipped output:
    nonspatial k=30 group 11 took Vascular at real=-0.0134, and spatial k=30 group 12 took
    Myeloid at real=-0.0045 across 833 cells.
    """
    # core genes are hot OUTSIDE group A, so A's contrast on them is negative
    X = _counts(n_cells=600, seed=3, hot=range(10), hot_rows=range(300, 600))
    groups = np.array(["A"] * 300 + ["B"] * 300)
    ann = annotate_groups(X, GENES, groups, {"Dim": GENES[:10]}, n_draws=40, seed=0)
    a = ann.set_index("group").loc["A"]
    assert a.real < 0
    assert a.lineage == ""
    assert a.reason == "no_core_beats_its_matched_null"


def test_two_cores_that_cannot_be_told_apart_award_nothing():
    """Clearing is not winning.

    Two cores over the SAME genes must clear identically, so no z-ranking between them is
    meaningful. The paired null catches that and refuses, rather than handing the whole group
    to whichever came out 0.06 ahead -- which is what flipped four groups between k=20 and
    k=30 on donor 6539.
    """
    X = _counts(n_cells=600, seed=4, hot=range(10), hot_rows=range(300))
    groups = np.array(["A"] * 300 + ["B"] * 300)
    ann = annotate_groups(X, GENES, groups,
                          {"CoreOne": GENES[:10], "CoreTwo": list(GENES[:10])},
                          n_draws=40, seed=0)
    a = ann.set_index("group").loc["A"]
    assert a.lineage == ""
    assert a.reason == "cores_cleared_but_are_not_separable"
    assert set(a.runners_up.split("|")) == {"CoreOne", "CoreTwo"}


def test_one_bright_gene_cannot_carry_a_core():
    """Per-gene SD scaling is what stops a 58-gene core being one gene's contrast.

    On donor 6539 XBP1 alone was 39% of the Immune_B core mean and CUZD1 57% of Exocrine's,
    so the core's contrast was largely that one gene's -- a structure no random set of the
    same size shares, which is why the null could not see it.
    """
    rng = np.random.default_rng(5)
    X = _counts(n_cells=600, seed=5)
    X[np.ix_(range(300), [0])] += rng.poisson(60.0, (300, 1))            # one very loud gene
    X[np.ix_(range(300), range(20, 40))] += rng.poisson(1.5, (300, 20))  # a coherent core
    # Genes matched on detection and magnitude that carry NO group difference -- what the
    # null is meant to draw from. Without them the only bright genes on this toy panel are
    # the cores themselves, the matched pool collapses onto them, and the null inherits the
    # very signal it exists to be compared against.
    X[np.ix_(range(150, 450), range(100, 160))] += rng.poisson(1.5, (300, 60))
    groups = np.array(["A"] * 300 + ["B"] * 300)
    ann = annotate_groups(X, GENES, groups,
                          {"OneLoudGene": GENES[0:20], "Coherent": GENES[20:40]},
                          n_draws=40, seed=0)
    a = ann.set_index("group").loc["A"]
    assert a.lineage == "Coherent", f"got {a.lineage!r} runners_up={a.runners_up!r}"


def test_composition_never_renders_an_absence_as_a_zero():
    """A lineage that won nothing used to have no row, and a missing row reads as 0%.

    That is how donor 6539's spatial arm was reported as `Immune_T` 0.0%, when Immune_T had
    in fact cleared its null and been outranked in six separate groups. Three distinct states
    -- awarded, cleared-but-outranked, never-cleared -- were collapsed into one number that
    named none of them.
    """
    import pandas as pd
    from run_donor import composition

    ann = pd.DataFrame({
        "group": [0, 1, 2, 3],
        "lineage": ["Exocrine", "", "Endocrine", ""],
        "reason": ["", "no_core_beats_its_matched_null", "", "too_few_cells"],
        "runners_up": ["Immune_B", "", "", ""],
    })
    groups = np.array([0] * 500 + [1] * 300 + [2] * 150 + [3] * 50)
    c = composition(ann, groups,
                    lineages=["Exocrine", "Endocrine", "Immune_B", "Immune_T"]
                    ).set_index("lineage")

    assert c.loc["Exocrine", "status"] == "awarded"
    # outranked is not the same as absent, and neither is a zero
    assert c.loc["Immune_B", "status"] == "cleared_but_outranked"
    assert np.isnan(c.loc["Immune_B", "frac"])
    assert c.loc["Immune_T", "status"] == "no_group_cleared"
    assert np.isnan(c.loc["Immune_T", "frac"])

    # the two refusals stay apart -- one is about pooling, the other about evidence
    unavail = [i for i in c.index if str(i).startswith("Unavailable:")]
    assert sorted(unavail) == ["Unavailable:no_core_beats_its_matched_null",
                               "Unavailable:too_few_cells"]

    # and nothing is redistributed: the measured shares still account for every cell
    measured = c[c.status.isin(["awarded", "unavailable"])]
    assert measured.frac.sum() == pytest.approx(1.0)
