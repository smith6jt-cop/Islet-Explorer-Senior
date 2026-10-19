import numpy as np

from analysis.cohesive.crossstudy import shares, state


def test_states_follow_the_cohort_level_and_unmeasured_stays_unmeasured():
    tcf1 = np.array([0.5, 0.5, 0.0, 0.0, np.nan])
    tox = np.array([0.0, 0.3, 0.3, 0.0, 0.5])
    assert list(state(tcf1, tox)) == ["TCF1+TOX-", "TCF1+TOX+", "TCF1-TOX+", "TCF1-TOX-", "unmeasured"]


def test_shares_use_the_published_denominator():
    # Figure 7 divides by T cells positive for TCF1 and/or TOX: double negatives and unmeasured leave it.
    out = shares(np.array(["TCF1+TOX-", "TCF1+TOX+", "TCF1+TOX+", "TCF1-TOX+", "TCF1-TOX-", "unmeasured"]))
    assert out["n"] == 4
    assert out["TCF1+TOX+"] == 50.0 and out["TCF1+TOX-"] == 25.0 and out["TCF1-TOX+"] == 25.0
