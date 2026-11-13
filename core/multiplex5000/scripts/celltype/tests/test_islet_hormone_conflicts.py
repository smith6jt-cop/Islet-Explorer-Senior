import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import run_islet as islet
import run_sameslide_protein as bridge


def test_only_two_positive_coleaders_go_to_native_review():
    values = np.array([[1, 0, 1], [1, 0, 1], [1, 1, 1], [.2, 0, .8],
                       [0, 0, 0], [np.nan, np.nan, np.nan]])
    ids = ["a", "b", "c", "d", "e", "f"]

    def resolve(groups):
        assert groups == {"a": ("INS", "SST"), "b": ("INS", "SST")}
        return {"a": {"winner": "SST", "status": "typed"},
                "b": {"winner": None, "status": "insufficient_area"}}

    result = islet.subtype_endocrine(values, ["INS", "GCG", "SST"], .1, ids, resolve)
    assert list(result.call) == ["Endo-delta", "Endo-unresolved", "Endo-unresolved",
                                 "Endo-delta", "Endo-hormone-negative", "Endo-unresolved"]
    assert result.loc[0, "runner_up"] == "Endo-beta"
    assert result.loc[1, "hormone_status"] == "insufficient_area"
    assert islet.subtype_endocrine(values[:1], ["INS", "GCG", "SST"], .1).call[0] == "Endo-unresolved"


def test_unresolved_positive_conflict_retains_endocrine_prior():
    values = np.array([[0, 0], [1, 1], [1, .2]])
    core = np.array([False, True, True])
    out = bridge._islet_core_cascade(values, ["SST", "INS"], core, [], .1, None,
                                     np.full(3, ""), set())
    assert list(out.index) == [1, 2]
    assert list(out.call) == ["Endo-unresolved", "Endo-delta"]
    assert list(out.reason) == ["hormone_conflict_unresolved", "hormone_argmax"]
    assert out.reason.isin(bridge.ENDOCRINE_REASONS).all()
