"""The two SOM configurations pass distinct, explicit inputs to the shared worker boundary."""
from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "arms"))

import arm_base
import pins
import som_kintsugi as ranked
import som_kintsugi_aligned as aligned


@pytest.mark.parametrize("ablate", [None, "INS"])
def test_aligned_then_ranked_preserve_inputs_and_do_not_change_each_other(tmp_path, monkeypatch, ablate):
    monkeypatch.setattr(pins, "OUT", tmp_path)
    monkeypatch.setattr(ranked, "SCRATCH_ROOT", tmp_path / "scratch" / ranked.ARM)
    monkeypatch.setattr(ranked, "DEFAULT_PARAMS", copy.deepcopy(ranked.DEFAULT_PARAMS))
    monkeypatch.setattr(aligned, "DEFAULT_PARAMS", copy.deepcopy(aligned.DEFAULT_PARAMS))
    ranked_params = copy.deepcopy(ranked.DEFAULT_PARAMS)
    ranked_scratch = ranked.SCRATCH_ROOT
    X = pd.DataFrame(np.arange(192, dtype=float).reshape(64, 3) + 100,
                     columns=["CD3e", "INS", "CD31"])
    X.loc[0, "CD3e"] = np.nan
    if ablate:
        X = X.drop(columns=ablate)
    ids = np.array([f"cell-{i}" for i in range(len(X))])

    class WorkerReached(Exception):
        pass

    # Run the real input preparation, stopping at the foreign-environment boundary.
    for method in (aligned, ranked):
        ctx = arm_base.ArmContext(
            donor_id="synthetic", scope="frame", ids=ids, X=X.copy(), markers=list(X.columns),
            context=pd.DataFrame({"object_id": ids}), ablated=ablate,
        )
        expected = X if method is aligned else arm_base.within_donor_rank(X)
        captured = []

        def worker(features, worker_ids, fits, scratch):
            captured.append(scratch)
            np.testing.assert_array_equal(worker_ids, ids[1:])
            pd.testing.assert_frame_equal(features, expected.iloc[1:].reset_index(drop=True))
            assert scratch.parent.parent == tmp_path / "scratch" / method.ARM
            raise WorkerReached

        monkeypatch.setattr(ranked, "run_worker", worker)
        with pytest.raises(WorkerReached):
            method.run(ctx)
        assert len(captured) == 1
        pd.testing.assert_frame_equal(ctx.X, X)
        assert ctx.markers == list(X.columns)
        result = method.DEFAULT_PARAMS["result"]
        assert result is ctx.extra["som_kintsugi"]
        assert result["identity_markers"] == list(X.columns)
        assert result["n_non_finite"] == 1
        assert result["non_finite_by_marker"] == {"CD3e": 1}
        assert ranked.ARM == "som-kintsugi"
        assert ranked.SCRATCH_ROOT == ranked_scratch
        assert ranked.arm_base is arm_base
        if method is aligned:
            assert ranked.DEFAULT_PARAMS == ranked_params


def test_shared_som_pipeline_refuses_reordered_identity_features(tmp_path):
    X = pd.DataFrame({"CD3e": [1.0, 2.0], "INS": [3.0, 4.0]})
    ctx = arm_base.ArmContext(
        donor_id="synthetic", scope="frame", ids=np.array(["a", "b"]), X=X,
        markers=list(X.columns), context=pd.DataFrame({"object_id": ["a", "b"]}),
    )
    for features in (X.iloc[::-1], X[["INS", "CD3e"]]):
        with pytest.raises(ValueError, match="rows and identity marker order"):
            ranked.run_features(ctx, features, params={}, scratch_root=tmp_path)
