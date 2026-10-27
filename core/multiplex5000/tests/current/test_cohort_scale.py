"""The old pancreas lane and Stage 06 use one packaged marker-scale method."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from multiplex5000 import cohort_scale


def test_one_scale_spans_donors_and_keeps_missing_measurements_missing():
    values = {"CD8": np.array([10.0, 20.0, 40.0, np.nan])}
    scale = cohort_scale.marker_scale(values, {"CD8": 10.0}, q=100.0)
    assert scale.loc["CD8", "floor"] == 10.0
    assert scale.loc["CD8", "span"] == 30.0
    fitted = cohort_scale.apply_scale(
        np.array([[10.0], [25.0], [40.0], [np.nan]]), ["CD8"], scale)
    assert fitted[:3, 0].tolist() == [0.0, 0.5, 1.0]
    assert np.isnan(fitted[3, 0])


def test_historical_lane_delegates_its_scale_to_the_package():
    source = Path(__file__).resolve().parents[2] / "scripts/celltype/score.py"
    spec = importlib.util.spec_from_file_location("celltype_score_scale_check", source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.marker_scale is cohort_scale.marker_scale
    assert module.apply_scale is cohort_scale.apply_scale


def test_cohort_scale_uses_the_block_and_is_independent_of_parquet_row_order(tmp_path):
    paths = {}
    for donor, offset in (("A", 0.0), ("B", 100.0)):
        path = tmp_path / f"{donor}.parquet"
        pd.DataFrame({
            "object_id": ["a", "b", "c", "d", "e"],
            "qc_analysis_eligible": [True, True, True, True, False],
            "in_alignment_block": [True, True, True, False, True],
            "CD8": [10 + offset, 20 + offset, 30 + offset, 50_000, 60_000],
            "CD57": [np.nan] * 5,
        }).to_parquet(path, index=False)
        paths[donor] = path
    first = cohort_scale.build_scale(paths, ["CD8", "CD57"], n_per_donor=2, seed=7)
    for path in paths.values():
        pd.read_parquet(path).sample(frac=1, random_state=123).to_parquet(path, index=False)
    shuffled = cohort_scale.build_scale(paths, ["CD8", "CD57"], n_per_donor=2, seed=7)
    pd.testing.assert_frame_equal(first, shuffled)
    pd.testing.assert_frame_equal(
        first, cohort_scale.build_scale(dict(reversed(list(paths.items()))),
                                        ["CD8", "CD57"], n_per_donor=2, seed=7))
    assert first.loc["CD8", "top"] < 200, "other-block and ineligible cells never set the top"
    assert np.isnan(first.loc["CD57", "top"]), "unavailable is not converted to zero"


def test_cohort_scale_treats_absent_marker_column_as_unavailable(tmp_path):
    path = tmp_path / "donor.parquet"
    pd.DataFrame({"object_id": ["a", "b"], "qc_analysis_eligible": [True] * 2,
                  "in_alignment_block": [True] * 2, "CD8": [10.0, 20.0]}).to_parquet(path)
    scale = cohort_scale.build_scale({"donor": path}, ["CD8", "CD38"])
    assert np.isnan(scale.loc["CD38", "top"])
    assert scale.loc["CD8", "top"] > 0
