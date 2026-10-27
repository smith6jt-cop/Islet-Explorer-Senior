"""The old pancreas lane and Stage 06 use one packaged marker-scale method."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

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
