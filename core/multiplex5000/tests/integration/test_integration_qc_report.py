"""Geometry-only reports must not resurrect stale identity results."""

from __future__ import annotations

import pandas as pd

from phenocycler.config import PipelineConfig
from phenocycler.integration.qc import concordance_summary, sanity_flags


def test_geometry_report_states_that_identity_sections_are_omitted():
    stale = {"donor": pd.DataFrame({"spearman_composition_ev": [0.99]})}
    report = concordance_summary(stale, identity_attached=False)
    assert "identity not attached — composition sections omitted" in report
    assert "0.990" not in report


def test_missing_mirror_measurement_is_not_reported_as_mirrored(tmp_path):
    qc = pd.DataFrame({
        "scale": [float("nan")],
        "max_disp_um": [float("nan")],
        "mirrored": [float("nan")],
    })
    out = sanity_flags(qc, PipelineConfig(data_dir=tmp_path))
    assert out.loc[0, "sanity_flags"] == ""
