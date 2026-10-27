"""Shared T-cell marker comparison in an installed MultiPlex5000 package."""

import numpy as np

import pandas as pd

from multiplex5000.tissues.t_cell import (
    clipped_t_tie_profile, compare_cd8_cd4, unbounded_scale,
)


def test_t_lineage_comparison_uses_preclip_units_and_reports_ties():
    cd8 = unbounded_scale(np.array([0.5, 0.4, np.nan, 1.0, 4.0]),
                          floor=1.0, span=10.0)
    cd4 = unbounded_scale(np.array([0.4, 0.5, 2.0, 1.0, 1.0]),
                          floor=1.0, span=10.0)
    states = compare_cd8_cd4(cd8, cd4,
                             np.array([True, True, True, True, False]))
    assert states.tolist() == ["cd8", "cd4", "unavailable",
                               "review_total_marker_profile", "not_t_cell"]


def test_clipped_t_tie_uses_full_preclip_identity_profile():
    markers = ["CD3e", "CD8", "CD4", "CD45", "CD68"]
    raw = pd.DataFrame({"CD3e": [2., 2., np.nan], "CD8": [11., 11., 11.],
                        "CD4": [12., 12., 12.], "CD45": [1., 1., 1.],
                        "CD68": [20., 3., np.nan]})
    clipped = raw.clip(0, 10).div(10)
    scale = pd.DataFrame({"floor": [0.] * len(markers),
                          "span": [10.] * len(markers)}, index=markers)
    tied, top, withdrawn = clipped_t_tie_profile(clipped, raw, scale, markers)
    assert tied.tolist() == [True, True, True]
    assert top.tolist() == ["CD68", "CD4", "CD4"]
    assert withdrawn.tolist() == [True, False, False]
