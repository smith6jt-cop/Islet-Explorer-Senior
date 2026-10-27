"""Shared T-cell marker comparison in an installed MultiPlex5000 package."""

import numpy as np

from multiplex5000.tissues.t_cell import compare_cd8_cd4, unbounded_scale


def test_t_lineage_comparison_uses_preclip_units_and_reports_ties():
    cd8 = unbounded_scale(np.array([0.5, 0.4, np.nan, 1.0, 4.0]),
                          floor=1.0, span=10.0)
    cd4 = unbounded_scale(np.array([0.4, 0.5, 2.0, 1.0, 1.0]),
                          floor=1.0, span=10.0)
    states = compare_cd8_cd4(cd8, cd4,
                             np.array([True, True, True, True, False]))
    assert states.tolist() == ["cd8", "cd4", "unavailable",
                               "review_total_marker_profile", "not_t_cell"]
