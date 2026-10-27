"""Compare CD8 and CD4 within an already established protein T-cell lineage.

The broad T-cell rule supplies the lineage. This comparison has no second
positivity cutoff. Values must be in the same marker-normalized units within
the physical section; compare before clipping so two values at a scale floor
or bright end are not made artificially equal. An exact tie needs review of
the whole marker profile, not an arbitrary subtype.
"""

from __future__ import annotations

import numpy as np

NOT_T_CELL = "not_t_cell"
UNAVAILABLE = "unavailable"
CD8 = "cd8"
CD4 = "cd4"
REVIEW_TIE = "review_total_marker_profile"


def unbounded_scale(values: np.ndarray, *, floor: float, span: float) -> np.ndarray:
    """Put one measured marker in background-to-bright-end units without clipping."""
    raw = np.asarray(values, dtype=np.float64)
    if raw.ndim != 1:
        raise ValueError("marker values must be one-dimensional")
    if not np.isfinite(floor) or not np.isfinite(span) or span <= 0:
        return np.full(len(raw), np.nan)
    return (raw - floor) / span


def compare_cd8_cd4(cd8: np.ndarray, cd4: np.ndarray,
                    t_cell: np.ndarray) -> np.ndarray:
    """Return one auditable comparison state per cell, with no 0.10 cutoff."""
    a, b = np.asarray(cd8, dtype=np.float64), np.asarray(cd4, dtype=np.float64)
    selected = np.asarray(t_cell)
    if (a.ndim != 1 or a.shape != b.shape or selected.shape != a.shape
            or selected.dtype != np.bool_):
        raise ValueError("CD8, CD4 and T-cell masks must align one-to-one")
    result = np.full(len(a), NOT_T_CELL, dtype=object)
    result[selected] = UNAVAILABLE
    measured = selected & np.isfinite(a) & np.isfinite(b)
    result[measured & (a > b)] = CD8
    result[measured & (b > a)] = CD4
    result[measured & (a == b)] = REVIEW_TIE
    return result
