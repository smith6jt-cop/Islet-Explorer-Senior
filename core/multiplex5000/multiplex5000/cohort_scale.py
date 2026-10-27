"""One cohort-wide marker scale for the pancreas protein typing vector.

UniFORM puts a marker in a common unit across donors. This scale then clips
that marker at its measured background and robust bright end before the
hierarchical protein method compares markers within a cell. It is fitted once
on the aligned cohort, never separately per donor.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SCALE_Q = 99.9


def marker_scale(values_by_marker: dict, floor_by_marker: dict | None = None,
                 q: float = SCALE_Q) -> pd.DataFrame:
    """Cohort-wide floor and top per marker, so every marker enters the angle comparably."""
    rows = []
    for marker, values in values_by_marker.items():
        data = np.asarray(values, dtype=float)
        data = data[np.isfinite(data)]
        floor = float(floor_by_marker.get(marker, 0.0)) if floor_by_marker else 0.0
        top = float(np.percentile(data, q)) if len(data) else np.nan
        rows.append({"marker": marker, "floor": floor, "top": top,
                     "span": top - floor if np.isfinite(top) else np.nan})
    return pd.DataFrame(rows).set_index("marker")


def apply_scale(values: np.ndarray, markers: list, scale: pd.DataFrame) -> np.ndarray:
    """(cells x markers) to [0, 1], clipped at the measured floor and robust top."""
    floor = scale.loc[markers, "floor"].to_numpy(float)
    span = scale.loc[markers, "span"].to_numpy(float)
    span = np.where(np.isfinite(span) & (span > 0), span, 1.0)
    return np.clip((np.asarray(values, dtype=np.float32) - floor) / span, 0.0, 1.0)
