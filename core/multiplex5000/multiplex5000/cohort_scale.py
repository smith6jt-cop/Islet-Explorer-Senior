"""One cohort-wide marker scale for the pancreas protein typing vector.

UniFORM puts a marker in a common unit across donors. This scale then clips
that marker at its measured background and robust bright end before the
hierarchical protein method compares markers within a cell. It is fitted once
on the aligned cohort, never separately per donor.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

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


def build_scale(aligned: dict[str, Path], markers: list[str], *, n_per_donor: int = 200_000,
                seed: int = 0) -> pd.DataFrame:
    """Fit one scale from QC-eligible cells in the recorded pancreas block.

    This is the selected October method's 512-bin background mode and p99.9
    bright end, with its 200,000-cell-per-donor ceiling. Sampling first orders
    by native ``object_id`` so a Parquet rewrite cannot choose different cells.
    A missing or refused marker remains unavailable and yields no usable span.
    """
    if not aligned or not markers or len(markers) != len(set(markers)):
        raise ValueError("cohort scale needs donor tables and unique marker names")
    if n_per_donor < 1:
        raise ValueError("n_per_donor must be positive")
    rng = np.random.default_rng(seed)
    pools: dict[str, list[np.ndarray]] = {marker: [] for marker in markers}
    columns = ["object_id", "qc_analysis_eligible", "in_alignment_block", *markers]
    for donor in sorted(aligned):
        path = aligned[donor]
        table = pq.ParquetFile(path).read(columns=columns).to_pandas()
        if table["object_id"].isna().any() or not table["object_id"].is_unique:
            raise ValueError(f"{donor}: aligned table has missing or duplicate object_id")
        sample = table.loc[table["qc_analysis_eligible"].to_numpy(bool)
                           & table["in_alignment_block"].to_numpy(bool)]
        if sample.empty:
            raise ValueError(f"{donor}: no eligible cells in the pancreas block")
        sample = sample.sort_values("object_id", kind="stable")
        if len(sample) > n_per_donor:
            sample = sample.iloc[rng.choice(len(sample), n_per_donor, replace=False)]
        for marker in markers:
            pools[marker].append(sample[marker].to_numpy("float64"))
    rows = []
    for marker, chunks in pools.items():
        combined = np.concatenate(chunks)
        positive = combined[np.isfinite(combined) & (combined > 0)]
        floor = None
        if len(positive) >= 1000:
            counts, edges = np.histogram(np.log1p(positive), bins=512)
            floor = float(np.expm1(edges[1:][int(np.argmax(counts[1:])) + 1]))
        rows.append(marker_scale({marker: combined},
                                 {marker: floor} if floor is not None else None))
    return pd.concat(rows)
