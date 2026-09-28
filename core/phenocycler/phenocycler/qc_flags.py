"""Reason-coded QC flags for the v3 QuPath export, beside geometry QC.

``geometry_qc.py`` decides three eligibility roles from a small, fixed geometry-metric table
(``GEOMETRY_METRIC_COLUMNS``) and consumes no marker intensity. The v3 export carries a great
deal more per cell -- shape statistics beyond that fixed set, per-compartment marker
Mean/Min/Max/Std.Dev., signed distance to the tissue boundary, and a Parent-derived
``cell_region`` that can read ``artifact`` -- and none of that belongs in the geometry-only
contract. This module is a sibling: it takes whatever columns the caller has already resolved
(REDSEA's own "prefer Nucleus/Cytoplasm/Membrane" choice is the caller's, not this module's) and
turns each of them into one independently switchable, reason-coded flag.

Every flag follows the same rule as the rest of this workflow: a measurement that is absent --
the column does not exist for this donor/table, or the cell's own value is not finite -- is
UNAVAILABLE (a nullable boolean, ``pandas.NA``), never ``False``. A cell that clears no fence is
``False``, not silently omitted. Nothing here derives a threshold from a desired frequency; the
Tukey fences and top-quantile fences are computed from the donor's own observed distribution, and
every fence actually used is carried on the result for audit (:meth:`QCFlagResult.to_frame`,
:func:`summarize_qc_flags`).

This module is not wired into any stage, ``build_layers.py`` node, or the Makefile. Calling it
changes no existing eligibility decision; it is additive capability only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .artifacts import ContractError

QC_FLAG_REASON_ORDER: tuple[str, ...] = (
    "tukey_outlier_log_cell_area",
    "tukey_outlier_solidity",
    "tukey_outlier_nucleus_cell_ratio",
    "dapi_top_quantile",
    "ecadherin_top_quantile",
    "channel_saturated",
    "edge_proximal",
    "artifact_region",
    "pan_channel_autofluorescence",
)

#: The nullable-boolean flag column for every reason above, in the same order.
QC_FLAG_COLUMNS: tuple[str, ...] = QC_FLAG_REASON_ORDER


@dataclass(frozen=True)
class QCFlagConfig:
    """Frozen, per-donor-fittable policy. Every flag can be switched off independently."""

    tukey_k: float = 3.0
    enable_tukey_log_cell_area: bool = True
    enable_tukey_solidity: bool = True
    enable_tukey_nucleus_cell_ratio: bool = True

    dapi_top_quantile: float = 0.99
    enable_dapi_top_quantile: bool = True
    ecadherin_top_quantile: float = 0.99
    enable_ecadherin_top_quantile: bool = True

    saturation_value: int = 65535
    enable_saturation: bool = True

    edge_distance_um: float = 10.0
    # +1 when a positive signed distance means "inside the tissue annotation"; -1 for the
    # opposite convention. Never guessed from the data -- the caller states it.
    inside_sign: int = 1
    enable_edge: bool = True

    artifact_region_value: str = "artifact"
    enable_artifact_region: bool = True

    autofluorescence_quantile: float = 0.99
    enable_autofluorescence: bool = True

    def __post_init__(self) -> None:
        if not (self.tukey_k > 0 and np.isfinite(self.tukey_k)):
            raise ContractError("tukey_k must be finite and > 0")
        for name in ("dapi_top_quantile", "ecadherin_top_quantile", "autofluorescence_quantile"):
            value = getattr(self, name)
            if not (0.0 < value < 1.0):
                raise ContractError(f"{name} must be in (0, 1)")
        if not (self.edge_distance_um > 0 and np.isfinite(self.edge_distance_um)):
            raise ContractError("edge_distance_um must be finite and > 0")
        if self.inside_sign not in (1, -1):
            raise ContractError("inside_sign must be 1 or -1")
        if not self.artifact_region_value:
            raise ContractError("artifact_region_value must be non-empty")


@dataclass(frozen=True)
class QCFlagResult:
    """Immutable per-donor result: one nullable-boolean Series per flag, plus the fences used."""

    donor_id: str
    object_id: tuple[str, ...]
    flags: Mapping[str, pd.Series]
    fences: Mapping[str, float]
    n_cells: int

    def to_frame(self) -> pd.DataFrame:
        """A new mutable convenience frame: one row per cell, one column per flag plus reasons."""
        frame = pd.DataFrame({"object_id": self.object_id})
        for name in QC_FLAG_COLUMNS:
            # `.values` keeps the pandas nullable BooleanArray; `.to_numpy()` would collapse it
            # to a plain numpy bool/object array and silently drop the "boolean" extension dtype.
            frame[name] = self.flags[name].values
        # `.fillna(False)` is only for this joined convenience string: an unavailable flag
        # contributes no reason code, but its own boolean column (above) still carries `NA`,
        # never `False`.
        reason_matrix = np.column_stack(
            [self.flags[name].fillna(False).to_numpy(dtype=bool) for name in QC_FLAG_REASON_ORDER]
        )
        reasons = [
            ";".join(name for name, hit in zip(QC_FLAG_REASON_ORDER, row, strict=True) if hit)
            for row in reason_matrix
        ]
        frame["qc_flag_reasons"] = reasons
        frame.insert(0, "donor_id", self.donor_id)
        return frame


def _nullable_bool(values: np.ndarray) -> pd.Series:
    return pd.Series(values, dtype="boolean")


def _tukey_flag(values: np.ndarray, *, k: float) -> tuple[pd.Series, float, float]:
    finite = np.isfinite(values)
    if finite.sum() < 2:
        return _nullable_bool(np.full(len(values), pd.NA, dtype=object)), float("nan"), float("nan")
    q1, q3 = np.nanpercentile(values[finite], [25, 75])
    iqr = q3 - q1
    low, high = q1 - k * iqr, q3 + k * iqr
    out = np.full(len(values), pd.NA, dtype=object)
    out[finite] = (values[finite] < low) | (values[finite] > high)
    return _nullable_bool(out), float(low), float(high)


def _top_quantile_flag(values: np.ndarray, *, quantile: float) -> tuple[pd.Series, float]:
    finite = np.isfinite(values)
    if finite.sum() < 2:
        return _nullable_bool(np.full(len(values), pd.NA, dtype=object)), float("nan")
    fence = float(np.nanquantile(values[finite], quantile))
    out = np.full(len(values), pd.NA, dtype=object)
    out[finite] = values[finite] > fence
    return _nullable_bool(out), fence


def _pan_channel_offset(channel_means: np.ndarray) -> np.ndarray:
    """Each cell's median rank across channels, in [0, 1]. Mirrors ``immune_subtype.per_cell_offset``.

    A cell that is bright on every channel at once (rather than on the one marker that would
    explain it) reads as pan-channel autofluorescence, not a positive call on any single marker.
    Ranks are taken per channel over this donor's finite values, so the statistic is invariant to
    multiplying any channel by a positive constant.
    """
    n, k = channel_means.shape
    if n == 0 or k == 0:
        return np.full(n, np.nan, dtype=np.float32)
    unit = np.full((n, k), np.nan, dtype=np.float32)
    for j in range(k):
        column = channel_means[:, j]
        finite = np.isfinite(column)
        if finite.sum() > 1:
            unit[finite, j] = (rankdata(column[finite]) / finite.sum()).astype(np.float32)
    any_finite = np.isfinite(unit).any(axis=1)
    out = np.full(n, np.nan, dtype=np.float32)
    out[any_finite] = np.nanmedian(unit[any_finite], axis=1)
    return out


def evaluate_qc_flags(
    table: pd.DataFrame,
    *,
    donor_id: str,
    config: QCFlagConfig,
    cell_area_column: str | None = "cell_area_um2",
    solidity_column: str | None = "cell_solidity",
    nucleus_cell_ratio_column: str | None = "nucleus_cell_area_ratio",
    dapi_column: str | None = None,
    ecadherin_column: str | None = None,
    channel_max_columns: Sequence[str] | None = None,
    channel_mean_columns: Sequence[str] | None = None,
    distance_to_tissue_column: str | None = None,
    region_column: str | None = "cell_region",
    object_id_column: str = "object_id",
) -> QCFlagResult:
    """Evaluate one donor's v3 export table against every switched-on flag.

    Every argument naming a column is optional: a column the caller does not name, or names but
    is absent from ``table``, makes that one flag UNAVAILABLE (nullable boolean, never ``False``)
    for every cell -- it does not raise, because a donor's v3 export legitimately may not carry
    every optional column (for example, no Max columns at all means saturation cannot be assessed
    for that donor, exactly as this workflow already treats a marker never acquired).
    """
    donor_id = str(donor_id).strip()
    if not donor_id:
        raise ContractError("donor_id is required")
    if object_id_column not in table.columns:
        raise ContractError(f"table is missing the object_id column {object_id_column!r}")
    if len(table) == 0:
        raise ContractError(f"donor {donor_id}: QC flag input table is empty")

    object_ids = tuple(str(v) for v in table[object_id_column].to_numpy())
    n = len(table)
    flags: dict[str, pd.Series] = {}
    fences: dict[str, float] = {}

    def _numeric(column: str | None) -> np.ndarray | None:
        if column is None or column not in table.columns:
            return None
        return pd.to_numeric(table[column], errors="coerce").to_numpy(dtype=float)

    def _unavailable() -> pd.Series:
        return _nullable_bool(np.full(n, pd.NA, dtype=object))

    # Tukey fences on log cell area, solidity, and nucleus/cell ratio.
    area = _numeric(cell_area_column)
    if config.enable_tukey_log_cell_area and area is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            log_area = np.where(area > 0, np.log(area), np.nan)
        (
            flags["tukey_outlier_log_cell_area"],
            fences["log_cell_area_low"],
            fences["log_cell_area_high"],
        ) = _tukey_flag(log_area, k=config.tukey_k)
    else:
        flags["tukey_outlier_log_cell_area"] = _unavailable()

    solidity = _numeric(solidity_column)
    if config.enable_tukey_solidity and solidity is not None:
        flags["tukey_outlier_solidity"], fences["solidity_low"], fences["solidity_high"] = (
            _tukey_flag(solidity, k=config.tukey_k)
        )
    else:
        flags["tukey_outlier_solidity"] = _unavailable()

    nc_ratio = _numeric(nucleus_cell_ratio_column)
    if config.enable_tukey_nucleus_cell_ratio and nc_ratio is not None:
        (
            flags["tukey_outlier_nucleus_cell_ratio"],
            fences["nucleus_cell_ratio_low"],
            fences["nucleus_cell_ratio_high"],
        ) = _tukey_flag(nc_ratio, k=config.tukey_k)
    else:
        flags["tukey_outlier_nucleus_cell_ratio"] = _unavailable()

    # Top-1% (parameterizable) DAPI and E-cadherin, per donor.
    dapi = _numeric(dapi_column)
    if config.enable_dapi_top_quantile and dapi is not None:
        flags["dapi_top_quantile"], fences["dapi_top_quantile"] = _top_quantile_flag(
            dapi, quantile=config.dapi_top_quantile
        )
    else:
        flags["dapi_top_quantile"] = _unavailable()

    ecad = _numeric(ecadherin_column)
    if config.enable_ecadherin_top_quantile and ecad is not None:
        flags["ecadherin_top_quantile"], fences["ecadherin_top_quantile"] = _top_quantile_flag(
            ecad, quantile=config.ecadherin_top_quantile
        )
    else:
        flags["ecadherin_top_quantile"] = _unavailable()

    # Saturation: any named Max column reads exactly the sensor ceiling.
    present_max_columns = [c for c in (channel_max_columns or ()) if c in table.columns]
    if config.enable_saturation and present_max_columns:
        max_values = (
            table[present_max_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        )
        saturated = np.any(max_values == config.saturation_value, axis=1)
        flags["channel_saturated"] = _nullable_bool(saturated)
    else:
        flags["channel_saturated"] = _unavailable()

    # Edge: within edge_distance_um of the tissue boundary, on the inside.
    distance = _numeric(distance_to_tissue_column)
    if config.enable_edge and distance is not None:
        signed = config.inside_sign * distance
        finite = np.isfinite(signed)
        out = np.full(n, pd.NA, dtype=object)
        out[finite] = (signed[finite] >= 0) & (signed[finite] <= config.edge_distance_um)
        flags["edge_proximal"] = _nullable_bool(out)
        fences["edge_distance_um"] = config.edge_distance_um
    else:
        flags["edge_proximal"] = _unavailable()

    # Artifact membership: the user-drawn fold/debris region, from Parent.
    if (
        config.enable_artifact_region
        and region_column is not None
        and region_column in table.columns
    ):
        region = table[region_column]
        out = np.where(region.isna(), pd.NA, region.astype(str) == config.artifact_region_value)
        flags["artifact_region"] = _nullable_bool(out)
    else:
        flags["artifact_region"] = _unavailable()

    # Pan-channel autofluorescence: median rank across every requested channel, above a
    # per-donor top-quantile fence -- the same statistic as immune_subtype.py's panel_offset,
    # applied here as a QC gate rather than a typing diagnostic.
    present_mean_columns = [c for c in (channel_mean_columns or ()) if c in table.columns]
    if config.enable_autofluorescence and len(present_mean_columns) >= 2:
        means = (
            table[present_mean_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        )
        offset = _pan_channel_offset(means)
        flags["pan_channel_autofluorescence"], fences["autofluorescence_offset"] = (
            _top_quantile_flag(offset, quantile=config.autofluorescence_quantile)
        )
    else:
        flags["pan_channel_autofluorescence"] = _unavailable()

    return QCFlagResult(
        donor_id=donor_id,
        object_id=object_ids,
        flags=flags,
        fences=fences,
        n_cells=n,
    )


def summarize_qc_flags(results: QCFlagResult | Sequence[QCFlagResult]) -> pd.DataFrame:
    """One row per donor: how many cells each flag caught, how many it could not assess, and the
    fence actually used. A flag with zero assessable cells is reported at zero, not omitted."""
    if isinstance(results, QCFlagResult):
        results = [results]
    rows = []
    for result in results:
        row: dict[str, object] = {"donor_id": result.donor_id, "n_cells": result.n_cells}
        for name in QC_FLAG_REASON_ORDER:
            series = result.flags[name]
            row[f"{name}_flagged"] = int((series == True).sum())  # noqa: E712 (nullable boolean)
            row[f"{name}_unavailable"] = int(series.isna().sum())
        for fence_name, fence_value in result.fences.items():
            row[f"fence_{fence_name}"] = fence_value
        rows.append(row)
    return pd.DataFrame(rows)


__all__ = [
    "QC_FLAG_COLUMNS",
    "QC_FLAG_REASON_ORDER",
    "QCFlagConfig",
    "QCFlagResult",
    "evaluate_qc_flags",
    "summarize_qc_flags",
]
