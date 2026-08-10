"""Production typing dispatcher.

The label-free production default is the validated RESTORE ordered residual
hierarchy.  The separately signed probabilistic path remains available only
when an immutable promoted model bundle is explicitly configured.
"""

from __future__ import annotations

import pandas as pd

from .hierarchical_typing import AssignmentThresholds, TypingRegistry
from .marker_registry import MarkerRegistry
from .ordered_residual_typing import (
    METHOD_VERSION as ORDERED_METHOD_VERSION,
    OrderedResidualPolicy,
    type_ordered_residual_cells,
)
from .typing_model_bundle import TypingModelBundle, TypingModelPromotionManifest
from .typing_stage import type_calibrated_cells as type_hierarchical_cells


METHOD_VERSION = f"production-default-{ORDERED_METHOD_VERSION}"


def type_production_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    ordered_policy: OrderedResidualPolicy | None = None,
    typing_registry: TypingRegistry | None = None,
    thresholds: AssignmentThresholds | None = None,
    model_bundle: TypingModelBundle | None = None,
    model_promotion_manifest: TypingModelPromotionManifest | None = None,
) -> pd.DataFrame:
    """Run ordered RESTORE typing unless a signed model is explicitly supplied."""

    if model_bundle is None:
        if model_promotion_manifest is not None:
            raise ValueError(
                "an ordered RESTORE run cannot configure a probabilistic "
                "promotion manifest"
            )
        if thresholds is not None:
            raise ValueError(
                "ordered RESTORE typing has no probabilistic thresholds"
            )
        return type_ordered_residual_cells(
            evidence_wide,
            marker_registry=marker_registry,
            policy=ordered_policy,
        )
    return type_hierarchical_cells(
        evidence_wide,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        thresholds=thresholds,
        model_bundle=model_bundle,
        model_promotion_manifest=model_promotion_manifest,
    )


__all__ = ["METHOD_VERSION", "type_production_cells"]
