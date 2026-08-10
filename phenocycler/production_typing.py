"""Production typing dispatcher.

The label-free production default is simultaneous broad typing over RESTORE's
tri-state evidence: every gate is evaluated against the same evidence, and a
cell supported by more than one gate is ``Ambiguous`` rather than awarded to
whichever gate happened to be listed first.  The separately signed
probabilistic path remains available only when an immutable promoted model
bundle is explicitly configured.
"""

from __future__ import annotations

import pandas as pd

from .hierarchical_typing import AssignmentThresholds, TypingRegistry
from .marker_registry import MarkerRegistry
from .simultaneous_typing import (
    METHOD_VERSION as SIMULTANEOUS_METHOD_VERSION,
    BroadGatePolicy,
    type_simultaneous_cells,
)
from .typing_model_bundle import TypingModelBundle, TypingModelPromotionManifest
from .typing_stage import type_calibrated_cells as type_hierarchical_cells


METHOD_VERSION = f"production-default-{SIMULTANEOUS_METHOD_VERSION}"


def type_production_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    broad_gate_policy: BroadGatePolicy | None = None,
    typing_registry: TypingRegistry | None = None,
    thresholds: AssignmentThresholds | None = None,
    model_bundle: TypingModelBundle | None = None,
    model_promotion_manifest: TypingModelPromotionManifest | None = None,
) -> pd.DataFrame:
    """Run simultaneous broad typing unless a signed model is explicitly supplied."""

    if model_bundle is None:
        if model_promotion_manifest is not None:
            raise ValueError(
                "a simultaneous RESTORE run cannot configure a probabilistic "
                "promotion manifest"
            )
        if thresholds is not None:
            raise ValueError(
                "simultaneous RESTORE typing has no probabilistic thresholds"
            )
        return type_simultaneous_cells(
            evidence_wide,
            marker_registry=marker_registry,
            policy=broad_gate_policy,
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
