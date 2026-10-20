"""Production simultaneous broad typing over RESTORE evidence.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  Simultaneous / hierarchical / probabilistic broad typing over
RESTORE evidence is on the never-revive list (root ``CLAUDE.md``), RESTORE itself was retired on
2026-08-14, and the ``type`` stage is blocked by ``NO_VERIFIED_TYPING``.  The typing that exists is
the subpanel classification under ``core/multiplex5000/scripts/typing/``.  Kept as a record.

The production path evaluates every broad gate against RESTORE's
tri-state evidence: every gate is evaluated against the same evidence, and a
cell supported by more than one gate is ``Ambiguous`` rather than awarded to
whichever gate happened to be listed first.
"""

from __future__ import annotations

import pandas as pd

from .artifacts import ContractError
from .marker_registry import NO_VERIFIED_TYPING_STATUS, MarkerRegistry, load_registry
from .simultaneous_typing import (
    METHOD_VERSION as SIMULTANEOUS_METHOD_VERSION,
    BroadGatePolicy,
    type_simultaneous_cells,
)
METHOD_VERSION = f"production-default-{SIMULTANEOUS_METHOD_VERSION}"


def type_production_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    broad_gate_policy: BroadGatePolicy | None = None,
) -> pd.DataFrame:
    """Run the sole production typing implementation."""

    registry = load_registry() if marker_registry is None else marker_registry
    if registry.production_status == NO_VERIFIED_TYPING_STATUS:
        raise ContractError(
            "production typing is blocked by NO_VERIFIED_TYPING: the marker "
            "registry declares no verified production pairs"
        )
    return type_simultaneous_cells(
        evidence_wide,
        marker_registry=registry,
        policy=broad_gate_policy,
    )


__all__ = ["METHOD_VERSION", "type_production_cells"]
