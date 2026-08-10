from __future__ import annotations

import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.marker_registry import load_registry
from phenocycler.ordered_residual_typing import (
    TYPING_MODE,
    load_ordered_residual_policy,
    type_ordered_residual_cells,
)
from phenocycler.restore_validation import ACCEPTED_PAIRS


def _evidence(rows: list[tuple[str, bool]]) -> pd.DataFrame:
    registry = load_registry()
    frame = pd.DataFrame(
        {
            "donor_id": ["D"] * len(rows),
            "object_id": [row[0] for row in rows],
            "qc_analysis_eligible": [row[1] for row in rows],
        }
    )
    for marker in registry.active_markers:
        frame[f"{marker.name}__state"] = "negative"
    return frame


def test_policy_is_the_exact_historical_gate_order_and_pair_target_set():
    registry = load_registry()
    policy = load_ordered_residual_policy(registry=registry)

    assert policy.compartment_order == (
        "Immune",
        "Endothelial",
        "Epithelial",
        "Neural",
        "Mesenchymal",
    )
    assert set(policy.markers) == {target for target, _ in ACCEPTED_PAIRS}
    assert len(policy.markers) == 24


def test_first_positive_gate_claims_the_residual_population():
    frame = _evidence(
        [
            ("immune_before_epithelial", True),
            ("endothelial_before_epithelial", True),
            ("epithelial_before_neural", True),
            ("mesenchymal", True),
            ("other", True),
        ]
    )
    frame.loc[0, ["CD3e__state", "E_cadherin__state"]] = "positive"
    frame.loc[1, ["CD31__state", "INS__state"]] = "positive"
    frame.loc[2, ["GCG__state", "B3TUBB__state"]] = "positive"
    frame.loc[3, "SMA__state"] = "positive"

    out = type_ordered_residual_cells(frame).set_index("object_id")

    assert out.loc["immune_before_epithelial", "broad_type"] == "Immune"
    assert out.loc["endothelial_before_epithelial", "broad_type"] == "Endothelial"
    assert out.loc["epithelial_before_neural", "broad_type"] == "Epithelial"
    assert out.loc["mesenchymal", "broad_type"] == "Mesenchymal"
    assert out.loc["other", "broad_type"] == "Other"
    assert out.loc["other", "assignment_reason"] == "all_negative"
    assert out["typing_mode"].eq(TYPING_MODE).all()


def test_unavailable_gate_blocks_lower_priority_but_positive_sibling_wins():
    frame = _evidence(
        [
            ("blocked_immune", True),
            ("immune_positive_wins", True),
            ("blocked_endothelial", True),
            ("geometry_excluded", False),
        ]
    )
    frame.loc[0, "CD3e__state"] = "unavailable"
    frame.loc[0, "E_cadherin__state"] = "positive"
    frame.loc[1, ["CD3e__state", "CD68__state"]] = [
        "unavailable",
        "positive",
    ]
    frame.loc[2, "CD31__state"] = "indeterminate"
    frame.loc[2, "E_cadherin__state"] = "positive"
    frame.loc[3, "CD3e__state"] = "positive"

    out = type_ordered_residual_cells(frame).set_index("object_id")

    assert out.loc["blocked_immune", "broad_type"] == "Unavailable"
    assert out.loc["blocked_immune", "blocked_by"] == "Immune"
    assert (
        out.loc["blocked_immune", "assignment_reason"]
        == "blocked_unavailable_gate"
    )
    assert out.loc["immune_positive_wins", "broad_type"] == "Immune"
    assert out.loc["blocked_endothelial", "broad_type"] == "Unavailable"
    assert out.loc["blocked_endothelial", "blocked_by"] == "Endothelial"
    assert out.loc["geometry_excluded", "broad_type"] == "Unavailable"
    assert (
        out.loc["geometry_excluded", "assignment_reason"]
        == "geometry_qc_ineligible"
    )


def test_missing_or_unknown_restore_state_fails_closed():
    frame = _evidence([("cell", True)])
    with pytest.raises(ContractError, match="missing columns"):
        type_ordered_residual_cells(frame.drop(columns="CD3e__state"))

    frame.loc[0, "CD3e__state"] = "maybe"
    with pytest.raises(ContractError, match="unknown cell states"):
        type_ordered_residual_cells(frame)
