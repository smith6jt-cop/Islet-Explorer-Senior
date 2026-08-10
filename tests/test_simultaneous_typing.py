import json

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.marker_registry import load_registry
from phenocycler.simultaneous_typing import (
    BroadGatePolicy,
    load_broad_gate_policy,
    type_simultaneous_cells,
)


def _policy():
    return load_broad_gate_policy(registry=load_registry())


def _evidence(rows):
    """rows: list of {marker: (state, normalized)} -> the wide evidence contract."""

    policy = _policy()
    frame = {
        "donor_id": [f"D{i}" for i in range(len(rows))],
        "object_id": [f"c{i}" for i in range(len(rows))],
        "qc_analysis_eligible": [row.pop("__eligible", True) for row in rows],
    }
    for marker in policy.markers:
        frame[f"{marker}__state"] = [
            row.get(marker, ("negative", 0.0))[0] for row in rows
        ]
        frame[f"{marker}__restore_normalized"] = [
            row.get(marker, ("negative", 0.0))[1] for row in rows
        ]
    return pd.DataFrame(frame)


def test_a_cell_supported_by_two_gates_is_ambiguous_not_awarded_to_one():
    """The defect this module exists to remove.

    CD3e (Immune) and E_cadherin (Epithelial) both positive is contradictory
    evidence about one cell. The ordered residual gave it to whichever gate was
    listed first; here it stays Ambiguous with both gates recorded.
    """

    out = type_simultaneous_cells(
        _evidence([{"CD3e": ("positive", 4.0), "E_cadherin": ("positive", 1.1)}])
    )
    assert out.loc[0, "broad_label"] == "Ambiguous"
    assert out.loc[0, "broad_assignment_status"] == "ambiguous"
    assert out.loc[0, "broad_supported_gates"] == "Epithelial|Immune"
    assert out.loc[0, "broad_supported_gate_n"] == 2
    # The support margin is reported so a later, separately versioned rule can
    # use it -- but no winner is chosen here.
    assert out.loc[0, "broad_support_margin"] == pytest.approx(2.9)
    assert out.loc[0, "broad_anchor_classes"] == ""


def test_a_single_supported_gate_is_an_anchor():
    out = type_simultaneous_cells(
        _evidence([{"CD3e": ("positive", 4.0)}])
    )
    assert out.loc[0, "broad_label"] == "Immune"
    assert out.loc[0, "broad_assignment_status"] == "anchor"
    assert out.loc[0, "broad_authoritative_markers"] == "CD3e"
    assert out.loc[0, "broad_supported_gate_n"] == 1


def test_no_support_with_an_unmeasured_marker_is_unavailable_not_other():
    """Absence of evidence is not evidence of absence."""

    out = type_simultaneous_cells(
        _evidence([{"CD3e": ("unavailable", float("nan"))}])
    )
    assert out.loc[0, "broad_label"] == "Unavailable"
    assert out.loc[0, "broad_reason"] == "unmeasured_gate_marker"
    assert "Immune" in out.loc[0, "blocked_by"]


def test_a_positive_sibling_beats_an_unmeasured_sibling_in_the_same_gate():
    out = type_simultaneous_cells(
        _evidence([{"CD3e": ("positive", 3.0), "CD68": ("unavailable", float("nan"))}])
    )
    assert out.loc[0, "broad_label"] == "Immune"
    assert out.loc[0, "broad_assignment_status"] == "anchor"


def test_everything_explicitly_negative_is_other():
    out = type_simultaneous_cells(_evidence([{}]))
    assert out.loc[0, "broad_label"] == "Other"
    assert out.loc[0, "broad_reason"] == "all_negative"


def test_geometry_ineligible_cells_are_unavailable():
    out = type_simultaneous_cells(
        _evidence([{"CD3e": ("positive", 9.0), "__eligible": False}])
    )
    assert out.loc[0, "broad_label"] == "Unavailable"
    assert out.loc[0, "broad_reason"] == "geometry_qc_ineligible"


def test_scores_that_were_never_estimated_are_nan_not_one():
    """The ordered lane stamped 1.0 into probability, margin and stability.

    Those are categorical gate outcomes, and the same rows recorded zero
    bootstrap replicates. They must be NaN.
    """

    out = type_simultaneous_cells(
        _evidence([{"CD3e": ("positive", 4.0)}, {}])
    )
    for column in (
        "broad_best_probability",
        "broad_margin",
        "broad_stability",
        "confidence",
    ):
        assert out[column].isna().all(), column
    assert (out["typing_broad_bootstrap_replicates"] == 0).all()


def test_assignments_are_invariant_to_gate_listing_order(tmp_path):
    """The property the ordered residual hierarchy does not have.

    Reordering the gates moved the cohort Immune fraction from 20.50% to 3.31%
    under the ordered residual. Here it must change nothing at all.
    """

    policy = _policy()
    rows = [
        {"CD3e": ("positive", 4.0), "E_cadherin": ("positive", 1.1)},
        {"CD31": ("positive", 2.0)},
        {"Vimentin": ("positive", 3.0), "CD31": ("positive", 2.5)},
        {},
        {"MPO": ("unavailable", float("nan"))},
    ]
    baseline = type_simultaneous_cells(_evidence([dict(r) for r in rows]))

    source = json.loads(policy.source_path.read_text())
    reversed_gates = {k: source["gates"][k] for k in reversed(list(source["gates"]))}
    assert list(reversed_gates) != list(source["gates"])
    source["gates"] = reversed_gates
    path = tmp_path / "reordered_gates.json"
    path.write_text(json.dumps(source))
    reordered = load_broad_gate_policy(path, registry=load_registry())

    assert reordered.fingerprint == policy.fingerprint
    other = type_simultaneous_cells(
        _evidence([dict(r) for r in rows]), policy=reordered
    )
    pd.testing.assert_frame_equal(baseline, other)


def test_a_gate_that_contradicts_the_registry_fails_the_contract(tmp_path):
    """SST is the real instance.

    ordered_residual_rules.json placed SST in the Epithelial broad gate while
    marker_registry.json declares uses=["subtype"] with no broad_compartment,
    and hierarchical_typing lists it in SUBTYPE_ONLY_MARKERS. The ordered
    policy validated only set equality, so nothing caught it.
    """

    source = json.loads(_policy().source_path.read_text())
    source["gates"]["Epithelial"] = [*source["gates"]["Epithelial"], "SST"]
    path = tmp_path / "sst_gates.json"
    path.write_text(json.dumps(source))

    with pytest.raises(ContractError, match="SST"):
        load_broad_gate_policy(path, registry=load_registry())


def test_a_marker_in_the_wrong_lineage_gate_fails_the_contract(tmp_path):
    source = json.loads(_policy().source_path.read_text())
    source["gates"]["Immune"] = [
        m for m in source["gates"]["Immune"] if m != "CD3e"
    ]
    source["gates"]["Epithelial"] = [*source["gates"]["Epithelial"], "CD3e"]
    path = tmp_path / "misplaced_gates.json"
    path.write_text(json.dumps(source))

    with pytest.raises(ContractError, match="broad_compartment='Immune'"):
        load_broad_gate_policy(path, registry=load_registry())


def test_every_broad_registry_marker_must_be_gated(tmp_path):
    source = json.loads(_policy().source_path.read_text())
    source["gates"]["Neural"] = ["B3TUBB"]
    source["gates"]["Mesenchymal"] = ["Vimentin"]  # drops SMA
    path = tmp_path / "incomplete_gates.json"
    path.write_text(json.dumps(source))

    with pytest.raises(ContractError, match="SMA"):
        load_broad_gate_policy(path, registry=load_registry())


def test_the_policy_fingerprint_is_content_addressed():
    policy = BroadGatePolicy.read_json()
    assert policy.fingerprint == BroadGatePolicy.read_json().fingerprint
    assert len(policy.fingerprint) == 64
