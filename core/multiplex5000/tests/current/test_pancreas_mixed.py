"""A mixed IO60 slide keeps pancreas and pLN typing evidence in their own tissue."""

import pandas as pd
import pytest

from multiplex5000.identity import check_cell_keys
from multiplex5000.tissues.pancreas_mixed import MixedTypingError, assemble


def _inputs():
    aligned = pd.DataFrame({
        "object_id": ["shared", "pancreas-unfit", "node", "node-unfit"],
        "qc_analysis_eligible": [True, False, True, False],
        "cell_region": ["tissue", "core", "pLN", "pLN"],
        "in_alignment_block": [True, True, True, True],
    })
    pancreas = pd.DataFrame({"object_id": ["shared"], "cell_class": ["Immune"],
                             "refined_class": ["T_cell"]})
    pln = pd.DataFrame({"object_id": ["node", "node-unfit"],
                        "candidate_phenotype": ["B_cell", None],
                        "candidate_state": ["unique_candidate", "unavailable"]})
    blocks = pd.DataFrame({"object_id": ["node", "node-unfit"],
                           "block": ["same-block", "same-block"]})
    return aligned, pancreas, pln, blocks


def test_mixed_slide_preserves_all_cells_without_cross_tissue_calls():
    aligned, pancreas, pln, blocks = _inputs()
    result = assemble(aligned, pancreas, pln, donor_id="6450", acquisition_id="scan:1",
                      segmentation_id="v3", pln_blocks=blocks).set_index("object_id")
    assert len(result) == len(aligned)
    check_cell_keys(result.reset_index())
    assert result.loc["shared", "cell_class"] == "Immune"
    assert result.loc["shared", "refined_class"] == "T_cell"
    assert pd.isna(result.loc["shared", "candidate_phenotype"])
    assert result.loc["node", "tissue"] == "pancreatic_lymph_node"
    assert pd.isna(result.loc["node", "cell_class"])
    assert result.loc["node", "candidate_phenotype"] == "B_cell"
    assert result.loc["node", "pln_block"] == "same-block"
    assert result.loc["node-unfit", "typing_state"] == "unavailable_qc"
    assert result.loc["pancreas-unfit", "typing_state"] == "unavailable_qc"


def test_mixed_slide_refuses_a_pancreas_call_on_pln_and_missing_node_cells():
    aligned, pancreas, pln, blocks = _inputs()
    wrong = pancreas.copy()
    wrong.loc[0, "object_id"] = "node"
    with pytest.raises(MixedTypingError, match="eligible non-pLN"):
        assemble(aligned, wrong, pln, donor_id="6450", acquisition_id="scan",
                 segmentation_id="v3", pln_blocks=blocks)
    with pytest.raises(MixedTypingError, match="every recorded pLN"):
        assemble(aligned, pancreas, pln.iloc[:1], donor_id="6450", acquisition_id="scan",
                 segmentation_id="v3", pln_blocks=blocks)


def test_mixed_slide_refuses_unknown_region_and_ineligible_candidate():
    aligned, pancreas, pln, blocks = _inputs()
    unknown = aligned.copy()
    unknown.loc[0, "cell_region"] = "unknown"
    with pytest.raises(MixedTypingError, match="unknown tissue route"):
        assemble(unknown, pancreas, pln, donor_id="6450", acquisition_id="scan",
                 segmentation_id="v3", pln_blocks=blocks)
    invalid = pln.copy()
    invalid.loc[1, "candidate_state"] = "unique_candidate"
    with pytest.raises(MixedTypingError, match="ineligible pLN"):
        assemble(aligned, pancreas, invalid, donor_id="6450", acquisition_id="scan",
                 segmentation_id="v3", pln_blocks=blocks)
