"""A same-section candidate join preserves exact IDs and refuses false pairings."""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiplex5000.integration.paired_typing import PairedTypingError, link_candidate_cells
from multiplex5000.integration.sameslide import match_by_centroid
from multiplex5000.integration.transform import Transform
from multiplex5000.same_slide_sources import SameSlidePair


def _fixture():
    pair = SameSlidePair(
        "spleen", "HDL090", "HDL090_PC67_0110618-S1",
        Path("/source/HDL090_PC67_0110618-S1.csv.gz"), Path("/source/0110618.qptiff"),
        "xenium:0110618:090",
        Path("/source/output-XETG00298__0110618__090__20260826__224731"),
        "same_section_after_phenocycler", "crop_local_um", "region_local_um")
    protein = pd.DataFrame({
        "cell_uid": [f"phenocycler:0110618:{pair.protein_entry}:digest:{i}"
                     for i in (1, 2, 3)], "object_id": ["1", "2", "3"],
        "donor_id": ["HDL090"] * 3, "acquisition_id": ["phenocycler:0110618"] * 3,
        "tissue": ["spleen"] * 3,
        "x_um": [100.0, 110.0, 120.0], "y_um": [0.0] * 3,
        "protein_eligible": [True, True, False]})
    protein_calls = pd.DataFrame({
        "cell_uid": protein["cell_uid"], "object_id": ["1", "2", "3"],
        "candidate_phenotype": ["T cell", None, None],
        "candidate_state": ["unique_candidate", "unavailable", "unavailable"]})
    rna = pd.DataFrame({
        "cell_uid": [f"{pair.xenium_acquisition_id}:cells:{i}" for i in (1, 2, 3)],
        "donor_id": ["HDL090"] * 3, "tissue": ["spleen"] * 3,
        "acquisition_id": ["xenium:0110618:090"] * 3,
        "x_um": [0.5, 10.5, 20.5], "y_um": [0.0] * 3,
        "typing_eligible": [True, True, False]})
    rna_calls = pd.DataFrame({
        "cell_uid": rna["cell_uid"],
        "typing_eligible": [True, True, False], "gene:CD3E": [1, -1, -1]})
    transform = Transform(
        affine=np.array([[1, 0, 100], [0, 1, 0]]),
        meta={"protein_entry": pair.protein_entry,
              "xenium_acquisition_id": pair.xenium_acquisition_id,
              "moving_frame": "region_local_um", "fixed_frame": "crop_local_um",
              "reviewed": True})
    return pair, protein, protein_calls, rna, rna_calls, transform


def test_candidate_links_use_distinct_native_ids_and_keep_missing_gene_state():
    pair, protein, protein_calls, rna, rna_calls, transform = _fixture()
    result = link_candidate_cells(pair, protein, protein_calls, rna, rna_calls, transform)
    assert result.links[["protein_cell_uid", "rna_cell_uid"]].values.tolist() == [
        [protein.loc[0, "cell_uid"], rna.loc[0, "cell_uid"]],
        [protein.loc[1, "cell_uid"], rna.loc[1, "cell_uid"]]]
    assert result.links["distance_um"].tolist() == pytest.approx([0.5, 0.5])
    assert result.links["gene:CD3E"].tolist() == [1, -1]
    assert result.links["protein_candidate_state"].tolist() == [
        "unique_candidate", "unavailable"]
    assert (result.protein_cells, result.rna_cells) == (3, 3)
    assert (result.protein_eligible, result.rna_eligible) == (2, 2)
    assert (result.unmatched_protein_eligible, result.unmatched_rna_eligible) == (0, 0)


@pytest.mark.parametrize("change", ["wrong_acquisition", "wrong_tissue", "wrong_section",
                                         "unreviewed", "serial", "lost_call", "wrong_object",
                                         "missing_eligibility"])
def test_candidate_links_refuse_unverified_or_inconsistent_sources(change):
    pair, protein, protein_calls, rna, rna_calls, transform = _fixture()
    if change == "wrong_acquisition":
        rna.loc[0, "acquisition_id"] = "xenium:0110618:082"
    elif change == "wrong_tissue":
        protein.loc[0, "tissue"] = "pancreas"
    elif change == "wrong_section":
        protein.loc[0, "cell_uid"] = protein.loc[0, "cell_uid"].replace("-S1", "-S2")
        protein_calls.loc[0, "cell_uid"] = protein.loc[0, "cell_uid"]
    elif change == "unreviewed":
        transform.meta["reviewed"] = False
    elif change == "serial":
        pair = replace(pair, relationship="serial")
    elif change == "lost_call":
        rna_calls = rna_calls.iloc[:-1]
    elif change == "wrong_object":
        protein_calls.loc[0, "object_id"] = "other"
    elif change == "missing_eligibility":
        rna["typing_eligible"] = pd.Series([None, True, False], dtype="boolean")
    with pytest.raises(PairedTypingError):
        link_candidate_cells(pair, protein, protein_calls, rna, rna_calls, transform)


def test_candidate_links_report_unmatched_cells_without_forcing_a_pair():
    pair, protein, protein_calls, rna, rna_calls, transform = _fixture()
    rna.loc[1, "x_um"] = 50.0
    result = link_candidate_cells(pair, protein, protein_calls, rna, rna_calls, transform)
    assert result.links[["protein_cell_uid", "rna_cell_uid"]].values.tolist() == [
        [protein.loc[0, "cell_uid"], rna.loc[0, "cell_uid"]]]
    assert (result.unmatched_protein_eligible, result.unmatched_rna_eligible) == (1, 1)


def test_nonrigid_registration_cannot_silently_load_as_affine_only(tmp_path):
    (tmp_path / "transform.json").write_text(json.dumps({
        "affine": [[1, 0, 0], [0, 1, 0]], "grid": None,
        "has_displacement": True, "meta": {}}))
    with pytest.raises(FileNotFoundError, match="declares a displacement"):
        Transform.load(tmp_path)


def _bound(value):
    return {"bound_um": value, "frame": "crop_local_um", "kind": "radial_displacement_bound",
            "evidence_ref": "synthetic-registration-review"}


@pytest.mark.parametrize("side", ["RNA", "protein"])
def test_equal_centroids_are_ambiguous_without_arbitrary_links(side):
    pair, p, pc, x, xc, tf = _fixture()
    if side == "RNA":
        x.loc[1, "x_um"] = x.loc[0, "x_um"]
    else:
        p.loc[1, "x_um"] = p.loc[0, "x_um"]
    result = link_candidate_cells(pair, p, pc, x, xc, tf)
    assert result.links.empty
    assert len(result.candidates[result.candidates.nominal_state.eq("ambiguous_nearest")]) == 2
    assert result.cell_dispositions.link_state.eq("ambiguous_nearest").sum() == 3
    assert result.unmatched_protein_eligible == result.unmatched_rna_eligible == 2
    assert set(result.candidates.segmentation_relation_state) == {"unavailable_without_overlap_evidence"}


def test_tied_candidate_outputs_are_invariant_to_independent_row_permutations():
    pair, p, pc, x, xc, tf = _fixture()
    p.loc[1, "x_um"] = 101.0  # RNA 1 is equidistant from both protein centroids.
    x.loc[1, "x_um"] = 50.0
    first = link_candidate_cells(pair, p, pc, x, xc, tf)
    second = link_candidate_cells(pair, p.iloc[[2, 0, 1]], pc.iloc[[1, 2, 0]],
                                  x.iloc[[1, 0, 2]], xc.iloc[[2, 1, 0]], tf)
    pd.testing.assert_frame_equal(first.candidates, second.candidates)
    pd.testing.assert_frame_equal(first.cell_dispositions, second.cell_dispositions)
    assert first.links.empty


def test_legacy_centroid_matcher_also_refuses_ties():
    moving = np.array([[0.0, 0.0]])
    fixed = np.array([[-1.0, 0.0], [1.0, 0.0]])
    assert match_by_centroid(moving, fixed).empty
    assert match_by_centroid(fixed, moving).empty


def test_missing_uncertainty_is_unavailable_not_zero():
    result = link_candidate_cells(*_fixture())
    assert result.metadata["registration_uncertainty"]["status"] == "unavailable"
    assert result.links.registration_uncertainty_um.isna().all()
    assert result.links.stable_under_provided_bound.isna().all()
    assert result.links.rna_identity_call_present.eq(False).all()


def test_provided_uncertainty_reveals_distance_boundary_sensitivity():
    pair, p, pc, x, xc, tf = _fixture()
    x.loc[0, "x_um"] = 4.9
    result = link_candidate_cells(pair, p, pc, x, xc, tf, registration_uncertainty=_bound(0.2))
    assert len(result.links) == 2  # The method's nominal 5-um cap is unchanged.
    row = result.links.loc[result.links.rna_cell_uid == x.loc[0, "cell_uid"]].iloc[0]
    assert row.distance_lower_um == pytest.approx(4.7)
    assert row.distance_upper_um == pytest.approx(5.1)
    assert not row.stable_under_provided_bound
    assert result.metadata["max_dist_um"] == 5.0


def test_provided_uncertainty_reveals_competing_nearest_candidate():
    pair, p, pc, x, xc, tf = _fixture()
    p.loc[1, "x_um"] = 101.1
    x.loc[1, "x_um"] = 50.0
    result = link_candidate_cells(pair, p, pc, x, xc, tf, registration_uncertainty=_bound(0.1))
    assert len(result.links) == 1
    assert result.links.iloc[0].distance_upper_um < 5
    assert not result.links.iloc[0].stable_under_provided_bound
    assert result.links.iloc[0].rna_candidate_count == 2


def test_uncertainty_candidates_outside_nominal_cap_do_not_become_links():
    pair, p, pc, x, xc, tf = _fixture()
    p.loc[1, "x_um"] = 200.0
    x.loc[0, "x_um"], x.loc[1, "x_um"] = 5.1, 200.0
    result = link_candidate_cells(pair, p, pc, x, xc, tf, registration_uncertainty=_bound(0.2))
    assert result.links.empty
    assert len(result.candidates) == 1
    assert result.candidates.iloc[0].nominal_state == "outside_nominal_cap"
    assert result.cell_dispositions.link_state.eq("only_uncertainty_candidates").sum() == 2


def test_provided_zero_and_small_uncertainty_can_preserve_unique_candidates():
    for bound in [0.0, 0.1]:
        result = link_candidate_cells(*_fixture(), registration_uncertainty=_bound(bound))
        assert result.links.stable_under_provided_bound.all()
        assert result.metadata["registration_uncertainty"]["status"] == "provided"


@pytest.mark.parametrize("field,value", [("bound_um", -1), ("bound_um", float("nan")),
                                        ("frame", "region_local_um"), ("evidence_ref", "")])
def test_uncertainty_requires_explicit_compatible_evidence(field, value):
    bound = _bound(0.2)
    bound[field] = value
    with pytest.raises(PairedTypingError, match="uncertainty"):
        link_candidate_cells(*_fixture(), registration_uncertainty=bound)


def test_actual_rna_calls_disagreements_and_unavailable_states_are_preserved():
    pair, p, pc, x, xc, tf = _fixture()
    xc = xc.drop(columns="gene:CD3E")
    xc["cell_class"] = ["B cell", "B cell", None]
    xc["call_state"] = ["resolved", "unavailable", "unavailable"]
    pc.loc[1, "candidate_phenotype"] = "T cell"
    mapping = {"tissue": "spleen", "level": "lineage", "protein_column": "candidate_phenotype",
        "rna_column": "cell_class", "protein": {"T cell": "T"}, "rna": {"B cell": "B"},
        "protein_state_column": "candidate_state", "protein_available_states": ["unique_candidate"],
        "rna_state_column": "call_state", "rna_available_states": ["resolved"]}
    result = link_candidate_cells(pair, p, pc, x, xc, tf, hierarchy_map=mapping)
    assert result.links.protein_candidate_phenotype.tolist() == ["T cell", "T cell"]
    assert result.links.rna_cell_class.tolist() == ["B cell", "B cell"]
    assert result.links.rna_call_state.tolist() == ["resolved", "unavailable"]
    assert result.links.identity_concordance.tolist() == ["disagree", "unavailable"]
    assert "cell_class" not in result.links  # No fused identity.


def test_comparison_requires_declared_tissue_hierarchy_mapping():
    pair, p, pc, x, xc, tf = _fixture()
    xc["cell_class"] = ["T cell", None, None]
    result = link_candidate_cells(pair, p, pc, x, xc, tf)
    assert result.links.identity_concordance.eq("unavailable").all()
    with pytest.raises(PairedTypingError, match="exact tissue"):
        link_candidate_cells(pair, p, pc, x, xc, tf, hierarchy_map={"tissue": "pancreas", "level": "lineage"})


@pytest.mark.parametrize("change", ["frame", "units", "segmentation", "call_frame", "call_acquisition"])
def test_mixed_cell_or_call_source_provenance_is_rejected(change):
    pair, p, pc, x, xc, tf = _fixture()
    if change == "frame":
        p["coordinate_frame"] = ["crop_local_um", "full_slide_um", "crop_local_um"]
    elif change == "units":
        x["coordinate_units"] = ["um", "pixels", "um"]
    elif change == "segmentation":
        x["segmentation_version"] = ["v1", "v2", "v1"]
    elif change == "call_frame":
        x["coordinate_frame"] = "region_local_um"
        xc["coordinate_frame"] = "full_slide_um"
    else:
        xc["acquisition_id"] = "wrong"
    with pytest.raises(PairedTypingError):
        link_candidate_cells(pair, p, pc, x, xc, tf)


def test_zero_eligible_cells_reconcile_without_fabricated_candidates():
    pair, p, pc, x, xc, tf = _fixture()
    x["typing_eligible"] = xc["typing_eligible"] = False
    result = link_candidate_cells(pair, p, pc, x, xc, tf)
    assert result.links.empty and result.candidates.empty
    assert result.rna_eligible == 0 and result.unmatched_protein_eligible == 2
    excluded = result.cell_dispositions.loc[result.cell_dispositions.assay == "rna"]
    assert excluded.candidate_count.isna().all()


def test_call_only_frame_metadata_is_still_validated():
    pair, p, pc, x, xc, tf = _fixture()
    xc["coordinate_frame"] = ["region_local_um", "full_slide_um", "region_local_um"]
    with pytest.raises(PairedTypingError, match="coordinate frame"):
        link_candidate_cells(pair, p, pc, x, xc, tf)


def test_independent_call_columns_cannot_be_silently_overwritten():
    pair, p, pc, x, xc, tf = _fixture()
    xc["candidate_count"] = [42, 42, 42]
    with pytest.raises(PairedTypingError, match="collide with link metadata"):
        link_candidate_cells(pair, p, pc, x, xc, tf)


def test_duplicate_cell_ids_and_cell_columns_are_rejected():
    pair, p, pc, x, xc, tf = _fixture()
    with pytest.raises(PairedTypingError, match="duplicate cell columns"):
        link_candidate_cells(pair, pd.concat([p, p[["x_um"]]], axis=1), pc, x, xc, tf)
    p.loc[1, "cell_uid"] = p.loc[0, "cell_uid"]
    with pytest.raises(PairedTypingError, match="unique"):
        link_candidate_cells(pair, p, pc, x, xc, tf)


def test_declared_frames_and_missing_segmentation_provenance_are_distinguished():
    pair, p, pc, x, xc, tf = _fixture()
    p["coordinate_frame"] = "crop_local_um"
    p["segmentation_id"] = "reviewed-protein-segmentation"
    p["segmentation_version"] = 1
    result = link_candidate_cells(pair, p, pc, x, xc, tf)
    assert result.metadata["coordinate_frame_evidence"] == {
        "protein": "cell_metadata_and_pair", "rna": "declared_pair_only"}
    assert result.metadata["segmentation_provenance"]["rna"] == {
        "segmentation_id": None, "segmentation_version": None}
    assert result.metadata["segmentation_provenance"]["protein"]["segmentation_version"] == 1
    json.dumps(result.metadata, allow_nan=False)


def test_unavailable_typed_eligibility_cannot_pass_nullable_equality_check():
    pair, p, pc, x, xc, tf = _fixture()
    xc["typing_eligible"] = pd.Series([pd.NA, True, False], dtype="boolean")
    with pytest.raises(PairedTypingError, match="typing_eligible changed"):
        link_candidate_cells(pair, p, pc, x, xc, tf)
