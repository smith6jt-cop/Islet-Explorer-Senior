"""A same-section candidate join preserves exact IDs and refuses false pairings."""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiplex5000.integration.paired_typing import PairedTypingError, link_candidate_cells
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
        "cell_uid": ["pc:1", "pc:2", "pc:3"], "object_id": ["1", "2", "3"],
        "donor_id": ["HDL090"] * 3, "acquisition_id": ["phenocycler:0110618"] * 3,
        "x_um": [100.0, 110.0, 120.0], "y_um": [0.0] * 3,
        "protein_eligible": [True, True, False]})
    protein_calls = pd.DataFrame({
        "cell_uid": ["pc:1", "pc:2", "pc:3"], "object_id": ["1", "2", "3"],
        "candidate_phenotype": ["T cell", None, None],
        "candidate_state": ["unique_candidate", "unavailable", "unavailable"]})
    rna = pd.DataFrame({
        "cell_uid": ["rna:1", "rna:2", "rna:3"], "donor_id": ["HDL090"] * 3,
        "acquisition_id": ["xenium:0110618:090"] * 3,
        "x_um": [0.5, 10.5, 20.5], "y_um": [0.0] * 3,
        "typing_eligible": [True, True, False]})
    rna_calls = pd.DataFrame({
        "cell_uid": ["rna:1", "rna:2", "rna:3"],
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
        ["pc:1", "rna:1"], ["pc:2", "rna:2"]]
    assert result.links["distance_um"].tolist() == pytest.approx([0.5, 0.5])
    assert result.links["gene:CD3E"].tolist() == [1, -1]
    assert result.links["protein_candidate_state"].tolist() == [
        "unique_candidate", "unavailable"]
    assert (result.protein_cells, result.rna_cells) == (3, 3)
    assert (result.protein_eligible, result.rna_eligible) == (2, 2)
    assert (result.unmatched_protein_eligible, result.unmatched_rna_eligible) == (0, 0)


@pytest.mark.parametrize("change", ["wrong_acquisition", "unreviewed", "serial",
                                         "lost_call", "wrong_object", "missing_eligibility"])
def test_candidate_links_refuse_unverified_or_inconsistent_sources(change):
    pair, protein, protein_calls, rna, rna_calls, transform = _fixture()
    if change == "wrong_acquisition":
        rna.loc[0, "acquisition_id"] = "xenium:0110618:082"
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
        ["pc:1", "rna:1"]]
    assert (result.unmatched_protein_eligible, result.unmatched_rna_eligible) == (1, 1)


def test_nonrigid_registration_cannot_silently_load_as_affine_only(tmp_path):
    (tmp_path / "transform.json").write_text(json.dumps({
        "affine": [[1, 0, 0], [0, 1, 0]], "grid": None,
        "has_displacement": True, "meta": {}}))
    with pytest.raises(FileNotFoundError, match="declares a displacement"):
        Transform.load(tmp_path)
