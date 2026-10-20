"""Ordinary endocrine uncertainty stays admitted; exact exclusions stay excluded."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))
from protein_analysis_admission import load_protein_exclusions, protein_analysis_mask
import build_endocrine_object_composition as composition


def script(name):
    path = REPO / "analysis/xenium/scripts" / (name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_only_explicit_single_status_and_keyed_exclusions_remove_uncertain_cells():
    cells = pd.DataFrame({
        "object_id": ["ordinary-zero", "ordinary-tie", "excluded-key", "unsupported-single", "valid-single", "ineligible"],
        "cell_class": ["Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)", "Alpha",
                       "Endocrine (ambiguous subtype)", "Beta", "Immune"],
        "endocrine_object_class": ["Islet", "EndocrineCluster", "Islet", "EndocrineSingle", "EndocrineSingle", ""],
        "qc_analysis_eligible": [True] * 6,
        "analysis_eligible": [True] * 5 + [False]})
    exclusions = pd.DataFrame({"donor_id": ["112", "115"], "object_id": ["excluded-key", "ordinary-zero"]})
    expected = ["ordinary-zero", "ordinary-tie", "valid-single"]
    assert cells.loc[protein_analysis_mask(cells, "0112", exclusions), "object_id"].tolist() == expected
    assert cells.loc[composition.composition_cell_analysis_mask(cells, "0112", exclusions), "object_id"].tolist() == expected


def test_exclusion_registry_without_a_matching_cell_key_is_not_silently_ignored():
    with pytest.raises(ValueError, match="no exact cell key"):
        protein_analysis_mask(pd.DataFrame({"cell_class": ["Alpha"], "object_id": ["a"]}),
                              "112", pd.DataFrame({"donor_id": ["112"], "wrong_key": ["a"]}))


def test_assignment_registry_selects_exact_typing_exclusions_not_all_singles_or_qc_failures(tmp_path):
    path = tmp_path / "donor_id=112/data.parquet"
    path.parent.mkdir()
    pd.DataFrame({"donor_id": ["112"] * 3, "object_id": ["excluded", "valid-single", "failed-QC"],
                  "analysis_eligible": [False, True, False], "qc_analysis_eligible": [True, True, False],
                  "analysis_exclusion_reason": ["single_reason", "", "geometry_reason"]}).to_parquet(path, index=False)
    assert load_protein_exclusions(tmp_path, "0112").object_id.tolist() == ["excluded"]
    with pytest.raises(ValueError, match="expected one"):
        load_protein_exclusions(tmp_path, "115")


def test_uncertain_endocrine_contributes_to_broad_denominator_without_inventing_subtype_fraction():
    cells = pd.DataFrame({"object_id": ["alpha", "uncertain", "immune", "removed"],
                          "endocrine_object_id": ["object"] * 4,
                          "cell_class": ["Alpha", "Endocrine (unresolved subtype)", "Immune", "Immune"],
                          "refined_class": ["Alpha", "Endocrine (unresolved subtype)", "Immune", "Immune"],
                          "is_endocrine": [True, True, False, False],
                          "qc_analysis_eligible": [True] * 4,
                          "INS": [0.] * 4, "GCG": [.8, 0., 0., 0.], "SST": [0.] * 4,
                          "immune_class": [None, None, "Imm-T-cell", "Imm-T-cell"],
                          "immune_subtype": [None] * 4})
    row = {"donor_id": "112", "endocrine_object_id": "object", "n_cells_in_object": 4,
           "endocrine_denominator": 2, "identity_unresolved_count": 1,
           "stromal_identity_available_in_donor": False,
           "available__immune__all_immune_per_endocrine": True}
    for group, (_, suffix) in composition.OUTLINE_STRUCTURES.items():
        row[f"available__{group}__cell_bodies_per_endocrine"] = True
    for label in composition.IMMUNE_LABELS:
        row["count_immune__" + label] = 0
    for label in composition.IMMUNE_SUBTYPES:
        row["count_immune_subtype__" + label] = np.nan
    for hormone in ("INS", "GCG", "SST"):
        row["supplementary_available__" + hormone] = True
    exclusions = pd.DataFrame({"donor_id": ["112"], "object_id": ["removed"]})
    out = composition.rebuild_analytic_object_composition(pd.DataFrame([row]), cells, "112", exclusions).iloc[0]
    assert out.n_original_cells_in_object == 4
    assert out.n_analysis_admitted_cells_in_object == 3
    assert out.endocrine_denominator == 2
    assert out.count_endocrine__alpha == 1
    assert out.count_endocrine__unresolved_subtype == 1
    assert out.feature__endocrine__unresolved_fraction == .5
    assert out.feature__immune__all_immune_per_endocrine == .5
    assert not out.endocrine_subtype_composition_complete
    assert np.isnan(out.feature__endocrine__alpha_fraction)
    assert np.isnan(out.feature__endocrine__beta_fraction)
    assert not out.available__endocrine__beta_fraction


def test_aggregate_marker_profiles_retain_uncertainty_only_from_an_admitted_source():
    paper = script("build_complete_typing_paper_figures")
    profiles = pd.DataFrame({"marker": ["INS"] * 3,
                             "cell_class": ["Beta", "Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)"],
                             "mean_calibrated_marker_signal": [.8, 0., .5],
                             "denominator_scope": ["typing_analysis_eligible_cells"] * 3})
    assert set(paper.identity_marker_matrix(profiles).columns) == set(profiles.cell_class)
    with pytest.raises(ValueError, match="Regenerate marker profiles"):
        paper.identity_marker_matrix(profiles.drop(columns="denominator_scope"))


def test_native_display_filters_exact_keys_before_aggregation_and_keeps_missing_signal(tmp_path, monkeypatch):
    display = script("prepare_native_typing_display")
    source = tmp_path / "source/donor_id=0112/data.parquet"
    source.parent.mkdir(parents=True)
    pd.DataFrame({"cell_id": ["ordinary", "excluded-single", "valid-single"],
                  "qc_analysis_eligible": [True] * 3,
                  "cell_class": ["Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)", "Beta"],
                  "X_centroid": [1., 2., 3.], "Y_centroid": [4., 5., 6.],
                  "disease_status": ["ND"] * 3, "typing_scaled_INS": [np.nan, 1., .8]}).to_parquet(source, index=False)
    registry = tmp_path / "registry.parquet"
    pd.DataFrame({"donor_id": ["112"], "cell_id": ["excluded-single"]}).to_parquet(registry, index=False)
    out = tmp_path / "display"
    monkeypatch.setattr(sys, "argv", ["display", "--source", str(tmp_path / "source"),
                      "--out", str(out), "--exclusions", str(registry)])
    display.main()
    sample = pd.read_parquet(out / "pancreas_native_centroid_display_sample.parquet")
    assert sample.cell_id.tolist() == ["ordinary", "valid-single"]
    assert sample.analysis_eligible.all()
    profiles = pd.read_csv(out / "pancreas_native_identity_marker_profiles.csv")
    uncertain = profiles.loc[profiles.cell_class.eq("Endocrine (unresolved subtype)")].iloc[0]
    assert uncertain.n_class_cells == 1
    assert uncertain.n_marker_measured == 0
    assert np.isnan(uncertain.mean_calibrated_marker_signal)
    assert profiles.denominator_scope.eq("typing_analysis_eligible_cells").all()


@pytest.mark.parametrize("excluded_source_qc", [True, False])
def test_current_assignment_materialization_computes_census_and_preserves_all_rows(tmp_path, monkeypatch,
                                                                                excluded_source_qc):
    material = script("materialize_protein_analysis_assignments")
    source = tmp_path / "source/donor_id=0112/data.parquet"
    source.parent.mkdir(parents=True)
    pd.DataFrame({"cell_id": ["ordinary", "excluded", "valid-single", "source-ineligible"],
                  "object_id": ["ordinary", "excluded", "valid-single", "source-ineligible"],
                  "donor_id": ["0112"] * 4,
                  "cell_class": ["Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)", "Beta", "Immune"],
                  "refined_class": ["Endocrine (unresolved subtype)", "Endocrine (ambiguous subtype)", "Beta", "Immune"],
                  "is_endocrine": [True, True, True, False],
                  "qc_analysis_eligible": [True, excluded_source_qc, True, True],
                  "analysis_eligible": [True, True, True, False],
                  "analysis_reason": [""] * 4, "analysis_exclusion_reason": ["", "", "", "source restriction"],
                  "excluded_reason": [""] * 4, "resolved_by": ["cd99", "single_hormone", "single_hormone", "supporting"],
                  "endocrine_subtype_status": ["no_usable_hormone", "tied_highest_hormones", "unique_highest_hormone", ""]}).to_parquet(source, index=False)
    registry = tmp_path / "registry.parquet"
    pd.DataFrame({"donor_id": ["112"], "object_id": ["excluded"],
                  "analysis_exclusion_reason": ["endocrine_single_tied_highest_hormones"]}).to_parquet(registry, index=False)
    out = tmp_path / "fresh_assignments"
    monkeypatch.setattr(sys, "argv", ["materialize", "--source", str(tmp_path / "source"),
                      "--out", str(out), "--policy", "current_generation",
                      "--mandatory-exclusions", str(registry)])
    material.main()
    frame = pd.read_parquet(out / "donor_id=112/data.parquet")
    assert frame.object_id.tolist() == ["ordinary", "excluded", "valid-single", "source-ineligible"]
    assert frame.analysis_eligible.tolist() == [True, False, True, False]
    assert frame.cell_class.tolist() == ["Endocrine (unresolved subtype)", "Excluded", "Beta", "Immune"]
    assert frame.endocrine_subtype_status.tolist() == ["no_usable_hormone", "tied_highest_hormones", "unique_highest_hormone", ""]
    assert frame.qc_analysis_eligible.tolist() == [True, excluded_source_qc, True, True]
    assert frame.mandatory_analysis_exclusion.tolist() == [False, True, False, False]
    assert frame.loc[1, "analysis_assignment_state"] == "mandatory_typing_excluded"
    assert frame.loc[1, "analysis_exclusion_reason"] == "endocrine_single_tied_highest_hormones"
    import json
    census = json.loads((out / "manifest.json").read_text())["census"]
    assert census["mandatory_exclusion_keys"] == 1
    assert census["mandatory_exclusions_QC_eligible"] == int(excluded_source_qc)
    assert census["mandatory_exclusions_already_QC_ineligible"] == int(not excluded_source_qc)
    # Reusing a fresh assignment as the registry must preserve mandatory keys
    # even after preparation has already set their current source QC to false.
    assert load_protein_exclusions(out, "0112").object_id.tolist() == ["excluded"]
    assert not (tmp_path / "historical_package").exists()
    with pytest.raises(FileExistsError, match="fresh output"):
        material.main()


def test_historical_assignment_materialization_requires_package_before_writing(tmp_path, monkeypatch, capsys):
    material = script("materialize_protein_analysis_assignments")
    out = tmp_path / "assignments"
    monkeypatch.setattr(sys, "argv", ["materialize", "--out", str(out)])
    with pytest.raises(SystemExit) as error:
        material.main()
    assert error.value.code == 2
    assert "--package is required for --policy historical_20261005" in capsys.readouterr().err
    assert not out.exists()
