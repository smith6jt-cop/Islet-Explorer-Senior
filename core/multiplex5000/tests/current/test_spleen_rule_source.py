"""The spleen pilot reuses the selected hierarchy and declares panel gaps."""

import csv
import hashlib
from pathlib import Path

import pytest
import yaml

from multiplex5000.integration import tissues
from multiplex5000.same_slide_sources import load_pairs
from multiplex5000.tissues.lymph_node import RuleSourceError
from multiplex5000.tissues.spleen import load_rule_source, normalized_markers, observed_markers

REPO = Path(__file__).resolve().parents[4]
MANIFEST = REPO / "core/panels/data/tissues/spleen/manifest.yml"
PAIRS = REPO / "core/multiplex5000/sites/hipergator/spleen_same_slide_20261026.yaml"


def test_spleen_profile_is_discoverable_without_another_tissue_list():
    assert tissues.normalize("splenic") == "spleen"
    assert tissues.for_roi("spl_1") == "spleen"
    assert tissues.rois("spleen") == ("spl",)
    assert tissues.panel_tissue_dir("spleen") is None
    assert not tissues.has_structures("spleen")


def test_pc67_panel_is_observed_on_all_eight_exports():
    markers = observed_markers(MANIFEST)
    assert len(markers) == 67
    assert {"CD21", "CD45RO", "CD141", "CD20", "CD3e", "CD38"} <= markers
    manifest = yaml.safe_load(MANIFEST.read_text())
    assert manifest["protein_panel"]["checked_exports"] == 8
    assert len(manifest["protein_panel"]["source_header_sha256"]) == 64


def test_spleen_rules_drop_pancreas_branches_and_reuse_one_source():
    source = load_rule_source(MANIFEST)
    classes = source.by_phenotype
    assert source.sha256 == yaml.safe_load(MANIFEST.read_text())["rule_table"]["source_sha256"]
    assert not {"Endocrine", "Beta_cell", "Ductal_Epithelial", "Epithelial",
                "Myofibroblast_CAF", "HEV", "Basement_membrane"} & classes.keys()
    assert {"FRC", "GC_B_cell", "T_cell", "B_cell", "CD21_positive_B_cell",
            "CD45RO_positive_T_cell", "CD141_positive_cDC", "SMA_positive_stromal"} <= classes.keys()
    assert "INS" not in source.inherited_markers("T_cell")
    assert "GCG" not in source.inherited_markers("B_cell")


def test_missing_markers_are_unavailable_not_negative():
    source = load_rule_source(MANIFEST)
    missing = source.missing_by_phenotype(observed_markers(MANIFEST))
    assert missing["T_cell"] == missing["CD21_positive_B_cell"] == ()
    assert missing["Schwann_cell"] == ("B3TUBB",)
    assert "Iba1" in missing["Macrophage"]


def test_scale_table_is_declared_as_partial_coverage_not_full_panel():
    manifest = yaml.safe_load(MANIFEST.read_text())
    scales = manifest["protein_normalization"]
    assert scales["pc67_entries"] == 8
    assert scales["markers_per_pc67_entry"] == 33
    assert len(observed_markers(MANIFEST)) == 67
    assert scales["role"] == "existing_classifier_scale_comparison"


def test_scale_audit_binds_all_eight_paired_sections_and_refuses_a_gap(tmp_path):
    entries = {pair.protein_entry for pair in load_pairs(PAIRS)}
    assert len(entries) == 8
    table = tmp_path / "scales.csv"
    columns = ["image", "group", "marker", "reference", "scale", "shift_log",
               "provenance", "reagent_changed"]

    def write(gap=False):
        with table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for image in sorted(entries):
                for marker in ("CD20", "CD21"):
                    if gap and image == sorted(entries)[0] and marker == "CD21":
                        continue
                    writer.writerow({"image": image, "group": "PC67", "marker": marker,
                                     "reference": image, "scale": 1.0, "shift_log": 0.0,
                                     "provenance": "automatic_registration",
                                     "reagent_changed": "False"})

    manifest = yaml.safe_load(MANIFEST.read_text())
    manifest["protein_panel"]["file"] = str(MANIFEST.parent / "pc67_export_header.tsv")
    manifest["protein_normalization"].update({"rows": 16, "markers_per_pc67_entry": 2})
    fixture = tmp_path / "manifest.yml"
    write()
    manifest["protein_normalization"]["sha256"] = hashlib.sha256(table.read_bytes()).hexdigest()
    fixture.write_text(yaml.safe_dump(manifest))
    assert normalized_markers(fixture, table_path=table, expected_entries=entries) == {
        "CD20", "CD21"}
    with pytest.raises(RuleSourceError, match="Stage 00 pairs"):
        normalized_markers(fixture, table_path=table, expected_entries=entries - {min(entries)})
    write(gap=True)
    manifest["protein_normalization"]["rows"] = 15
    manifest["protein_normalization"]["sha256"] = hashlib.sha256(table.read_bytes()).hexdigest()
    fixture.write_text(yaml.safe_dump(manifest))
    with pytest.raises(RuleSourceError, match="same normalized marker set"):
        normalized_markers(fixture, table_path=table, expected_entries=entries)


def test_unselected_source_digest_is_refused(tmp_path):
    manifest = yaml.safe_load(MANIFEST.read_text())
    manifest["rule_table"]["source_sha256"] = "0" * 64
    manifest["rule_table"]["source_manifest"] = str(
        MANIFEST.parent.parent / "pancreatic_lymph_node/manifest.yml")
    wrong = tmp_path / "manifest.yml"
    wrong.write_text(yaml.safe_dump(manifest))
    with pytest.raises(RuleSourceError):
        load_rule_source(wrong)


def test_changed_pc67_header_is_refused(tmp_path):
    manifest = yaml.safe_load(MANIFEST.read_text())
    changed_header = tmp_path / "changed.tsv"
    changed_header.write_bytes((MANIFEST.parent / "pc67_export_header.tsv").read_bytes() + b"\n")
    manifest["protein_panel"]["file"] = changed_header.name
    changed = tmp_path / "manifest.yml"
    changed.write_text(yaml.safe_dump(manifest))
    with pytest.raises(RuleSourceError, match="header differs"):
        observed_markers(changed)
