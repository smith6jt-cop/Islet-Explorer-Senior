"""The spleen pilot reuses the selected hierarchy and declares panel gaps."""

from pathlib import Path

import pytest
import yaml

from multiplex5000.integration import tissues
from multiplex5000.tissues.lymph_node import RuleSourceError
from multiplex5000.tissues.spleen import load_rule_source, observed_markers

REPO = Path(__file__).resolve().parents[4]
MANIFEST = REPO / "core/panels/data/tissues/spleen/manifest.yml"


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
