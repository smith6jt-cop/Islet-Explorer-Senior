"""The selected pLN protein reference remains exact and exposes unavailable calls."""

from pathlib import Path

import pytest
import yaml

from multiplex5000.tissues.lymph_node import RuleSourceError, load_rule_source

REPO = Path(__file__).resolve().parents[4]
MANIFEST = REPO / "core/panels/data/tissues/pancreatic_lymph_node/manifest.yml"
PANEL = REPO / "core/multiplex5000/profiles/platforms/phenocycler.yaml"


def test_full_selected_sheet_and_explicit_pln_adaptation():
    source = load_rule_source(MANIFEST)
    assert len(source.rules) == 69
    assert len(source.omitted) == 16
    assert {m for _, m in source.omitted} == {"INS", "GCG", "SST", "Ker8-18"}
    assert {"FRC", "HEV", "GC_B_cell"} <= set(source.by_phenotype)
    assert "INS" not in source.inherited_markers("T_cell")


def test_missing_markers_are_reported_for_the_actual_pln_panel():
    source = load_rule_source(MANIFEST)
    panel = yaml.safe_load(PANEL.read_text())["panel_names"]["reference"]["pancreatic_lymph_node"]
    gaps = source.missing_by_phenotype(panel)
    assert sum(bool(missing) for missing in gaps.values()) == 14
    assert gaps["T_cell"] == gaps["FRC"] == gaps["HEV"] == ()
    assert gaps["GC_B_cell"] == gaps["Plasma_cell"] == ("CD38",)
    assert gaps["Beta_cell"] == ("CD99", "GCG", "INS", "SST")
    # The sheet's spelling and the dedicated panel's spelling are an explicit alias.
    assert gaps["PD-L1_expressing"] == ()
    assert source.missing_by_phenotype(["PD-L1"])["PD-L1_expressing"] == ()


def test_a_changed_source_cannot_silently_change_the_rules(tmp_path):
    changed = tmp_path / "changed.csv"
    changed.write_bytes((MANIFEST.parent / "phenotype_rules.csv").read_bytes() + b"\n")
    with pytest.raises(RuleSourceError, match="digest differs"):
        load_rule_source(MANIFEST, table_path=changed)
