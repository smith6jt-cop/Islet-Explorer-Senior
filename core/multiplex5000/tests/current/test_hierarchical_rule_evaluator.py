"""One rule evaluator serves the selected pLN source and its spleen adaptation."""

from pathlib import Path

import pytest

from multiplex5000.tissues.evaluate import evaluate_rules
from multiplex5000.tissues.lymph_node import RuleSourceError
from multiplex5000.tissues.lymph_node import load_rule_source as load_pln
from multiplex5000.tissues.spleen import load_rule_source as load_spleen

REPO = Path(__file__).resolve().parents[4]
PLN = REPO / "core/panels/data/tissues/pancreatic_lymph_node/manifest.yml"
SPLEEN = REPO / "core/panels/data/tissues/spleen/manifest.yml"


def _negative_panel(source):
    return {marker: False for rule in source.rules for marker in rule.markers}


def test_parent_and_any_positive_rule_are_required_for_a_t_cell():
    source = load_pln(PLN)
    calls = _negative_panel(source)
    assert evaluate_rules(source, calls).states["T_cell"] == "not_matched"
    calls["CD3e"] = True
    result = evaluate_rules(source, calls)
    assert result.states["Immune"] == "matched"
    assert result.states["T_cell"] == "matched"
    assert result.states["B_cell"] == "not_matched"


def test_missing_inherited_marker_is_unavailable_not_a_negative():
    source = load_pln(PLN)
    calls = _negative_panel(source)
    calls["CD20"] = True
    calls["CD38"] = None
    result = evaluate_rules(source, calls)
    assert result.states["B_cell"] == "matched"
    assert result.states["GC_B_cell"] == "unavailable"
    assert result.missing_markers["GC_B_cell"] == ("CD38",)
    assert result.missing_markers["T_cell"] == ()


def test_spleen_adaptation_uses_the_same_evaluator():
    source = load_spleen(SPLEEN)
    calls = _negative_panel(source)
    calls.update({"CD20": True, "CD21": True})
    result = evaluate_rules(source, calls)
    assert result.states["B_cell"] == "matched"
    assert result.states["CD21_positive_B_cell"] == "matched"
    assert "Beta_cell" not in result.states


def test_alias_is_a_single_observation_and_conflicts_are_refused():
    source = load_pln(PLN)
    calls = _negative_panel(source)
    calls.pop("PD-L1")
    calls["PDL1"] = True
    assert evaluate_rules(source, calls).states["PD-L1_expressing"] == "matched"
    calls["PD-L1"] = False
    with pytest.raises(RuleSourceError, match="disagree"):
        evaluate_rules(source, calls)


def test_non_boolean_intensity_cannot_be_mistaken_for_a_marker_call():
    source = load_pln(PLN)
    with pytest.raises(RuleSourceError, match="true, false or unavailable"):
        evaluate_rules(source, {"CD3e": 0.5})
