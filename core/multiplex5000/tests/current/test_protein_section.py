"""One section-local candidate engine serves the selected pLN and spleen rules."""

from pathlib import Path

import numpy as np
import pytest

from multiplex5000.tissues.lymph_node import Rule, RuleSource
from multiplex5000.tissues.lymph_node import load_rule_source as load_pln
from multiplex5000.tissues.protein_section import (ProteinSectionError,
                                                   evaluate_protein_section)
from multiplex5000.tissues.spleen import load_rule_source as load_spleen

REPO = Path(__file__).resolve().parents[4]
PANELS = REPO / "core/panels/data/tissues"


@pytest.mark.parametrize("loader,manifest", [
    (load_pln, PANELS / "pancreatic_lymph_node/manifest.yml"),
    (load_spleen, PANELS / "spleen/manifest.yml"),
])
def test_both_lymphoid_tissues_use_section_local_evidence(loader, manifest):
    source = loader(manifest)
    required = {marker for rule in source.rules for marker in rule.markers}
    measured = (required - {"CD38", "PD-L1"}) | {"PDL1"}
    values = np.r_[np.linspace(5, 6, 1199), 100.0]
    admitted = np.ones(len(values), dtype=bool)
    admitted[0] = False
    section = evaluate_protein_section(
        source, "section-1", {marker: values for marker in measured}, admitted,
        {marker: 1.0 for marker in measured})
    assert section.markers["PDL1"].cutoff is not None
    assert "PD-L1" not in section.missing_markers_by_rule["PD-L1_expressing"]
    assert section.markers["PDL1"].states[0] == -1
    assert all(states[0] == -1 for states in section.rules.states.values())
    if loader is load_pln:
        assert "CD38" in section.missing_markers_by_rule["Naive_B_cell"]
        assert np.all(section.rules.states["Naive_B_cell"] == -1)


def test_unscaled_marker_is_unavailable_and_nucleus_rule_applies_to_empty_rule():
    source = RuleSource("synthetic", (Rule("all", "Empty", ()),
                                    Rule("Empty", "CD20_positive", (("CD20", "pos"),))),
                        (), {})
    values = np.r_[np.linspace(5, 6, 1199), 100.0]
    admitted = np.ones(len(values), dtype=bool)
    admitted[0] = False
    section = evaluate_protein_section(source, "section-1", {"CD20": values}, admitted,
                                       {"CD20": None})
    assert np.all(section.markers["CD20"].states == -1)
    assert section.rules.states["Empty"].tolist() == [-1] + [1] * 1199
    assert np.all(section.rules.states["CD20_positive"] == -1)
    with pytest.raises(ProteinSectionError, match="absent"):
        evaluate_protein_section(source, "section-1", {"CD20": values}, admitted,
                                 {"CD38": 1.0})
