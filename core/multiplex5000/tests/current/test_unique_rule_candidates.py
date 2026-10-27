"""The hierarchy can yield a review candidate without forcing overlapping leaves."""

from __future__ import annotations

import numpy as np
import pytest

from multiplex5000.tissues.evaluate import (
    BatchRuleEvaluation,
    evaluate_rules_batch,
    resolve_unique_matched_rules,
)
from multiplex5000.tissues.lymph_node import Rule, RuleSource, RuleSourceError


def _source() -> RuleSource:
    return RuleSource("synthetic", (
        Rule("all", "Immune", (("CD45", "pos"),)),
        Rule("Immune", "T_cell", (("CD3", "pos"),)),
        Rule("T_cell", "CD8_T", (("CD8", "pos"),)),
        Rule("Immune", "B_cell", (("CD20", "pos"),)),
        Rule("all", "Stromal", (("CD31", "pos"),)),
        Rule("any", "Proliferating", (("Ki67", "pos"),)),
    ), (), {})


def test_unique_deepest_rule_and_ambiguous_branches():
    source = _source()
    calls = {
        "CD45": np.array([1, 1, 0, -1, 1, 1, 1, 0], dtype=np.int8),
        "CD3": np.array([1, 1, 0, -1, 1, 1, 1, 0], dtype=np.int8),
        "CD8": np.array([1, 0, 0, -1, -1, 1, 1, 0], dtype=np.int8),
        "CD20": np.array([0, 1, 0, -1, 0, 0, 0, 0], dtype=np.int8),
        "CD31": np.array([0, 0, 0, -1, 0, 1, 0, 0], dtype=np.int8),
        "Ki67": np.array([0, 0, 0, -1, 0, 0, 1, -1], dtype=np.int8),
    }
    resolved = resolve_unique_matched_rules(source, evaluate_rules_batch(source, calls))
    assert resolved.phenotype.tolist() == ["CD8_T", "", "", "", "T_cell", "", "CD8_T", ""]
    assert resolved.status.tolist() == ["unique_candidate", "ambiguous", "no_match",
                                        "unavailable", "unique_candidate", "ambiguous",
                                        "unique_candidate", "no_match"]
    assert resolved.maximal_matches.tolist() == [1, 2, 0, 0, 1, 2, 1, 0]
    assert resolved.unavailable_lineage_rules[4] == 1  # CD8 subtype unavailable; T confirmed.
    assert resolved.unavailable_lineage_rules[7] == 0  # State gap cannot erase lineage no-match.


def test_resolver_refuses_a_matched_child_without_its_parent():
    source = _source()
    bad = BatchRuleEvaluation({name: np.array([1 if name == "T_cell" else 0], dtype=np.int8)
                               for name in source.by_phenotype}, {})
    with pytest.raises(RuleSourceError, match="matched child has no matched parent"):
        resolve_unique_matched_rules(source, bad)


def test_tissue_candidate_roots_exclude_a_pancreas_rule_without_erasing_evidence():
    source = RuleSource("synthetic", (
        Rule("all", "Immune", (("CD45", "pos"),)),
        Rule("all", "Endocrine", (("INS", "pos"),)),
    ), (), {}, candidate_lineage_roots=("Immune",))
    evidence = evaluate_rules_batch(source, {
        "CD45": np.array([0, 1], dtype=np.int8),
        "INS": np.array([-1, 1], dtype=np.int8),
    })
    assert evidence.states["Endocrine"].tolist() == [-1, 1]
    resolved = resolve_unique_matched_rules(source, evidence)
    assert resolved.status.tolist() == ["no_match", "unique_candidate"]
    assert resolved.phenotype.tolist() == ["", "Immune"]
    assert resolved.unavailable_lineage_rules.tolist() == [0, 0]
