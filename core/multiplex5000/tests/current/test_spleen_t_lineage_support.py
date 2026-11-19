"""A marker-qualified label cannot hide a double-negative spleen T call."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from multiplex5000.positivity import hierarchy as H
from multiplex5000.positivity import registry as RG
from multiplex5000.stages import spleen_typing as stage
from multiplex5000.tissues import spleen
from multiplex5000.tissues.lymph_node import Rule, RuleSourceError

MANIFEST = Path(__file__).resolve().parents[4] / "core/panels/data/tissues/spleen/manifest.yml"


def _name(positive, unavailable=()):
    rules = spleen.load_rule_source(MANIFEST)
    reg = RG.load("phenocycler-pc67", "spleen")
    names = H.sheet_names(rules, reg.table.index)
    markers = list(names.values())
    calls = {m: np.array([-1 if m in unavailable else int(m in positive)], np.int8) for m in markers}
    policy = spleen.naming_policy(MANIFEST)
    evaluated = H.evaluate_profiles(rules, calls, spleen.structure_only_branches(MANIFEST),
                                    defer_exclusions=policy["defer_exclusions"])
    p = np.column_stack([calls[m] == 1 for m in markers])
    fold, cover = np.ones(p.shape), np.full(p.shape, .2)
    fold[0, markers.index("CD4")], cover[0, markers.index("CD4")] = 4, .8
    evidence = {"columns": markers, "typing": [names[m] for m in reg.typing_markers() if m in names],
                "P": p, "C": p.copy(), "F": fold, "COV": cover, "interior_valid": np.ones(1, bool)}
    out = H.name_cells(rules, evaluated, evidence, reg.const("CANARY_MARGIN"), reg.const("CANARY_MIN_DEPTH"),
                       reg.const("CANARY_TIE_INTERIOR_LEAD"), require_specific=policy["require_specific"],
                       inherit_depth=policy["inherit_depth"])
    return out.iloc[0], calls, evaluated


@pytest.mark.parametrize("unavailable", [(), ("CD4",)])
def test_cd45ro_cannot_name_a_t_cell_without_cd4_or_cd8_support(unavailable):
    row, calls, evaluated = _name({"CD3e", "CD45RO"}, unavailable)
    assert row.refined_class == ""
    assert calls["CD3e"][0] == calls["CD45RO"][0] == 1
    assert calls["CD8"][0] == 0
    assert evaluated.states["CD45RO_positive_T_cell"][0] == (-1 if unavailable else 0)


def test_treg_profile_with_cd45ro_keeps_its_specificity():
    row, _, _ = _name({"CD3e", "CD4", "FOXP3", "CD45RO"})
    assert row.refined_class == "Treg" and row.cell_class == "T_cell"
    rules = spleen.load_rule_source(MANIFEST).by_phenotype
    assert dict(rules["T_cell"].constraints)["CD4"] == "anypos"
    assert dict(rules["T_cell"].constraints)["CD8"] == "anypos"
    assert dict(rules["T_cell"].constraints)["CD3e"] == "pos"
    assert rules["CD45RO_positive_T_cell"].constraints == (("CD45RO", "pos"),)


def test_semantic_check_sees_descendants_and_unavailable_support():
    rules = spleen.load_rule_source(MANIFEST)
    refined = np.array(["CD45RO_positive_T_cell", "Activated_Treg", "Memory_CD8",
                        "CD45RO_positive_T_cell", "Endothelial_Sinusoidal_CD8+", ""], object)
    calls = {"CD3e": np.array([1, 0, 1, 1, 0, 1]),
             "CD4": np.array([0, 1, 0, -1, 0, 0]), "CD8": np.array([0, 0, 1, 0, 1, 0])}
    assert stage._unsupported_t_cells(rules, refined, calls).tolist() == [True, True, False, True, False, False]


@pytest.mark.parametrize("unsupported", [1, None])
def test_semantic_gate_blocks_bad_calls_or_a_missing_support_audit(tmp_path, unsupported):
    per = {"cells": 1, "named_below_level": 0, "forbidden_final_classes": 0,
           "unsupported_t_cells": unsupported, "refined_class": {"CD45RO_positive_T_cell": 1},
           "status": {}, "cell_class": {}, "outside_regions": 0}
    summary = {"markers": {"required_stopped": [], "stop": []}, "typed": True,
               "sections": ["s"], "per_section": {"s": per}, "not_typed": {}}
    gate = stage.gate(SimpleNamespace(staging=tmp_path), summary)
    check = next(c for c in gate.checks if c.name == "T-cell calls without CD3e and CD4-or-CD8 support")
    assert check.blocking and not check.passed


def test_support_adaptation_cannot_weaken_an_existing_or_group():
    rules = [Rule("all", "X", (("A", "anypos"), ("B", "anypos")))]
    with pytest.raises(RuleSourceError, match="cannot replace or merge"):
        spleen._require_anypos(rules, {"X": ["CD4", "CD8"]}, frozenset({"A", "B", "CD4", "CD8"}))
