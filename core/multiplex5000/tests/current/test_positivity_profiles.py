"""Regression cases from the spleen QuPath review, 2026-11-02."""
from pathlib import Path

import numpy as np

from multiplex5000.positivity import hierarchy as H
from multiplex5000.tissues import spleen
from multiplex5000.tissues.lymph_node import Rule, RuleSource

RULES = RuleSource("0" * 64, (
    Rule("all", "Immune", (("CD45", "pos"),)),
    Rule("Immune", "T_cell", (("CD3e", "pos"), ("CD20", "allneg"))),
    Rule("Immune", "B_cell", (("CD20", "pos"), ("CD3e", "allneg"))),
    Rule("T_cell", "CD4_T", (("CD4", "pos"), ("CD8", "allneg"), ("FOXP3", "allneg"))),
    Rule("T_cell", "CD8_T", (("CD8", "pos"), ("CD4", "allneg"))),
    Rule("T_cell", "Treg", (("CD4", "pos"), ("FOXP3", "pos"), ("CD8", "allneg"))),
), (), {}, "all", "any", ("Immune",))
MARKERS = ("CD45", "CD3e", "CD20", "CD4", "CD8", "FOXP3")
POLICY = {"Immune": ["CD20", "CD3e", "CD4", "CD8"]}


def classify(positive, *, stronger=None, unavailable=(), inherit_depth=True, raw_pixels=None):
    calls = {m: np.array([-1 if m in unavailable else int(m in positive)], np.int8) for m in MARKERS}
    states = H.evaluate_profiles(RULES, calls, defer_exclusions=POLICY)
    p = np.column_stack([calls[m] == 1 for m in MARKERS])
    f = np.ones(p.shape)
    cover = np.full(p.shape, .2)
    if stronger:
        f[0, MARKERS.index(stronger)] = 4
        cover[0, MARKERS.index(stronger)] = .8
    ev = {"columns": MARKERS, "typing": ["CD3e", "CD20", "CD4", "CD8"],
          "P": p, "C": p.copy(), "F": f, "COV": cover, "interior_valid": np.ones(1, bool)}
    if raw_pixels is not None:
        ev["pixels"] = raw_pixels(MARKERS, [{m: np.repeat(f[0, j], 10) for j, m in enumerate(MARKERS)}])
    result = H.name_cells(RULES, states, ev, 1, 2, 0, require_specific=["Immune", "T_cell"],
                          inherit_depth=inherit_depth)
    return result.iloc[0], calls, states


def test_cd3_cd4_profile_reaches_subtyping_despite_a_competing_b_marker():
    row, calls, states = classify({"CD45", "CD3e", "CD4", "CD20"})
    assert row.refined_class == "CD4_T" and row.cell_class == "T_cell"
    assert states.states["CD4_T"][0] == 1
    assert calls["CD20"][0] == 1
    assert dict(RULES.by_phenotype["T_cell"].constraints)["CD20"] == "allneg"


def test_cd8_is_compared_with_cd4_instead_of_both_subtypes_being_vetoed(raw_conflict_pixels):
    row, _, _ = classify({"CD45", "CD3e", "CD20", "CD4", "CD8"}, stronger="CD8", raw_pixels=raw_conflict_pixels)
    assert row.refined_class == "CD8_T" and row.status == "tie_settled"
    without_parent, _, _ = classify({"CD45", "CD3e", "CD20", "CD4", "CD8"},
                                    stronger="CD8", inherit_depth=False)
    assert without_parent.refined_class == ""


def test_supported_foxp3_cd4_profile_is_not_lost_to_cd8_exclusion(raw_conflict_pixels):
    row, _, _ = classify({"CD45", "CD3e", "CD20", "CD4", "CD8", "FOXP3"}, stronger="CD4", raw_pixels=raw_conflict_pixels)
    assert row.refined_class == "Treg"


def test_cd3_is_still_required_for_a_t_cell():
    row, _, states = classify({"CD45", "CD8", "FOXP3"})
    assert states.states["T_cell"][0] == 0
    assert row.refined_class == ""


def test_missing_competing_measurement_remains_unavailable():
    row, calls, states = classify({"CD45", "CD3e", "CD4"}, unavailable={"CD20"})
    assert calls["CD20"][0] == -1 and states.states["CD4_T"][0] == -1
    assert row.refined_class == "" and row.status == "unavailable"


def test_absence_of_cd4_cd8_does_not_name_a_t_subtype_or_generic_t_cell():
    row, _, _ = classify({"CD45", "CD3e"})
    assert row.refined_class == "" and row.status == "no_supported_subtype"


def test_unresolved_b_vs_t_is_not_served_as_immune():
    row, _, _ = classify({"CD45", "CD3e", "CD20"})
    assert row.refined_class == "" and row.status == "lineage_ambiguous"


def test_spleen_policy_removes_dn_and_keeps_foxp3_specificity():
    manifest = Path(__file__).resolve().parents[4] / "core/panels/data/tissues/spleen/manifest.yml"
    rules = spleen.load_rule_source(manifest).by_phenotype
    assert "DN_T_cell" not in rules
    assert dict(rules["CD4_T_helper"].constraints)["FOXP3"] == "allneg"
    policy = spleen.naming_policy(manifest)
    assert policy["require_specific"] == ["Immune", "T_cell"]
    assert "FOXP3" not in policy["defer_exclusions"]["Immune"]
    assert policy["inherit_depth"] is True
