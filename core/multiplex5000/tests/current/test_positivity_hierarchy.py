"""positivity.hierarchy: per-image calls named by a tissue's recorded hierarchy (spleen, maintainer 2026-11-01)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiplex5000.positivity import hierarchy as H
from multiplex5000.positivity import registry as RG
from multiplex5000.tissues import spleen
from multiplex5000.tissues.lymph_node import Rule, RuleSource

REPO = Path(__file__).resolve().parents[4]
MANIFEST = REPO / "core/panels/data/tissues/spleen/manifest.yml"

RULES = RuleSource("0" * 64, (
    Rule("all", "Immune", (("CD45", "pos"),)),
    Rule("Immune", "T_cell", (("CD3e", "pos"), ("CD20", "allneg"))),
    Rule("T_cell", "CD4_T", (("CD4", "pos"), ("CD8", "allneg"))),
    Rule("T_cell", "CD8_T", (("CD8", "pos"), ("CD4", "allneg"))),
    Rule("CD8_T", "CD8_a", (("GZMB", "pos"),)),
    Rule("CD8_T", "CD8_b", (("CD57", "pos"),)),
    Rule("Immune", "B_cell", (("CD20", "pos"), ("CD3e", "allneg"))),
    Rule("all", "Stromal", (("Vim", "pos"),)),
    Rule("all", "Endothelial", (("CD31", "pos"), ("CD3e", "allneg"))),
    Rule("Endothelial", "Endo_CD8", (("CD8", "pos"),)),
    Rule("any", "Prolif", (("Ki67", "pos"),)),
), (), {}, "all", "any", ("Immune", "Stromal", "Endothelial"))
MARKERS = ["CD45", "CD3e", "CD20", "CD4", "CD8", "GZMB", "CD57", "Vim", "CD31", "Ki67"]


def cells(*positives):
    """One row of tri-state calls per cell from the markers positive in it (all others measured negative)."""
    return {m: np.array([1 if m in p else 0 for p in positives], np.int8) for m in MARKERS}


def evidence(n, canary=(), fold=None, cover=None, typing=MARKERS):
    C = np.zeros((n, len(MARKERS)), bool)
    for i, m in canary:
        C[i, MARKERS.index(m)] = True
    F = np.zeros((n, len(MARKERS)))
    COV = np.zeros((n, len(MARKERS)))
    for i, m, v in fold or ():
        F[i, MARKERS.index(m)] = v
    for i, m, v in cover or ():
        COV[i, MARKERS.index(m)] = v
    return {"columns": MARKERS, "C": C, "F": F, "COV": COV, "interior_valid": np.ones(n, bool), "typing": list(typing)}


def name(calls, ev, structure_only=("Endothelial",), in_structure=None, min_depth=2):
    e = H.evaluate(RULES, calls, structure_only, in_structure)
    ev = dict(ev)
    ev["P"] = np.column_stack([calls[m] == 1 for m in MARKERS])
    return H.name_cells(RULES, e, ev, margin=1, min_depth=min_depth, lead=0.15)


def test_a_unique_match_names_the_cell_and_its_broad_step():
    out = name(cells({"CD45", "CD3e", "CD4", "Ki67"}, {"Vim"}), evidence(2))
    assert list(out["refined_class"]) == ["CD4_T", "Stromal"]
    assert list(out["cell_class"]) == ["T_cell", "Stromal"]     # Immune's children have subtypes: a compartment
    assert list(out["status"]) == ["unique", "unique"]
    assert out.loc[0, "states"] == "Prolif"


def test_a_structure_only_branch_is_unavailable_without_its_objects_never_negative():
    calls = cells({"CD31", "CD8"}, {"CD31"}, {"CD45", "CD3e", "CD8"})
    out = name(calls, evidence(3))
    assert list(out["status"][:2]) == ["unavailable", "unavailable"]
    assert out.loc[2, "refined_class"] == "CD8_T"               # a CD3e+ CD8+ cell is unaffected
    inside = name(calls, evidence(3), in_structure={"Endothelial": np.array([True, False, False])})
    assert inside.loc[0, "refined_class"] == "Endo_CD8" and inside.loc[0, "cell_class"] == "Endothelial"
    assert inside.loc[1, "status"] in ("no_match", "unavailable") and inside.loc[1, "refined_class"] == ""


def test_a_lineage_tie_is_settled_by_the_canary_and_otherwise_stays_ambiguous():
    calls = cells({"CD45", "CD3e", "CD4", "Vim"}, {"CD45", "CD3e", "CD4", "Vim"})
    ev = evidence(2, canary=[(0, "CD45"), (0, "CD3e"), (0, "CD4"),     # cell 0: only the T side canary-positive
                             (1, "CD45"), (1, "Vim")])                 # cell 1: depth 1 each, below CANARY_MIN_DEPTH
    out = name(calls, ev)
    assert (out.loc[0, "refined_class"], out.loc[0, "status"], out.loc[0, "tie_basis"]) == ("CD4_T", "tie_settled", "canary")
    assert (out.loc[1, "refined_class"], out.loc[1, "status"]) == ("", "lineage_ambiguous")
    assert out.loc[1, "cell_class"] == ""


def test_a_plain_only_tie_must_resolve_too():
    calls = cells({"CD45", "CD3e", "CD4", "Vim"})
    out = name(calls, evidence(1))          # nothing canary-positive: plain depth 3 vs 1, margin 1 -> resolved
    assert (out.loc[0, "refined_class"], out.loc[0, "tie_basis"]) == ("CD4_T", "plain")
    out = name(cells({"CD45", "Vim"}), evidence(1))     # plain depth 1 vs 1 and no settling lead: stays ambiguous
    assert out.loc[0, "status"] == "lineage_ambiguous"


def test_a_subtype_tie_uses_brightness_and_interior_coverage_else_keeps_the_shared_lineage():
    calls = cells(*[{"CD45", "CD3e", "CD8", "GZMB", "CD57"}] * 3)
    canary = [(i, m) for i in range(3) for m in ("GZMB", "CD57")]
    ev = evidence(3, canary=canary,
                  fold=[(0, "GZMB", 2.0), (0, "CD57", 1.0), (1, "GZMB", 1.0), (1, "CD57", 2.0),
                        (2, "GZMB", 2.0), (2, "CD57", 1.0)],
                  cover=[(0, "GZMB", 0.8), (0, "CD57", 0.3), (1, "GZMB", 0.3), (1, "CD57", 0.8),
                         (2, "GZMB", 0.5), (2, "CD57", 0.45)])            # brighter but no clear interior lead
    out = name(calls, ev, min_depth=1)
    assert list(out["refined_class"]) == ["CD8_a", "CD8_b", "CD8_T"]
    assert list(out["status"]) == ["tie_settled", "tie_settled", "subtype_ambiguous"]
    assert set(out["cell_class"]) == {"T_cell"}
    # markers without a typing role detect nothing: the tie keeps the lineage
    out = name(calls, evidence(3, canary=canary, typing=[m for m in MARKERS if m not in ("GZMB", "CD57")]), min_depth=1)
    assert set(out["refined_class"]) == {"CD8_T"} and set(out["status"]) == {"subtype_ambiguous"}


def test_calls_are_unavailable_where_a_value_or_level_is_missing():
    class Reg:
        def column(self, m, meas):
            return f"{m}_{meas}"
    d = pd.DataFrame({"A_whole_cell": [0.0, 10.0, np.nan], "B_whole_cell": [5.0, 5.0, 5.0]})
    levels = pd.DataFrame({"measurement": ["whole_cell", "whole_cell"], "level_log": [np.log1p(1.0), np.nan]},
                          index=["A", "B"])
    calls = H.marker_calls(d, Reg(), levels, ["A", "B"])
    assert calls["A"].tolist() == [0, 1, -1]
    assert calls["B"].tolist() == [-1, -1, -1]


def test_registry_names_reach_the_sheet_and_collisions_are_refused():
    source = spleen.load_rule_source(MANIFEST)
    reg = RG.load("phenocycler-pc67", "spleen")
    names = H.sheet_names(source, reg.table.index)
    assert names["Pan_Cytokeratin"] == "Pan-Cytokeratin" and names["PD_L1"] == "PD-L1"
    assert names["HLA_DR"] == "HLA-DR" and names["CD45RO"] == "CD45RO" and "Keratin_14" not in names
    sheet = {m for r in source.rules for m in r.markers}
    missing = sorted(m for m in sheet if m not in names.values())
    assert missing == ["B3TUBB", "Keratin 5"]           # every sheet marker but the two PC67 lacks
    with pytest.raises(H.HierarchyError, match="cannot be told apart"):
        H.sheet_names(RuleSource("0" * 64, (Rule("all", "X", (("PD-L1", "pos"), ("PD_L1", "allneg"))),),
                                 (), {}, "all", "any", ("X",)), ["PD_L1"])


def test_spleen_littoral_cells_wait_for_vessels_and_cdcs_are_myeloid():
    source = spleen.load_rule_source(MANIFEST)
    reg = RG.load("phenocycler-pc67", "spleen")
    names = H.sheet_names(source, reg.table.index)
    sheet = list(names.values())

    def call(*positive):
        return {m: np.array([1 if m in positive else 0], np.int8) for m in sheet}

    branches = spleen.structure_only_branches(MANIFEST)
    littoral = H.evaluate(source, call("CD8", "CD31", "Caveolin"), branches)
    assert littoral.states["Endothelial_Sinusoidal_CD8+"][0] == -1
    assert littoral.states["CD8_T_cytotoxic"][0] == 0 and littoral.states["T_cell"][0] == 0
    ev = {"columns": sheet, "P": np.zeros((1, len(sheet)), bool), "C": np.zeros((1, len(sheet)), bool),
          "F": np.zeros((1, len(sheet))), "COV": np.zeros((1, len(sheet))), "interior_valid": np.ones(1, bool),
          "typing": sheet}
    out = H.name_cells(source, littoral, ev, 1, 2, 0.15)
    assert (out.loc[0, "refined_class"], out.loc[0, "status"]) == ("", "unavailable")
    cdc = H.name_cells(source, H.evaluate(source, call("CD11c", "HLA-DR"), branches), ev, 1, 2, 0.15)
    assert (cdc.loc[0, "refined_class"], cdc.loc[0, "cell_class"]) == ("cDC", "Myeloid")
    treg = H.name_cells(source, H.evaluate(source, call("CD45", "CD3e", "CD4", "FOXP3"), branches), ev, 1, 2, 0.15)
    assert (treg.loc[0, "refined_class"], treg.loc[0, "cell_class"]) == ("Treg", "T_cell")
