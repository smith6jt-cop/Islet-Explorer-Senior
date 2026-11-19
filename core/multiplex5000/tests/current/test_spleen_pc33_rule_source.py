"""The smaller assay preserves observed evidence and reports unmeasured clauses."""
from pathlib import Path

import numpy as np
import pytest
import yaml

from multiplex5000.positivity import hierarchy as H
from multiplex5000.tissues import spleen
from multiplex5000.tissues.lymph_node import RuleSourceError

MANIFEST = Path(__file__).resolve().parents[4] / "core/panels/data/tissues/spleen/manifest.yml"
PANEL = "spleen_pc33"


def test_pc33_inventory_and_pc67_default_are_distinct():
    pc33 = spleen.observed_markers(MANIFEST, panel_id=PANEL)
    assert len(pc33) == 33
    assert {"CD3e", "CD20", "CD4", "CD8", "FOXP3", "HLA-DR", "PDL1"} <= pc33
    assert not {"CD79a", "CD11b", "Pan-Cytokeratin", "EpCAM", "Granzyme B", "CD56"} & pc33
    assert len(spleen.observed_markers(MANIFEST)) == 67
    assert spleen.load_rule_source(MANIFEST) == spleen.load_rule_source(MANIFEST, panel_id="spleen_pc67")
    assert not spleen.panel_adaptation_audit(MANIFEST)["omitted_clauses"]


def test_pc33_projection_preserves_every_measured_clause_and_defining_positive():
    original = spleen.load_rule_source(MANIFEST)
    projected = spleen.load_rule_source(MANIFEST, panel_id=PANEL)
    observed = spleen.observed_markers(MANIFEST, panel_id=PANEL)
    audit = spleen.panel_adaptation_audit(MANIFEST, panel_id=PANEL)
    removed = {(item["phenotype"], item["marker"]): item["token"]
               for item in audit["omitted_clauses"]}
    actual = {}
    assert original.by_phenotype.keys() == projected.by_phenotype.keys()
    for name, before in original.by_phenotype.items():
        after = projected.by_phenotype[name]
        assert before.parent == after.parent
        assert set(after.constraints) <= set(before.constraints)
        for marker, token in set(before.constraints) - set(after.constraints):
            assert marker not in observed and original.panel_aliases.get(marker) not in observed
            assert token not in {"pos", "allpos"}
            if token == "anypos":
                assert any(t == "anypos" and (m in observed or original.panel_aliases.get(m) in observed)
                           for m, t in after.constraints)
            actual[name, marker] = token
    assert actual == removed
    assert {("T_cell", "CD79a"), ("Immune", "Pan-Cytokeratin"),
            ("B_cell", "CD79a")} <= removed.keys()
    assert set(removed) <= set(projected.omitted)


def test_pc33_has_specific_measured_profiles_but_no_invented_positive_evidence():
    source = spleen.load_rule_source(MANIFEST, panel_id=PANEL)
    missing = source.missing_by_phenotype(spleen.observed_markers(MANIFEST, panel_id=PANEL))
    for name in ("CD4_T_helper", "CD8_T_cytotoxic", "Treg", "Tfh", "Activated_CD4",
                 "Memory_CD4", "Exhausted_CD8", "Senescent_CD8", "GC_B_cell",
                 "CD21_positive_B_cell", "Macrophage", "cDC", "CD141_positive_cDC",
                 "Fibroblast", "FRC", "SMA_positive_stromal", "Epithelial"):
        assert missing[name] == (), (name, missing[name])
    unsupported = {
        "NK_cell": "CD56", "Neutrophil": "MPO", "Monocyte": "CD11b",
        "Plasma_cell": "CD79a", "Plasmablast": "CD79a", "Naive_B_cell": "Bcl-2",
        "Memory_B_cell": "Bcl-2", "M1_macrophage": "iNOS", "M2_macrophage": "CD206",
        "DC_SIGN_DC": "CD209", "Effector_CD8": "Granzyme B", "CD8_CTL": "Granzyme B",
        "Memory_CD8": "Bcl-2", "Activated_Treg": "CD39", "Exhausted_Treg": "LAG3",
        "CD4_exhausted": "LAG3", "Pericyte": "Caveolin", "Smooth_muscle": "Caveolin",
        "Endothelial": "Caveolin", "Progenitor_Epithelial": "SOX2",
    }
    for name, marker in unsupported.items():
        assert marker in missing[name]
    assert "DN_T_cell" not in source.by_phenotype
    assert spleen.structure_only_branches(MANIFEST, panel_id=PANEL) == {"Endothelial"}


def _evaluate(positive, unavailable=()):
    source = spleen.load_rule_source(MANIFEST, panel_id=PANEL)
    names = H.sheet_names(source, spleen.observed_markers(MANIFEST, panel_id=PANEL))
    calls = {marker: np.array([-1 if marker in unavailable else int(marker in positive)], np.int8)
             for marker in names.values()}
    policy = spleen.naming_policy(MANIFEST, panel_id=PANEL)
    result = H.evaluate_profiles(source, calls, spleen.structure_only_branches(MANIFEST, panel_id=PANEL),
                                 defer_exclusions=policy["defer_exclusions"])
    return result, calls


@pytest.mark.parametrize("positive, expected", [
    ({"CD3e", "CD4"}, "CD4_T_helper"),
    ({"CD3e", "CD8"}, "CD8_T_cytotoxic"),
    ({"CD3e", "CD4", "FOXP3"}, "Treg"),
    ({"CD11c", "HLA-DR", "CD141"}, "CD141_positive_cDC"),
    ({"CD20", "CD38", "Ki67"}, "GC_B_cell"),
])
def test_pc33_observed_profiles_match_without_manufacturing_missing_calls(positive, expected):
    result, calls = _evaluate(positive)
    assert result.states[expected][0] == 1
    assert "CD79a" not in calls and "Pan-Cytokeratin" not in calls


@pytest.mark.parametrize("positive, unavailable", [
    ({"CD3e", "CD45RO"}, ()),
    ({"CD3e", "CD45RO"}, ("CD4",)),
    ({"CD4", "CD45RO"}, ()),
])
def test_pc33_t_ancestor_never_matches_without_measured_cd3_and_cd4_or_cd8(positive, unavailable):
    result, _ = _evaluate(positive, unavailable)
    assert result.states["T_cell"][0] != 1
    assert result.states["CD45RO_positive_T_cell"][0] != 1


def test_acquired_measurement_failure_remains_unavailable():
    result, _ = _evaluate({"CD3e", "CD4"}, unavailable=("FOXP3",))
    assert result.states["CD4_T_helper"][0] == result.states["Treg"][0] == -1
    result, _ = _evaluate({"CD11c", "CD66"})
    assert result.states["Neutrophil"][0] == -1  # MPO was not acquired.


@pytest.mark.parametrize("field, phenotype, markers, match", [
    ("omit_absent_allneg", "CD4_T_helper", ["FOXP3"], "measures it"),
    ("omit_absent_allneg", "NK_cell", ["CD56"], "only an allneg"),
    ("omit_absent_anypos", "Endothelial", ["Caveolin"], "empty its OR group"),
])
def test_panel_adaptation_refuses_observed_exclusions_and_defining_positive_removal(
        tmp_path, field, phenotype, markers, match):
    manifest = yaml.safe_load(MANIFEST.read_text())
    manifest["rule_table"]["source_manifest"] = str(MANIFEST.parent.parent / "pancreatic_lymph_node/manifest.yml")
    manifest["protein_panels"][PANEL]["file"] = str(MANIFEST.parent / "pc33_export_header.tsv")
    manifest["rule_table"]["panel_adaptations"][PANEL][field][phenotype] = markers
    fixture = tmp_path / "manifest.yml"
    fixture.write_text(yaml.safe_dump(manifest))
    with pytest.raises(RuleSourceError, match=match):
        spleen.load_rule_source(fixture, panel_id=PANEL)


def test_unknown_panel_is_not_silently_treated_as_pc67():
    with pytest.raises(RuleSourceError, match="unknown spleen protein panel"):
        spleen.load_rule_source(MANIFEST, panel_id="unknown")
