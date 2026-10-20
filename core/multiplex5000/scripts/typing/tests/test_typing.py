"""Invariants of the typing lane.

These pin the decisions that were hard to get right, so a later change cannot quietly undo them.

    ../../../../.venv/bin/python -m pytest -q core/multiplex5000/scripts/typing/tests
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve()
LANE = HERE.parents[1]
REPO = HERE.parents[5]
sys.path.insert(0, str(LANE))

import compartments as C  # noqa: E402
import lineages as L  # noqa: E402


def _registry_entries() -> list[dict]:
    reg = json.loads((REPO / "core/multiplex5000/multiplex5000/marker_registry.json").read_text())
    return reg["markers"] if isinstance(reg, dict) and "markers" in reg else reg


def registry_markers() -> list[str]:
    return [m["name"] for m in _registry_entries()]


# The metabolic panel carried by the two pLNpanc slides (0080870/6533, 0081366/112), added
# to the registry on 2026-09-16 so REDSEA could run on them. The maintainer rules the same
# day that these are deferred: measured, never typed, and so not required to declare a
# typing compartment.
#
# Listed BY NAME on purpose. `calibration_status == "deferred"` does not identify them --
# 45 of 72 registry markers are deferred and 32 of those are placed. Neither does `uses`:
# 42 markers are exclusively `deferred_*` and 30 of those are placed. Exempting by either
# field would quietly drop ~42 markers from the compartment check instead of these 12.
DEFERRED_METABOLIC_PANEL = frozenset({
    "ASCT2", "ATP5A", "CD141", "CD21", "CS", "G6PD",
    "GLUT1", "HK1", "Histone_H3_(p_Ser28)", "IDH2", "LDHA", "SDHA",
})


def typing_markers() -> list[str]:
    """Registry markers that typing may use, so must be placed in a compartment."""
    return [m["name"] for m in _registry_entries()
            if m["name"] not in DEFERRED_METABOLIC_PANEL]


# --------------------------------------------------------------------- the panel

def test_only_iapp_and_dapi_are_held_back():
    """The maintainer's rule: IAPP and DAPI are the only channels unusable for typing."""
    assert set(C.NOT_FOR_TYPING) == {"IAPP", "DAPI"}


def test_every_typing_marker_has_a_compartment():
    """A channel typing can use must be placed deliberately, never default to the cell mean."""
    unplaced = [m for m in typing_markers()
                if m not in C.NOT_FOR_TYPING and m not in C._ASSIGNMENT]
    assert unplaced == [], f"markers with no declared compartment: {unplaced}"


def test_the_metabolic_exemption_covers_exactly_the_twelve_it_names():
    """The exemption is a named maintainer ruling, and it may not widen on its own.

    Every exempted marker must still be in the registry, still declare only `deferred_*`
    uses, and still be the ONLY thing the compartment test lets through unplaced. If a 13th
    channel arrives without a compartment this fails and someone has to decide about it,
    which is the behaviour the original test was protecting.
    """
    entries = {m["name"]: m for m in _registry_entries()}
    for name in DEFERRED_METABOLIC_PANEL:
        assert name in entries, f"{name} is exempted but not in the registry"
        uses = entries[name]["uses"]
        assert uses and all(str(u).startswith("deferred_") for u in uses), (
            f"{name} is exempted from placement but declares a live use: {uses}"
        )

    unplaced = {m for m in registry_markers()
                if m not in C.NOT_FOR_TYPING and m not in C._ASSIGNMENT}
    assert unplaced == set(DEFERRED_METABOLIC_PANEL), (
        "the set of unplaced markers drifted from the recorded exemption: "
        f"unexpected={sorted(unplaced - DEFERRED_METABOLIC_PANEL)}, "
        f"now placed={sorted(DEFERRED_METABOLIC_PANEL - unplaced)}"
    )


def test_no_marker_is_in_two_compartments():
    groups = (C.NUCLEUS, C.CYTOPLASM, C.MEMBRANE)
    seen: dict[str, int] = {}
    for i, g in enumerate(groups):
        for m in g:
            assert m not in seen, f"{m} declared in two compartments"
            seen[m] = i


def test_nuclear_proteins_are_read_in_the_nucleus():
    """FOXP3 and TP63 are transcription factors; a whole-cell mean dilutes them over cytoplasm
    they cannot occupy, which is what coupled them to DAPI."""
    for m in ("FOXP3", "TP63"):
        assert C.compartment_for(m) == "Nucleus"


def test_hormones_are_read_in_the_cytoplasm():
    """Insulin, glucagon and somatostatin live in secretory granules."""
    for m in ("INS", "GCG", "SST"):
        assert C.compartment_for(m) == "Cytoplasm"


def test_surface_markers_are_read_in_the_membrane_band():
    for m in ("CD3e", "CD4", "CD8", "CD20", "CD31", "CD68", "E_cadherin"):
        assert C.compartment_for(m) == "Membrane"


def test_unknown_marker_raises_rather_than_defaulting():
    with pytest.raises(KeyError):
        C.compartment_for("NotAMarker")


# --------------------------------------------------------------------- the label space

def test_no_lineage_uses_a_held_back_channel():
    bad = [m for m in L.identity_markers() if m in C.NOT_FOR_TYPING]
    assert bad == [], f"lineages must not use {bad}"


def test_lineages_use_the_markers_an_earlier_lane_dropped():
    """B3TUBB, CD34, CD56, CD66 and CD11b are usable channels; dropping B3TUBB also deleted the
    Neural lineage and folded granulocytes into myeloid."""
    used = set(L.identity_markers())
    for m in ("B3TUBB", "CD34", "CD56", "CD66", "CD11b"):
        assert m in used, f"{m} is usable and must not be dropped"
    assert "Neural" in L.LINEAGE_MARKERS
    assert "Granulocyte" in L.LINEAGE_MARKERS


def test_no_state_marker_defines_a_lineage():
    """State markers say what a cell is doing, not what it is.  A Ki67-high T cell is a T cell."""
    reg = json.loads((REPO / "core/multiplex5000/multiplex5000/marker_registry.json").read_text())
    markers = reg["markers"] if isinstance(reg, dict) and "markers" in reg else reg
    state = {m["name"] for m in markers if m.get("kind") == "state"}
    leaked = state & set(L.identity_markers())
    assert leaked == set(), f"state markers used for identity: {sorted(leaked)}"


def test_subtypes_are_resolved_within_a_declared_parent():
    for parent in L.SUBTYPE_MARKERS:
        assert parent in L.LINEAGE_MARKERS, f"{parent} subtypes have no parent lineage"


def test_exocrine_is_a_single_label():
    """The 60-plex carries no acinar-specific and no ductal-specific marker, so the exocrine
    compartment cannot be split and must not pretend to be."""
    assert "Exocrine" in L.LINEAGE_MARKERS
    assert "Exocrine" not in L.SUBTYPE_MARKERS
    names = {n.lower() for n in L.LINEAGE_MARKERS}
    names |= {s.lower() for subs in L.SUBTYPE_MARKERS.values() for s in subs}
    assert not (names & {"acinar", "ductal", "duct"})




def _raw_like(Z, labels, cluster, cols, base=0.02, lifted=0.70):  # lifted>MIN_CONF so most cells count positive
    """A raw-confidence matrix matching a planted Z: baseline everywhere, raised where planted."""
    X = np.full(Z.shape, base)
    if cluster is not None:
        X[np.ix_(labels == cluster, cols)] = lifted
    return X

# --------------------------------------------------------------------- the statistic

# --------------------------------------------------------------------- the export contract

def test_export_columns_match_the_groovy_importer():
    """The importer refuses a CSV missing any required column, so the two lists must agree."""
    import export_qupath as E
    groovy = (REPO / "core/multiplex5000/scripts/groovy/import_cell_assignments.groovy").read_text()
    block = groovy.split("def requiredColumns = [", 1)[1].split("]", 1)[0]
    required = [t.strip().strip('"') for t in block.split(",") if t.strip()]
    missing = [c for c in required if c not in E.COLUMNS]
    assert missing == [], f"export is missing columns the importer requires: {missing}"


def test_uncertainty_never_receives_a_biological_class():
    """``ambiguous`` and ``unavailable`` get their own QuPath classes; the candidate is carried
    separately in best_broad_type so a reviewer can see it without it being asserted."""
    import export_qupath as E
    assert "best_broad_type" in E.COLUMNS
    groovy = (REPO / "core/multiplex5000/scripts/groovy/import_cell_assignments.groovy").read_text()
    assert 'className = "Phenocycler - Ambiguous"' in groovy
    assert 'className = "Phenocycler - Unavailable"' in groovy


def test_every_declared_lineage_has_a_qupath_colour():
    """A lineage with no colour renders the same grey as every other, which makes the QuPath
    review -- the only place these calls can be judged -- unreadable."""
    groovy = (REPO / "core/multiplex5000/scripts/groovy/import_cell_assignments.groovy").read_text()
    palette = groovy.split("def broadColors = [", 1)[1].split("]\n", 1)[0]
    for lineage in L.LINEAGE_MARKERS:
        assert f'"{lineage}"' in palette, f"{lineage} has no colour in the importer"


# --------------------------------------------------------------------- positivity

def test_anchor_threshold_is_declared_once_and_only_labels_evidence():
    import export_qupath as E
    assert E.ANCHOR_CONF == 0.5
    src = Path(E.__file__).read_text()
    assert src.count("ANCHOR_CONF") == 2, "the anchor boundary must be used in exactly one place"


# --------------------------------------------------------------------- multiplicity

# --------------------------------------------------------------------- the merge stage

def test_merge_ties_break_toward_more_clusters():
    """Equally stable counts must resolve to the finer one.  Measured on a synthetic set with one
    rare group, breaking toward the coarser count merged it away (purity 1.00 -> 0.11)."""
    import pandas as pd
    import metacluster as MC
    sweep = pd.DataFrame({"n_final": [5, 6, 7], "stability_ari": [1.0, 1.0, 0.9],
                          "smallest_merged_cluster": [39, 5, 5], "singletons": [0, 0, 0]})
    assert MC.choose(sweep, 4, 12) == 6


def test_merge_keeps_a_rare_population_separate():
    """The whole reason for clustering fine and then merging is that a 0.4%-abundance population
    disappears if you cluster straight to a few dozen."""
    import metacluster as MC
    rng = np.random.default_rng(0)
    true = np.concatenate([np.repeat(np.arange(5), 39), np.repeat(5, 5)])
    centres = rng.normal(size=(6, 34)) * 2
    P = centres[true] + rng.normal(size=(200, 34)) * 0.15
    n = MC.choose(MC.sweep(P, 4, 12, draws=10), 4, 12)
    res = MC.rare_population_survives(np.arange(200), MC.merge(P, n), true == 5)
    assert res["completeness"] == 1.0
    assert res["purity"] > 0.9, f"rare population was merged away: {res}"


def test_merge_weights_clusters_equally_not_by_size():
    """The merge runs on the profile matrix, one row per cluster, so a 44-cell cluster carries the
    same weight as a 40,000-cell one.  That is what carries a rare type through."""
    import inspect
    import metacluster as MC
    src = inspect.getsource(MC.merge)
    assert "pdist(P[ok]" in src, "merge must operate on profiles, not on cells"



def test_subtypes_within_a_parent_are_mutually_exclusive():
    """They are resolved by argmax, so a subtype that is a SUBSET of another steals its cells.
    Treg as a peer of CD4_T produced 11,463 Tregs -- a third of every T cell on donor 6442."""
    assert "Treg" not in L.SUBTYPE_MARKERS["T_NK"]
    assert "M2_Macrophage" not in L.SUBTYPE_MARKERS["Myeloid"]
    assert L.SUBTYPE_MODIFIERS["CD4_T"][0] == "FOXP3_high"
    assert L.SUBTYPE_MODIFIERS["Macrophage"][0] == "CD163_CD206_high"


def test_modifiers_are_named_for_the_measurement_not_the_inference():
    """A modifier flags the upper quartile of its parent, so its size is 25% by construction.
    Naming that quartile 'Treg' would read as a biological call the panel cannot support --
    Nimbus scores FOXP3 above its own boundary in 90.4% of cells."""
    for parent, (name, ms) in L.SUBTYPE_MODIFIERS.items():
        assert any(m in name for m in ms), f"{name!r} does not name the markers it measures"
        assert name not in {"Treg", "M2", "M2_Macrophage"}, f"{name!r} asserts an inference"


def test_every_modifier_refines_a_declared_subtype():
    declared = {s for subs in L.SUBTYPE_MARKERS.values() for s in subs}
    for parent in L.SUBTYPE_MODIFIERS:
        assert parent in declared, f"{parent} is not a subtype anything can refine"


def test_islet_composition_records_the_size_artefact():
    """Boundary bleed grows with islet size (+0.52 with double-positive fraction, -0.39 with
    INS-dominant fraction on 6442).  T1D islets are smaller, so a naive cross-donor comparison
    would report them as more beta-rich purely from the artefact."""
    import inspect
    import validate as V
    src = inspect.getsource(V.islet_composition)
    assert "0.52" in src and "condition on" in src, (
        "the size dependence and its consequence for the disease comparison must stay recorded")


def test_islet_grouping_excludes_cells_outside_any_islet():
    """islet_num carries an empty string, not NaN, for cells outside an islet; notna() alone
    grouped 2,533 non-islet endocrine cells into one pseudo-islet."""
    import inspect
    import validate as V
    assert '!= ""' in inspect.getsource(V.islet_composition)


# --------------------------------------------------------------- per-cell resolution

def test_unmeasured_objects_are_omitted_from_the_export_by_default():
    """Objects with no nucleus are filtered out before typing, so writing them back as an
    Unavailable class would CHANGE their classification in QuPath rather than leave them alone."""
    import inspect
    import export_qupath as E
    src = inspect.getsource(E.build)
    assert 'out["assignment_status"] != "unavailable"' in src
    assert "include_unavailable" in inspect.signature(E.build).parameters
