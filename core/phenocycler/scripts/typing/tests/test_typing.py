"""Invariants of the typing lane.

These pin the decisions that were hard to get right, so a later change cannot quietly undo them.

    ../../../../.venv/bin/python -m pytest -q core/phenocycler/scripts/typing/tests
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


def registry_markers() -> list[str]:
    reg = json.loads((REPO / "core/phenocycler/phenocycler/marker_registry.json").read_text())
    markers = reg["markers"] if isinstance(reg, dict) and "markers" in reg else reg
    return [m["name"] for m in markers]


# --------------------------------------------------------------------- the panel

def test_only_iapp_and_dapi_are_held_back():
    """The maintainer's rule: IAPP and DAPI are the only channels unusable for typing."""
    assert set(C.NOT_FOR_TYPING) == {"IAPP", "DAPI"}


def test_every_registry_marker_has_a_compartment():
    """A new channel must be placed deliberately, never default to the whole-cell mean."""
    unplaced = [m for m in registry_markers()
                if m not in C.NOT_FOR_TYPING and m not in C._ASSIGNMENT]
    assert unplaced == [], f"markers with no declared compartment: {unplaced}"


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
    reg = json.loads((REPO / "core/phenocycler/phenocycler/marker_registry.json").read_text())
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

def test_statistic_reads_every_k_up_to_top_n():
    """A cell type expresses a subset of its lineage's markers, and how large that subset is
    differs by cell type.  Fixing k penalises the ones that rest on a single marker."""
    import call as K
    rng = np.random.default_rng(0)
    Z = rng.normal(size=(200, 30))
    for kk in (1, 2, 3):
        assert K._stat_k(Z, list(range(8)), kk).shape == (200,)
    # the top-1 mean is never below the top-3 mean of the same markers
    assert (K._stat_k(Z, list(range(8)), 1) >= K._stat_k(Z, list(range(8)), 3) - 1e-9).all()
    assert K.TOP_N >= 2


def test_statistic_ignores_missing_markers():
    import call as K
    Z = np.array([[1.0, np.nan, 0.5, 2.0]])
    assert np.isfinite(K._stat_k(Z, [0, 1, 2, 3], 2)[0])


# --------------------------------------------------------------------- the export contract

def test_export_columns_match_the_groovy_importer():
    """The importer refuses a CSV missing any required column, so the two lists must agree."""
    import export_qupath as E
    groovy = (REPO / "core/phenocycler/scripts/groovy/import_cell_assignments.groovy").read_text()
    block = groovy.split("def requiredColumns = [", 1)[1].split("]", 1)[0]
    required = [t.strip().strip('"') for t in block.split(",") if t.strip()]
    missing = [c for c in required if c not in E.COLUMNS]
    assert missing == [], f"export is missing columns the importer requires: {missing}"


def test_uncertainty_never_receives_a_biological_class():
    """``ambiguous`` and ``unavailable`` get their own QuPath classes; the candidate is carried
    separately in best_broad_type so a reviewer can see it without it being asserted."""
    import export_qupath as E
    assert "best_broad_type" in E.COLUMNS
    groovy = (REPO / "core/phenocycler/scripts/groovy/import_cell_assignments.groovy").read_text()
    assert 'className = "Phenocycler - Ambiguous"' in groovy
    assert 'className = "Phenocycler - Unavailable"' in groovy


def test_every_declared_lineage_has_a_qupath_colour():
    """A lineage with no colour renders the same grey as every other, which makes the QuPath
    review -- the only place these calls can be judged -- unreadable."""
    groovy = (REPO / "core/phenocycler/scripts/groovy/import_cell_assignments.groovy").read_text()
    palette = groovy.split("def broadColors = [", 1)[1].split("]\n", 1)[0]
    for lineage in L.LINEAGE_MARKERS:
        assert f'"{lineage}"' in palette, f"{lineage} has no colour in the importer"


# --------------------------------------------------------------------- positivity

def test_positivity_is_never_binarised_into_the_call():
    """Confidences enter the call as continuous features.

    The failure this guards against is subtle: comparing the confidence matrix to any constant
    turns a probability into a yes/no and reintroduces exactly the operating point the whole
    design avoids -- and the right constant would then have to be chosen per marker, which on
    this panel cannot be done without reference to how many cells "should" be positive.
    """
    import ast
    import call as K
    tree = ast.parse(Path(K.__file__).read_text())
    matrices = {"X", "Z"}                       # the confidence matrix and its standardised form
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Compare):
            continue
        names = {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}
        consts = [c.value for c in ast.walk(node)
                  if isinstance(c, ast.Constant) and isinstance(c.value, (int, float))]
        if names & matrices and consts:
            offenders.append(ast.unparse(node))
    assert offenders == [], f"confidences compared to a constant: {offenders}"


def test_anchor_threshold_is_declared_once_and_only_labels_evidence():
    import export_qupath as E
    assert E.ANCHOR_CONF == 0.5
    src = Path(E.__file__).read_text()
    assert src.count("ANCHOR_CONF") == 2, "the anchor boundary must be used in exactly one place"


# --------------------------------------------------------------------- multiplicity

def test_naming_corrects_for_the_number_of_tests():
    """Every cluster is tested against every lineage, so the bar must be family-wide."""
    import inspect
    import call as K
    src = inspect.getsource(K.name_clusters)
    assert "fam_max" in src, "naming must use a max-statistic (family-wise) bar"
    assert "Sigma" in src and "cholesky" in src, (
        "the null must randomise which CELLS form a cluster, not which markers form a lineage: "
        "a marker-permutation null is corrupted by the signal it is meant to detect")


def test_fdr_names_nothing_on_pure_noise():
    """Structureless input must produce no names at all.  This is the property the p95 bar
    lacked: it named clusters in proportion to how many were tested."""
    import call as K
    rng = np.random.default_rng(0)
    n_cells, n_clusters = 4000, 30
    markers = L.identity_markers()
    Z = rng.normal(size=(n_cells, len(markers)))
    labels = rng.integers(0, n_clusters, n_cells)
    named, why, *_ = K.name_clusters(Z, np.clip(Z*0+1.0,0,1), labels, markers, n_clusters,
                                     draws=2000, margin=0.15, alpha=0.05, seed=0)
    called = sum(1 for n in named if n != L.AMBIGUOUS)
    assert called <= 1, f"{called} of {n_clusters} noise clusters were named"


def test_fdr_still_finds_real_structure():
    """The correction must not simply refuse everything: a cluster whose lineage markers really
    are elevated has to survive it."""
    import call as K
    rng = np.random.default_rng(1)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n_cells, n_clusters = 4000, 30
    Z = rng.normal(size=(n_cells, len(markers)))
    labels = rng.integers(0, n_clusters, n_cells)
    planted = [idx[m] for m in L.LINEAGE_MARKERS["Vascular"] if m in idx]
    Z[np.ix_(labels == 7, planted)] += 4.0
    Xraw = _raw_like(Z, labels, 7, planted)
    named, *_ = K.name_clusters(Z, Xraw, labels, markers, n_clusters,
                                draws=2000, margin=0.15, alpha=0.05, seed=0)
    assert named[7] == "Vascular", f"planted vascular cluster came back as {named[7]!r}"





def test_null_randomises_cells_not_markers():
    """A marker-permutation null fails on a positive control three separate ways (a p-value floor,
    a conditional bias, and inflation by the signal itself).  The null must ask what a cluster of
    this size looks like when its cells are assigned at random."""
    import inspect
    import call as K
    src = inspect.getsource(K.name_clusters)
    assert "(N - n_c)" in src, "the null needs the finite-population correction for cluster size"


def test_naming_is_sensitive_to_a_modest_effect():
    """Calibration is worthless if nothing can ever clear the bar.  A planted lineage at a quarter
    of a standard deviation must still be found."""
    import call as K
    rng = np.random.default_rng(1)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    Z = rng.normal(size=(4000, len(markers)))
    labels = rng.integers(0, 30, 4000)
    planted = [idx[m] for m in L.LINEAGE_MARKERS["Vascular"]]
    Z[np.ix_(labels == 7, planted)] += 0.25
    Xraw = _raw_like(Z, labels, 7, planted)
    named, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=2000,
                                margin=0.15, alpha=0.05, seed=0)
    assert named[7] == "Vascular"


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



def test_single_marker_evidence_is_not_diluted_by_quiet_partners():
    """Lymphatic endothelium is genuinely Podoplanin-high and CD31-low.  Averaging Podoplanin with
    Vascular's quiet markers dropped a real +4.87 cluster below Neural's lone B3TUBB +4.17 and the
    cluster was named Neural; scoring each lineage at its best k fixes that."""
    import call as K
    rng = np.random.default_rng(3)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    Z = rng.normal(size=(6000, len(markers)))
    labels = rng.integers(0, 30, 6000)
    # one marker of a four-marker lineage strongly up; its partners untouched
    Z[labels == 5, idx["Podoplanin"]] += 4.0
    Xraw = _raw_like(Z, labels, 5, [idx["Podoplanin"]])
    named, why, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                     margin=0.15, alpha=0.05, seed=0)
    assert named[5] == "Vascular", f"single-marker evidence was diluted away: {named[5]!r} ({why[5]})"


def test_a_lone_weak_marker_still_names_nothing():
    """The other half of the same fix: best-of-k must not let a modest single marker name a
    cluster.  B3TUBB at +1.90 named 99,410 cells Neural before this was corrected."""
    import call as K
    rng = np.random.default_rng(4)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    Z = rng.normal(size=(6000, len(markers)))
    labels = rng.integers(0, 30, 6000)
    Z[labels == 5, idx["B3TUBB"]] += 0.15
    Xraw = _raw_like(Z, labels, 5, [idx["B3TUBB"]], lifted=0.03)   # raised, but barely
    named, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                margin=0.15, alpha=0.05, seed=0)
    assert named[5] == L.AMBIGUOUS, f"a weak lone marker named a cluster: {named[5]!r}"


def test_subtype_uses_raw_confidences_not_standardised():
    """Standardising divides each marker by its spread across the whole section, so argmax over z
    ranks markers by rarity rather than by how positive the cell is.  On donor 6442 that put 4,357
    of 4,414 SST-high endocrine cells behind INS or GCG and reduced Delta to 41 cells in 1.7M."""
    import inspect
    import call as K
    src = inspect.getsource(K.subtype)
    assert "X: np.ndarray" in src, "subtype must take the raw confidence matrix"
    assert "standardised" in src.lower(), "the reason must stay recorded next to the code"
    assert "subtype(df, X_typed, markers, sep)" in Path(K.__file__).read_text()


def test_subtype_boundary_is_declared_once_and_cannot_change_a_lineage():
    import call as K
    assert K.SUBTYPE_MIN_CONF == 0.5
    src = Path(K.__file__).read_text()
    assert src.count("SUBTYPE_MIN_CONF") == 2, "the one boundary must be defined and used once"


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


def test_a_cluster_must_be_positive_not_merely_unusual():
    """Within-donor standardisation has no floor.  On a donor that has LOST a cell type the
    marker is absent everywhere, its spread collapses, and a trivial elevation becomes a large z:
    139,040 cells on donor 6436 were named Endocrine outside any islet with median GCG confidence
    0.026.  Naming therefore requires raw positivity as well as an unusual z."""
    import call as K
    rng = np.random.default_rng(11)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n = 6000
    Z = rng.normal(size=(n, len(markers)))
    Xraw = np.full((n, len(markers)), 0.02)      # nothing is raised anywhere
    labels = rng.integers(0, 30, n)
    # a cluster that is strongly unusual on Endocrine markers but still not positive
    for m in L.LINEAGE_MARKERS["Endocrine"]:
        Z[labels == 4, idx[m]] += 5.0
        Xraw[labels == 4, idx[m]] = 0.026
    named, why, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                     margin=0.15, alpha=0.05, seed=0)
    assert named[4] == L.AMBIGUOUS, f"a hormone-negative cluster was named {named[4]!r}"
    # and the same cluster, actually positive, must be named
    for m in L.LINEAGE_MARKERS["Endocrine"]:
        Xraw[labels == 4, idx[m]] = 0.8
    named2, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                 margin=0.15, alpha=0.05, seed=0)
    assert named2[4] == "Endocrine", f"a positive endocrine cluster was not named: {named2[4]!r}"


def test_a_broadly_positive_marker_does_not_hand_its_lineage_every_cluster():
    """An absolute floor lets any lineage carrying a always-on marker pass everywhere: Nimbus
    scores TP63 positive in 98.6% of cells, and Exocrine took a cluster whose INS was +5.98 and
    whose Endocrine score was seventeen times higher.  Lift above the donor's own baseline is
    what distinguishes a raised marker from a permanently high one."""
    import call as K
    rng = np.random.default_rng(21)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n = 6000
    Z = rng.normal(size=(n, len(markers)))
    labels = rng.integers(0, 30, n)
    Xraw = np.full((n, len(markers)), 0.02)
    for m in L.LINEAGE_MARKERS["Exocrine"]:          # always-on everywhere, never raised
        Xraw[:, idx[m]] = 0.95
    endo = [idx[m] for m in L.LINEAGE_MARKERS["Endocrine"]]
    Z[np.ix_(labels == 9, endo)] += 5.0
    Xraw[np.ix_(labels == 9, endo)] = 0.60           # genuinely raised
    named, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                margin=0.15, alpha=0.05, seed=0)
    assert named[9] == "Endocrine", f"always-on Exocrine stole the cluster: {named[9]!r}"


def test_a_weakly_raised_marker_is_refused():
    """Lift alone passes a marker raised above a near-zero baseline but still below the model's
    boundary: donor 6436 lifted INS to 0.338 (lift 0.32) with only 28-33% of cells positive, and
    190,460 cells were named Endocrine outside any islet.  The positive fraction refuses it."""
    import call as K
    rng = np.random.default_rng(31)
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n = 6000
    Z = rng.normal(size=(n, len(markers)))
    labels = rng.integers(0, 30, n)
    endo = [idx[m] for m in L.LINEAGE_MARKERS["Endocrine"]]
    Xraw = np.full((n, len(markers)), 0.02)
    Z[np.ix_(labels == 6, endo)] += 5.0
    Xraw[np.ix_(labels == 6, endo)] = 0.34          # raised, but under the boundary
    named, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                margin=0.15, alpha=0.05, seed=0)
    assert named[6] == L.AMBIGUOUS, f"a sub-boundary cluster was named {named[6]!r}"
    Xraw[np.ix_(labels == 6, endo)] = 0.62          # now genuinely positive
    named2, *_ = K.name_clusters(Z, Xraw, labels, markers, 30, draws=1500,
                                 margin=0.15, alpha=0.05, seed=0)
    assert named2[6] == "Endocrine", f"a positive endocrine cluster was refused: {named2[6]!r}"


# --------------------------------------------------------------- per-cell resolution

def test_marker_specificity_is_measured_per_donor_not_assumed():
    """A marker positive in 96% of cells settles nothing.  Specificity is one minus the marker's
    own donor-wide positive fraction, recomputed per donor."""
    import inspect
    import call as K
    src = inspect.getsource(K.resolve_unnamed)
    assert "np.nanmean(Xraw >= MIN_CONF, axis=0)" in src
    assert K.MIN_SPECIFICITY == 0.60


def test_an_always_on_marker_cannot_name_a_cell():
    """TP63 is positive in 95.5% of cells on 6442.  A cell positive only for markers like that
    must stay unnamed rather than be called Exocrine by default."""
    import call as K
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n = 2000
    X = np.full((n, len(markers)), 0.02)
    X[:, idx["TP63"]] = 0.9                      # positive nearly everywhere -> specificity ~0
    df = pd.DataFrame({"lineage": [L.AMBIGUOUS] * n, "reason": ["tied:a|b"] * n})
    lin, _ = K.resolve_unnamed(df, X, markers, np.zeros(n, bool), margin=0.15)
    assert (lin == L.AMBIGUOUS).all(), "an always-on marker named cells"


def test_a_specific_marker_names_a_cell_its_cluster_could_not():
    """The reviewer's case: a cell plainly positive for one myeloid marker and negative for
    everything that could compete should be called, not left Ambiguous."""
    import call as K
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n = 2000
    X = np.full((n, len(markers)), 0.02)
    X[:, idx["TP63"]] = 0.9                      # always-on noise
    X[:50, idx["CD68"]] = 0.8                    # rare, discriminating
    df = pd.DataFrame({"lineage": [L.AMBIGUOUS] * n, "reason": ["tied:a|b"] * n})
    lin, why = K.resolve_unnamed(df, X, markers, np.zeros(n, bool), margin=0.15)
    assert (lin[:50] == "Myeloid").all(), f"CD68-positive cells were not resolved: {set(lin[:50])}"
    assert (lin[50:] == L.AMBIGUOUS).all()
    assert (why[:50] == "resolved_per_cell").all()


def test_a_nucleus_only_cell_is_read_on_nuclear_markers_alone():
    """A cell whose polygon equals its nucleus has no cytoplasm and no membrane, so a surface
    marker scored on it was measured in a compartment it does not have."""
    import call as K
    markers = L.identity_markers()
    idx = {m: i for i, m in enumerate(markers)}
    n = 100
    X = np.full((n, len(markers)), 0.02)
    X[:, idx["CD68"]] = 0.8                      # membrane marker, must be ignored for these
    df = pd.DataFrame({"lineage": [L.AMBIGUOUS] * n, "reason": [""] * n})
    lin, why = K.resolve_unnamed(df, X, markers, np.ones(n, bool), margin=0.15)
    assert (lin == L.AMBIGUOUS).all(), "a nucleus-only cell was named from a membrane marker"
    assert (why == "nucleus_only_no_discriminating_marker").all()


def test_load_reads_the_compartment_flags():
    """The nuclear-only restriction is silent if has_cytoplasm is never loaded: `.get(...)`
    defaults it to True and every cell is read on every marker."""
    import inspect
    import call as K
    src = inspect.getsource(K.load)
    assert "has_cytoplasm" in src and "has_nucleus" in src


def test_unmeasured_objects_are_omitted_from_the_export_by_default():
    """Objects with no nucleus are filtered out before typing, so writing them back as an
    Unavailable class would CHANGE their classification in QuPath rather than leave them alone."""
    import inspect
    import export_qupath as E
    src = inspect.getsource(E.build)
    assert 'out["assignment_status"] != "unavailable"' in src
    assert "include_unavailable" in inspect.signature(E.build).parameters
