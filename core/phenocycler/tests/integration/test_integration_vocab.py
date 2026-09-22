"""Vocabulary harmonisation: totality of the crosswalk, and the off-panel marker trap."""

from __future__ import annotations

import pytest

from phenocycler import load_config
from phenocycler.integration.phenospace import PhenotypeSpaceError
from phenocycler.integration.vocab import (
    COMMON_LINEAGES, GATE_MARKERS, NON_CALLS,
    PROTEIN_GENE_CANDIDATES, RESOLVABLE, SURROGATE_PANEL, XENIUM_TO_COMMON, check_panel_roles_agreement, coverage_report, crosswalk_frame,
    load_panel_taxonomy, pheno_to_common, protein_gene_pairs, resolvable_for, usable_anchors,
    xenium_to_common,
)

#: The Xenium pancreas broad vocabulary (scripts/tissue_markers/pancreas.py BROAD_TO_FINE).
XENIUM_BROAD = ("Endocrine", "Exocrine", "Stromal", "Vascular",
                "Immune_T", "Immune_B", "Myeloid", "Neural")

#: Markers verified absent from the deployed 5,081-gene panel.
KNOWN_OFF_PANEL = ("INS", "GCG", "SST", "Vimentin", "HLA_DR", "HLA_A",
                   "Iba1", "Caveolin", "M2Gal3")


@pytest.fixture(scope="module")
def cfg():
    return load_config()


@pytest.fixture(scope="module")
def tax(cfg):
    try:
        return load_panel_taxonomy(cfg)
    except Exception as exc:  # pragma: no cover - submodule not checked out
        pytest.skip(f"XeniumPanelExplorer submodule unavailable: {exc}")


# --------------------------------------------------------------------------- #
# Totality — nothing may fall through silently
# --------------------------------------------------------------------------- #

def test_every_declared_phenocycler_label_resolves(legacy_space):
    """The producer declaration, rather than core config, owns mapping totality."""
    for label in legacy_space.rules:
        common, resolution = legacy_space.common_for(label)
        if resolution == "mapped":
            assert common in COMMON_LINEAGES


def test_every_xenium_lineage_maps():
    for lineage in XENIUM_BROAD:
        assert lineage in XENIUM_TO_COMMON, lineage
        assert XENIUM_TO_COMMON[lineage] in COMMON_LINEAGES


def test_immune_splits_three_ways(legacy_space):
    """PhenoCycler collapses T/B/myeloid; the split is recovered from the per-marker gates."""
    assert pheno_to_common(legacy_space, "Immune", "T") == "T_NK"
    assert pheno_to_common(legacy_space, "Immune", "B") == "B_Plasma"
    assert pheno_to_common(legacy_space, "Immune", "Myeloid") == "Myeloid"
    # Multi-gate-positive cells go to the least-specific class rather than being resolved
    # by column order.
    assert pheno_to_common(legacy_space, "Immune", "Mixed") == "Myeloid"
    assert pheno_to_common(legacy_space, "Immune", "") == "Myeloid"


def test_fibroblast_and_muscle_collapse_to_stromal(legacy_space):
    """Xenium cannot split them: ACTA2 is not in the 5K reference and 12_pericyte
    contributes one SMC-hinted gene."""
    assert (pheno_to_common(legacy_space, "Fibroblast")
            == pheno_to_common(legacy_space, "Muscle") == "Stromal")
    assert XENIUM_TO_COMMON["Stromal"] == "Stromal"


def test_xenium_reaches_granulocyte_only_through_a_fine_label():
    """Xenium folds `07_immune_granulocyte` into `Myeloid` at broad level, whatever the
    protein side does — so a fine label is the only route."""
    assert RESOLVABLE["Granulocyte"]["xenium"] == "fine_only"
    assert xenium_to_common("Myeloid", "Mast") == "Granulocyte"
    assert xenium_to_common("Myeloid", "Myeloid") == "Myeloid"


def test_phenocycler_granulocyte_route_is_release_declared(legacy_space, v4_space):
    """Neither synthetic release declares a granulocyte call, so none is advertised."""
    with pytest.raises(PhenotypeSpaceError, match="undeclared"):
        pheno_to_common(legacy_space, "Neutrophil")
    assert resolvable_for(space=legacy_space)["Granulocyte"]["phenocycler"] == "no"
    assert resolvable_for(space=v4_space)["Granulocyte"]["phenocycler"] == "no"


def test_epithelial_is_a_default_sink_and_is_flagged(legacy_space):
    assert pheno_to_common(legacy_space, "Epithelial") == "Exocrine"
    assert RESOLVABLE["Exocrine"]["phenocycler"] == "default_sink"


def test_unknown_producer_labels_raise_instead_of_falling_through(legacy_space):
    with pytest.raises(PhenotypeSpaceError, match="undeclared"):
        pheno_to_common(legacy_space, "SomethingNew")
    assert xenium_to_common("SomethingNew") == "Unassigned"
    assert pheno_to_common(legacy_space, "") == "Unassigned"


def test_every_common_lineage_has_a_resolvability_entry():
    for lineage in COMMON_LINEAGES:
        assert lineage in RESOLVABLE, lineage
        assert set(RESOLVABLE[lineage]) == {"phenocycler", "xenium"}


# --------------------------------------------------------------------------- #
# Panel taxonomy (needs the pinned submodule)
# --------------------------------------------------------------------------- #

def test_panel_roles_drift_guard(tax):
    """Upstream renaming or adding a lineage must be caught, not silently Unassigned."""
    assert check_panel_roles_agreement(tax) == []


def test_identity_cores_loaded(tax):
    assert tax.panel_to_lineage, "no broad_mask panels found"
    for stem in tax.panel_to_lineage:
        assert tax.core_genes.get(stem), f"{stem} has no identity_core genes"
    # The endocrine core is TF/secretory machinery — the hormones are not on the panel.
    endo = set(tax.lineage_genes("Endocrine"))
    assert {"PDX1", "ISL1", "NEUROD1", "ABCC8"} <= endo
    assert not ({"INS", "GCG", "SST"} & endo)


def test_granulocyte_surrogate_genes_available(tax):
    """Granulocytes have no Xenium broad class; the hinted genes are the only route."""
    genes = tax.genes_by_hint("07_immune_granulocyte", "neutrophil")
    assert {"CSF3R", "CXCR1", "CXCR2", "FCGR3B", "MPO"} & set(genes) or genes
    assert len(genes) >= 5


# --------------------------------------------------------------------------- #
# Protein <-> gene crosswalk
# --------------------------------------------------------------------------- #

def test_pairs_are_validated_against_the_actual_panel(tax):
    """A gene absent from the run must not become a phantom anchor."""
    pairs = protein_gene_pairs(["PECAM1", "CD3E"], tax=tax)
    on = {p.protein for p in pairs if p.on_panel}
    assert on == {"CD31", "CD3e"}
    assert all(not p.on_panel for p in pairs if p.protein not in on)


def test_off_panel_markers_route_to_surrogate_panels(tax):
    """INS/GCG/SST/Vimentin define PhenoCycler's Endocrine and Fibroblast classes and are
    absent from the Xenium panel — so those comparisons must go through surrogate panels."""
    panel = sorted({g for gs in tax.core_genes.values() for g in gs})
    pairs = {p.protein: p for p in protein_gene_pairs(panel, tax=tax)}
    for marker in ("INS", "GCG", "SST", "Vimentin"):
        p = pairs[marker]
        assert not p.on_panel
        assert not p.usable_anchor, "an off-panel marker must never count as an anchor"
        assert p.surrogate_panel, f"{marker} has no surrogate panel"
        assert p.surrogate_genes, f"{marker} surrogate panel resolved to no genes"


def test_surrogate_genes_are_filtered_to_the_run_panel(tax):
    """A surrogate is only useful if its genes are actually present in the run.

    Filtering here rather than at use-time is what stops a surrogate panel from looking
    available when none of its genes were measured.
    """
    pairs = {p.protein: p for p in protein_gene_pairs(["PECAM1"], tax=tax)}
    assert pairs["INS"].surrogate_panel == "01_pancreas_endocrine"
    assert pairs["INS"].surrogate_genes == []


def test_known_off_panel_set(tax):
    """Guards the empirical finding against silent drift in the candidate table."""
    real_panel = sorted({g for gs in tax.core_genes.values() for g in gs})
    # Add the markers we know ARE on the deployed panel via the custom-100 add-on.
    real_panel += ["KRT19", "ACTA2", "MS4A1", "CD163", "CD68", "CD99", "TUBB3", "MPO",
                   "PECAM1", "CD3E", "CD3G"]
    pairs = {p.protein: p for p in protein_gene_pairs(real_panel, tax=tax)}
    for marker in KNOWN_OFF_PANEL:
        assert not pairs[marker].on_panel, f"{marker} unexpectedly resolved on-panel"


def test_gate_markers_are_flagged():
    for marker in GATE_MARKERS:
        assert marker in PROTEIN_GENE_CANDIDATES, marker
    pairs = protein_gene_pairs(["PECAM1"])
    assert {p.protein for p in pairs if p.is_gate} == set(GATE_MARKERS)


def test_usable_anchors_requires_a_real_gene(tax):
    pairs = protein_gene_pairs(["PECAM1", "MS4A1", "CD163"], tax=tax)
    anchors = usable_anchors(pairs)
    assert {a.protein for a in anchors} == {"CD31", "CD20", "CD163"}
    assert all(a.genes for a in anchors)


def test_type_only_anchors_exclude_state_markers(tax):
    pairs = protein_gene_pairs(["PECAM1", "MKI67"], tax=tax)
    assert {a.protein for a in usable_anchors(pairs)} == {"CD31", "Ki67"}
    assert {a.protein for a in usable_anchors(pairs, type_only=True)} == {"CD31"}


def test_surrogate_panels_reference_real_cores(tax):
    for marker, stem in SURROGATE_PANEL.items():
        if stem:
            assert stem in tax.core_genes, f"{marker} -> unknown identity core {stem}"


def test_crosswalk_frame_and_report(tax, legacy_space):
    pairs = protein_gene_pairs(["PECAM1", "CD3E", "MS4A1"], tax=tax)
    df = crosswalk_frame(legacy_space, pairs, tax=tax)
    assert set(df["kind"]) == {"lineage", "noncall", "marker"}
    assert len(df[df["kind"] == "lineage"]) == len(COMMON_LINEAGES)
    # One row per recorded non-call, outside the lineage block, and never a lineage's is-immune.
    assert len(df[df["kind"] == "noncall"]) == len(NON_CALLS)
    assert set(df.loc[df["kind"] == "noncall", "is_immune"]) == {"not_applicable"}
    report = coverage_report(legacy_space, pairs, tax)
    assert "OFF-PANEL LINEAGE GATES" in report
    assert "default_sink" in report


def test_swapped_space_changes_crosswalk_without_mutating_state(tax, legacy_space, v4_space):
    pairs = protein_gene_pairs(["PECAM1"], tax=tax)
    legacy = crosswalk_frame(legacy_space, pairs, tax=tax)
    v4 = crosswalk_frame(v4_space, pairs, tax=tax)
    legacy_exocrine = legacy.loc[legacy["common"] == "Exocrine", "phenocycler"].iloc[0]
    v4_exocrine = v4.loc[v4["common"] == "Exocrine", "phenocycler"].iloc[0]
    assert legacy_exocrine == "Epithelial"
    assert v4_exocrine == ""
    assert pheno_to_common(legacy_space, "Epithelial") == "Exocrine"


def test_crosswalk_covers_exactly_the_declared_labels(tax, legacy_space):
    pairs = protein_gene_pairs([], tax=tax)
    frame = crosswalk_frame(legacy_space, pairs, tax=tax)
    listed = set()
    for value in frame.loc[frame["kind"] == "lineage", "phenocycler"]:
        listed.update(part for part in str(value).split("|") if part)
    expected = {
        label for label in legacy_space.rules
        if legacy_space.common_for(label)[1] == "mapped"
    }
    assert listed == expected


def test_mpo_only_cells_reach_myeloid_in_legacy_declaration(legacy_space):
    """Broad types only: with MPO inside `Immune`, a neutrophil has no subclass and maps to
    Myeloid — which is exactly where Xenium puts a granulocyte at broad level."""
    assert pheno_to_common(legacy_space, "Immune", "") == "Myeloid"
    assert pheno_to_common(legacy_space, "Immune", "Myeloid") == "Myeloid"


def test_granulocyte_never_enters_a_comparable_set():
    """Xenium cannot call it at broad level, so it is not evidence either way."""
    from phenocycler.integration.donor import COMPARABLE, comparable_for
    from phenocycler.integration.tissues import TISSUES

    assert "Granulocyte" not in COMPARABLE
    for tissue in TISSUES:
        assert "Granulocyte" not in comparable_for(tissue), tissue


# --------------------------------------------------------------------------- #
# Is-immune — recorded, total, and scoped to the tissues the ruling covers
# --------------------------------------------------------------------------- #

def test_immune_set_is_exactly_the_four_immune_lineages():
    from phenocycler.integration.vocab import IMMUNE_LINEAGES, immune_status

    assert set(IMMUNE_LINEAGES) == {"T_NK", "B_Plasma", "Myeloid", "Granulocyte"}
    assert {lin for lin in COMMON_LINEAGES if immune_status(lin) == "yes"} == set(IMMUNE_LINEAGES)


def test_every_common_lineage_is_classified_and_unassigned_is_neither():
    """Totality, and the one bucket that is not answerable is not answered."""
    from phenocycler.integration.vocab import immune_status

    values = {lin: immune_status(lin) for lin in COMMON_LINEAGES}
    assert set(values.values()) <= {"yes", "no", "unresolvable"}
    assert values["Unassigned"] == "unresolvable"
    assert {lin for lin, v in values.items() if v == "no"} == {
        "Endocrine", "Exocrine", "Stromal", "Vascular", "Neural"}


def test_an_unclassified_lineage_raises_instead_of_defaulting_to_not_immune():
    from phenocycler.integration.vocab import immune_status

    with pytest.raises(KeyError, match="no recorded is-immune classification"):
        immune_status("Tuft")


def test_immune_scope_covers_pancreas_and_lymph_node_and_cites_the_ruling():
    """Maintainer ruling 2026-09-17: immune calls from pancreas apply to pLN."""
    from phenocycler.integration.vocab import (
        IMMUNE_TISSUE_SCOPE, immune_note, immune_tissue_scope,
    )

    assert set(IMMUNE_TISSUE_SCOPE) == {"pancreas", "pancreatic_lymph_node"}
    assert immune_tissue_scope("yes") == "pancreas|pancreatic_lymph_node"
    assert immune_tissue_scope("no") == ""
    note = immune_note("yes")
    assert "2026-09-17" in note and "pancreatic_lymph_node" in note
    assert immune_note("unresolvable").startswith("no lineage call")


def test_crosswalk_frame_records_is_immune_per_lineage(tax, legacy_space):
    pairs = protein_gene_pairs(["PECAM1", "CD3E", "MS4A1"], tax=tax)
    df = crosswalk_frame(legacy_space, pairs, tax=tax)
    lineages = df[df["kind"] == "lineage"].set_index("common")
    assert set(lineages.loc[lineages["is_immune"] == "yes"].index) == {
        "T_NK", "B_Plasma", "Myeloid", "Granulocyte"}
    assert lineages.loc["T_NK", "immune_tissue_scope"] == "pancreas|pancreatic_lymph_node"
    assert lineages.loc["Endocrine", "immune_tissue_scope"] == ""
    assert lineages.loc["Unassigned", "is_immune"] == "unresolvable"
    # A marker row is not a lineage, and must not read as "not immune".
    assert set(df.loc[df["kind"] == "marker", "is_immune"]) == {"not_applicable"}


def test_shipped_crosswalk_csv_records_the_immune_mapping(cfg):
    """The file a builder actually reads, not just the in-memory frame."""
    import csv as _csv

    path = cfg.vocab_crosswalk_csv
    if not path.exists():  # pragma: no cover - crosswalk not built in this checkout
        pytest.skip(f"{path} has not been built")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [r for r in _csv.DictReader(handle) if r["kind"] == "lineage"]
    assert rows, f"{path} carries no lineage rows"
    immune = {r["common"] for r in rows if r["is_immune"] == "yes"}
    assert immune == {"T_NK", "B_Plasma", "Myeloid", "Granulocyte"}
    scopes = {r["immune_tissue_scope"] for r in rows if r["is_immune"] == "yes"}
    assert scopes == {"pancreas|pancreatic_lymph_node"}
    assert all("2026-09-17" in r["immune_note"] for r in rows if r["is_immune"] == "yes")
