"""Thymus marker definitions for the Xenium pipeline.

Hand-curated against the panel-availability data in
`data/processed/panel_audit_thymus.csv`. Each scored panel has ≥3
panel-present genes (CLAUDE.md decision #12 — no single-gene panels).
Cell-type vocabulary follows Park et al. 2020 Science thymic atlas /
Bautista 2021 Cell Reports (AIRE+ mTEC) / standard immunology nomenclature.

NOTABLE PANEL ABSENCES (do NOT promote into a panel — they're absent from
the 5K hMulti_100g panel for these two thymus samples):
  CD3D, COL1A1, COL3A1, DCN, LUM, MYH11, TAGLN, CNN1, VWF, GNLY,
  IVL, LOR, FLG, KRT5, KRT10, KRT15, HLA-DRA, C1QA, C1QB, LYZ,
  TPSAB1, CPA3, TPSB2, S100B, PLP1, TYROBP.

Substitutions used in this module are documented inline where they
deviate from "canonical" marker textbook expectations. Verifier:
`scripts/tissue_markers/verify_panel_coverage.py --tissue thymus`.
"""
import numpy as np
import pandas as pd

from ._shared import SHARED_MARKERS

# -----------------------------------------------------------------------------
# Score-genes argmax panels (cell 10 of 02_phenotyping.ipynb)
# -----------------------------------------------------------------------------
# Notes on substitutions:
#   - CD3D absent → use CD3E (33%/25%) + CD3G + TRAC + TRBC2 throughout.
#   - COL1A1/COL3A1/DCN/LUM absent → fibroblast uses COL5A1/COL5A2/
#     COL6A1 (10x stronger detection than COL1) + CDH11 + THY1 + POSTN
#     + FAP + PDGFRA.
#   - MYH11/TAGLN/CNN1 absent → smooth-muscle uses ACTA2 + MCAM +
#     PDGFRB + RGS5. (RGS5 is also pericyte, intentional.)
#   - VWF absent → endothelial uses PECAM1 + CDH5 + CLDN5 + PLVAP + KDR
#     + ENG + CD34 + TIE1 + ERG + CD93.
#   - HLA-DRA absent → DC uses ITGAX + CD1C + CD83 + BATF3 + IRF8 +
#     FLT3 + LAMP3 + CLEC9A + CLEC10A + CD1A.
#   - C1QA/B/LYZ absent → myeloid uses CD68 + CD163 + CSF1R + MARCO +
#     CD14 + ITGAM + CD86 + MRC1 + STAB1.
#   - GNLY absent → NK uses NKG7 + PRF1 + GZMA + KLRD1 + KLRB1 +
#     FCGR3A + NCR1 + KLRF1 + EOMES + TBX21.
#   - TPSAB1/CPA3/TPSB2 absent → mast restricted to KIT + CMA1 + HDC +
#     MS4A2 + FCER1A + HRH4 (all <0.5% detection — mast genuinely rare;
#     panel-coverage gate may downgrade most calls).
#   - S100B/PLP1 absent → Schwann uses MPZ + PMP22 + NES + NGFR +
#     S100A1.
#   - KRT5/CLDN3/CLDN4 absent → mTEC uses KRT14 + AIRE + FEZF2 + CD83 +
#     DSG3 (DSG3 = desmoglein-3, mTEC-medullary).
#   - IVL/LOR/FLG/KRT10 absent → Hassall's corpuscles cannot be scored
#     via argmax; identified post-hoc via apply_rare_celltype_rules()
#     using KRT1 + TGM1/TGM3 + DSG1.

# Non-overlapping celltype panels — every gene appears in only one
# panel. Sourced from XeniumPanelExplorer thymus subpanel CSVs
# (https://github.com/smith6jt-cop/XeniumPanelExplorer) crossed with
# detection ≥0.5% in at least one thymus sample (relaxed to ≥0.1%
# for the rare types Mast / Hassalls / Lymphatic_EC / Schwann).
# Detection floors come from data/processed/panel_audit_thymus.csv.
#
# Thymocyte SUBTYPE final calls come from apply_thymocyte_subtyping()
# (per-cell CD4/CD8/CD3/FOXP3 direct detection), independent of
# MARKER_DICT — so SP4/SP8/DN/DP score panels can be fully
# non-overlapping without breaking subtyping.
MARKER_DICT = {
    # Thymocyte developmental stages — no shared canonical markers
    "DN_thymocyte":  ["CD34", "NOTCH1", "HES1", "IL7R", "BCL11B", "CD7"],
    # Stem/progenitor + Notch-pathway anchors. CD34 owned here.
    "DP_thymocyte":  ["CD8B", "CCR9", "LEF1", "TCF7"],
    # DP-specific TFs (LEF1/TCF7 peak in DP); CD8B owned here.
    # CD4 / CD8A NOT in this panel — assigned to SP4 / SP8 respectively.
    "SP4_thymocyte": ["CD4", "GATA3", "CCR4", "CD40LG", "ICOS",
                       "KLF2", "CCR7", "FOXP1"],
    # CD4 owned here; mature helper-T markers.
    "SP8_thymocyte": ["CD8A", "RUNX3", "EOMES", "CRTAM"],
    # CD8A owned here; CD8-lineage TF + late effector.
    # NKG7/PRF1/GZMA moved out to NK (more canonical for NK biology).
    "Treg":          ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT",
                       "TNFRSF18"],
    # CCR4 moved to SP4 (helper-T context).

    # Thymic epithelial cells
    "cTEC":          ["PSMB11", "CTSL", "PRSS16", "CCL25", "LY75",
                       "FOXN1", "CXCL12", "EPCAM"],
    "mTEC":          ["KRT14", "AIRE", "FEZF2", "DSG3"],
    # CD83 moved to DC (DC-activation marker, more specific there).
    # Hassall's + AIRE_mTEC + pDC handled post-hoc via
    # apply_rare_celltype_rules() — see below.

    # Other immune
    "B_cells":       ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5"],
    "Plasma_cells":  list(SHARED_MARKERS["Plasma_cells"]),  # MZB1,TNFRSF17,PRDM1 (W2)
    "Myeloid":       ["CD68", "CD163", "CSF1R", "MARCO", "CD14",
                       "ITGAM", "MRC1", "STAB1", "CD86"],
    "DC":            ["ITGAX", "CD1C", "CLEC9A", "CLEC10A", "BATF3",
                       "IRF8", "FLT3", "LAMP3", "FSCN1", "IRF4", "CD83"],
    # CD1A dropped (0.12% — too low); + FSCN1, IRF4, CD83 (DC-activation).
    "NK":            ["KLRK1", "NKG7", "PRF1", "KLRB1", "KLRF1",
                       "NCR3", "KLRD1", "GZMA", "FCGR3A"],
    # NKG7 / PRF1 / GZMA moved IN from SP8 (canonical NK cytotoxic).
    # NCR1 / EOMES / TBX21 dropped (very low detection or moved to SP8).
    "Mast":          ["KIT", "CMA1", "HDC", "MS4A2", "FCER1A",
                       "ENPP3", "HRH4"],
    # ENPP3 added (mast/basophil-specific; 5.3% in 0041435).

    # Stromal / vascular / neural
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR",
                       "ENG", "TIE1", "ERG", "CD93"],
    # CD34 moved to DN_thymocyte (stem/progenitor anchor); 9 strong
    # endothelial markers remain.
    "Lymphatic_EC":  ["LYVE1", "PROX1", "PDPN", "CCL21", "FLT4"],
    "Fibroblast":    ["COL5A1", "COL5A2", "COL6A1", "CDH11", "THY1",
                       "FBLN2", "POSTN", "FAP", "PDGFRA"],
    "Smooth_muscle": ["ACTA2", "MCAM", "PDGFRB", "RGS5"],
    "Schwann":       ["MPZ", "PMP22", "NES", "NGFR", "S100A1",
                       "GLDN", "PRPH"],
    # + GLDN, PRPH from XPE shared neural_glial subpanel.

    # Generic
    "Proliferating": ["MKI67", "TOP2A", "PCNA", "BIRC5", "MCM5",
                       "MCM3", "CDC20"],
}

# Lineage groupings used by notebook 02 cell 15 to compute the
# lineage-margin form of celltype_confidence. The thymocyte subtypes
# share the bulk of their markers (CD3E/CD3G/TRAC/TRBC2) so the old
# top-2-margin metric reported ~zero confidence for every clear
# thymocyte. Aggregating to lineages (Thymocyte vs everything else)
# restores the meaningful margin: high for any clear thymocyte; small
# only at lineage-boundary cells. The subtype call itself is preserved
# in predicted_celltype (full all-types argmax) and the within-lineage
# margin is preserved as subtype_confidence.
LINEAGE_GROUPS = {
    "Thymocyte":  {"DN_thymocyte", "DP_thymocyte", "SP4_thymocyte",
                   "SP8_thymocyte", "Treg"},
    "TEC":        {"cTEC", "mTEC"},
    "B_lymphoid": {"B_cells", "Plasma_cells"},
    "NK_lineage": {"NK"},
    "Myeloid":    {"Myeloid", "DC", "Mast"},
    "Vascular":   {"Endothelial", "Lymphatic_EC"},
    "Stromal":    {"Fibroblast", "Smooth_muscle", "Schwann"},
    "Cycling":    {"Proliferating"},
}

# Low-confidence types — rare and/or thin panel detection.
LOW_CONFIDENCE_TYPES = {
    "Mast", "Lymphatic_EC", "Schwann", "Treg", "Smooth_muscle",
}

INDETERMINATE_BUCKETS = {
    "Treg":          "Indeterminate_thymocyte",
    "Mast":          "Indeterminate_immune",
    "Lymphatic_EC":  "Indeterminate_vascular",
    "Schwann":       "Indeterminate_neural",
    "Smooth_muscle": "Indeterminate_stromal",
}

# -----------------------------------------------------------------------------
# Anti-DP_thymocyte override (analog of pancreas anti-acinar)
#   DP thymocytes dominate pediatric thymus (~50–70% of cells).
#   Cells argmax'd as DP but expressing strong markers of a non-thymocyte
#   type get reassigned.
# -----------------------------------------------------------------------------

# Unset (2026-05-29): thymocyte subtyping is now the Tier-2 detection-gated argmax
# (DETECTION_SUBTYPE_LINEAGES['Thymocyte'] below), NOT a dominant-type default. A
# non-None value would make the Thymocyte lineage take tier2_fine's dominant branch
# (which runs BEFORE detect_cfg), silently ignoring the gated argmax.
ANTI_OVERRIDE_DOMINANT_TYPE = None

# When True, the legacy per-cell boolean apply_thymocyte_subtyping() (CD4/CD8/CD3/
# FOXP3 co-expression) overwrites the Tier-2 thymocyte call inside
# apply_rare_celltype_rules. Set False (2026-05-29, user choice): the Tier-2
# detection-gated {DN,DP,SP4,SP8,Treg}_score argmax IS the call. Flip True to run
# the gated-vs-boolean side-by-side (W5 verification).
THYMOCYTE_BOOLEAN_SUBTYPING = False

# Retained for the legacy notebook anti-override path (unused in the hierarchical
# path now that ANTI_OVERRIDE_DOMINANT_TYPE is None).
_LEGACY_ANTI_OVERRIDE_DOMINANT_TYPE = "DP_thymocyte"

ANTI_OVERRIDE_MARKERS = {
    "cTEC":          ["PSMB11", "CCL25", "PRSS16", "FOXN1"],
    "mTEC":          ["KRT14", "AIRE", "FEZF2"],
    "Myeloid":       ["CD68", "CD163", "CSF1R", "MARCO"],
    "DC":            ["ITGAX", "CLEC9A", "CLEC10A"],
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP"],
    "Fibroblast":    ["COL5A1", "COL6A1", "CDH11", "THY1"],
    "Smooth_muscle": ["ACTA2", "MCAM", "PDGFRB"],
    "B_cells":       ["MS4A1", "CD79A", "CD19"],
    "NK":            ["NKG7", "PRF1", "KLRD1"],
    "Mast":          ["KIT", "CMA1"],
}

# -----------------------------------------------------------------------------
# Lineage canary panels (cell 22 + scripts/verify_upper_tail_doublets.py +
# scripts/lineage_phenotype.py)
# -----------------------------------------------------------------------------

LINEAGE_CANARY_PANELS = {
    "DP_thymocyte":  ["CD4", "CD8A", "CD3E", "TRAC"],
    "SP_thymocyte":  ["CD3E", "TRAC", "TRBC2", "CCR7"],
    "DN_thymocyte":  ["CD7", "CD44", "CD34"],
    "cTEC":          ["PSMB11", "CCL25", "PRSS16", "FOXN1"],
    "mTEC":          ["KRT14", "AIRE", "FEZF2"],
    "Hassalls":      ["KRT1", "TGM1", "DSG1"],
    "B_cells":       ["MS4A1", "CD79A", "CD19"],
    "Myeloid":       ["CD68", "CD163", "MARCO"],
    "DC":            ["ITGAX", "CLEC9A"],
    "pDC":           ["LILRA4", "IRF7"],
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5"],
    "Fibroblast":    ["COL5A1", "COL6A1", "CDH11"],
    "NK":            ["NKG7", "KLRD1", "PRF1"],
    "Mast":          ["KIT", "CMA1", "HDC"],
}

LINEAGE_EXCLUSIVE_GROUP = {
    "DP_thymocyte":  "Thymocyte",
    "SP_thymocyte":  "Thymocyte",
    "DN_thymocyte":  "Thymocyte",
    "cTEC":          "TEC",
    "mTEC":          "TEC",
    "Hassalls":      "TEC",
    "B_cells":       "Lymphoid_non_T",
    "NK":            "Lymphoid_non_T",
    "Myeloid":       "Myeloid",
    "DC":            "Myeloid",
    "pDC":           "Myeloid",
    "Endothelial":   "Endothelial",
    "Fibroblast":    "Stromal",
    "Smooth_muscle": "Stromal",
    "Mast":          "Mast",
}

# -----------------------------------------------------------------------------
# Ligand-receptor lineage map (cell 14 of 03_spatial_analysis.ipynb)
# -----------------------------------------------------------------------------

LINEAGE_MAP = {
    "DP_thymocyte":  "Thymocyte",
    "SP4_thymocyte": "Thymocyte",
    "SP8_thymocyte": "Thymocyte",
    "DN_thymocyte":  "Thymocyte",
    "Treg":          "Thymocyte",
    "cTEC":          "TEC",
    "mTEC":          "TEC",
    "Hassalls":      "TEC",
    "AIRE_mTEC":     "TEC",
    "B_cells":       "Other_Immune",
    "Plasma_cells":  "Other_Immune",
    "NK":            "Other_Immune",
    "Mast":          "Other_Immune",
    "Myeloid":       "Myeloid",
    "DC":            "Myeloid",
    "pDC":           "Myeloid",
    "Endothelial":   "Vascular",
    "Lymphatic_EC":  "Vascular",
    "Fibroblast":    "Stromal",
    "Smooth_muscle": "Stromal",
    "Schwann":       "Neural",
    "Proliferating": "Cycling",
}

# -----------------------------------------------------------------------------
# Consensus-annotation 3-panel plot dispatch (notebook 02 cell 20).
# Thymus vocab: specific celltypes shown in panel 2; panel 3 drops the
# bulk cortical thymocytes (DP + DN) so TECs, mature SP4/SP8/Treg,
# AIRE_mTEC, Hassalls, and rare immune/stromal types are visible.
# -----------------------------------------------------------------------------

SPECIFIC_KEEP = {
    # Thymocyte subtypes
    "DP_thymocyte", "SP4_thymocyte", "SP8_thymocyte", "DN_thymocyte", "Treg",
    # TEC subtypes
    "cTEC", "mTEC", "AIRE_mTEC", "Hassalls",
    # Immune
    "B_cells", "Plasma_cells", "NK", "Mast", "Myeloid", "DC", "pDC",
    # Stromal / vascular / neural / cycling
    "Endothelial", "Lymphatic_EC",
    "Fibroblast", "Smooth_muscle",
    "Schwann", "Proliferating",
}

PALETTE = {
    # Thymocytes
    "DP_thymocyte":  "#ff7f0e",
    "SP4_thymocyte": "#d62728",
    "SP8_thymocyte": "#bcbd22",
    "DN_thymocyte":  "#8c564b",
    "Treg":          "#17becf",
    "Thymocyte":     "#ff7f0e",
    # TECs
    "cTEC":          "#1f77b4",
    "mTEC":          "#9467bd",
    "AIRE_mTEC":     "#e83e8c",
    "Hassalls":      "#fdd0a2",
    "TEC":           "#1f77b4",
    # Immune
    "B_cells":       "#87ceeb",
    "Plasma_cells":  "#5e3c99",
    "NK":            "#1f3a93",
    "Mast":          "#fa8072",
    "Other_Immune":  "#87ceeb",
    "Myeloid":       "#4682b4",
    "DC":            "#6a3d9a",
    "pDC":           "#7fbf3f",
    # Stromal / vascular / neural / cycling
    "Endothelial":   "#e377c2",
    "Lymphatic_EC":  "#f7b6d2",
    "Vascular":      "#e377c2",
    "Fibroblast":    "#2ca02c",
    "Smooth_muscle": "#aec7e8",
    "Stromal":       "#2ca02c",
    "Schwann":       "#808000",
    "Neural":        "#808000",
    "Proliferating": "#7f7f7f",
    "Cycling":       "#7f7f7f",
}

THIRD_PANEL = {
    "exclude_types": {"DP_thymocyte", "DN_thymocyte"},
    "title":         "Specific celltype, no bulk thymocyte (DP/DN)",
}

# -----------------------------------------------------------------------------
# Lineage-subset panels for predicted_celltypes.png / scvi_vs_markers.png.
# Each subset is rendered as its own UMAP with high-contrast colors and a
# light-grey background of out-of-subset cells (scanpy's `groups=` mechanism).
# Lineage names refer to LINEAGE_MAP values above.
# -----------------------------------------------------------------------------

LINEAGE_SUBSETS = [
    {"title": "Thymocytes",   "lineages": ["Thymocyte"]},
    {"title": "TEC subtypes", "lineages": ["TEC"]},
    {"title": "Immune / Myeloid / Stromal / Vascular / Neural / Cycling",
     "lineages": ["Other_Immune", "Myeloid", "Stromal", "Vascular",
                  "Neural", "Cycling"]},
]

# -----------------------------------------------------------------------------
# Panel-coverage panels for the refine_rare_celltypes panel-coverage gate
# -----------------------------------------------------------------------------

# Mirror of MARKER_DICT — kept in sync so the panel-coverage gate in
# refine_rare_celltypes.py evaluates the same gene set that score_genes
# scored.
PANEL_COVERAGE_PANELS = {
    "DN_thymocyte":  ["CD34", "NOTCH1", "HES1", "IL7R", "BCL11B", "CD7"],
    "DP_thymocyte":  ["CD8B", "CCR9", "LEF1", "TCF7"],
    "SP4_thymocyte": ["CD4", "GATA3", "CCR4", "CD40LG", "ICOS",
                       "KLF2", "CCR7", "FOXP1"],
    "SP8_thymocyte": ["CD8A", "RUNX3", "EOMES", "CRTAM"],
    "Treg":          ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT",
                       "TNFRSF18"],
    "cTEC":          ["PSMB11", "CTSL", "PRSS16", "CCL25", "LY75",
                       "FOXN1", "CXCL12", "EPCAM"],
    "mTEC":          ["KRT14", "AIRE", "FEZF2", "DSG3"],
    "B_cells":       ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5"],
    "Myeloid":       ["CD68", "CD163", "CSF1R", "MARCO", "CD14",
                       "ITGAM", "MRC1", "STAB1", "CD86"],
    "DC":            ["ITGAX", "CD1C", "CLEC9A", "CLEC10A", "BATF3",
                       "IRF8", "FLT3", "LAMP3", "FSCN1", "IRF4", "CD83"],
    "NK":            ["KLRK1", "NKG7", "PRF1", "KLRB1", "KLRF1",
                       "NCR3", "KLRD1", "GZMA", "FCGR3A"],
    "Mast":          ["KIT", "CMA1", "HDC", "MS4A2", "FCER1A",
                       "ENPP3", "HRH4"],
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR",
                       "ENG", "TIE1", "ERG", "CD93"],
    "Lymphatic_EC":  ["LYVE1", "PROX1", "PDPN", "CCL21", "FLT4"],
    "Fibroblast":    ["COL5A1", "COL5A2", "COL6A1", "CDH11", "THY1",
                       "FBLN2", "POSTN", "FAP", "PDGFRA"],
    "Smooth_muscle": ["ACTA2", "MCAM", "PDGFRB", "RGS5"],
    "Schwann":       ["MPZ", "PMP22", "NES", "NGFR", "S100A1",
                       "GLDN", "PRPH"],
}

COMMON_CLUSTERING_TYPES = {
    "DP_thymocyte", "SP4_thymocyte", "SP8_thymocyte", "DN_thymocyte",
    "cTEC", "mTEC", "Myeloid", "Fibroblast", "Endothelial",
}

# Map celltype labels (as written in adata.obs['celltype']) to canary
# lineage names for verify_upper_tail_doublets.
LABEL_TO_CANARY = {
    "DN_thymocyte":  "DN_thymocyte",
    "DP_thymocyte":  "DP_thymocyte",
    "SP4_thymocyte": "SP_thymocyte",
    "SP8_thymocyte": "SP_thymocyte",
    "Treg":          "SP_thymocyte",
    "cTEC":          "cTEC",
    "mTEC":          "mTEC",
    "AIRE_mTEC":     "mTEC",
    "Hassalls":      "Hassalls",
    "B_cells":       "B_cells",
    "Myeloid":       "Myeloid",
    "DC":            "DC",
    "pDC":           "pDC",
    "NK":            "NK",
    "Mast":          "Mast",
    "Endothelial":   "Endothelial",
    "Fibroblast":    "Fibroblast",
}

# Candidate host cell-types for rule-based rare-type refinement.
_TEC_LIKE = {"mTEC", "cTEC",
             "Indeterminate_mTEC", "Indeterminate_cTEC",
             "Indeterminate_TEC"}
_DC_LIKE = {"DC", "Myeloid",
             "Indeterminate_DC", "Indeterminate_Myeloid",
             "Indeterminate_immune"}

# -----------------------------------------------------------------------------
# Rare cell-type rules (AIRE+ mTEC, Hassall's corpuscle, pDC)
# -----------------------------------------------------------------------------

def _get_lognorm_col(adata, gene):
    if gene not in adata.var_names:
        return None
    i = adata.var_names.get_loc(gene)
    L = adata.layers.get("lognorm", adata.X)
    col = L[:, i]
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def _vec(adata, gene):
    """Lognorm vector for `gene`, or zeros if not on the panel."""
    v = _get_lognorm_col(adata, gene)
    return v if v is not None else np.zeros(adata.n_obs)


_THYMOCYTE_LIKE = {
    "DN_thymocyte", "DP_thymocyte", "SP4_thymocyte", "SP8_thymocyte", "Treg",
    "Indeterminate_thymocyte",
    # cells routed here by margin / coverage gates
    "Indeterminate_DN_thymocyte",
    "Indeterminate_DP_thymocyte",
    "Indeterminate_SP4_thymocyte",
    "Indeterminate_SP8_thymocyte",
    "Indeterminate_DN_thymocyte_lowcoverage",
    "Indeterminate_DP_thymocyte_lowcoverage",
    "Indeterminate_SP4_thymocyte_lowcoverage",
    "Indeterminate_SP8_thymocyte_lowcoverage",
    "Indeterminate_Treg",
    "Indeterminate_Treg_lowcoverage",
}


def apply_thymocyte_subtyping(adata, log_detect, log_high, log_present=None):
    """Re-classify thymocyte-like cells into DN / DP / SP4 / SP8 / Treg
    using direct CD4 / CD8A / CD3E / FOXP3 per-cell detection.

    Why this exists: score_genes argmax with overlapping multi-marker
    panels (DP panel includes CD4 and CD8A, which also seed the SP4 and
    SP8 panels) underweights the canonical DP definition (CD4+ AND CD8+).
    Combined with the cluster-consensus step in cell 20 of
    02_phenotyping.ipynb, true DP cells get cluster-collapsed into the
    SP-dominant majority. This rule operates ONLY on cells already
    identified as thymocyte-family (via _THYMOCYTE_LIKE), so it does
    not hijack genuine non-thymocyte cells with stray CD signal.

    Detection threshold: `log_present` (default log1p(1) ≈ 0.69 = raw
    count ≥ 1). In FFPE Xenium any non-zero RNA count for CD4/CD8/CD3
    is biologically meaningful (the genes are not housekeeping; positive
    cells have raw count ≥ 1 with high specificity). The earlier
    `log_detect` (log1p(3) ≈ 1.39) bar was inherited from pancreas-tuned
    lineage thresholds calibrated for high-expression markers (AMY1A,
    INS) where raw counts run in the tens.

    Rule per cell:
      Treg: FOXP3+ AND CD4+
      DP:   CD4+ AND CD8A+
      SP4:  CD4+ AND CD8A- AND CD3E+
      SP8:  CD8A+ AND CD4- AND CD3E+
      DN:   CD4- AND CD8A- AND CD3E+ AND (CD7+ OR CD44+)
      else: leave label as-is (probably Indeterminate downstream).
    """
    if log_present is None:
        log_present = float(np.log1p(1.0))
    current = adata.obs["celltype"].astype(str).values
    host = pd.Series(current).isin(_THYMOCYTE_LIKE).values
    if host.sum() == 0:
        return {}, current

    cd4    = _vec(adata, "CD4")
    cd8a   = _vec(adata, "CD8A")
    cd8b   = _vec(adata, "CD8B")
    cd3e   = _vec(adata, "CD3E")
    foxp3  = _vec(adata, "FOXP3")
    cd7    = _vec(adata, "CD7")
    cd44   = _vec(adata, "CD44")

    cd4p  = cd4   >= log_present
    cd8p  = (cd8a >= log_present) | (cd8b >= log_present)
    cd3p  = cd3e  >= log_present
    foxp3p = foxp3 >= log_present

    new_lab = current.copy()
    # Treg first (highest priority within thymocyte family)
    is_treg = host & foxp3p & cd4p
    new_lab[is_treg] = "Treg"
    # DP next (canonical CD4+CD8+ co-expression)
    is_dp = host & cd4p & cd8p & ~is_treg
    new_lab[is_dp] = "DP_thymocyte"
    # SP4 (CD4+ exclusive)
    is_sp4 = host & cd4p & ~cd8p & cd3p & ~is_treg & ~is_dp
    new_lab[is_sp4] = "SP4_thymocyte"
    # SP8 (CD8+ exclusive)
    is_sp8 = host & cd8p & ~cd4p & cd3p & ~is_treg & ~is_dp
    new_lab[is_sp8] = "SP8_thymocyte"
    # DN (CD4-CD8- but T-lineage)
    is_dn = host & ~cd4p & ~cd8p & cd3p & ((cd7 >= log_detect) | (cd44 >= log_detect))
    new_lab[is_dn] = "DN_thymocyte"

    counts = {
        "Treg":          int(is_treg.sum()),
        "DP_thymocyte":  int(is_dp.sum()),
        "SP4_thymocyte": int(is_sp4.sum()),
        "SP8_thymocyte": int(is_sp8.sum()),
        "DN_thymocyte":  int(is_dn.sum()),
    }
    return counts, new_lab


def apply_rare_celltype_rules(adata, log_detect, log_high, log_present=None):
    """Multi-criteria rare-cell-type identification for thymus.

    Step 0 — thymocyte subtype reclassification (DN/DP/SP4/SP8/Treg)
             based on direct CD4/CD8/CD3/FOXP3 per-cell detection.
             Fixes the cluster-consensus collapse of DP cells into
             SP-dominant majority that affects score_genes argmax in
             FFPE Xenium thymus data.

    Step 1 — rare-cell rules (same multi-marker-co-expression principle
    as pancreas Gamma/Epsilon):

      AIRE_mTEC: AIRE ≥ LOG_HIGH AND KRT14 ≥ LOG_DETECT
                 AND no strong thymocyte (CD3E/CD4/CD8A) signal
                 AND currently TEC-like.
      Hassalls:  KRT1 ≥ LOG_HIGH AND (TGM1 | TGM3 | DSG1) ≥ LOG_DETECT
                 AND no strong cTEC (PSMB11/CCL25) signal
                 AND currently TEC-like.
      pDC:       LILRA4 ≥ LOG_HIGH AND ≥ 2 of (IRF7, IRF8, TCF4,
                 CLEC4C) ≥ LOG_DETECT
                 AND CD3E < LOG_DETECT AND currently DC/myeloid-like.

    If the obs column `predicted_celltype` is present (left by stage 02
    cell 9), the thymocyte subtyping uses celltype = predicted_celltype
    as input rather than the cluster-consensus celltype, so SP4/DP/DN
    diversity that the cluster step collapsed is recovered.
    """
    # Bypass cluster-consensus collapse: prefer per-cell predicted_celltype
    # if available (carries the original score_genes argmax + anti-X-override
    # diversity, which the cluster-consensus step in nb 02 cell 20
    # flattened into a per-cluster majority).
    if log_present is None:
        log_present = float(np.log1p(1.0))

    if "predicted_celltype" in adata.obs.columns:
        adata.obs["celltype"] = adata.obs["predicted_celltype"].astype(str).values
        print("    (using predicted_celltype as base for refinement; "
              "cluster-consensus celltype was collapsed)")

    if THYMOCYTE_BOOLEAN_SUBTYPING:
        thy_counts, post_thy_labels = apply_thymocyte_subtyping(
            adata, log_detect, log_high, log_present=log_present
        )
        if thy_counts:
            adata.obs["celltype"] = pd.Categorical(post_thy_labels)
            if "celltype_lineage" in adata.obs.columns:
                adata.obs["celltype_lineage"] = pd.Categorical(post_thy_labels)
            print("    thymocyte subtyping (boolean CD4/CD8/CD3/FOXP3):")
            for k, v in thy_counts.items():
                print(f"      {k}: {v:,}")
    else:
        print("    thymocyte subtyping: Tier-2 detection-gated argmax retained "
              "(boolean refiner OFF; THYMOCYTE_BOOLEAN_SUBTYPING=False)")

    current = adata.obs["celltype"].astype(str).values
    out_labels = current.copy()
    counts = {}

    # ---- AIRE+ mTEC ----
    aire  = _vec(adata, "AIRE")
    krt14 = _vec(adata, "KRT14")
    cd3e  = _vec(adata, "CD3E")
    cd4   = _vec(adata, "CD4")
    cd8a  = _vec(adata, "CD8A")

    # 2026-05-20 thresholds: use log_present (= log1p(1) ≈ 0.69) for the
    # "high marker present" floor — `log_high` (= log1p(5) ≈ 1.79) yielded
    # 0 cells because thymus rare-type markers (AIRE 0.6-1%, LILRA4 0.3-0.4%)
    # rarely accumulate enough raw counts to hit lognorm ≥ 1.79 in FFPE.
    # Exclusion gates widened from `log_detect` → `log_high`: only exclude
    # cells with STRONG thymocyte signal, not trace neighbor contamination.
    adata.obs["aire_mtec_aire_high"] = (aire >= log_present).astype(np.int8)
    aire_high  = aire >= log_present
    mtec_axis  = krt14 >= log_present
    not_thymo  = (cd3e < log_high) & (cd4 < log_high) & (cd8a < log_high)
    tec_host   = pd.Series(current).isin(_TEC_LIKE).values
    is_aire_mtec = aire_high & mtec_axis & not_thymo & tec_host
    out_labels[is_aire_mtec] = "AIRE_mTEC"
    counts["AIRE_mTEC"] = int(is_aire_mtec.sum())

    # ---- Hassall's corpuscle ----
    krt1   = _vec(adata, "KRT1")
    tgm1   = _vec(adata, "TGM1")
    # TGM3 dropped 2026-05-20: absent from 0041534 panel.
    dsg1   = _vec(adata, "DSG1")
    psmb11 = _vec(adata, "PSMB11")
    ccl25  = _vec(adata, "CCL25")

    krt1_high     = krt1 >= log_present
    cornified_any = (
        (tgm1 >= log_present)
        | (dsg1 >= log_present)
    )
    not_ctec      = (psmb11 < log_high) & (ccl25 < log_high)
    is_hassalls   = (krt1_high & cornified_any & not_ctec & tec_host
                     & ~is_aire_mtec)
    out_labels[is_hassalls] = "Hassalls"
    counts["Hassalls"] = int(is_hassalls.sum())

    # ---- pDC ----
    lilra4 = _vec(adata, "LILRA4")
    irf7   = _vec(adata, "IRF7")
    irf8   = _vec(adata, "IRF8")
    clec4c = _vec(adata, "CLEC4C")
    tcf4   = _vec(adata, "TCF4")

    pdc_high   = lilra4 >= log_present
    pdc_detect = (
        (irf7 >= log_present).astype(int)
        + (irf8 >= log_present).astype(int)
        + (clec4c >= log_present).astype(int)
        + (tcf4 >= log_present).astype(int)
    )
    pdc_not_t  = cd3e < log_high
    dc_host    = pd.Series(current).isin(_DC_LIKE).values
    is_pdc     = pdc_high & (pdc_detect >= 1) & pdc_not_t & dc_host
    out_labels[is_pdc] = "pDC"
    counts["pDC"] = int(is_pdc.sum())

    adata.obs["celltype"] = pd.Categorical(out_labels)
    if "celltype_lineage" in adata.obs.columns:
        cl = adata.obs["celltype_lineage"].astype(str).values.copy()
        cl[is_aire_mtec] = "AIRE_mTEC"
        cl[is_hassalls] = "Hassalls"
        cl[is_pdc] = "pDC"
        adata.obs["celltype_lineage"] = pd.Categorical(cl)

    return counts


# =============================================================================
# Two-tier hierarchical phenotyping — subpanel-driven broad lineages
# (added 2026-05-26; consumed by scripts/hierarchical_phenotype.py +
#  scripts/finalize_panel.py::annotate_subpanel_masks). APPENDED block.
# Mirrors the pancreas structures. Thymus specifics: dominant type is
# DP_thymocyte; the FINAL thymocyte subtyping (DN/DP/SP4/SP8/Treg) + the rare
# types (AIRE_mTEC/Hassalls/pDC) come from apply_rare_celltype_rules() above
# (CD4/CD8/CD3/FOXP3 boolean), which the orchestrator calls after Tier 2.
# =============================================================================
from pathlib import Path as _Path

# Vendored cell-type-category subpanels → broad lineage. Pathway/program panels
# (10_TCR_VDJ, 11_selection, 17-20 signaling, 21-49 HALLMARK/REACTOME/KEGG, 99*)
# and 17_cell_cycle (proliferation = cross-lineage program) are intentionally
# absent. data/reference/subpanels/{thymus,shared}/, pinned @ 79c5ab9.
SUBPANEL_ALLOWLIST = {
    "thymus/01_thymic_epithelial_cortical.csv":  "TEC",
    "thymus/02_thymic_epithelial_medullary.csv": "TEC",
    "thymus/03_hassall_keratinized.csv":         "TEC",
    "thymus/04_thymic_specialized_TEC.csv":      "TEC",
    "thymus/05_thymic_myoid.csv":                "Stromal",   # myoid (muscle-like) → stromal
    "thymus/06_thymocyte_DN.csv":                "Thymocyte",
    "thymus/07_thymocyte_DP.csv":                "Thymocyte",
    "thymus/08_thymocyte_SP.csv":                "Thymocyte",
    "thymus/09_thymic_treg.csv":                 "Thymocyte",
    "thymus/12_thymic_DC_myeloid.csv":           "Myeloid",
    "thymus/13_thymic_B_plasma.csv":             "B_lymphoid",
    "thymus/14_thymic_NK_ILC.csv":               "NK",
    "thymus/15_thymic_stroma_fibroblast.csv":    "Stromal",
    "thymus/16_thymic_endothelial_pericyte.csv": "Vascular",
    "shared/15_neural_glial.csv":                "Neural",
}

# Embedding-union-only (no per-lineage score panel): general epithelial markers.
SUBPANEL_EMBEDDING_EXTRA = ["shared/16_epithelial_general.csv"]

# Broad lineage → fine MARKER_DICT types (the Tier-2 universe). Consistent with
# LINEAGE_GROUPS / LINEAGE_MAP. Proliferating excluded (cycling flag); Hassalls /
# AIRE_mTEC / pDC are post-hoc rule-based (apply_rare_celltype_rules), not here.
BROAD_TO_FINE = {
    "Thymocyte":  {"DN_thymocyte", "DP_thymocyte", "SP4_thymocyte",
                   "SP8_thymocyte", "Treg"},
    "TEC":        {"cTEC", "mTEC"},
    "Myeloid":    {"Myeloid", "DC", "Mast"},
    "B_lymphoid": {"B_cells", "Plasma_cells"},
    "NK":         {"NK"},
    "Vascular":   {"Endothelial", "Lymphatic_EC"},
    "Stromal":    {"Fibroblast", "Smooth_muscle"},
    "Neural":     {"Schwann"},
}

# Detection-gated fine subtyping for every multi-fine broad lineage (mirrors the
# pancreas Endocrine pattern; routed through hierarchical_phenotype.tier2_fine's
# detect_cfg branch, which bypasses sub-clustering). gate = lineage marker union;
# residual = least-specific fine type (assigned only to cells with no specific
# signal). gate_markers built from MARKER_DICT so they never drift.
# Thymocyte is INCLUDED (user choice 2026-05-29): the {DN,DP,SP4,SP8,Treg}_score
# argmax replaces the sub-cluster path AND the boolean apply_thymocyte_subtyping
# (THYMOCYTE_BOOLEAN_SUBTYPING=False). Requires ANTI_OVERRIDE_DOMINANT_TYPE=None
# so the Thymocyte lineage reaches detect_cfg (not the dominant branch).
DETECTION_SUBTYPE_LINEAGES = {
    "Thymocyte": {
        "residual": "DP_thymocyte",
        "subtypes": ["DN_thymocyte", "DP_thymocyte", "SP4_thymocyte",
                     "SP8_thymocyte", "Treg"],
        "gate_markers": sorted(
            set(MARKER_DICT["DN_thymocyte"]) | set(MARKER_DICT["DP_thymocyte"])
            | set(MARKER_DICT["SP4_thymocyte"]) | set(MARKER_DICT["SP8_thymocyte"])
            | set(MARKER_DICT["Treg"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    "TEC": {
        "residual": "cTEC",
        "subtypes": ["cTEC", "mTEC"],
        "gate_markers": sorted(set(MARKER_DICT["cTEC"]) | set(MARKER_DICT["mTEC"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    "Vascular": {
        "residual": "Endothelial",
        "subtypes": ["Endothelial", "Lymphatic_EC"],
        "gate_markers": sorted(set(MARKER_DICT["Endothelial"])
                               | set(MARKER_DICT["Lymphatic_EC"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    "Stromal": {
        "residual": "Fibroblast",
        "subtypes": ["Fibroblast", "Smooth_muscle"],
        "gate_markers": sorted(set(MARKER_DICT["Fibroblast"])
                               | set(MARKER_DICT["Smooth_muscle"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    "Myeloid": {
        "residual": "Myeloid",
        "subtypes": ["Myeloid", "DC", "Mast"],
        "gate_markers": sorted(set(MARKER_DICT["Myeloid"]) | set(MARKER_DICT["DC"])
                               | set(MARKER_DICT["Mast"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    "B_lymphoid": {
        "residual": "B_cells",
        "subtypes": ["B_cells", "Plasma_cells"],
        "gate_markers": sorted(set(MARKER_DICT["B_cells"])
                               | set(MARKER_DICT["Plasma_cells"])),
        "gate_thr": float(np.log1p(1.0)),
    },
}

# Collapse broad vocab to the LINEAGE_MAP groups for renderers wanting the
# coarse coloring of celltype_broad.
BROAD_TO_LINEAGEMAP = {
    "Thymocyte": "Thymocyte", "TEC": "TEC", "Myeloid": "Myeloid",
    "B_lymphoid": "Other_Immune", "NK": "Other_Immune", "Vascular": "Vascular",
    "Stromal": "Stromal", "Neural": "Neural",
}

BROAD_PALETTE = {
    "Thymocyte": "#ff7f0e", "TEC": "#1f77b4", "Myeloid": "#4682b4",
    "B_lymphoid": "#87ceeb", "NK": "#1f3a93", "Vascular": "#e377c2",
    "Stromal": "#2ca02c", "Neural": "#808000",
}

N_MIN_SUBCLUSTER = 2000

# Proliferation is a cross-lineage PROGRAM, not a lineage (e.g. a cycling DP
# thymocyte is DP_thymocyte, cycling — not "Cycling"). Drives the orthogonal
# `cycling` flag + `celltype_cycling` composite, never a standalone label.
CYCLING_MARKERS = ["MKI67", "TOP2A", "PCNA"]


def load_subpanel_genes(reference_dir, allowlist=None, embedding_extra=None):
    """Read vendored subpanel CSVs → ``{broad_lineage: set(gene)}`` (thymus).
    Only the ``gene`` column is read; files mapped to the same lineage are
    unioned. ``embedding_extra`` files merge under ``'_embedding_extra'``.
    Raises ``FileNotFoundError`` if any allowlisted file is missing."""
    allowlist = SUBPANEL_ALLOWLIST if allowlist is None else allowlist
    embedding_extra = (SUBPANEL_EMBEDDING_EXTRA if embedding_extra is None
                       else embedding_extra)
    root = _Path(reference_dir)
    out = {}
    for rel, lineage in allowlist.items():
        df = pd.read_csv(root / rel)
        genes = {g for g in df["gene"].dropna().astype(str) if g}
        out.setdefault(lineage, set()).update(genes)
    extra = set()
    for rel in (embedding_extra or []):
        df = pd.read_csv(root / rel)
        extra |= {g for g in df["gene"].dropna().astype(str) if g}
    out["_embedding_extra"] = extra
    return out


# =============================================================================
# scANVI-classifier redesign (2026-06-08) — tissue tag + curated SEED markers
# =============================================================================
TISSUE = "thymus"

# IDENTITY-CORE rewire (2026-06-08 session 2): the broad CALL is now driven STRICTLY
# by the upstream ontology-curated identity_core. SEED_MARKERS is NO LONGER consumed
# by the broad call (kept for verify_panel_coverage + the legacy notebook path +
# documentation). AIRE/FEZF2 remain promoted into the TEC identity mask by
# finalize_panel (upstream Rec 3).
SEED_MARKERS = {
    "Thymocyte":  ["CD3E", "CD3G", "TRAC", "TRBC2", "CD7", "CD8A", "CD8B", "CD4",
                   "BCL11B", "LEF1", "TCF7"],
    "TEC":        ["EPCAM", "KRT14", "PSMB11", "PRSS16", "CCL25", "FOXN1",
                   "AIRE", "FEZF2", "CTSL", "LY75"],
    "Myeloid":    ["CD68", "CD163", "CSF1R", "MARCO", "CD14", "ITGAM", "ITGAX", "CD1C"],
    "B_lymphoid": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5", "MZB1",
                   "TNFRSF17"],
    "NK":         ["NKG7", "KLRD1", "KLRF1", "NCR3", "KLRK1", "PRF1", "GZMA", "FCGR3A"],
    "Vascular":   ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR", "ENG", "TIE1", "ERG",
                   "CD93", "LYVE1", "PROX1", "PDPN"],
    "Stromal":    ["COL5A1", "COL5A2", "COL6A1", "CDH11", "THY1", "FAP", "PDGFRA",
                   "ACTA2", "MCAM", "PDGFRB", "RGS5"],
    # Schwann-specific only: NES (progenitor/endothelial) + S100A1 (broad) dropped as
    # SEEDS — too non-specific to anchor supervision (kept in MARKER_DICT for scoring).
    "Neural":     ["MPZ", "PMP22", "NGFR", "GLDN", "PRPH"],
}

# identity_core `subtype_hint` → FINE cell-type label (2026-06-08 session 2). Fine
# subtyping prefers identity_core genes whose hint maps here (≥ IC_FINE_MIN of them)
# over MARKER_DICT (sanctioned thin-fine fallback). Thymocyte (DN/DP/CD4→SP4/CD8→SP8/
# Treg), TEC (cTEC/mTEC), B (plasma), myeloid (DC), vascular (lymphatic) are well
# supported by the cores; specialized-TEC hints (microfold/ionocyte/tuft/keratinocyte/
# neuroendocrine) and pan-lineage hints are left unmapped (those fine types use
# MARKER_DICT / the rare rules).
IC_HINT_FINE = {
    "DN":        "DN_thymocyte",
    "DP":        "DP_thymocyte",
    "CD4":       "SP4_thymocyte",
    "CD8":       "SP8_thymocyte",
    "Treg":      "Treg",
    "cTEC":      "cTEC",
    "mTEC":      "mTEC",
    "plasma":    "Plasma_cells",
    "DC":        "DC",
    "lymphatic": "Lymphatic_EC",
}
