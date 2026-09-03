"""Pancreas marker definitions for the Xenium pipeline.

Lifted verbatim from the pre-refactor inline sources:
  - notebooks/02_phenotyping.ipynb cell 10 (MARKER_DICT, LOW_CONFIDENCE_TYPES)
  - notebooks/02_phenotyping.ipynb cell 15 (ANTI_ACINAR_MARKERS)
  - notebooks/02_phenotyping.ipynb cell 22 (LINEAGE_CANARY_PANELS,
                                              LINEAGE_EXCLUSIVE_GROUP)
  - notebooks/03_spatial_analysis.ipynb cell 14 (LINEAGE_MAP)
  - scripts/verify_upper_tail_doublets.py (CANARY_PANELS, EXCLUSIVE_GROUP)
  - scripts/refine_rare_celltypes.py (PAN_ENDOCRINE, GAMMA/EPSILON rules)

DO NOT EDIT inline values here without also updating the pancreas oracle
test in Verification §3 of the plan, or you'll silently break the
pancreas pipeline regression.
"""
import numpy as np
import pandas as pd

from ._shared import SHARED_MARKERS

# -----------------------------------------------------------------------------
# Score-genes argmax panels (cell 10 of 02_phenotyping.ipynb)
# Pan-tissue types come from _shared.SHARED_MARKERS (single source of truth);
# pancreas-specific islet/exocrine/stromal panels + per-key overrides below.
# All genes verified on the 5k panel 2026-05-29; build_available_markers
# fail-loud asserts re-check per sample (decision #12).
# -----------------------------------------------------------------------------

# Pancreas-specific fine-type panels. Off-panel genes removed (absent from the
# 5k panel): Alpha -TM4SF4,-TTR (+IRX1); Delta -RBP4,-BHLHE41 (+POU2F2);
# Acinar -AMY2A,-AMY2B (had left AMY1A only — the 1-gene bug); Ductal -KRT7
# (+TSPAN8); Stellate -DCN,-COL1A1 (+THY1,+CDH11). Acinar is panel-limited to
# AMY1A+CUZD1 (2 genes); its FINE call is the detection-based Exocrine dominant
# branch (not score-argmax), so 2 is acceptable; W5 (validate_markers.py) mines
# a data-driven 3rd from rank_genes_groups.
_PANCREAS_SPECIFIC = {
    "Beta":      ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "ABCC8",
                  "KCNJ11", "GLP1R", "PCSK1"],
    "Alpha":     ["ARX", "MAFB", "GC", "IRX1"],
    "Delta":     ["HHEX", "LEPR", "PCSK1", "POU2F2"],   # PCSK1 Beta-shared (W5 flags)
    # Epsilon + Gamma intentionally absent — identified post-hoc via
    # rule-based multi-criteria in apply_rare_celltype_rules() below.
    "Endocrine": ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV", "PAX6", "SCG5",
                  "SCG2", "SCGN"],
    "Acinar":    ["AMY1A", "CUZD1"],
    "Ductal":    ["KRT19", "SOX9", "HNF1B", "MUC1", "CFTR", "TSPAN8"],
    "Stellate":  ["ACTA2", "PDGFRB", "RGS5", "PDGFRA", "THY1", "CDH11"],
}

# Pan-tissue types pancreas adopts from _shared (Lymphatic_EC not expected in
# pancreas). DC/Mast/NK added for the T1D insulitis focus (user choice
# 2026-05-29); routed via BROAD_TO_FINE + DETECTION_SUBTYPE_LINEAGES below.
_PANCREAS_SHARED_KEYS = ("T_cells", "B_cells", "Plasma_cells", "NK",
                         "Myeloid", "DC", "Mast", "Endothelial", "Schwann",
                         "Proliferating")

MARKER_DICT = {
    **{k: list(SHARED_MARKERS[k]) for k in _PANCREAS_SHARED_KEYS},
    **_PANCREAS_SPECIFIC,
    # pancreas Endothelial override: shared 6 + arterial/capillary markers
    "Endothelial": list(SHARED_MARKERS["Endothelial"]) + ["ERG", "TIE1", "CD93"],
}

# Lineage groupings used by notebook 02 cell 15 to compute the
# lineage-margin form of celltype_confidence. Aggregating sub-types that
# share markers (e.g. islet endocrine sub-types share CHGA family; Acinar
# and Ductal share epithelial program) makes the margin metric report
# "how clearly does this cell belong to the lineage" rather than "how
# clearly is it Alpha vs Beta". Epsilon and Gamma are intentionally
# omitted from the Endocrine set because they're identified post-hoc via
# scripts/refine_rare_celltypes.py rules, not via score_genes argmax.
LINEAGE_GROUPS = {
    "Endocrine":  {"Alpha", "Beta", "Delta", "Endocrine"},
    "Exocrine":   {"Acinar", "Ductal"},
    "Stromal":    {"Stellate"},
    "Vascular":   {"Endothelial"},
    "Lymphoid":   {"T_cells", "B_cells", "Plasma_cells", "NK"},
    "Myeloid":    {"Myeloid", "DC", "Mast"},
    "Neural":     {"Schwann"},
    "Cycling":    {"Proliferating"},
}

LOW_CONFIDENCE_TYPES = {"Alpha", "Delta", "Acinar"}

# Maps a LOW_CONFIDENCE_TYPES member to the indeterminate bucket name.
# Cell 15 originally hard-coded: endocrine for Alpha/Delta/Epsilon,
# else exocrine.
INDETERMINATE_BUCKETS = {
    "Alpha":   "Indeterminate_endocrine",
    "Delta":   "Indeterminate_endocrine",
    "Epsilon": "Indeterminate_endocrine",
    "Acinar":  "Indeterminate_exocrine",
}

# -----------------------------------------------------------------------------
# Anti-acinar override (cell 15 of 02_phenotyping.ipynb)
# -----------------------------------------------------------------------------

ANTI_OVERRIDE_DOMINANT_TYPE = "Acinar"

ANTI_OVERRIDE_MARKERS = {
    "Beta":        ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8"],
    "Endocrine":   ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV"],
    "Ductal":      ["KRT19", "CFTR", "MUC1"],
    "Stellate":    ["ACTA2", "PDGFRB", "RGS5"],
    "Endothelial": ["PECAM1", "CDH5", "CLDN5", "PLVAP"],
    "T_cells":     ["CD3D", "CD3E", "CD8A", "PTPRC"],
    "B_cells":     ["MS4A1", "CD79A"],
    "Myeloid":     ["CD68", "CD163", "MARCO"],
    "Schwann":     ["MPZ", "SOX10"],
}

# -----------------------------------------------------------------------------
# Lineage refinement / doublet detection canary panels (cell 22 +
# scripts/verify_upper_tail_doublets.py + scripts/lineage_phenotype.py)
# -----------------------------------------------------------------------------

LINEAGE_CANARY_PANELS = {
    "Acinar":        ["AMY1A", "CUZD1"],
    "Ductal":        ["KRT19", "CFTR", "MUC1"],
    "Beta":          ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8"],
    "Alpha":         ["ARX", "MAFB"],
    "Endocrine_pan": ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV"],
    "Stellate":      ["ACTA2", "PDGFRB", "RGS5"],
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP"],
    "T_cells":       ["CD3D", "CD3E", "CD8A", "PTPRC"],
    "B_cells":       ["MS4A1", "CD79A"],
    "Myeloid":       ["CD68", "CD163", "MARCO"],
    "Schwann":       ["MPZ", "SOX10"],
}

LINEAGE_EXCLUSIVE_GROUP = {
    "Acinar":        "Exocrine_epi",
    "Ductal":        "Exocrine_epi",
    "Beta":          "Endocrine_islet",
    "Alpha":         "Endocrine_islet",
    "Endocrine_pan": "Endocrine_islet",
    "Stellate":      "Stromal",
    "Endothelial":   "Endothelial",
    "T_cells":       "Lymphoid",
    "B_cells":       "Lymphoid",
    "Myeloid":       "Myeloid",
    "Schwann":       "Neural",
}

# -----------------------------------------------------------------------------
# Ligand-receptor lineage map (cell 14 of 03_spatial_analysis.ipynb)
# -----------------------------------------------------------------------------

LINEAGE_MAP = {
    "Alpha":         "Endocrine",
    "Beta":          "Endocrine",
    "Delta":         "Endocrine",
    "Gamma":         "Endocrine",
    "Epsilon":       "Endocrine",
    "Endocrine":     "Endocrine",
    "Endocrine_pan": "Endocrine",
    "Acinar":        "Exocrine",
    "Ductal":        "Exocrine",
    "T_cells":       "Immune",
    "B_cells":       "Immune",
    "Plasma_cells":  "Immune",
    "NK":            "Immune",
    "Myeloid":       "Immune",
    "DC":            "Immune",
    "Mast":          "Immune",
    "Stellate":      "Stromal",
    "Endothelial":   "Vascular",
    "Schwann":       "Neural",
    "Proliferating": "Cycling",
}

# -----------------------------------------------------------------------------
# Consensus-annotation 3-panel plot dispatch (notebook 02 cell 20).
# Pancreas vocab: specific celltypes shown in panel 2; panel 3 drops the
# exocrine bulk (Acinar+Ductal) so rare endocrine + immune + stromal types
# are visible.
# -----------------------------------------------------------------------------

SPECIFIC_KEEP = {
    "Alpha", "Beta", "Delta", "Gamma", "Epsilon",
    "Acinar", "Ductal",
    "T_cells", "B_cells", "Plasma_cells", "NK", "Myeloid", "DC", "Mast",
    "Stellate", "Endothelial", "Schwann", "Proliferating",
}

PALETTE = {
    "Alpha":         "#d62728",
    "Beta":          "#ff7f0e",
    "Delta":         "#8c564b",
    "Gamma":         "#fa8072",
    "Epsilon":       "#e377c2",
    "Endocrine":     "#b22222",
    "Acinar":        "#7fbf3f",
    "Ductal":        "#17becf",
    "Exocrine":      "#2ca02c",
    "T_cells":       "#1f3a93",
    "B_cells":       "#87ceeb",
    "Plasma_cells":  "#5e3c99",
    "NK":            "#2c7fb8",
    "Myeloid":       "#4682b4",
    "DC":            "#6a3d9a",
    "Mast":          "#fa8072",
    "Immune":        "#1f77b4",
    "Stellate":      "#9467bd",
    "Stromal":       "#6a3d9a",
    "Endothelial":   "#e83e8c",
    "Vascular":      "#e83e8c",
    "Schwann":       "#808000",
    "Neural":        "#808000",
    "Proliferating": "#7f7f7f",
    "Cycling":       "#7f7f7f",
}

THIRD_PANEL = {
    "exclude_types": {"Acinar", "Ductal"},
    "title":         "Specific celltype, no exocrine",
}

# -----------------------------------------------------------------------------
# Lineage-subset panels for predicted_celltypes.png / scvi_vs_markers.png.
# Each subset is rendered as its own UMAP with high-contrast colors and a
# light-grey background of out-of-subset cells (scanpy's `groups=` mechanism).
# Lineage names refer to LINEAGE_MAP values above.
# -----------------------------------------------------------------------------

LINEAGE_SUBSETS = [
    {"title": "Endocrine",                     "lineages": ["Endocrine"]},
    {"title": "Exocrine",                      "lineages": ["Exocrine"]},
    {"title": "Immune / Stromal / Vascular / Neural / Cycling",
     "lineages": ["Immune", "Stromal", "Vascular", "Neural", "Cycling"]},
]

# predicted_celltypes.png / scvi_vs_markers.png colored-panel mode:
#   "major_lineage"   — collapse fine cell types through LINEAGE_MAP into the
#                       ~7 major lineages and color a single UMAP panel. The
#                       de-crowded view the user asked for (2026-05-26): the
#                       fine 14-15-type palette was too scattered to read the
#                       acinar-bulk vs islet vs immune/stroma separation.
#   "lineage_subsets" — the older 3-panel LINEAGE_SUBSETS split (one high-
#                       contrast subset per panel over a grey background).
# Pancreas uses "major_lineage". Tissues that omit this attribute (e.g. thymus,
# where the thymocyte DN/DP/SP detail must stay visible) default to
# "lineage_subsets" in the renderers.
PREDICTED_VIEW = "major_lineage"

# -----------------------------------------------------------------------------
# Panel-coverage panels (refine_rare_celltypes.py lines 280-297)
# -----------------------------------------------------------------------------

PAN_ENDOCRINE_PANEL = ["CHGA", "INSM1", "ISL1", "NEUROD1", "FEV", "SCG2", "SCGN"]

PANEL_COVERAGE_PANELS = {
    "Beta":          ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "DLK1"],
    "Alpha":         ["ARX", "MAFB", "GC", "TM4SF4", "TTR"],
    "Delta":         ["HHEX", "PCSK1", "RBP4", "BHLHE41"],
    "Acinar":        ["CTRB1", "CPA1", "CPA2", "CELA1", "PRSS1", "PRSS2",
                       "AMY2A", "AMY2B", "REG1A"],
    "Ductal":        ["KRT19", "KRT7", "CFTR", "MUC1", "SOX9"],
    "Stellate":      ["DCN", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "COL3A1"],
    "Endothelial":   ["PECAM1", "CDH5", "VWF", "CLDN5", "FLT1", "KDR"],
    "T_cells":       ["CD3E", "CD3G", "CD2", "CD5", "PTPRC"],
    "B_cells":       ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
    "Myeloid":       ["CD68", "CD163", "CD14", "ITGAM", "FCGR3A"],
    "Schwann":       ["S100B", "MPZ", "PMP22", "PLP1", "NRXN1"],
    "Endocrine":     list(PAN_ENDOCRINE_PANEL),
    "Endocrine_pan": list(PAN_ENDOCRINE_PANEL),
}

COMMON_CLUSTERING_TYPES = {"Beta", "Alpha", "Acinar", "Ductal", "Stellate",
                            "Schwann"}

# Map celltype labels (as written in adata.obs['celltype']) to their canary
# lineage in LINEAGE_CANARY_PANELS. Used by
# verify_upper_tail_doublets.py to check whether each cell's assigned
# celltype's canary panel was detected.
LABEL_TO_CANARY = {
    "Acinar":        "Acinar",
    "Ductal":        "Ductal",
    "Beta":          "Beta",
    "Alpha":         "Alpha",
    "Delta":         "Endocrine_pan",
    "Endocrine":     "Endocrine_pan",
    "Epsilon":       "Endocrine_pan",
    "Gamma":         "Endocrine_pan",
    "Stellate":      "Stellate",
    "Endothelial":   "Endothelial",
    "T_cells":       "T_cells",
    "B_cells":       "B_cells",
    "Myeloid":       "Myeloid",
    "Schwann":       "Schwann",
}

# Celltype labels (or Indeterminate variants) that
# apply_rare_celltype_rules() considers as candidate Gamma / Epsilon hosts.
_ENDO_LIKE = {"Beta", "Alpha", "Delta", "Endocrine", "Endocrine_pan",
              "Indeterminate_Beta", "Indeterminate_Alpha", "Indeterminate_Delta",
              "Indeterminate_Endocrine", "Indeterminate_endocrine"}

# -----------------------------------------------------------------------------
# Rare cell-type rules (Gamma / Epsilon)
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


def apply_rare_celltype_rules(adata, log_detect, log_high):
    """Multi-criteria Gamma + Epsilon identification on the LOGNORM layer.

    Logic copied verbatim from refine_rare_celltypes.py apply_rare_celltype_rules().
    Returns dict with the count per rare type.
    """
    pan_layers = []
    for g in PAN_ENDOCRINE_PANEL:
        v = _get_lognorm_col(adata, g)
        if v is not None:
            pan_layers.append(v >= log_detect)
    if not pan_layers:
        return {"Gamma": 0, "Epsilon": 0}
    pan_count = np.stack(pan_layers, axis=1).sum(axis=1)
    adata.obs["pan_endocrine_marker_count"] = pan_count.astype(np.int16)

    iapp = _get_lognorm_col(adata, "IAPP")
    mafa = _get_lognorm_col(adata, "MAFA")
    if iapp is None:
        iapp = np.zeros(adata.n_obs)
    if mafa is None:
        mafa = np.zeros(adata.n_obs)
    not_beta = (iapp < log_detect) & (mafa < log_detect)

    hhex = _get_lognorm_col(adata, "HHEX")
    if hhex is None:
        hhex = np.zeros(adata.n_obs)
    not_delta = (hhex < log_detect)

    aqp3 = _get_lognorm_col(adata, "AQP3")
    etv1 = _get_lognorm_col(adata, "ETV1")
    ghrl = _get_lognorm_col(adata, "GHRL")
    if aqp3 is None:
        aqp3 = np.zeros(adata.n_obs)
    if etv1 is None:
        etv1 = np.zeros(adata.n_obs)
    if ghrl is None:
        ghrl = np.zeros(adata.n_obs)

    gamma_markers_count = ((aqp3 >= log_detect).astype(int)
                            + (etv1 >= log_detect).astype(int))
    adata.obs["gamma_marker_count"] = gamma_markers_count.astype(np.int16)

    not_gamma_for_epsilon = (aqp3 < log_high) | (etv1 < log_detect)

    epsilon_markers_count = (ghrl >= log_detect).astype(int)
    adata.obs["epsilon_marker_count"] = epsilon_markers_count.astype(np.int16)
    not_epsilon_for_gamma = (ghrl < log_detect)

    is_gamma = (
        (pan_count >= 2)
        & ((gamma_markers_count >= 2) | (aqp3 >= log_high))
        & not_beta
        & not_delta
        & not_epsilon_for_gamma
    )
    is_epsilon = (
        (pan_count >= 2)
        & (ghrl >= log_detect)
        & not_beta
        & not_delta
        & not_gamma_for_epsilon
    )

    current = adata.obs["celltype"].astype(str).values
    endo_mask = pd.Series(current).isin(_ENDO_LIKE).values
    is_gamma &= endo_mask
    is_epsilon &= endo_mask
    overlap = is_gamma & is_epsilon
    if overlap.any():
        is_epsilon[overlap & (aqp3[overlap] >= log_high)] = False
        is_gamma[overlap & ~(aqp3[overlap] >= log_high)] = False

    new_lab = current.copy()
    new_lab[is_gamma] = "Gamma"
    new_lab[is_epsilon] = "Epsilon"
    adata.obs["celltype"] = pd.Categorical(new_lab)

    if "celltype_lineage" in adata.obs.columns:
        cl = adata.obs["celltype_lineage"].astype(str).values.copy()
        cl[is_gamma] = "Gamma"
        cl[is_epsilon] = "Epsilon"
        adata.obs["celltype_lineage"] = pd.Categorical(cl)

    return {"Gamma": int(is_gamma.sum()), "Epsilon": int(is_epsilon.sum())}


# =============================================================================
# Two-tier hierarchical phenotyping — subpanel-driven broad lineages
# (added 2026-05-26; consumed by scripts/hierarchical_phenotype.py and
#  scripts/finalize_panel.py::annotate_subpanel_masks). APPENDED block — the
#  structures above are unchanged.
# =============================================================================
from pathlib import Path as _Path

# Vendored cell-type-category subpanels → broad lineage. Each file is a curated
# gene pool for one biological compartment (data/reference/subpanels/, pinned to
# XeniumPanelExplorer @ 79c5ab9; see that dir's SOURCE.txt). These are EMBEDDING
# FEATURE-SELECTION sets (large, overlapping, include housekeeping genes), NOT
# tight score panels — fine subtyping still uses MARKER_DICT (decision #12).
# Pathway/program panels are intentionally absent (never ingested).
SUBPANEL_ALLOWLIST = {
    "pancreas/01_pancreas_endocrine.csv":     "Endocrine",
    "pancreas/02_pancreas_exocrine.csv":      "Exocrine",
    "pancreas/13_fibroblast_stellate.csv":    "Stromal",
    "shared/03_immune_T_cell.csv":            "Immune_T",
    "shared/04_immune_B_plasma.csv":          "Immune_B",
    "shared/05_immune_NK_ILC.csv":            "Immune_T",  # NK folded into T (no NK fine type in pancreas MARKER_DICT)
    "shared/06_immune_myeloid.csv":           "Myeloid",
    "shared/07_immune_granulocyte.csv":       "Myeloid",   # granulocytes folded into Myeloid for pancreas
    "shared/11_endothelial_vascular.csv":     "Vascular",
    "shared/12_pericyte_smooth_muscle.csv":   "Stromal",   # pericyte/SMC join stellate under Stromal
    "shared/15_neural_glial.csv":             "Neural",
    # NOTE: 17_cell_cycle_proliferation is intentionally NOT here — proliferation
    # is a cross-lineage PROGRAM, not a lineage. Including cell-cycle genes in the
    # embedding union would also make cells cluster by cycle phase, not lineage.
    # Cycling is carried as an orthogonal flag (see CYCLING_MARKERS below).
}

# Files whose genes join the Tier-1 embedding union but are NOT used as a
# per-lineage broad-score panel: general epithelial markers (KRT*/EPCAM) are
# shared by ductal + acinar + endocrine, so scoring them as one lineage would
# bias the broad-score argmax. They still help the embedding separate epithelium
# from non-epithelium.
SUBPANEL_EMBEDDING_EXTRA = ["shared/16_epithelial_general.csv"]

# Broad lineage → fine MARKER_DICT types eligible inside it (the Tier-2 universe).
# Consistent with LINEAGE_GROUPS / LINEAGE_MAP above. Gamma/Epsilon are NOT
# listed — they are added post-hoc inside Endocrine via apply_rare_celltype_rules.
BROAD_TO_FINE = {
    "Endocrine": {"Beta", "Alpha", "Delta", "Endocrine"},
    "Exocrine":  {"Acinar", "Ductal"},
    "Stromal":   {"Stellate"},
    "Vascular":  {"Endothelial"},
    "Immune_T":  {"T_cells", "NK"},          # NK split via DETECTION_SUBTYPE_LINEAGES
    "Immune_B":  {"B_cells", "Plasma_cells"},
    "Myeloid":   {"Myeloid", "DC", "Mast"},
    "Neural":    {"Schwann"},
    # No "Cycling" — proliferation is a cross-lineage program (see CYCLING_MARKERS).
}

# Detection-gated fine subtyping for broad lineages that pair a GENERIC residual
# fine type with specific subtypes. Endocrine is the only such lineage: cluster-
# mean argmax over {Beta,Alpha,Delta,Endocrine} let the broad pan-endocrine
# `Endocrine_score` (CHGA/INSM1/ISL1…, high in EVERY islet cell) out-compete the
# specific Beta/Delta panels and erased β/δ in the 2026-05-26 pilot (decision #12
# dilution, at cluster level). Fix (mirrors Exocrine/thymocyte detection logic):
# assign each endocrine cell by argmax of the SPECIFIC {subtype}_score ONLY
# (generic `Endocrine` dropped from the competition); a cell with NO specific
# marker detected (gate, lognorm >= log1p(1) — log1p(3) is too strict for islet
# TFs, cf. thymocyte subtyping) stays the generic `residual`. Validated to give
# islet-plausible β:α:δ (0041323 ~45:46:9; 0041326 ~16:71:13, advanced-T1D
# β-depletion). Consumed by hierarchical_phenotype.tier2_fine.
DETECTION_SUBTYPE_LINEAGES = {
    "Endocrine": {
        "residual": "Endocrine",
        "subtypes": ["Beta", "Alpha", "Delta"],
        "gate_markers": ["IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8",
                         "ARX", "MAFB", "HHEX"],
        "gate_thr": float(np.log1p(1.0)),
    },
    # B vs plasma: both specific, both compete; residual B_cells (cells with B-
    # lineage signal but no plasma marker). gate = union of both fine panels.
    "Immune_B": {
        "residual": "B_cells",
        "subtypes": ["B_cells", "Plasma_cells"],
        "gate_markers": sorted(set(MARKER_DICT["B_cells"])
                               | set(MARKER_DICT["Plasma_cells"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    # NK split out of the Immune_T broad lineage (panel-marginal: NK markers
    # overlap CD8 effector — W5 confusion check audits NK-vs-CD8). residual T_cells.
    "Immune_T": {
        "residual": "T_cells",
        "subtypes": ["T_cells", "NK"],
        "gate_markers": sorted(set(MARKER_DICT["T_cells"])
                               | set(MARKER_DICT["NK"])),
        "gate_thr": float(np.log1p(1.0)),
    },
    # DC + Mast split out of the Myeloid broad lineage. residual Myeloid (the
    # least-specific macrophage/monocyte bulk).
    "Myeloid": {
        "residual": "Myeloid",
        "subtypes": ["Myeloid", "DC", "Mast"],
        "gate_markers": sorted(set(MARKER_DICT["Myeloid"]) | set(MARKER_DICT["DC"])
                               | set(MARKER_DICT["Mast"])),
        "gate_thr": float(np.log1p(1.0)),
    },
}

# Collapse the split-immune broad vocabulary back to the LINEAGE_MAP 7-group view
# for any renderer that wants the old coarse coloring of celltype_broad.
BROAD_TO_LINEAGEMAP = {
    "Endocrine": "Endocrine", "Exocrine": "Exocrine", "Stromal": "Stromal",
    "Vascular": "Vascular", "Immune_T": "Immune", "Immune_B": "Immune",
    "Myeloid": "Immune", "Neural": "Neural",
}

# Palette for celltype_broad labels (Immune_T/Immune_B reuse the Immune family).
BROAD_PALETTE = {
    "Endocrine": "#b22222", "Exocrine": "#2ca02c", "Stromal": "#6a3d9a",
    "Vascular": "#e83e8c", "Immune_T": "#1f3a93", "Immune_B": "#87ceeb",
    "Myeloid": "#4682b4", "Neural": "#808000",
}

# Lineages with fewer than this many cells skip the Tier-2 sub-embedding and fall
# back to per-cell score-argmax (leiden on a few hundred cells is unstable).
N_MIN_SUBCLUSTER = 2000

# Proliferation is a cross-lineage PROGRAM, not a lineage: a cycling ductal cell
# is Ductal (cycling), not "Cycling". These markers drive an orthogonal `cycling`
# flag + a `celltype_cycling` composite (e.g. "Ductal_cycling") — proliferation is
# never a standalone lineage or celltype.
CYCLING_MARKERS = ["MKI67", "TOP2A", "PCNA"]


def load_subpanel_genes(reference_dir, allowlist=None, embedding_extra=None):
    """Read vendored subpanel CSVs → ``{broad_lineage: set(gene)}``.

    Only the ``gene`` column is read; files mapped to the same lineage are
    unioned. ``embedding_extra`` files (no per-lineage role) are merged under the
    special key ``'_embedding_extra'`` (callers add these to the embedding union
    but build NO per-lineage mask for them). ``reference_dir`` is the directory
    holding the ``pancreas/`` + ``shared/`` subfolders; allowlist keys are
    relative to it. Raises ``FileNotFoundError`` if any allowlisted file is
    missing (reproducibility: the vendored snapshot must be complete).
    """
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
TISSUE = "pancreas"

# IDENTITY-CORE rewire (2026-06-08 session 2): the broad CALL is now driven STRICTLY
# by the upstream ontology-curated identity_core (scripts/tissue_markers/_identity_core
# + finalize_panel's var['subpanel_<L>'] = identity_core) — absolute mean lognorm,
# raw-count anchors, scANVI clean-anchor+Unknown+propagate, expression-masked predict.
# SEED_MARKERS is NO LONGER consumed by the broad call (kept here only for
# verify_panel_coverage + the legacy notebook path + as documentation of the curated
# surrogates; canonical hormones INS/GCG/SST/PPY + acinar zymogens are OFF-panel).
SEED_MARKERS = {
    "Endocrine": ["CHGA", "INSM1", "ISL1", "NEUROD1", "PAX6", "SCG2", "SCGN",
                  "IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "ARX", "MAFB", "HHEX"],
    "Exocrine":  ["AMY1A", "CUZD1", "KRT19", "SOX9", "HNF1B", "CFTR", "MUC1", "TSPAN8"],
    "Immune_T":  list(SHARED_MARKERS["T_cells"]),
    "Immune_B":  list(SHARED_MARKERS["B_cells"]) + list(SHARED_MARKERS["Plasma_cells"]),
    "Myeloid":   ["CD68", "CD163", "CSF1R", "MARCO", "CD14", "ITGAM", "ITGAX", "CD1C"],
    "Stromal":   ["ACTA2", "PDGFRB", "RGS5", "PDGFRA", "THY1", "CDH11"],
    "Vascular":  ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR", "ENG", "ERG", "TIE1", "CD93"],
    # Schwann-specific only: NES (neural-progenitor/endothelial) dropped as a SEED —
    # too broad to anchor supervision (it stays in MARKER_DICT for scoring).
    "Neural":    ["MPZ", "PMP22", "NGFR", "SOX10"],
}

# identity_core `subtype_hint` → FINE cell-type label (2026-06-08 session 2). The
# fine tier (hierarchical_phenotype.tier2_fine) prefers the identity_core genes whose
# subtype_hint maps here (≥ IC_FINE_MIN of them) over MARKER_DICT — the latter is the
# sanctioned thin-fine fallback. Only specific, well-supported endocrine/exocrine/
# plasma/DC hints are mapped; pan-lineage hints (monocyte/macrophage/CD4/CD8/neuron/
# pericyte/…) are left unmapped so those fine types keep their tuned MARKER_DICT panel.
IC_HINT_FINE = {
    "beta":   "Beta",
    "alpha":  "Alpha",
    "delta":  "Delta",
    "ductal": "Ductal",
    "acinar": "Acinar",
    "plasma": "Plasma_cells",
    "DC":     "DC",
}
