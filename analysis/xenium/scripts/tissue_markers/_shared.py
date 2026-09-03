"""Pan-tissue cell-type marker panels shared across tissues (pancreas + thymus).

Single source of truth for cell types that are biologically identical across
tissues. Tissue modules (pancreas.py, thymus.py) import SHARED_MARKERS and build
their MARKER_DICT as ``{**shared_subset, **tissue_specific}``, overriding
individual keys where a tissue needs a richer / different panel (e.g. thymus keeps
its tuned Endothelial/Myeloid/DC/NK/Schwann sets; pancreas Endothelial adds
ERG/TIE1/CD93).

Every gene here is verified PRESENT on the Xenium 5k panel for the pancreas
sample (0041323) checked 2026-05-29, and each panel has >= 3 genes
(CLAUDE.md decision #12). Build-time asserts in
hierarchical_phenotype.build_available_markers re-check on-panel membership
per sample and fail loud on a thin (<2-gene) survivor set — so a future panel
swap that drops one of these is caught, not silently absorbed.

OFF-PANEL genes intentionally excluded (do NOT add them back):
  GNLY, NCAM1                 (NK)
  VWF, FLT1                   (Endothelial)
  LYZ, C1QA, C1QB             (Myeloid)
  TPSAB1, TPSB2, CPA3         (Mast tryptases)
  S100B, PLP1, MBP            (Schwann)
  CD3D                        (T cells)
  XBP1, DERL3, FKBP11,
  TENT5C, POU2AF1             (Plasma — pan-secretory / pan-B, not plasma-specific:
                               they fire in secretory acinar and drove over-calls)

Import (relative — these modules always load as
``scripts.tissue_markers.<tissue>`` via importlib, so the package is set):
    from ._shared import SHARED_MARKERS
"""

SHARED_MARKERS = {
    # ---- Lymphoid ----
    # Pan-T (CD3D off-panel; FOXP3 reserved for Treg in thymus). CD4/CD8A/CD8B
    # included as T-lineage markers (pancreas has no thymocyte split).
    "T_cells":       ["CD3E", "CD3G", "CD2", "CD5", "CD7", "PTPRC", "TRAC",
                      "TRBC2", "CD8A", "CD8B", "CD4"],
    "B_cells":       ["MS4A1", "CD79A", "CD79B", "CD19", "CD22"],
    # Plasma — tight & plasma-specific only (no pan-secretory genes; see header).
    "Plasma_cells":  ["MZB1", "TNFRSF17", "PRDM1"],
    # NK — GNLY off-panel; set proven on-panel via the thymus NK panel.
    "NK":            ["NKG7", "KLRD1", "KLRF1", "NCR3", "KLRK1", "PRF1",
                      "GZMA", "KLRB1", "FCGR3A"],
    # ---- Myeloid / DC / Mast ----
    "Myeloid":       ["CD68", "CD163", "CSF1R", "MARCO", "CD14", "ITGAM"],
    "DC":            ["ITGAX", "CD1C", "CLEC9A", "CLEC10A", "BATF3", "IRF8",
                      "FLT3", "LAMP3"],
    "Mast":          ["KIT", "CMA1", "HDC", "MS4A2", "FCER1A"],
    # ---- Vascular ----
    "Endothelial":   ["PECAM1", "CDH5", "CLDN5", "PLVAP", "KDR", "ENG"],
    "Lymphatic_EC":  ["LYVE1", "PROX1", "PDPN", "CCL21", "FLT4"],
    # ---- Neural / cycling ----
    "Schwann":       ["MPZ", "PMP22", "NES", "NGFR", "SOX10"],
    "Proliferating": ["MKI67", "TOP2A", "PCNA"],
}
