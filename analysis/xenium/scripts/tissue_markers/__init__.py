"""Tissue-specific marker definitions for the Xenium pipeline.

Each submodule (pancreas, thymus) exports a uniform interface:

    MARKER_DICT                   dict[str, list[str]]
                                  cell-type → score_genes panel (≥3 genes per panel)

    LOW_CONFIDENCE_TYPES          set[str]
                                  celltypes to downgrade to Indeterminate_*
                                  when celltype_confidence < median-MAD.

    INDETERMINATE_BUCKETS         dict[str, str]
                                  Maps a LOW_CONFIDENCE_TYPES member to the
                                  bucket name written when downgraded.

    ANTI_OVERRIDE_DOMINANT_TYPE   str
                                  The cell-type whose argmax-mislabels are
                                  being corrected (Acinar / DP_thymocyte).

    ANTI_OVERRIDE_MARKERS         dict[str, list[str]]
                                  alt cell-type → specific marker panel.
                                  For each <dominant> cell, if any alt's
                                  per-cell mean exceeds its threshold, reassign.

    LINEAGE_CANARY_PANELS         dict[str, list[str]]
                                  lineage → small specific marker panel for
                                  per-cell depth counting (cell 22 lineage
                                  refinement + scripts/lineage_phenotype.py).

    LINEAGE_EXCLUSIVE_GROUP       dict[str, str]
                                  Maps each LINEAGE_CANARY_PANELS lineage to
                                  its mutually-exclusive group name; cells
                                  with ≥3 groups detected are tagged
                                  doublet_suspected.

    LINEAGE_MAP                   dict[str, str]
                                  celltype → broad-lineage mapping for the
                                  stage 03 ligrec cluster_key='celltype_general'
                                  aggregation.

    PANEL_COVERAGE_PANELS         dict[str, list[str]]
                                  cell-type → full marker panel used for the
                                  panel-coverage gate in refine_rare_celltypes.

    COMMON_CLUSTERING_TYPES       set[str]
                                  Types that get the spatial-coherence gate
                                  in refine_rare_celltypes (rare types are
                                  exempted to avoid losing genuine isolated
                                  cells).

    apply_rare_celltype_rules(adata, log_detect, log_high) -> dict
                                  Applies tissue-specific multi-criteria
                                  rules for rare cell-type identification.
                                  Returns counts per rare type. Updates
                                  adata.obs['celltype'] and 'celltype_lineage'
                                  in place.
"""
