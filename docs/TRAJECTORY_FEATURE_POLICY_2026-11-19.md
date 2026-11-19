# Islet composition graph inputs: 2026-11-19

The islet-composition producer selects graph inputs separately from displayed
phenotypes. QC/ambiguous labels (including unresolved and conflict labels)
do not enter either graph. Endocrine classes enter the core graph only;
exocrine classes enter the peri-islet graph only. Other supported named
classes may enter either region. When multiple measured macrophage marker
profiles identify Gal3+ or Gal3−, their counts are pooled by sign before the
transformation. Missing data are not imputed, and an aggregate count supplied
alongside its component counts is refused rather than double-counted.

The existing Core mode remains core-only. The Core + Peri mode additionally
includes eligible peri-islet classes, including exocrine. The graph feature
names are recorded in H5AD metadata; raw class counts and displayed per-islet
proportions are not removed.

This is a feature policy, not a claim that all pipelines use the same fit.
The Senior composition method retains its additive log ratio with pseudocount
0.5, Euclidean distance, root policy, cohort and validation gates. CODEX uses a
different log1p/cosine calculation. Legacy input tables can contain only broad
immune classes without measured Gal3 subtypes; those subtypes cannot be
invented from the available columns.

The separate core-only ND-reference method retains its modality-specific
features and model; it excludes QC-like features and archived split macrophage
features and admits the already-pooled Gal3+/- fields only when measured.

**Deployment boundary:** the current Senior app's adata_ins_root.h5ad
pseudotime is rebuilt from a separate within-donor single-cell marker
trajectory by scripts/rebuild_adata_ins_root.py; this change to the pooled
composition producer does not replace that ordering. A scientific migration
requires a separately built, validated candidate and the existing checks. The
lymph configuration describes a different tissue and marker trajectory; no
pancreas endocrine features are inferred for it.
