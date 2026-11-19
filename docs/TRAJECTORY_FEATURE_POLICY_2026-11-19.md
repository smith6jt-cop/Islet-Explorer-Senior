# Islet composition graph inputs: 2026-11-19

The islet-composition producer selects graph inputs separately from displayed
phenotypes. QC/ambiguous labels (including unresolved and conflict labels)
do not enter either graph. Endocrine classes enter the core graph only;
exocrine classes enter the peri-islet graph only. Other supported named
classes may enter either region. When multiple measured macrophage marker
profiles identify Gal3+ or Gal3−, their counts are pooled by sign before the
transformation. Missing data are not imputed, and an aggregate count supplied
alongside its component counts is refused rather than double-counted.
Senior's coarse Immune class is a marker-supported lineage, not the CODEX
unresolved immune residual, and is retained when it is the only immune class
in the table. If broad Immune and named subtype columns are both supplied,
their overlapping coverage is unknown from column names alone; the producer
refuses to build until a non-overlapping measured schema is supplied.

The composition graph has one column per measured class. Endocrine counts come
from core only and exocrine counts from peri only. All other retained counts are
summed across regions before log1p(count / max(core endocrine cells, 1)).
The accepted CODEX table provides 14 features. Senior uses its own measured
classes, cohort, Euclidean metric, root choice and validation gates. The H5AD
records feature names in trajectory_pooled_inputs and the fitted matrix in
X_composition_pooled. Both composition pseudotime fields use the pooled fit.

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
