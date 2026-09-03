# Data layout

Project-specific reference. Generic Xenium data-loading docs live in the [10x
Xenium documentation](https://www.10xgenomics.com/support/in-situ-gene-expression).

## What this project actually has

Two human pancreas FFPE samples on the **hAtlas v1.1 + 100-gene custom T1D
add-on** panel (Xenium 5K base + custom probes). Both samples are FFPE,
probe-based, chemistry v2.

| Sample | Cells (post-QC) | Donor sex | Notes |
|---|---:|---|---|
| 0041323 | ~1,145,295 | female | (KDM5D ~0%) |
| 0041326 | ~1,168,473 | male | (KDM5D ~10%, more advanced T1D progression by IAPP loss) |

`donor_sex` is auto-derived from KDM5D + DDX3Y detection — see
`scripts/finalize_panel.py`.

## Directory structure

```
data/
├── raw/
│   ├── 0041323.zarr/        # spatialdata zarr ingested via notebook 00
│   └── 0041326.zarr/
├── processed/
│   ├── 0041323/
│   │   ├── 0041323_phenotyped.h5ad         # PRIMARY output of stage 02
│   │   ├── 0041323_spatial_analysis.h5ad   # PRIMARY output of stage 03
│   │   ├── 0041323_immune_phenotyped.h5ad  # gated immune subset (notebook 07)
│   │   └── *.csv                           # per-stage tables
│   ├── 0041326/                            # same layout
│   ├── panel_audit.csv                     # candidate panel exclusions
│   ├── islet_*.csv                         # insulitis headline tables
│   ├── immune_proximity_summary.csv
│   └── legacy/                             # archived pre-rebuild artifacts
└── xenium_explorer/
    ├── {sample}/analysis.zarr.zip          # drop next to experiment.xenium
    ├── {sample}_cell_groups.csv            # Cell → Add cell categorization
    └── {sample}_color_palette.json
```

## h5ad layer / obs / obsm contracts

`*_phenotyped.h5ad` (stage 02 + hierarchical merge output) MUST contain:

- `layers['counts']` — raw integer counts
- `layers['lognorm']` — `log1p(normalize_total)` (sample-comparable scale)
- `obsm['spatial']` — (n_cells, 2) xy in microns
- `obsm['X_scvi']` — flat-scheme scVI latent (30 dims; preserved across merge)
- `obsm['X_umap']` — flat-scheme scVI-derived UMAP (preserved across merge)
- `obsm['X_scvi_broad']` — broad-tier scANVI latent from hierarchical
  phenotyping (added by `scripts/hierarchical_phenotype.py`); since the
  2026-06-08 redesign this is the scANVI semi-supervised CLASSIFIER latent,
  trained on `var['identity_panel']`
- `obsm['X_umap_broad']` — diagnostic broad UMAP (subsample, in `uns['umap_broad']`)
- `obs['celltype']` — **hierarchical** final fine label (Indeterminate-free;
  downstream consumers MUST use this). Since 2026-06-08 the broad lineage is the
  scANVI `predict()` argmax (not a marker argmax); fine subtype is detection-gated
- `obs['celltype_flat']` — the prior flat-scheme `celltype` label
  preserved by the in-place hierarchical merge; for audit/comparison only
- `obs['celltype_broad']` — hierarchical broad lineage (Endocrine /
  Exocrine / Stromal / Vascular / Immune_* / Myeloid / Neural for
  pancreas; Thymocyte / TEC / Myeloid / B_lymphoid / NK / etc. for
  thymus; see `scripts/tissue_markers/{tissue}.py::BROAD_TO_FINE`)
- `obs['celltype_lineage']` — FINE lineage label (set by
  `hierarchical_phenotype.py`, = `celltype`); never carries broad labels
- `obs['leiden_scvi']` — flat-scheme leiden cluster ID (preserved)
- `obs['leiden_broad']`, `obs['leiden_sub_<L>']` — hierarchical broad +
  per-lineage sub-clusterings
- `obs['broad_confidence']` (= `broad_confidence_prob`) — scANVI **max posterior
  ∈ [0,1]** for the broad call (redesign 2026-06-08; was a score-margin). Reported,
  never gated. `obs['broad_marker_margin']` — top1−top2 posterior margin (figure).
  `obs['subtype_confidence']` — fine score-margin. `obs['broad_lowsignal']` — cell
  was seeded 'Unknown' (placed by the latent, not by a marker seed).
- `obs['state_<group>_score']` — Tier-3 state/program module scores (IFN_response,
  Cytokine_signaling, Matrix_remodeling, Proliferation, ER_stress, the pathway
  panels, pancreas `38_insulin_biology_T1D`, …) from `scripts/state_panels.py`.
  Orthogonal — NEVER used for the identity call (decision #17 GUARD); for
  QC/annotation/downstream only.
- `obs['donor_sex']` — `'male'` / `'female'` (added by `scripts/finalize_panel.py`)
- `var['panel_for_scoring']`, `var['panel_for_embedding']` — panel-audit boolean masks
- `var['identity_panel']` — the scANVI broad-classifier FEATURE set (trimmed
  cell-type-panel union ∩ embedding; pure-state/pathway genes excluded — the GUARD).
  `var['is_state_only']` — audit: gene only in state/fine_tier panels.
  `var['identity_promoted']` — gene promoted into an identity mask (thymus AIRE/FEZF2).
  `var['subpanel_<L>']` — per-lineage identity mask, trimmed of cross-lineage genes
  (`scripts/tissue_markers/_identity_trim.py`). `var['celltype_panel']` — untrimmed
  union (bookkeeping). All stamped by `scripts/finalize_panel.py` (decisions #16/#17).
- `var['celltype_panel']`, `var['subpanel_<L>']` — hierarchical
  subpanel masks (stamped by `finalize_panel.py`)

Pre-merge backups at `data/processed/{s}/{s}_phenotyped.prepilot.bak.h5ad`
hold the flat-only schema and can be removed once the refreshed figures
are confirmed.

`*_spatial_analysis.h5ad` (stage 03 output) adds:
- `obs['spatial_niche']` — KMeans cluster ID
- `uns['nhood_enrichment']`, `uns['co_occurrence']`, `uns['ligrec']` — squidpy outputs

## Spatial-method benchmark outputs (`data/processed/benchmark/`)

Produced by `scripts/benchmark/` (BANKSY/UTAG/KODAMA/MISTy + niche baseline on
`0041323`+`0041435`). Per `(sample, variant)` — `variant` ∈
`{crop25k, crop50k, thin50k, crop100k, crop200k, full}`, headline `crop50k`:

```
benchmark/{sample}/
├── {s}_bench_{v}_input.h5ad        # shared Python input: X=lognorm on var['embedding_panel'],
│                                   #   layers counts+lognorm, obsm spatial+X_scvi
├── {s}_bench_{v}_{expr.mtx,coords.csv,meta.csv,barcodes.csv,genes.csv}  # identical R bridge
├── {s}_bench_{v}_misty_markers.txt # tissue MARKER_DICT union (MISTy targets)
├── {s}_bench_{v}_cropinfo.json     # n_cells, gene_panel, composition (reported, not gated)
├── results/{method}_{v}_labels.csv # cell_id,domain (banksy/utag/kodama/baseline)
├── results/{method}_{v}_metrics.json  # wall_s, peak_rss_mb, n_domains, status (profile_wrap)
├── results/misty_{v}/ + misty_{v}_*.csv  # MISTy per-target results + tidy improvements/importances
├── {s}_benchmark_scorecard_{v}.csv # master scorecard (rows=methods)
└── {s}_benchmark_report_{v}.md     # human-readable report + bias notes
```
Figures: `figures/{sample}/benchmark/*.png`. All identical-input across methods
(parity-asserted). REPORT-ONLY (Tier-C): the scorecard tabulates objective
metrics; no proportion-plausibility verdict. See `scripts/benchmark/README.md`.

## Headline CSV outputs

| File | Owner script | Schema |
|---|---|---|
| `panel_audit.csv` | `scripts/audit_panel_exclusions.py` | gene-by-gene flagging w/ exclude category |
| `islet_infiltration_per100endo.csv` | `scripts/insulitis_analysis.py` | per-islet rows; per-phenotype zone counts, rates, grades; composition |
| `islet_insulitis_grades.csv` | `scripts/insulitis_analysis.py` | long format: per-(sample × subtype) grade prevalence |
| `islet_insulitis_thresholds.csv` | `scripts/insulitis_analysis.py` | per-(sample × size_class × subtype) rotation-null thresholds |
| `islet_insulitis_regression.csv` | `scripts/insulitis_analysis.py` | per-(sample × subtype) odds ratios, no p-values (n=2) |
| `immune_proximity_summary.csv` | `scripts/insulitis_analysis.py` | density enrichment per immune subtype × sample |

`scripts/insulitis_analysis.py` is the **single writer** of the islet/immune
headline CSVs; it consumes `{s}_immune_phenotyped.h5ad` (written by notebook 07).
Both run as the pancreas `final` stage of `run_pipeline.sh`.

## Xenium Explorer round-trip

`utils/xenium_explorer_export.py::export_for_xenium_explorer(adata,
group_keys=[...], output_dir=..., sample_name=...)` writes both formats and a
color palette JSON. The CLI wrapper is
`scripts/export_xenium_explorer_groups.py`.

## Storage footprint per sample

- raw zarr: ~5–10 GB
- `_phenotyped.h5ad`: ~60 GB
- `_spatial_analysis.h5ad`: ~10 GB
- `_immune_phenotyped.h5ad`: ~1 GB

## Data privacy

Human pancreas tissue from de-identified donors (consent on file). Standard
HIPAA/IRB practices apply; no PHI in any committed file.
