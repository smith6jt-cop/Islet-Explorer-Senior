# Thymus pipeline — full execution walkthrough

> **2026-06-09 — driver & post-processing superseded.** The thymus arm now runs through
> the single consolidated driver:
> `STAGES="01 02 phenotype 03 final" TISSUE=thymus bash scripts/run_pipeline.sh`.
> The old `run_thymus_post_processing.sh` (finalize → refine_rare → lineage →
> resolve_indeterminate → zonation) has been **REMOVED**: its middle three steps are
> SUBSUMED by `hierarchical_phenotype.py` (the `phenotype` stage writes an
> Indeterminate-free `celltype` + `celltype_lineage` and scores state programs), and
> zonation is the thymus `final` stage. Re-running that old chain on identity_core labels
> would corrupt them. The **Phase B/C commands below are historical**; use `run_pipeline.sh`.
> The scVI / zonation **parameter values remain accurate**.

**Scope:** the thymus arm of the Xenium analysis (samples **0041435** and **0041534**),
documented *as actually executed* and verified against the live notebooks/scripts on
**2026-05-25**. Where this document and the older notes in `CLAUDE.md` / `HANDOFF.md`
disagree, **this document is authoritative for thymus** (CLAUDE.md is pancreas-centric and
has drifted on several parameters — see "Authoritative values" below).

| Sample | Donor sex | Cells (post-QC) | Notes |
|---|---|---:|---|
| 0041435 | female | ~542,055 | medulla-rich section |
| 0041534 | male | ~302,056 | cortex-dense / partially involuted section |

Both: 10x Xenium 5k `hMulti_100g` panel, FFPE, probe-based, chemistry v2. Tissue dispatch is
driven by `XENIUM_TISSUE=thymus`, which makes cell 1 of every notebook import
`scripts/tissue_markers/thymus.py` as `_markers`.

## Authoritative values (verified against live code; supersede conflicting older notes)
| Item | Older note said | Live value |
|---|---|---|
| scVI `n_latent` | 10 | **30** (`n_layers=2`, `gene_likelihood='nb'`, `dropout_rate=0.1`) |
| scVI `max_epochs` / patience | 100 / 200 | **400 / 40** |
| scVI-graph `n_neighbors` | 30 | **10** (thymus UMAP retune folded into nb02) |
| stage-01 & stage-02 UMAP | `min_dist=0.5` | **`min_dist=0.05, spread=2.5`** |
| zonation lobules (0041435 / 0041534) | 6/1 or 25/25 | **29 / 16** (live `medulla_zones.csv`) |
| post-processing | manual / unclear | **`scripts/run_thymus_post_processing.sh`** (5 ordered steps) |

---

## Update — 2026-05-27 — hierarchical phenotyping promoted to canonical

Both thymus `{s}_phenotyped.h5ad` files have been merged in-place with the
two-tier hierarchical (subpanel-expression) categorization from
`scripts/hierarchical_phenotype.py --tissue thymus --output inplace
--reuse-scvi --dl-workers 0` (commit `6c6c5de`). The labels downstream
consumers should use are now:

- `obs['celltype']` — **hierarchical** FINE label (17/17 thymus types
  present; Indeterminate-free; built by per-cell subpanel-expression
  argmax + kNN denoise → broad lineage, then per-lineage sub-PCA +
  sub-leiden + tight-`MARKER_DICT` argmax → fine subtype). See CLAUDE.md
  decision #16 for the algorithm.
- `obs['celltype_broad']` — broad lineage (Thymocyte, TEC, Myeloid,
  B_lymphoid, NK, etc.; see `scripts/tissue_markers/thymus.py::
  BROAD_TO_FINE`).
- `obs['celltype_flat']` — the prior flat-scheme `celltype` (the
  Phase B / Phase C output described below). Preserved by the in-place
  merge for audit; **not** the recommended label for new analyses.
- `obsm['X_scvi_broad']` / `obsm['X_umap_broad']` — broad-tier scVI
  latent + UMAP from the hierarchical step.
- `obsm['X_scvi']` / `obsm['X_umap']` / `obs['leiden_scvi']` — the
  flat-scheme embeddings are preserved unchanged.

The flat Phase B + Phase C flow described below remains accurate for how
the flat-scheme label (`celltype_flat`) and the flat scVI embedding
(`X_scvi`/`X_umap`) were produced — that's the input the hierarchical
step refines. The thymus-specific rules in Phase C (DN/DP/SP4/SP8/Treg
boolean subtyping, AIRE_mTEC / Hassalls / pDC rule-based detection) are
**not** re-applied after a hierarchical merge — they are encoded in the
`MARKER_DICT` + sub-leiden argmax of `hierarchical_phenotype.py` via
`scripts/tissue_markers/thymus.py`. Post-processing
`resolve_indeterminate.py` is no longer needed (hierarchical emits
Indeterminate-free `celltype`); `lineage_phenotype.py` and the panel
mask step still run.

Pre-merge backup: `data/processed/{s}/{s}_phenotyped.prepilot.bak.h5ad`.

---

## Phase A — Ingest (`scripts/ingest_thymus.py`)
Xenium output bundle → `data/raw/{sample}.zarr` (spatialdata), loading cell + nucleus
boundaries and transcripts. One-time per sample. **Logic:** decouples the slow Xenium→zarr
conversion from analysis so stages 01–03 start from a fast random-access store.

## Phase B — Per-sample stages via `run_local_pipeline.sh` (TISSUE=thymus)
Launch:
```bash
TISSUE=thymus STAGES="01 02 03" bash scripts/run_local_pipeline.sh
```
The driver sets thymus defaults: `SAMPLES="0041435 0041534"`, **`XENIUM_MIN_COUNTS=20`,
`XENIUM_MIN_GENES=10`** (vs pancreas 50/20 — thymus median is ~107–124 counts/cell; the
pancreas floor would cut ~20–25% of legitimate small thymocytes). It forwards
`XENIUM_TISSUE=thymus` to every `nbconvert`.

### Stage 01 — `01_preprocessing_v2.ipynb` → `{sample}_preprocessed.h5ad`
1. **Ingest zarr** → AnnData; **drop control probes** (regex `NegControlProbe_|UnassignedCodeword_|NegControlCodeword_|antisense_|BLANK_`). *Logic:* controls must not contaminate QC/normalization.
2. **QC filter** — lower bounds absolute (`min_genes=10`, `min_counts=20`, a noise/segmentation floor); **upper bound is DENSITY-based**: drop cells with `total_counts/cell_area > q98`, not absolute counts. *Logic:* in FFPE Xenium `cell_area` varies 5–10× across cell types, so an absolute counts cap would predominantly remove large *healthy* cells; density catches genuine doublets/merges. `min_cells=100` on genes.
3. **Panel masks (set inline from `panel_audit_thymus.csv`):** `panel_for_scoring` (drops tissue-restricted genes) and `panel_for_embedding` (also drops 6 sex-chromosome genes). *Logic:* genes biologically implausible in thymus, and Y-chromosome genes, must not drive the embedding.
4. **donor_sex inference:** KDM5D/DDX3Y detected (≥`log1p(1)`) in ≥2% of cells → `male`, else `female` (→ 0041435 female, 0041534 male). *Logic:* feeds the scVI categorical covariate so the M/F pair is not confounded.
5. **counts layer** → **HVG** (`seurat_v3`, `n_top_genes=2000`, **metadata only — never used to subset features**) → **`normalize_total` + `log1p`** → **`lognorm` layer**.
6. **GPU scale** (`zero_center=False, max_value=10`) + **GPU PCA** `n_comps=50` on the `panel_for_embedding` subset. 50 PCs are computed; **30 are used downstream**. Sanity check: `Spearman(PC1, KDM5D/DDX3Y)` must be < 0.1. (No gene-space scree figure — that diagnostic was a z-scored-noise-floor artifact and is replaced by the stage-02 `scvi_latent_diagnostic.png`, which diagnoses the latent phenotyping actually uses.)
7. **GPU neighbors** `n_neighbors=10, n_pcs=30, cosine` → **CPU UMAP** `min_dist=0.05, spread=2.5` → **CPU leiden** (igraph) `resolution=1.5` → `leiden_1.5`. *Logic:* `leiden_1.5` is a diagnostic/QC clustering only; final cell types come from stage-02's scVI clusters. (UMAP/leiden stay on CPU per CLAUDE.md decision #4.)
8. **Spatial:** Delaunay graph (full sample); a 200k subsample for graph-stats. Centrality and Moran's I are skipped here (decorative / infeasible at 1.2M cells); co-occurrence runs on the subsample.

Figures: `qc_*`, `highly_variable_genes.png`, `umap.png`,
`spatial_scatter.png`, `co_occurrence_and_spatial.png`.

### Stage 02 — `02_phenotyping.ipynb` → `{sample}_phenotyped.h5ad`
Loads the preprocessed h5ad. Marker vocabulary is `_markers.MARKER_DICT` (21 thymus types:
thymocyte DN/DP/SP4/SP8/Treg; cTEC/mTEC; B_cells/Myeloid/DC/NK/Mast; Endothelial/Lymphatic_EC/
Fibroblast/Smooth_muscle/Schwann/Proliferating). AIRE_mTEC, Hassalls and pDC are **not** scored
here — they are rule-based, added in post-processing. Panel-aware substitutions handle genes
absent from this panel (CD3D→CD3E/CD3G/TRAC/TRBC2; VWF→PECAM1/CDH5/…; COL1A1→COL5A1/COL6A1/…);
every scored panel has ≥3 panel-present genes (decision #12).

1. **Score** each type with `sc.tl.score_genes` on `lognorm`.
2. **Argmax → `predicted_celltype`.** Legacy top-2 cell-type margin saved as `celltype_confidence_margin_legacy`.
3. **Anti-`DP_thymocyte` override** (the thymus analog of the pancreas anti-acinar override; DP is ~50–70% of thymus). A DP-called cell whose per-cell lognorm mean over a curated non-thymocyte alt panel (cTEC/mTEC/Myeloid/DC/Endothelial/Fibroblast/Smooth_muscle/B/NK/Mast) exceeds that type's baseline is reassigned. **T-lineage protection:** cells with `CD3E` or `TRAC` ≥ `log1p(1)` are *immune* to the override — prevents FFPE neighbor-contamination from flipping confirmed thymocytes.
4. **Lineage-margin confidence:** `celltype_confidence` is recomputed as the top−second **lineage** score (lineages from `_markers.LINEAGE_GROUPS`). *Logic:* DP/SP4/SP8 share CD3E/CD3G/TRAC/TRBC2, so the per-celltype top-2 margin was ≈0 for every clear thymocyte. The within-lineage margin is preserved as `subtype_confidence`.
5. **Indeterminate flagging** for `LOW_CONFIDENCE_TYPES` (Mast, Lymphatic_EC, Schwann, Treg, Smooth_muscle) where `confidence < median − 1·MAD`; routed into `_markers.INDETERMINATE_BUCKETS`. (Resolved later — never a final label, decision #15.)
6. **scVI** (notebook section 7): `setup_anndata(layer='counts', continuous_covariate_keys=['cell_area'], categorical_covariate_keys=['donor_sex'])`; `SCVI(n_layers=2, n_latent=30, gene_likelihood='nb', dropout_rate=0.1)`; `train(max_epochs=400, early_stopping=True, early_stopping_patience=40)` → `X_scvi`. *Logic:* NB likelihood models library size via a size factor (so `total_counts` is deliberately **not** a covariate, decision #1); `cell_area` absorbs segmentation-size effects; `donor_sex` absorbs the M/F batch.
7. **Embedding on scVI:** GPU neighbors `use_rep='X_scvi', n_neighbors=10, cosine` (`key_added='scvi'`) → CPU UMAP `min_dist=0.05, spread=2.5` → leiden `resolution=1.0` → `leiden_scvi`. A PCA-based fallback (`leiden_pca`, `n_pcs=30`) mirrors this. The aggressive UMAP params are the thymus retune, now baked into the notebook; `scripts/retune_umap_thymus.py` remains only to backfill pre-2026-05-22 h5ads. Deferred figures `predicted_celltypes.png` / `scvi_vs_markers.png` render here with Indeterminate cells hidden. `scvi_latent_diagnostic.png` (scVI-latent scree + per-dim usage + library-size leak gate) also renders here — the diagnostic for the representation phenotyping actually uses; re-render via `scripts/render_scvi_latent_diagnostic.py`.
8. **Cluster consensus** (section 8): `manual_celltype` = the majority `predicted_celltype` within each `leiden_scvi` cluster. **Thymocyte-preservation guard:** for any cell whose `predicted_celltype` ∈ {DN/DP/SP4/SP8/Treg}, restore the per-cell label — cluster-majority would otherwise collapse DP cells into SP-dominant cortical clusters. The 3-panel consensus UMAP uses tissue-dispatched `LINEAGE_MAP` / `SPECIFIC_KEEP` / `PALETTE` / `THIRD_PANEL` (thymus panel-3 drops the bulk DP+DN so rare types are visible).
9. **Finalize** (sections 9–10): `celltype = manual_celltype`; `celltype_lineage` via `LINEAGE_CANARY_PANELS` composite argmax. Save `{sample}_phenotyped.h5ad` and `{sample}_celltype_assignments.csv` (a stage-02 snapshot of `celltype, celltype_confidence, celltype_lineage` — later overwritten post-resolve).

### Stage 03 — `03_spatial_analysis.ipynb` → `{sample}_spatial_analysis.h5ad`
The tissue-shared squidpy stage, run on the phenotyped labels:
`spatial_neighbors` (Delaunay, k=10); `nhood_enrichment` (`n_perms=100, n_jobs=1`);
co-occurrence (150k subsample); Moran's I (`rsc` GPU, 200k-subsample fallback);
**spatial niches** (cuML kNN k=50 + cuML KMeans `n_clusters=12` + majority-vote smoothing);
**ligrec & Ripley grouped by `_markers.LINEAGE_MAP`** (7–8 thymus lineages —
Thymocyte/TEC/Myeloid/Other_Immune/Vascular/Stromal/Neural/Cycling — because rare types have
too few cells for stable per-type statistics), 100k subsample, `n_perms=50`.
Outputs: `*_neighborhood_enrichment.csv`, `*_morans_i.csv`, `*_spatial_niches.csv`,
`*_ligrec_lineage_{means,pvalues,metadata}.csv` + figures.

## Phase C — Post-processing via `run_thymus_post_processing.sh` (5 ordered steps)
Run after stages 01/02/03 finish; each step writes back into the phenotyped (and spatial)
h5ads.
1. **`finalize_panel.py --audit-csv panel_audit_thymus.csv`** — (re)applies the `panel_for_scoring` / `panel_for_embedding` var-masks and `donor_sex` (decisions #13/#14).
2. **`refine_rare_celltypes.py --tissue thymus`** — the *authoritative* cell-type refinement. It resets the working `celltype` to `predicted_celltype` (bypassing the cluster-consensus collapse), then:
   - **`apply_thymocyte_subtyping`** — the FINAL DN/DP/SP4/SP8/Treg calls, made by **per-cell boolean detection** at `log1p(1)`:
     `Treg` = FOXP3+ & CD4+ · `DP` = CD4+ & CD8+ · `SP4` = CD4+ & CD8− & CD3E+ ·
     `SP8` = CD8+ & CD4− & CD3E+ · `DN` = CD4− & CD8− & CD3E+ & (CD7|CD44).
     *Logic:* score_genes argmax underweights the canonical CD4+CD8+ DP definition because CD4/CD8 also seed the SP panels.
   - **`apply_rare_celltype_rules`** — rule-based rare types (not argmax):
     **AIRE_mTEC** (AIRE & KRT14 present, no strong T signal, currently TEC-like);
     **Hassalls** (KRT1 & (TGM1|DSG1), not cTEC, TEC-like);
     **pDC** (LILRA4 & ≥1 of IRF7/IRF8/CLEC4C/TCF4, not T, DC/myeloid-like).
     *Logic:* thresholds relaxed to `log1p(1)` because these markers (AIRE ~0.6–1%, LILRA4 ~0.3–0.4%) rarely exceed `log1p(5)` in FFPE.
3. **`lineage_phenotype.py --tissue thymus`** — per-cell canary-panel lineage depth (`LINEAGE_CANARY_PANELS`); sets `celltype_lineage` and flags `doublet_suspected` when ≥3 mutually-exclusive lineages co-detect. *Logic:* post-density-filter doublet refinement; `celltype` is never overwritten.
4. **`resolve_indeterminate.py --target both --write --no-downstream`** — propagates real labels onto every `Indeterminate*` cell via a three-tier fallback (label-suffix → `predicted_celltype` → scVI-latent kNN). *Logic:* Indeterminate must never be a final label (#15).
5. **`thymic_zonation_analysis.py`** — the stage-07 zonation (below).

*Historical / backfill only (not in the orchestrator; the logic now lives in nb02):*
`regenerate_celltype_confidence.py` (lineage-margin confidence) and `retune_umap_thymus.py`
(aggressive UMAP) — retained to fix pre-2026-05-22 h5ads.

## Stage 07 — Thymic zonation (`thymic_zonation_analysis.py`)
Standalone; reads `{sample}_phenotyped.h5ad` (uses `celltype_lineage` if present, else
`celltype`). Replaces the old DBSCAN "bullseye" method, which fragmented sparse
marker-positive seeds into dozens of tiny concentric rings that did not match thymic anatomy.

**kNN-density method:**
1. **Source classes by label** — cortex = {cTEC, DP_thymocyte}; medulla = {mTEC, AIRE_mTEC, Hassalls, SP4_thymocyte, SP8_thymocyte, Treg, DC, pDC}; everything else neutral (dilutes both denominators equally).
2. **Spatial kNN** `k=50` (cuML brute, sklearn fallback) on `obsm['spatial']` (~150 µm neighborhood).
3. **Two-pass density** — per-cell `cortex_frac` / `medulla_frac` over the k neighbors, then averaged once more over the same graph. *Logic:* denoises salt-and-pepper without smearing across lobule borders.
4. **Signed dominance zone** — `score = medulla_frac − cortex_frac`: `> 0.05` → medulla; `< −0.10` → cortex; between → cortico_medullary_junction; `max(frac) < 0.15` → cortex (off-tissue default). *Logic:* the signed score is symmetric and robust to the two samples' very different medulla:cortex ratios; an absolute per-class threshold would systematically under-call medulla in the cortex-dense 0041534.
5. **Connected components** on the medulla-restricted kNN graph; components `< 200` cells are demoted to peri_medullary, leaving true lobule-scale regions. **Live result: 29 lobules (0041435), 16 lobules (0041534).**

**Outputs** — `data/processed/medulla_zones.csv` (per-lobule centroid/size/extent),
`thymic_zone_composition.csv` (zone × celltype fractions),
`thymocyte_maturation_gradient.csv` (DN→DP→SP→Treg per zone),
`AIRE_neighborhood.csv` (cell composition within 50 µm of each AIRE_mTEC).
Figures in `figures/{sample}/07_thymic_zonation/`: `medulla_spatial_overview`,
`medulla_lobules_overview`, `zone_composition_heatmap`, `thymocyte_gradient_per_zone`,
`aire_neighborhood_composition`.

CLI tunables (defaults shown):
```bash
python scripts/thymic_zonation_analysis.py --samples 0041435 0041534 \
  --knn 50 --score-medulla 0.05 --score-cortex -0.10 \
  --presence-floor 0.15 --smooth-passes 2 --min-medulla-cells 200
```

## Phase D — Sample QC: `diagnose_0041534_zonation.py`
0041534's medulla/cortex appear spatially "inverted" (medulla on the outer rim). Three tests
(marker-vs-label, spatial coherence, zone-vs-marker) return a **verdict of real biology**, not
a labeling error: cortical involution/hypoplasia in the older male donor — cTEC marker
expression is globally depleted (~½ of 0041435), and residual DP thymocytes pack into small
cortical remnants (blue core) ringed by the surviving mTEC/Treg signal (red rim). Emits
`figures/0041534/07_thymic_zonation/diagnose_inversion.png`.

---

## Parameter quick-reference (thymus)
QC `min_counts=20 / min_genes=10 / density-q98 / min_cells=100` · HVG 2000 (metadata only) ·
PCA 50 computed / 30 used · neighbors `k=10, cosine` · stage-01 `leiden_1.5` ·
scVI `n_latent=30, gene_likelihood='nb', cell_area + donor_sex, max_epochs=400, patience=40` ·
scVI UMAP `min_dist=0.05, spread=2.5`; `leiden_scvi` resolution 1.0 ·
thymocyte subtyping @ `log1p(1)` · rare-type "present" gates @ `log1p(1)`, exclusion gates @ `log1p(5)` ·
niches KMeans `n_clusters=12` · ligrec/Ripley 100k subsample, `n_perms=50`, grouped by `LINEAGE_MAP` ·
zonation `k=50`, `score>0.05 / <−0.10`, min 200 cells/lobule.

## Known caveats (live)
- **0041534 library-size leak (G2):** scVI-UMAP Spearman ≈ (+0.11 on X, −0.34 on Y) — the Y-axis is slightly over the 0.3 target; non-blocking, because cell types come from marker panels, not from `leiden_scvi` directly.
- **0041534 "inversion" is real biology** (see Phase D), not an artifact.
- **CLAUDE.md is pancreas-centric** and has historically drifted on thymus-specific
  parameters (this file remains the thymus reference of record). The 2026-05-28 doc
  sweep refreshed CLAUDE.md decision #16 to the current subpanel-expression
  hierarchical scheme that now applies to both tissues; older drift on stage-01/02
  UMAP `min_dist=0.05/spread=2.5` was already corrected in the 2026-05-25 sweep.

## How to verify any claim here
```bash
# scVI / embedding params (ground truth in the notebook)
grep -nE "n_latent|max_epochs|min_dist|spread|n_neighbors|resolution=" notebooks/02_phenotyping.ipynb
# post-processing order
sed -n '30,64p' scripts/run_thymus_post_processing.sh
# thymus markers, rules, lineage maps
less scripts/tissue_markers/thymus.py
# live zonation lobule counts per sample
awk -F, 'NR>1{print $1}' data/processed/medulla_zones.csv | sort | uniq -c
```
