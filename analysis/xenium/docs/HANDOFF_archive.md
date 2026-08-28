# HANDOFF archive — sessions pre-2026-05-27

This file holds session notes that were prepended to `HANDOFF.md` before
the 2026-05-27 hierarchical-phenotyping promotion. Split off 2026-05-28
to keep the live handoff under 40 KB.

For current pipeline state and how to resume, see [`../HANDOFF.md`](../HANDOFF.md)
at the repo root. That file holds the 2026-05-27 promotion / validation /
broad-tier-rationale sessions, the durable "Current state / Resume /
Verify / Failure" reference, and the "Files modified today (2026-04-25)"
note.

What's in this archive:
- The 2026-05-26 broad-tier two-tier pilot, PCA-scree replacement, and
  predicted_celltypes major-lineage rewire session notes.
- The 2026-05-25 UMAP overhaul / HVG honesty / PCA-elbow session notes.
- The full 2026-04-27 → 2026-05-22 historical tail (stage 07 immune,
  Epsilon mislabel + Acinar hotfix, niche-pair ligrec, thymus 5-fix +
  zonation rewrite + UMAP retune, full rebuild + inversion diagnostic).

These sections are not loaded by Claude on session start; they exist
only for historical / forensic lookup.

---

## 2026-05-26 — two-tier hierarchical phenotyping (subpanel-driven) — PILOT

User flagged that the scVI UMAPs fail to separate cell types into well-formed
clusters *before* subtyping. Root cause is architectural, not a bug: the flat
scheme embeds + clusters + argmaxes over the full undifferentiated panel, so
broad lineages never separate first (and the flat argmax is what produced the
historical Epsilon-150K / AMY1A-overspill mislabels). Diagnosis + approved plan:
`/home/smith6jt/.claude/plans/vivid-hugging-wind.md`.

Built this session (pancreas pilot; thymus is a follow-on):
- `data/reference/subpanels/` — vendored cell-type-category subpanels pinned to
  XeniumPanelExplorer **@ 79c5ab9** (see `SOURCE.txt`). 12 lineage files + 1
  embedding-only (`16_epithelial_general`). Pathway/program panels excluded.
- `scripts/tissue_markers/pancreas.py` — APPENDED `SUBPANEL_ALLOWLIST`,
  `BROAD_TO_FINE`, `BROAD_TO_LINEAGEMAP`, `BROAD_PALETTE`, `N_MIN_SUBCLUSTER`,
  `load_subpanel_genes()`. Existing structures untouched.
- `scripts/finalize_panel.py` — `annotate_subpanel_masks()` writes
  `var['subpanel_<L>']` + `var['celltype_panel']` (union ∩ `panel_for_embedding`);
  new `--tissue` / `--subpanel-dir` / `--no-subpanels`.
- `scripts/hierarchical_phenotype.py` — **NEW.** Tier-1 scVI on `celltype_panel`
  → `leiden_broad` (res 0.3) → cluster-mean broad-score argmax → `celltype_broad`.
  Tier-2 per-lineage sub-PCA+leiden (res 0.5) → tight-`MARKER_DICT` fine argmax;
  within-Exocrine Acinar→Ductal refinement; Endocrine `apply_rare_celltype_rules`
  (Gamma/Epsilon). **Subsumes `fix_acinar_argmax.py` + the rare pass of
  `refine_rare_celltypes.py`** → emits a final, Indeterminate-free `celltype`, so
  refine/resolve are NOT run after it.
- `verify_lineage_columns.py` (+new cols), `slim_phenotyped.py` + notebook 03
  `_drop_obsm` (drop `X_scvi_broad`/`X_umap_broad`).

Contract preserved: canonical `obs['celltype']` stays FINE; **`obs['celltype_lineage']`
stays FINE** (insulitis + stage-07 `isin` it against fine endocrine names — never
put broad labels there). New cols: `celltype_broad`, `leiden_broad`,
`leiden_sub_<L>`, `celltype_fine`, `celltype_hier`, `broad_confidence`,
`subtype_confidence`; new obsm `X_scvi_broad`/`X_umap_broad`. Old flat
`X_scvi`/`leiden_scvi`/`X_umap` preserved; old label saved as `celltype_flat`.

Run the pilot (non-destructive **slim sidecar** by default — review before the
60 GB in-place merge + downstream rerun):
```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
source scripts/env_rapids.sh
python scripts/finalize_panel.py --tissue pancreas         # stamp subpanel masks
python scripts/hierarchical_phenotype.py --tissue pancreas --samples 0041323 0041326
```
→ writes `data/processed/{s}/{s}_hierarchical.h5ad` +
`figures/{s}/02_phenotyping/hierarchical_umap.png`. Validation (leak gate,
ARI, Acinar≈69% in 0041323, Epsilon<0.05%, 0041326 immune/ductal shrink) prints
at the end of each sample. To commit: re-run with `--output inplace`.

## 2026-05-26 — PCA scree replaced by scVI-latent diagnostic

User flagged that the `pca_variance_ratio.png` scree (a) diagnosed a
representation phenotyping never uses, and (b) showed alarmingly low cumulative
variance (thymus 0041534: PC1 0.55%, cum-at-30 5.2%). Both are real: phenotyping
runs on the **scVI latent** (`X_scvi` → cosine kNN → `leiden_scvi`/`X_umap`) +
`score_genes`, never `X_pca`; and the gene-space scree ran on z-scored ~5k genes,
so the total-variance denominator was dominated by the unit-variance noise floor
of the (mostly undetected) panel genes. In low-detection samples (0041534: median
90 counts/cell, 93% of genes detected in <5% of cells) that produced a near-flat
spectrum that *looked* like low intrinsic dimensionality but was a scaling
artifact. Proof: a PCA on the 30-dim `X_scvi` reaches 90% variance in **13 dims**
for 0041534 (26 for pancreas 0041323), cumulative → 1.0 for both.

**Change.** Removed `plot_pca_variance_ratio` + `scripts/render_pca_variance_plot.py`
+ the stage-01 cell-19 call + the four `figures/*/01_preprocessing/pca_variance_ratio.png`.
Added `utils.analysis_utils.plot_scvi_latent_diagnostic` (2-panel: scree of PCA on
`X_scvi` with effective-dim marker; per-latent-dim std bars + the library-size
Spearman leak gate) → `figures/{sample}/02_phenotyping/scvi_latent_diagnostic.png`,
wired into stage-02 cell 17 and a new `scripts/render_scvi_latent_diagnostic.py`.
Stage-01 `rsc.pp.pca` itself is kept (still feeds `X_pca`/stage-01 neighbors); only
the figure is gone. Rendered for 0041323 + 0041534. Note: 0041534's leak gate
reports `max|rho|=0.331` (borderline WARN ≥ 0.3) — a pre-existing condition the
diagnostic now surfaces, not introduced by this change. The historical entries
below referencing `plot_pca_variance_ratio`/`pca_variance_ratio.png` are superseded.

## 2026-05-26 — predicted_celltypes → major-lineage + stage-01 UMAP retune

Two user-requested changes. Pancreas only (0041323, 0041326); thymus untouched.

**1. `predicted_celltypes.png` / `scvi_vs_markers.png` collapsed to 7 major
lineages.** The fine 14–15-type palette was too scattered to read islet/exocrine/
immune separation. New shared helper `utils.analysis_utils.render_major_lineage_panel`
maps `celltype` (script) / `predicted_celltype` (notebook) through `LINEAGE_MAP`
to the 7 groups (Endocrine/Exocrine/Immune/Stromal/Vascular/Neural/Cycling) and
colors ONE panel with `PALETTE`. Dispatched per-tissue by a new
`scripts/tissue_markers/pancreas.py::PREDICTED_VIEW = "major_lineage"` flag;
tissues that omit it (thymus) fall back to the old 3-panel `LINEAGE_SUBSETS`
split (thymocyte DN/DP/SP detail must stay visible). Both the notebook
(cells 17/18) and `scripts/rerender_stage02_celltype_figs.py` branch on the flag,
so a fresh run and a re-render agree. Re-rendered both pancreas samples.

**2. Stage-01 preprocessed UMAP retuned (de-blobbed).** The on-disk
`_preprocessed.h5ad` (written 2026-04-28) was still the OLD blob — its
`uns['umap']['params']` decoded to `min_dist=0.5, spread=1.0` and neighbors used
`n_neighbors=15`; the 2026-05-25 retune had only touched the stage-02/03 h5ads.
New `scripts/retune_stage01_umap.py` reuses the existing (sex-excluded) `X_pca`
and rebuilds neighbors→UMAP→**leiden_1.5** to match notebook 01 cell 18 EXACTLY
(`n_neighbors=10, n_pcs=30, min_dist=0.05, spread=2.5, leiden res=1.5`). leiden
is recomputed because the neighbor graph changes (the required "downstream" step);
old `X_umap`/`leiden_1.5` preserved under `X_umap_blob`/`leiden_1.5_blob`, h5ads
backed up to `*.bak.preUmapRetune`. Both umap.png re-rendered. Result:

| sample | a/b (was 0.583/1.334 blob) | n_neighbors | leiden_1.5 | leak ρ (was blob) |
|---|---|---|---|---|
| 0041323 | 0.4206 / 0.8109 | 15→10 | 21 clusters | (+0.163, −0.027)  ← was (+0.219, −0.037) |
| 0041326 | 0.4206 / 0.8109 | 15→10 | 25 clusters | (+0.193, +0.349)  ← was (+0.235, +0.342) |

**Leak-gate note (important).** The stage-01 leak gate uses the SAME logic as
`retune_umap_thymus.py` — "don't worsen vs the blob baseline by >0.05 on either
axis + absolute ceiling 0.40" — NOT a fixed `<0.30` absolute. Reason: 0041326's
stage-01 PCA embedding already carried a Y-axis library-size correlation of
+0.342 in the blob (the PCA retains library-size structure that the scVI latent
— which drives the actual cell calls — does not). A fixed `<0.30` gate would have
aborted the de-blob (Y=0.349) even though it does NOT worsen the blob and is
under 0.40, leaving the worse blob in place. This matches the thymus 0041534
precedent (accepted at ~0.35). **0041326's stage-01 diagnostic UMAP has ~0.35
Y-leak by nature of its PCA; this is unchanged by the retune and does NOT affect
cell-type calls (those come from the clean scVI latent, leak ~0.14).**

**Downstream scope.** Stage-01's `leiden_1.5`/`X_umap` do NOT feed stage-02 cell
typing (stage 02 cell-types on the independent scVI latent). The `leiden_1.5`
column carried into the phenotyped/spatial h5ads is from the pre-retune stage-01
graph and is a stage-02 marker-discovery DIAGNOSTIC only — left as-is; stages
02/03 were NOT re-run and their cell calls are unaffected.

Files: new `scripts/retune_stage01_umap.py`; `utils/analysis_utils.py`
(+`render_major_lineage_panel`) + `utils/__init__.py` export;
`scripts/tissue_markers/pancreas.py` (+`PREDICTED_VIEW`);
`scripts/rerender_stage02_celltype_figs.py` + `notebooks/02_phenotyping.ipynb`
cells 17/18 (branch on the flag).

## 2026-05-25 — UMAP overhaul: pancreas retune + adaptive points + lineage-split panels

Three user-requested UMAP changes, applied to all four samples.

**1. Pancreas got the thymus UMAP params.** Ran the (now general-purpose)
`scripts/retune_umap_thymus.py 0041323 0041326`: rebuilt neighbors+UMAP from the
existing `X_scvi`/`X_pca` with `n_neighbors=10, min_dist=0.05, spread=2.5` (no
scVI retrain, no cell-call change). Library-size leak gate PASSED — leak ρ
0041323 (+0.080,−0.123), 0041326 (−0.055,−0.109), both improved on the X-axis;
old embeddings preserved under `X_umap*_blob`, h5ads backed up to
`*.bak.preUmapRetune`, UMAP keys mirrored into the spatial h5ads. ~55 min.
NOTE: like the thymus retune, this touched only the phenotyped/spatial h5ads
(stage-02/03). The stage-01 *preprocessed* UMAP embedding was NOT retuned — if
the stage-01 diagnostic UMAP should also use the new params, re-run stage 01.

**2. Adaptive point size everywhere** via new `utils.analysis_utils.umap_dot_size`
(`k=9e6, lo=10, hi=20`): full-population panels → 10, lineage-subset/rare-type
panels → up to 20, so panels look equally dense regardless of cell count. Shared
by the notebooks AND the re-render scripts (single knob for the whole pipeline).

**3. Lineage-split the remaining crowded celltype UMAPs.** `final_celltypes_overview.png`
now shows a de-crowded general-lineage panel (7–8 groups) instead of a ~15-25-type
crush; `predicted_celltypes.png` / `scvi_vs_markers.png` split into 3 lineage-subset
panels (already done in the script, now also adaptive); `cluster_consensus_annotation.png`
keeps its 3-panel split. New shared helpers: `render_celltype_subset_panels`,
`lineage_subsets_for`, `high_contrast_palette` (all in utils, exported).

**Ownership de-duplicated:** `rerender_stage02_celltype_figs.py` is now the SOLE
renderer of `final_celltypes_overview.png` (was also in `regen_umap_overlays.py`
at a different size — that script is now QC-overlay-only). New
`scripts/regen_stage01_umap.py` re-renders `umap.png` (the only UMAP with no
prior script) with adaptive size from a backed read.

**Notebook↔script parity:** nb01 cell 21 and nb02 cells 17/18/22 now route their
plotting through the same `utils` helpers (scVI-training and celltype_lineage
canary code preserved verbatim by slicing the cell source). A fresh notebook run
and a post-hoc re-render now produce the same figures.

Re-rendered all 4 samples (`rerender_stage02_celltype_figs.py`,
`regen_umap_overlays.py`, `regen_stage01_umap.py`) — celltype CSVs Indeterminate-free.

## 2026-05-25 — HVG plot honesty + doc/figure consistency sweep

Completes the HVG sibling of the PCA-elbow fix below and clears the
docs-vs-code drift the user flagged. **Labeling/docs/figures only — no
embeddings, h5ads, or cell calls changed.**
- HVG plot now states in-figure that highly-variable genes are QC/reference
  ONLY and do NOT subset features (PCA/scVI use `panel_for_embedding`). The
  disclaimer lives in the shared helper
  `utils.analysis_utils.plot_highly_variable_genes`, so the notebook and
  `scripts/render_hvg_plot.py` render an identical figure.
- `notebooks/01_preprocessing_v2.ipynb` cell 40 rewritten to CALL that helper
  (was hand-rolled inline with 5 labels + no disclaimer; now mirrors cell 19's
  PCA pattern — so they truly share one implementation).
- `utils/__init__.py` now also lazy-exports `plot_pca_variance_ratio` (symmetry).
- Doc fixes: CLAUDE.md decision #2 + HANDOFF decision #3 now say PCA/scVI run on
  `panel_for_embedding`, not "all genes" / "the HVG subset"; CLAUDE.md & HANDOFF
  UMAP `min_dist` corrected 0.5 → 0.05 (spread=2.5) to match both notebooks;
  CLAUDE.md raw cell-count table labeled "(raw, pre-QC)".
- HVG figures re-rendered for all four samples via `scripts/render_hvg_plot.py`.
- OPEN: both notebooks apply `min_dist=0.05, spread=2.5` unconditionally (no
  tissue branch). If the pancreas UMAP was meant to keep a larger `min_dist`,
  that is a code change to confirm — flagged, not silently changed.

## 2026-05-25 — PCA "elbow" is diagnostic-only (no embedding change)

A user flagged that thymus sample **0041534** showed "elbow at PC 5" on
`figures/0041534/01_preprocessing/pca_variance_ratio.png` and worried only
5 PCs were in use. Investigation: the elbow is **cosmetic** — computed in
stage-01 cell 19 and drawn on the scree plot, but it never fed
`n_comps`/`n_pcs`. The pipeline uses a **fixed 30 PCs** (stage-01
`rsc.pp.neighbors(n_pcs=30)`; the stage-02 PCA fallback is also 30), and
stage-02's final cell types come from the **30-dim scVI latent**, not PCA.
0041534's elbow lands lower (PC5 vs PC7–8 elsewhere) only because its scree
is genuinely flatter (5.2% cumulative variance at PC30 vs 10–12.7% for the
pancreas samples) — expected for homogeneous thymus; it does not move the
fixed cutoff.

**Changed this session (labeling/docs only — embeddings untouched):**
- New shared helper `utils.analysis_utils.plot_pca_variance_ratio`: the scree
  plot now foregrounds a solid "n_pcs used downstream = 30" line and demotes
  the elbow to a faint "diagnostic only" marker; the title no longer reads
  "elbow at PC N". The elbow off-by-one is corrected (`argmax + 2`, was `+2+1`).
- Stage-01 cell 19 rewritten to call the helper (mirrors the cell-40 HVG pattern).
- New `scripts/render_pca_variance_plot.py` backfills the figure from a saved
  `_preprocessed.h5ad` (backed-mode read of `uns['pca']`) — **no pipeline
  re-run needed**; re-rendered for all four samples.
- CLAUDE.md stage-01 block corrected (n_comps=50 computed / 30 used downstream).

No `.h5ad` outputs were modified; only the four `pca_variance_ratio.png`
figures were refreshed.

## 2026-04-27 — Stage 07 immune+islet proximity (added) + dual-seed islet ID

### What was added

1. **`notebooks/07_immune_islet_proximity.ipynb`** — pulls all immune
   cells + the indeterminate-with-immune-marker subset from both
   phenotyped h5ads, retrains scVI immune-only with `batch_key='sample'`
   to remove sample-of-origin variation, scores immune subtypes
   (T_helper/cytotoxic/reg/exhausted, NK, Macro_M1/M2/resident,
   Monocyte, DC, B_pan, B_plasma), identifies islets via DBSCAN on
   endocrine seed coordinates, computes per-cell distance-to-islet
   via cuML NN, bins (intra <25 μm, peri 25–50 μm, proximal 50–200 μm,
   distal >200 μm), and produces contingency χ² + per-islet
   infiltration scores normalized per 100 endocrine cells.

2. **Sample-comparable thresholds (lognorm).** Replaced all
   `raw_count >= N` thresholds with `lognorm >= log1p(N)` because
   0041326 has ~10% lower median total_counts and 298 panel genes
   drop by min_cells=100 (vs. 14 in 0041323). A static raw threshold
   was systematically more stringent in 0041326. Touched files:
   `scripts/lineage_phenotype.py`, `scripts/verify_upper_tail_doublets.py`,
   `scripts/show_canary_panels.py` (DETECTION RULE text), notebook 02
   cell 22 (lineage block), notebook 07 cells 5 (immune-marker
   recovery for indeterminate cells) + 15 (islet seed).

3. **Dual islet seeding.** Two parallel seed strategies:
   - **Strict hormone-marker seed** (`scripts`-internal, in notebook 07
     cell 15): `lognorm >= log1p(5)` of `[CHGA, INSM1, ISL1, NEUROD1,
     FEV]`. Output: `data/processed/islets_dbscan.csv` (110 islets in
     0041323, **2 islets in 0041326** — biological asymmetry, see below).
   - **Endo-label seed** (`scripts/islet_proximity_endolabel.py`):
     `celltype_lineage in {Endocrine, Beta, Alpha, Delta, Endocrine_pan}`
     (Epsilon excluded — 150k count in 0041326 is artefactual).
     Output: `data/processed/islets_dbscan_endolabel.csv` (265 islets in
     0041323, 318 islets in 0041326 — comparable across samples).

### The 0041326 islet asymmetry — what it means

With the **strict hormone-marker seed**, 0041326 yields only 2 islets
because hormone markers (CHGA/INSM1/ISL1/NEUROD1/FEV) are genuinely
under-expressed across this section. Per-100-endocrine numbers from the
strict CSV for 0041326 are **statistically meaningless** (every immune
cell in the section is NN-assigned to one of those 2 seeds).

With the **endo-label seed**, 0041326 recovers 318 islets — comparable
to 0041323's 265. The endo-label seed uses score-genes-argmax cell-type
assignments that don't require any single canonical marker to be highly
expressed, so a cell with merely above-cluster-average islet-program
expression still seeds. **Both samples now have islet-rate parity**
(7.30 vs 9.93 islets per 1000 endocrine-labeled cells).

### Headline biological findings (use endo-label CSVs for cross-sample comparison)

- **0041323 (less T1D progression?):** Most islets show NO immune
  infiltration; the small fraction that do are dominated by T_cytotoxic
  and T_reg. Median per-100-endocrine for every immune subtype = 0.
  Only DC shows a meaningful islet-zone enrichment (OR=2.58).

- **0041326 (more T1D progression?):** Pervasive low-level immune
  infiltration. Median per-100-endocrine: Monocyte=16, T_reg=12,
  Macro_M2=9, T_helper=4, B_plasma=3. Top-5 most-infiltrated islets
  carry 150–223 immune cells each, dominated by Macro_M2 (17–36),
  T_reg (22–36), and T_cytotoxic (8–18). Mild OR enrichments:
  B_plasma=1.36, T_cytotoxic=1.30, NK=1.15.

This asymmetry is consistent with 0041326 being further along the T1D
disease course — broader immune infiltration corresponds to greater
β-cell hormone marker depletion, which is exactly why the strict-seed
approach failed for that sample.

### CSVs produced

| File | Purpose |
|---|---|
| `data/processed/islets_dbscan.csv` | strict hormone-marker DBSCAN islets (use for high-stringency 0041323-only analyses) |
| `data/processed/islets_dbscan_endolabel.csv` | endo-label DBSCAN islets (use for cross-sample comparisons) |
| `data/processed/immune_proximity_summary.csv` | strict-seed contingency (subtype × bin × sample) |
| `data/processed/immune_proximity_summary_endolabel.csv` | endo-label contingency |
| `data/processed/islet_infiltration_per100endo.csv` | strict-seed per-islet immune counts + per100endo rates |
| `data/processed/islet_infiltration_per100endo_endolabel.csv` | endo-label per-islet immune counts + per100endo rates |

### Helper scripts

- `scripts/islet_proximity_endolabel.py` — generates the endo-label CSVs.
- `scripts/compare_islet_seeds.py` — side-by-side strict vs endo-label
  numerical summary (run after the above).
- `scripts/summarize_07.py` — dumps the strict-seed notebook 07
  summary (immune subtype distribution, per-islet stats, OR table for
  0041323).

### 2026-04-27 update — replaced strict-seed proximity tables with insulitis grades

The strict hormone-marker seed gave 0041326 only 2 islets, which made
its `immune_proximity_summary.csv` row uninterpretable (every immune
cell got binned to 'distal' by construction). User flagged the table
as "absolute garbage" and was right.

Replaced the headline tables with `scripts/insulitis_analysis.py` which:
- Uses the endo-label seed (`celltype_lineage` ∈ Endocrine/Beta/Alpha/Delta/
  Endocrine_pan) as the islet basis. Islet counts: 0041323=265, 0041326=318.
- Assigns each islet an `insulitis_grade` per Campbell-Thompson 2013
  nPOD criteria adapted to total CD45+ pool: `no_insulitis` (<3 immune
  within 50 μm), `peri_insulitis` (3-5), `insulitis` (≥6).
- Computes per-subtype `density_enrichment` = (zone fraction) /
  (distal fraction). >1 means concentrated near islets; <1 depleted.
  Sample-comparable, doesn't depend on per-islet sparsity.

Headline files now:
| File | Content |
|---|---|
| `data/processed/immune_proximity_summary.csv` | density enrichment per subtype × sample |
| `data/processed/islet_infiltration_per100endo.csv` | per-islet immune counts + insulitis_grade |
| `data/processed/islet_insulitis_grades.csv` | per-sample insulitis prevalence |
| `data/processed/islets_dbscan.csv` | endo-label seed islets (was strict-seed) |
| `*_strict.csv` | original strict-seed tables, archived |

Findings (**SUPERSEDED — see 2026-04-28 entry below**): 0041326 = 99.4%
insulitis (316/318 islets), 0041323 = 16.6% insulitis. Those numbers came
from an earlier immune-gating version (mean 42 immune cells/islet) that
no longer matches the current `*_immune_phenotyped.h5ad` (mean ~10
cells/islet). On the current immune set, the CT-3/6 absolute grade is
60.8% / 66.7% — much closer between samples — and the per-phenotype
rotation-null grades (the new headline) are 0–4% by subtype.

### 2026-04-28 update — per-phenotype dynamic insulitis grading

The single CT-3/6 absolute grade above flattens 0041326 to "everything is
insulitis". User asked for a size-, composition-, and immune-type-aware
replacement. `scripts/insulitis_analysis.py` rewritten to:

1. **Composition per islet** from per-cell `Beta_score` / `Alpha_score` /
   `Delta_score` (already present in `_phenotyped.h5ad` from stage 02 cell 13).
   Argmax with margin ≥ 0.1 else `Endo_unresolved`. `composition_class` ∈
   {Beta_rich, Alpha_rich, Delta_rich, Mixed} per 60% threshold.
2. **Per-phenotype rotation null** (1000 rotations of immune xy around the
   clustered-seed centroid) for each of 12 immune subtypes
   (`B_pan`, `B_plasma`, `DC`, `Macro_M1`, `Macro_M2`, `Macro_resident`,
   `Monocyte`, `NK`, `T_cytotoxic`, `T_exhausted`, `T_helper`, `T_reg`).
   Acinar/Schwann/Endothelial leakage dropped from immune count.
3. **Per-stratum thresholds**: `(sample × size_class)` with 5 size bins
   (`tiny ≤20`, `small 21-50`, `med 51-100`, `large 101-200`, `xl >200`),
   pooled per-100-endocrine null values, 95th / 99th percentile →
   `thresh_peri` / `thresh_insulitis`.
4. **Per-phenotype grade** per islet — independent {no, peri, insulitis}
   grade for each subtype, plus the CT-3/6 absolute grade preserved on
   `total_immune`. **Singular `insulitis_grade` column dropped.**
5. **Descriptive logistic regression** per (sample × subtype) with odds
   ratios + 95% CI, NO p-values (n=2 caveat).

#### New / changed CSVs (single writer = `insulitis_grade_absolute.py`)

| File | Schema |
|---|---|
| `islet_infiltration_per100endo.csv` | per-islet, 53 cols incl. `<S>_zone`, `<S>_per100endo`, `<S>_grade` for each of 12 subtypes, plus `composition_class`, `size_class`, `total_immune`, `total_per100endo`, `insulitis_grade_absolute` |
| `islet_insulitis_grades.csv` | long format: `sample, subtype, no_insulitis, peri_insulitis, insulitis, total_islets, pct_insulitis, pct_peri_or_insulitis` — 26 rows (2 samples × 13 = 12 subtypes + `absolute_total`) |
| `islet_insulitis_thresholds.csv` | NEW — `sample, size_class, subtype, n_islets, n_null_samples, thresh_peri, thresh_insulitis` (120 rows) |
| `islet_insulitis_regression.csv` | NEW — `sample, subtype, term, odds_ratio, ci_low, ci_high` (no p-values) |
| `*_preweighted.csv` | one-shot archive of pre-rewrite headline CSVs |

#### Findings (current run, 2026-04-28)

| Subtype | 0041323 % insulitis | 0041326 % insulitis |
|---|---:|---:|
| absolute_total (CT 3/6) | 60.8% (161/265) | 66.7% (212/318) |
| NK | 2.3% | 3.8% |
| Monocyte | 2.3% | 3.1% |
| T_helper | 2.3% | 2.2% |
| T_cytotoxic | 0.4% | 2.2% |
| DC | 1.5% | 1.6% |
| (others) | 0–1% | 0–1% |

The "99.4% / 16.6%" headlines from the 2026-04-27 entry above are **stale**
— they were computed from an earlier immune-gating version (mean 42
immune cells/islet) that no longer matches the current
`*_immune_phenotyped.h5ad` (mean ~10 cells/islet). The CT-3/6 absolute
grade on the current immune set is 60.8% / 66.7% — much closer between
samples than the prior headline suggested.

Per-phenotype prevalences are intentionally conservative — the
rotation-null 99th-percentile threshold flags only islets that are
statistically elevated *for that subtype × size_class stratum* relative
to the same sample's immune-density background. To loosen, drop the
threshold percentile in `insulitis_analysis.py` (currently 95/99 for
peri/insulitis; try 90/95 for a more permissive cut).

#### Composition findings

- 0041323: 64 Beta_rich islets (24%), 200 Mixed (76%), 1 Alpha_rich.
- 0041326: 6 Beta_rich (2%), 5 Alpha_rich (2%), 307 Mixed (96%) — β-cell
  identity loss consistent with more-advanced T1D progression.

#### Rotation-null caveat

Rotation around the clustered-seed centroid means rotated immune cells
can land outside the original tissue convex hull. Those positions
contribute zero to per-islet zone counts, biasing the null *low* and
making real-vs-null gaps look slightly larger than reality. A
tissue-mask-aware permutation would fix this; v1 does not implement it.

#### Duplicate-writer reconciliation

`scripts/rerun_immune_pipeline.py` step 6 now early-exits with a
redirect message; `scripts/build_notebook_07.py` writes its derivatives
to `*_legacy_notebook07.csv` so headline CSVs aren't silently
overwritten with the old schema.

#### Known follow-up

- Some per-(sample × subtype) regressions hit `LinAlgError: Singular
  matrix` because of perfect separation in small composition strata
  (e.g. 0041326 has only 5 Alpha_rich islets). Skipped silently; not a
  blocker.
- Per-phenotype thresholds may be too strict if the user's question is
  "where is *any* T-cell signal elevated"; loosen percentile or fold
  T_cytotoxic + T_exhausted into a single `T_effector` aggregate column
  if that question becomes important.

### 2026-04-28 (later) — Epsilon mislabel + full reanalysis pending

#### What's broken in the current saved h5ads

`scripts/audit_panel_exclusions.py` + manual review surfaced two structural
issues with the current `_phenotyped.h5ad` files:

1. **150,230 cells in 0041326 mislabeled `Epsilon`.** Cause: the stage 02
   cell 10 marker dict had `'Epsilon': ['GHRL']` — a 1-gene "panel".
   Under detection asymmetry (0041326 has lower per-cell detection),
   the score_genes argmax overweights single-gene panels relative to
   multi-gene panels. Result: ~13% of all cells in 0041326 got
   argmax-mislabeled Epsilon. 0041323 has 0 Epsilon (correct).
2. **Donor sex differs between samples.** KDM5D = 10.3% in 0041326 vs
   0.010% in 0041323; DDX3Y = 7.7% vs 0.017%. **0041326 is male,
   0041323 is female.** Y-chromosome variance is leaking into PCA, scVI,
   UMAP, leiden, and likely contaminating cell-type calls beyond
   Epsilon.

#### Resolution path — full reanalysis with the locked plan at `/home/smith6jt/.claude/plans/full-reanalysis-plan.md`

User-confirmed parameters:
- Margin gate: **0.2** (strict)
- Panel coverage threshold: **30%**
- Spatial coherence: **disabled for rare cell types** (only common
  clustering types get the check)
- scVI sex covariate: **YES**
- Cross-sample prevalence ratio: **reported only, not a gate** (n=2)
- Rare cell types (Epsilon, Gamma): identified via **rule-based
  multi-criteria definitions** (not score_genes argmax) — see
  `scripts/refine_rare_celltypes.py`

Helper scripts written 2026-04-28 (ready to use):
- `scripts/audit_panel_exclusions.py` — produces `data/processed/panel_audit.csv` (106 candidate exclusions)
- `scripts/finalize_panel.py` — adds `var['panel_for_scoring']`, `var['panel_for_embedding']`, `obs['donor_sex']`
- `scripts/refine_rare_celltypes.py` — post-stage-02 margin gate + panel coverage + rule-based Gamma/Epsilon + spatial coherence
- `scripts/verify_rare_celltypes.py` — per-(sample × type) verification suite

Notebook edits **still pending** (must be done before pipeline kickoff):
- `notebooks/02_phenotyping.ipynb` cell 10: drop `'Epsilon': ['GHRL']`; mask
  marker panels via `var['panel_for_scoring']`
- `notebooks/02_phenotyping.ipynb` cell 17: add `donor_sex` to scVI
  `categorical_covariate_keys`; restrict gene set to `panel_for_embedding`
- `notebooks/01_preprocessing_v2.ipynb` cell 18: HVG / PCA on
  `panel_for_embedding` only

Wall-time after kickoff: ~3.5–4 h on the B200 GPU node for stages 01 + 02 + 03,
followed by `refine_rare_celltypes.py` (~5 min), `insulitis_analysis.py`
(~30 s), `figs_07_immune.py` (~5 min), and `verify_rare_celltypes.py` (~10 min).

### 2026-04-29 — Reanalysis EXECUTED + Acinar hotfix

Full pipeline run: 2026-04-28 18:52 → 2026-04-29 01:47 (~6h 55m). Stage timings:
- 01/0041323 = 47:46, 01/0041326 = 1h 1m 47s
- 02/0041323 = 1h 35m 57s, 02/0041326 = 2h 6m 8s
- 03/0041323 = 16m 55s, 03/0041326 = 16m 24s
- regate_immune ≈ 5 min, rerun_immune_pipeline ≈ 4 min
- fix_acinar_argmax ≈ 11 min, refine_rare_celltypes ≈ 8 min
- post chain (insulitis + figs + verify) ≈ 17 min total

**Critical mid-run hotfix**: my notebook 02 cell 10 patch had Acinar markers
as `['AMY1A', 'AMY2A', 'AMY2B']` — but only AMY1A + CUZD1 are actually on the
panel. The ≥2-gene rule then dropped Acinar, and ~700K acinar cells argmax-
distributed across other types (B_cells ballooned to 464K). `scripts/fix_acinar_argmax.py`
(NEW, post-stage-02 hotfix) re-scored Acinar with `[AMY1A, CUZD1]`, re-argmaxed
across all 13 score columns, re-applied the anti-acinar override, and updated
both phenotyped and spatial h5ads in-place. Backups at `*.bak.preAcinarFix`.
The patch script + refine_rare_celltypes panel definitions both now use the
correct Acinar panel for future runs.

**Final celltype verification (post-Acinar-fix + refine):**

| Type | 0041323 (female) | 0041326 (male, T1D-progressed) | dotplot_pass |
|---|---:|---:|---|
| Acinar | 786,797 (68.7%) | 179,540 (15.4%) | (not in rare-verify) |
| Indeterminate (margin <0.2) | 147,900 | 614,974 | n/a |
| Indeterminate_*_lowcoverage | ~179K | ~317K | n/a |
| Beta | 178 (0.02%) | 17 | ✅ both |
| Alpha | 4,560 (0.40%) | 11,640 (1.00%) | ✅ both |
| Delta | 514 | 1,149 | ✅ both |
| Gamma (rule-based) | 0 | 2 | ✅ sample 2 |
| **Epsilon (rule-based)** | **9** | **3** | ✅ — was 150,230 with old bug, **99.99% reduction** |
| Schwann | 752 | 877 | ✅ both |
| Endothelial | 813 | 1,567 | ✅ both |

The high Indeterminate counts in 0041326 are an artefact of the user-locked
strict 0.2 margin gate × that sample's inherently lower per-cell detection.
Lowering the margin to 0.1 would recover most. Decision deferred to user.

**Headline insulitis (current run, post-Acinar-fix):**

| Subtype | 0041323 % insulitis | 0041326 % insulitis |
|---|---:|---:|
| absolute_total (CT 3/6) | **64.9%** (146/225 islets) | **78.7%** (266/338) |
| T_helper | 2.7% | 4.7% |
| NK | 2.7% | 3.0% |
| Monocyte | 3.1% | 2.1% |
| DC | 0.9% | 2.7% |
| B_plasma | 0.0% | 1.8% |
| T_cytotoxic | 0.0% | 1.5% |
| Macro_resident | 0.0% | 1.5% |
| T_reg | 0.4% | 1.2% |
| T_exhausted | 0.0% | 0.9% |

0041326 has elevated T_helper, T_cytotoxic, T_exhausted, T_reg, B_plasma,
DC, Macro_resident — consistent with active T1D auto-immune insulitis.
0041323 shows mostly innate (Monocyte) + some T_helper / NK.

**Cross-sample prevalence ratios (n=2, informational):**
- Beta: 10.3× lower in 0041326 → β-cell loss in T1D ✓ biological signal
- Alpha: 2.5× higher in 0041326 → α-cell hyperplasia (classic T1D)
- Delta: 2.2× higher in 0041326

**All headline outputs**: see DATA_README.md for the canonical table; CSVs
+ figures all timestamped 2026-04-29 01:30–01:47.

**Backups**: `*.bak.preAcinarFix` for both phenotyped + spatial h5ads. Total
~36 GB of backup data; can be cleaned later if disk pressure.

### 2026-05-04 — Stage 02 + Stage 03 figures regenerated post-resolve_indeterminate

`scripts/resolve_indeterminate.py` ran 2026-04-30 17:30 and chained the
stage-07 regenerators (`insulitis_analysis.py`, `figs_07_immune.py`,
`verify_rare_celltypes.py`), but did **not** touch stage 02 or stage 03
figures. Those still showed the pre-resolution `Indeterminate_*_lowcoverage`
legends until today.

**`scripts/refresh_figs_post_resolve.py`** (NEW) — single driver that
shells out to the existing renderers and adds the three gaps
(`spatial_neighbors`+`nhood_enrichment` recompute + CSV write,
`sq.gr.ripley` recompute, `spatial_connectivity` scatter). Idempotent.
Run end-to-end:

```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
source scripts/env_rapids.sh
python scripts/refresh_figs_post_resolve.py \
    2>&1 | tee logs/refresh_figs_post_resolve_$(date +%Y%m%d_%H%M%S).log
```

Wall time: **17.9 min** for both samples (much faster than the 25–35 min
estimate; `nhood_enrichment` was 6–8 s rather than the 5–10 min predicted
because squidpy's permutation kernel is fast on the resolved 14–15-category
problem).

**Regenerated figures** (mtimes 2026-05-04 15:20–15:37):

| Stage | File |
|---|---|
| 02 | `cluster_consensus_annotation.png`, `final_celltypes_overview.png`, `celltype_scores_umap.png`, `predicted_celltypes.png`, `scvi_vs_markers.png` |
| 03 | `co_occurrence.png`, `neighborhood_enrichment.png` (+ recomputed `*_neighborhood_enrichment.csv`), `niche_celltype_matrix.png`, `ripley_statistics.png`, `spatial_connectivity.png`, `spatial_domains.png` |

Final celltype categories visible in the legends (no Indeterminate
anywhere):

- 0041323 (14 types): `Acinar, Alpha, B_cells, Beta, Delta, Ductal, Endocrine, Endothelial, Epsilon, Myeloid, Proliferating, Schwann, Stellate, T_cells`
- 0041326 (15 types): same as above + `Gamma`

`figures/07_immune_islet/` was already current (post-2026-04-30 17:43,
chained inside the resolve_indeterminate run); ligand-receptor dotplots
were regenerated separately on 2026-05-04 14:54 by
`scripts/regen_ligrec_lineage.py`.

### 2026-05-05 — Niche-pair (neighborhood-restricted) ligand-receptor

**Status: scripts committed and pilot-validated; full sweep deferred to a
fatter node.** The active session has only 4 cores / 16 GB RAM and the
spatial h5ad loaded over Lustre stalled at 26 MB / 4 min when re-tried;
the OOM-prone `adata[mask].copy()` ordering was also fixed in the second
script revision. Re-run on a node with ≥ 64 GB RAM. Plan +
parameters preserved at `/home/smith6jt/.claude/plans/niche-pair-ligrec-plan.md`.

`scripts/regen_ligrec_lineage.py` runs squidpy ligrec globally on a 100k
cell subsample with `cluster_key='celltype_general'`. That mixes signal
across the whole tissue: a strong "Immune→Endocrine" hit can be driven by
immune cells in one tissue region and endocrine cells in a spatially
distant region — biologically implausible since L-R signaling needs
physical contact / short-range diffusion.

`scripts/niche_pair_ligrec.py` (NEW) restricts each ligrec run to the
cells in a pair of spatially adjacent niches, anchoring the inferred
interactions to a specific microenvironment.

**Adjacency rule.** Two niches A,B are adjacent when the count of
cross-niche edges in the squidpy spatial graph (k=10 Delaunay, already
in `adata.obsp['spatial_connectivities']`) is ≥ 2% of *each* niche's
self-edge count. Using the min fraction (rather than total cross-edges)
keeps small but tightly interlocked niches in scope and prevents huge
niches from being adjacent to everything by sheer size. Threshold lives
in `ADJACENCY_FRAC_MIN` and yielded 17 pairs / sample on the current run.

**Per pair.** Subset adata to cells in both niches, drop Indeterminate
and lineages with < 50 cells in the subset, subsample to 30k cells,
copy `lognorm` into `.X`, and call `sq.gr.ligrec(cluster_key='celltype_general',
n_perms=50, seed=0)`. Means / pvalues / metadata go to per-pair CSVs.
A dotplot of the top 25 significant L-R (means ≥ 0.3, p < 0.001) is
rendered with the same `_compact_ligrec_legends` styling as the global
script.

**Outputs**

```
data/processed/{sample}/niche_pair_ligrec/
    pair_{A}_vs_{B}_means.csv
    pair_{A}_vs_{B}_pvalues.csv
    pair_{A}_vs_{B}_metadata.csv
    niche_lineage_composition.csv
    niche_pair_summary.csv     # one row per pair: cross_edges, frac_a/b,
                                 n_cells, n_significant_lr, top 5 L-R
figures/{sample}/03_spatial_analysis/niche_pair_ligrec/
    niche_adjacency_overview.png   # spatial scatter + adjacency heatmap
    ligrec_niche_{A}_vs_{B}.png    # one per significant pair
```

**Run**

```bash
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
python scripts/niche_pair_ligrec.py            # both samples, all pairs
python scripts/niche_pair_ligrec.py --pilot    # first pair / sample only
python scripts/niche_pair_ligrec.py --render-only   # remake plots from CSVs
python scripts/niche_pair_ligrec.py --force         # ignore CSV cache
```

Wall time: ~30 s of ligrec per pair (50k → 30k subsample, 7 lineages,
~5k L-R candidates). 17 pairs × 2 samples ≈ 12 min compute + ~75 s of
h5ad load per sample; total ~15 min.

**Pilot finding (0041323, niche 1 ↔ niche 3 — "Endo/Exo islet boundary
↔ Endo islet core"):** 205 significant L-R pairs. Top hits by p-value:
ICAM1→PECAM1 (Immune→Vascular adhesion), GC→NR3C1 (Endocrine→Cycling
glucocorticoid signaling), THBS1→NOS3 (Stromal→Neural), ANPEP→NOS3
(Exocrine→Neural). Pilot CSVs + figure live under
`data/processed/0041323/niche_pair_ligrec/pair_1_vs_3_*.csv` and
`figures/0041323/03_spatial_analysis/niche_pair_ligrec/ligrec_niche_1_vs_3.png`;
the niche-adjacency overview (spatial scatter + adjacency heatmap with
the 17 pairs above the 2 % threshold marked) is at
`figures/0041323/03_spatial_analysis/niche_pair_ligrec/niche_adjacency_overview.png`.

**Companion aggregator** `scripts/niche_pair_ligrec_aggregate.py` reads
all `pair_*_means.csv` / `pair_*_pvalues.csv` and produces a per-sample
recurrence table (`niche_pair_ligrec_recurrence.csv`) plus a cross-
sample merge (`data/processed/niche_pair_ligrec_cross_sample.csv`).
Useful to separate constitutive interactions (significant in many
adjacent-niche contexts) from microenvironment-specific ones (significant
in only one or two pairs). Run after the main script.

### 2026-05-21 — Thymus pipeline: 5-fix bundle landed + verified

Two thymus samples (`0041435`, `0041534`) were added 2026-05-19; raw
zarrs at `data/raw/{S}.zarr`. The prior pipeline run produced
biologically wrong output (DP_thymocyte 1.6%/1.0%, AIRE=Hassalls=pDC=0,
~20–27% cell loss at QC, undersized DBSCAN). A 5-fix root-cause plan
(`/home/smith6jt/.claude/plans/the-last-session-was-fuzzy-starfish.md`)
was authored 2026-05-20, applied to the working tree, and committed
in `278f885` + `f9f9b4c` (2026-05-21).

**The five fixes (commit `278f885`):**

1. **Tissue-aware QC defaults** (`scripts/run_local_pipeline.sh`
   lines 41–58): thymus exports `XENIUM_MIN_COUNTS=20` / `MIN_GENES=10`
   (vs pancreas 50/20). Recovers ~95% of cells vs ~75% prior — thymus
   median is 107/124 counts/cell and the pancreas floor was cutting
   genuine small thymocytes.
2. **DP_thymocyte panel trim** (`scripts/tissue_markers/thymus.py`
   lines 55–62, 205–210): drop RAG1 (0.17%/0.27% detection), RAG2
   (0.09%/0.19%), PTCRA (0.26%/0.19%) — three near-zero markers that
   diluted the DP score_genes mean to ~(shared)/10 vs SP4/SP8's
   ~(shared)/8. Worked-example: a CD4+CD8 double-positive cell now
   scores DP=1.79 > SP=1.44 instead of DP=1.20 < SP=1.375.
3. **Relaxed rare-type rules** (same file, lines 305–368, 430–495):
   added `log_present = log1p(1) ≈ 0.693` threshold; AIRE_mTEC /
   Hassalls / pDC "high" gates now use `log_present` (was log_high =
   log1p(5) ≈ 1.792). Non-thymocyte exclusion widened from `< log_detect`
   to `< log_high`. TGM3 removed from Hassalls cornified panel (absent
   in 0041534's panel).
4. **Medulla DBSCAN tuning** (`scripts/thymic_zonation_analysis.py`
   lines 63, 68, 73, 121): `min_samples` 15→5; `marker_high` seed
   expanded to `AIRE | KRT14 | PSMB11 | CTSL | FOXN1` (cortico-
   medullary junction); lineage gate adds `cTEC`.
5. **Notebook guards** (`notebooks/02_phenotyping.ipynb`):
   - Cell 15: T-lineage protection — never flip cells with CD3E or
     TRAC ≥ log1p(1) via anti-DP override.
   - Cell 20: cluster-consensus guard preserves per-cell
     `predicted_celltype` for `_THY_SUBTYPES = {DN_thymocyte,
     DP_thymocyte, SP4_thymocyte, SP8_thymocyte, Treg}` so leiden-
     majority doesn't collapse the cortex back to one thymocyte type.

`scripts/refine_rare_celltypes.py` passes `log_present` through to
`apply_rare_celltype_rules`.

**Headline numbers (current on-disk artifacts, 2026-05-20):**

| Metric | 0041435 (f) | 0041534 (m) |
|---|---:|---:|
| Cells | 542,055 | 302,056 |
| Retention (% of raw) | 93.5% | 91.3% |
| DP_thymocyte | 118,015 (21.8%) | 99,231 (32.9%) |
| SP4 / SP8 / DN | 31K / 42K / 48K | 11K / 19K / 16K |
| Treg | 21,760 | 10,116 |
| AIRE_mTEC | 176 | 24 |
| Hassalls | 94 | 24 |
| pDC | 3 | 3 |
| mTEC / cTEC | 41,785 / 29,212 | 17,332 / 12,551 |
| Spatial niches | 12 | 12 |
| Library-size leak ρ (UMAP_scvi) | (0.22, 0.15) ✓ | (0.35, 0.10) ⚠ |
| Medulla lobules (kNN-CC, 2026-05-21) | 27 | 22 |
| AIRE+ neighborhood centers | 176 | 24 |

**Technical verification gates (G1–G7) — 6 of 7 pass:**

- G1 (retention ≥ 50%), G3 (spatial_niche present), G4 (zero null
  celltype), G5 (vocabulary OK, 21 labels each), G6 (rare-rule prereq
  genes in panel), G7 (4 zonation CSVs present): **all PASS** for both
  samples.
- G2 (library-size leak < 0.3): PASS for 0041435 (max ρ=0.216), **fail
  borderline for 0041534 (max ρ=0.353** on UMAP_scvi X-axis; Y-axis
  is 0.098). scVI under-converged on this smaller sample. Per CLAUDE.md
  fallback advice: bump `max_epochs=400` + `early_stopping_patience=40`
  in stage 02 cell 17 and rerun stage 02 only for 0041534. **Tracked
  but not blocking** — most cell-type calls are made from the marker
  panels, not the leiden_scvi clusters directly.

**Known residual issues (not biological blockers):**

- **Very low pDC counts (3 / 3).** Fix 3's relaxed thresholds couldn't
  recover more — LILRA4 / IRF7 / IRF8 / CLEC4C / TCF4 detection in
  this 5K panel + FFPE chemistry may be below the floor for pDC RNA
  signatures regardless of threshold.
- **0041534 G2 borderline scVI leak (ρ=0.353 X-axis).** Acceptable
  given cell-type calls are panel-driven, not embedding-driven. Can
  be tightened by bumping `max_epochs=400` in stage-02 cell 17.

**Files of record:**

```
data/processed/0041435/0041435_phenotyped.h5ad           1.2 GB
data/processed/0041435/0041435_spatial_analysis.h5ad     2.3 GB
data/processed/0041534/0041534_phenotyped.h5ad           596 MB
data/processed/0041534/0041534_spatial_analysis.h5ad     1.1 GB
data/processed/medulla_zones.csv                         49 rows (27+22 lobules; renamed 2026-05-21)
data/processed/thymic_zone_composition.csv               160 rows
data/processed/thymocyte_maturation_gradient.csv         40 rows
data/processed/AIRE_neighborhood.csv                     200 rows (176+24 centers)
data/processed/panel_audit_thymus.csv                    panel-exclusion audit
figures/{0041435,0041534}/{01_preprocessing,02_phenotyping,03_spatial_analysis,07_thymic_zonation}/
```

Forensic backups (untracked): `notebooks/02_phenotyping.ipynb.bak.thymus_fixes`,
`notebooks/03_spatial_analysis.ipynb.bak.ligrec_emptycat`.

**Resume work — quick commands:**

```bash
cd /blue/maigan/smith6jt/Xenium_Analysis
source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis

# Re-run stage 02 only for 0041534 to tighten G2 leak (~75 min on B200):
STAGES="02" TISSUE=thymus SAMPLES="0041534" \
    nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_thymus_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Re-run zonation (kNN-density algorithm; tune via CLI flags):
python scripts/thymic_zonation_analysis.py --samples 0041435 0041534
# Tighter medulla calling:   --score-medulla 0.10
# Looser medulla calling:    --score-medulla 0.00
# Larger smoothing radius:   --knn 100
# Single-pass (more detail): --smooth-passes 1
```

### 2026-05-21 (later) — Thymic zonation rewrite (DBSCAN → kNN density)

The DBSCAN-on-sparse-marker-seeds approach used through 2026-05-20
produced small "bullseye" zones (concentric distance bins around tiny
seed centroids) — anatomically wrong. Thymic medulla is broad,
contiguous *central* regions of each lobule (hundreds of μm wide),
not point-cluster bullseyes.

Replaced with cell-type-density zonation (commit `d889820`,
`scripts/thymic_zonation_analysis.py`; old script preserved at
`scripts/thymic_zonation_analysis.py.bak.bullseye`):

1. **Source classes** from post-5-fix labels:
   - cortex_source = `{cTEC, DP_thymocyte}`
   - medulla_source = `{mTEC, AIRE_mTEC, Hassalls, SP4_thymocyte,
     SP8_thymocyte, Treg, DC, pDC}`
2. **cuML brute kNN** (k=50, sklearn fallback) on `obsm['spatial']`.
3. **Two-pass density**: per-cell source-class fraction, then averaged
   over the same kNN graph for smoothing.
4. **Signed dominance score** `medulla_frac − cortex_frac`:
   - `> 0.05` → medulla
   - `< −0.10` → cortex
   - otherwise on-tissue → cortico_medullary_junction
   - `max(med, cor) < 0.15` → low-signal default to cortex
5. **`scipy.sparse.csgraph.connected_components`** labels medullary CCs
   on the kNN graph restricted to medulla cells; CCs < 200 cells are
   demoted to peri_medullary.

**Output schema changes:**

- `medulla_dbscan.csv` → renamed `medulla_zones.csv`. Columns:
  `sample, medulla_id, centroid_x, centroid_y, n_cells, x_extent_um,
  y_extent_um`. Rows are now lobule-scale CCs (sizes 200–91,800
  cells), not DBSCAN seed clusters.
- `thymic_zone_composition.csv` schema unchanged; zone values now
  `medulla / cortico_medullary_junction / peri_medullary / cortex`.
- `thymocyte_maturation_gradient.csv` schema unchanged.
- `AIRE_neighborhood.csv` unchanged.
- NEW figure: `figures/{sample}/07_thymic_zonation/medulla_lobules_overview.png`
  — scatter colored by connected-component id, grey for non-medulla.

**Results:**

| Metric | 0041435 | 0041534 |
|---|---:|---:|
| Lobules (CCs ≥ 200 cells) | 27 | 22 |
| Medulla cells total | 234,523 (43%) | 11,909 (4%) |
| Median lobule size (cells) | 1,188 | 392 |
| Max lobule size (cells / extent) | 91,800 / 4,532 μm | 2,202 / 1,305 μm |
| Median lobule extent | 393 μm | 405 μm |

0041534 is a cortex-dominant section by source-cell count (37% cortex
source vs 22% medulla source), so its medullary regions are
smaller / scattered relative to the medulla-rich 0041435 — consistent
with the male donor's broader phenotypic profile.

**Composition self-consistency (cortex → medulla enrichment):**

- 0041435: cTEC 11.5% → 2.0%; DP_thymocyte 35.2 → 12.9%; mTEC 2.0 →
  13.4%; Treg 2.1 → 5.2%.
- 0041534: cTEC ≈ ≈; DP_thymocyte 40.2 → 11.0%; mTEC 4.2 → 10.8%;
  Treg 2.5 → 7.7%.

**Thymocyte maturation gradient (cortex → medulla):**

- 0041435: DP 56.4 → 33.6%, SP4 10.3 → 13.8%, Treg 3.4 → 13.5%.
- 0041534: DP 69.3 → 31.5%, SP4 6.3 → 11.7%, Treg 4.3 → 22.1%.

DP dropping and Treg rising from cortex to medulla is the canonical
thymic maturation signature.

**Tunable CLI flags:**

- `--knn` (default 50) — k for kNN density
- `--score-medulla` (default 0.05) — score above this → medulla
- `--score-cortex` (default −0.10) — score below this → cortex
- `--presence-floor` (default 0.15) — min max(med, cor) to score
- `--smooth-passes` (default 2) — 1 = local only, 2 = + averaging
- `--min-medulla-cells` (default 200) — minimum CC size to keep

If lobules look over-merged (e.g., 91K-cell central lobule in 0041435
might be one large contiguous medullary swath OR a merging artifact),
re-run with `--score-medulla 0.10` to tighten the gate. If the section
looks under-segmented (medullary cells exist but lobules don't form),
try `--score-medulla 0.00` or `--smooth-passes 1`.

### 2026-05-21 (later) — Thymus UMAP retune (aggressive separation)

The post-5-fix thymus phenotyping ran through but the resulting UMAPs
(`figures/{0041435,0041534}/02_phenotyping/predicted_celltypes.png` and
siblings) rendered as single continuous blobs — 21 cell-type labels
overlapping each other with no visible inter-cluster whitespace. scVI
was actually separating types in latent space (the celltype_confidence
overlay showed clear local concentrations) but the
`min_dist=0.5/spread=1.0/n_neighbors=30` UMAP parameters from notebook
02 cell 17 (and 01 cell 18) collapsed those neighborhoods into one
visible mass.

**Fix (user-selected aggressive params):**

| Param | Old | New |
|---|---|---|
| `n_neighbors` (scVI graph) | 30 | 10 |
| `n_neighbors` (PCA graph) | 15 | 10 |
| `min_dist` | 0.5 | 0.05 |
| `spread` | 1.0 | 2.5 |
| `metric` | cosine | cosine (unchanged) |
| `random_state` | 0 | 0 (unchanged) |

**Scope: thymus only (0041435, 0041534).** Pancreas h5ads + figures NOT
re-processed in this pass. Once the notebooks below are rerun for
pancreas, those samples will pick up the new params too.

**Implementation (commit pending):**

- `scripts/retune_umap_thymus.py` (NEW). For each thymus sample:
  rebuilds `obsm['X_scvi']` and `obsm['X_pca']` neighbors on GPU via
  `rsc.pp.neighbors` (n_neighbors=10), recomputes UMAP on CPU via
  `sc.tl.umap` (min_dist=0.05, spread=2.5), preserves the prior
  embeddings under `obsm['X_umap_blob']`, `obsm['X_umap_scvi_blob']`,
  `obsm['X_umap_pca_blob']`, and mirrors the new keys into the matching
  `_spatial_analysis.h5ad`. One-shot safety copy at
  `*.h5ad.bak.preUmapRetune`. Library-size leak gate: aborts and
  restores from blob if `|spearman(X_umap[:,k], total_counts)|` exceeds
  old by > 0.05 on either axis or hits 0.40 absolute ceiling.
- `notebooks/02_phenotyping.ipynb` cell 17 (scVI + PCA branches) and
  `notebooks/01_preprocessing_v2.ipynb` cell 18 — params updated in
  place to the new values; print message updated to match.
- `scripts/refresh_figs_post_resolve.py 0041435 0041534` — invoked to
  re-render all stage 02 + 03 UMAP-based figures.

**Headline results (2026-05-21):**

| Sample | Old leak ρ (X_umap) | New leak ρ | Old range | New range |
|---|---|---|---|---|
| 0041435 | (+0.216, +0.148) | (+0.203, +0.149) | x [+1.2, +14.6] · y [-4.0, +8.6] | x [-7.2, +20.6] · y [-10.1, +14.7] |
| 0041534 | (+0.353, +0.098) | (+0.326, +0.117) | x [+1.2, +12.3] · y [-4.4, +7.5] | x [-6.0, +16.4] · y [-10.2, +12.9] |

Both samples pass G2 (leak stayed within tolerance; 0041534's X-axis
ρ actually improved from 0.353 → 0.326 because the smaller neighborhood
size (k=10) lets scVI noise dominate over residual library-size
structure). UMAP X-range roughly doubled in each sample — spread=2.5
took effect.

**Wall time:** 19.7 min for the retune script (both samples; GPU
neighbors ~3 min + CPU UMAP ~7-12 min each). Figure rerender via
`refresh_figs_post_resolve.py` follows separately (~15-18 min).

**How to revert** if the aggressive embeddings are unwanted later:

```python
import scanpy as sc
for s in ("0041435", "0041534"):
    a = sc.read_h5ad(f"data/processed/{s}/{s}_phenotyped.h5ad")
    for blob, dst in (("X_umap_blob","X_umap"),
                      ("X_umap_scvi_blob","X_umap_scvi"),
                      ("X_umap_pca_blob","X_umap_pca")):
        if blob in a.obsm:
            a.obsm[dst] = a.obsm[blob]
    a.write_h5ad(f"data/processed/{s}/{s}_phenotyped.h5ad")
    # mirror to spatial h5ad as well, then rerun refresh_figs_post_resolve
```

The `*.bak.preUmapRetune` files at `data/processed/{sample}/` are an
additional unconditional rollback path (full h5ad copies from before
the retune).

**Caveat for next pancreas rerun.** Notebooks 01 cell 18 + 02 cell 17
now contain the aggressive params, so the next `STAGES="01 02"` run on
0041323/0041326 will produce aggressively separated pancreas UMAPs too.
If that's not wanted, revert the four numerical changes in those cells
before kicking the pipeline off (n_neighbors 10→30/15; min_dist 0.05→0.5;
spread 2.5→1.0). For pancreas the moderate (min_dist=0.1, spread=1.5,
n_neighbors=15) values may be more appropriate given pancreatic
continuous gradients (e.g. Alpha/Beta/Delta hormonal continuum).

**Known risk under aggressive params.** The DN→DP→SP4/SP8→Treg
thymocyte maturation continuum may visually fragment into separate
islands rather than render as a smooth arc, since `min_dist=0.05`
tightens local cluster packing. If the user wants to *see* the
maturation gradient at a glance, re-run with moderate params on the
saved scVI latent — no full pipeline rerun needed.

### 2026-05-22 — Full thymus rebuild from scratch + scVI epoch bump + inversion diagnostic

**Trigger.** Apparent backwards medulla→cortex layer arrangement in
`figures/0041534/07_thymic_zonation/medulla_spatial_overview.png`:
medulla cells on the *outer* rim of each lobule, cortex cells in the
*core*, opposite of 0041435 (correct anatomy). User asked whether this
was real biology or a cell-ID artifact analogous to the Epsilon=150K
argmax mislabel fixed 2026-04-29.

**Diagnostic verdict (pre-rebuild and post-rebuild, identical):**
**Real biology — not a cell-ID artifact.**
`scripts/diagnose_0041534_zonation.py` runs three tests on
`_phenotyped.h5ad` for both samples, with 0041435 as a known-correct
internal control:
1. **T1 marker-vs-label**: 0041534 mTEC pass-rate 0.95 (own panel >
   other panel in 95% of mTEC-labeled cells); AIRE_mTEC + Hassalls
   1.00; cTEC degraded to 0.61 (vs 0.89 in 0041435 reference).
2. **T2 spatial coherence**: 0041534 mTEC same-label kNN median 0.10
   (vs 0.15 reference); cTEC 0.05 (vs 0.10). Modestly lower but not
   collapsed.
3. **T3 zone-vs-marker**: 0041534 recomputed medulla zone has mTEC mean
   lognorm 0.028 > cTEC 0.008 (correct direction, matches 0041435's
   medulla zone which has mTEC 0.038 > cTEC 0.013). The "inversion" is
   what the raw RNA says, not a label flip.

The most likely biological mechanism is **cortical involution /
hypoplasia** in 0041534 (the male donor): cTEC marker expression is
globally depleted (own_mean 0.081 vs 0041435's 0.192, ~half), residual
DP_thymocytes pack into small remnant cortical islands giving the
"blue-core" appearance, while the small remaining mTEC/Treg signal
forms the boundaries around those islands giving the "red-rim"
appearance.

**Rebuild details.** Pipeline restarted from raw zarrs (which were
intact at `data/raw/0041435.zarr` + `0041534.zarr`; processed h5ads
were missing). Steps:

1. Bumped scVI in `notebooks/02_phenotyping.ipynb` cell 17:
   `max_epochs 200 → 400`, `early_stopping_patience 20 → 40`. Target:
   tighten 0041534's borderline G2 leak (was ρ=0.353 X-axis).
2. Regenerated `data/processed/panel_audit_thymus.csv` (missing) via
   `compute_panel_detection.py --samples 0041435 0041534 --out-dir /tmp`
   then `audit_panel_exclusions.py --tissue thymus --samples 0041435 0041534
   --output data/processed/panel_audit_thymus.csv --staging-dir /tmp`.
3. `STAGES="01 02 03" TISSUE=thymus SAMPLES="0041435 0041534"
   bash scripts/run_local_pipeline.sh` — **2h 3m wall-clock total**.
   Per-stage timings (faster than HANDOFF estimates because scVI
   early-stopped well before the 400-epoch cap):
   - 01/0041435 = 15m, 01/0041534 = 8m
   - 02/0041435 = 42m, 02/0041534 = 23m
   - 03/0041435 = 16.5m, 03/0041534 = 18.5m
4. `scripts/refine_rare_celltypes.py --tissue thymus --samples 0041435 0041534`
5. `scripts/resolve_indeterminate.py --samples 0041435 0041534 --target both --write --no-downstream`
6. `scripts/thymic_zonation_analysis.py --samples 0041435 0041534`
7. `scripts/diagnose_0041534_zonation.py` — confirmed same verdict.

**Library-size leak (G2) on rebuilt h5ads:**

| Sample | UMAP_scvi leak ρ | Notes |
|---|---|---|
| 0041435 | (+0.203, −0.144) | ✓ both axes under 0.3 |
| 0041534 | (+0.112, −0.344) | X greatly improved (was +0.353); Y now slightly over 0.3 — same magnitude, rotated. Cell-type calls are panel-driven, not embedding-driven, so this is non-blocking. |

**Final celltype tally after refine + resolve_indeterminate
(rebuilt 2026-05-22):**

| Type | 0041435 | 0041534 |
|---|---:|---:|
| DP_thymocyte | 168,007 → ~similar | 99,057 |
| SP8_thymocyte | 34,196 | 18,630 |
| mTEC | 57,837 (refined) | 17,236 |
| Myeloid | 40,866 | 16,585 |
| B_cells | 43,549 | 16,394 |
| DN_thymocyte | 47,979 | 16,370 |
| Endothelial | — | 14,569 |
| cTEC | 28,438 → ~12,570 | 12,570 |
| SP4_thymocyte | — | 11,424 |
| Treg | 21,760 → ~similar | 10,124 |
| AIRE_mTEC | 176 | 24 |
| Hassalls | 94 | 24 |
| pDC | 3 | 3 |

Zonation: 25 medullary lobules (≥200 cells) in both samples now (was
27 / 22 in the prior run). 0041534 medulla cells = 13,139 vs 0041435's
234,867 — cortex-dominant section confirmed; consistent with the
inversion-by-involution biology.

**Artifacts produced (2026-05-22):**

- `data/processed/0041435/0041435_{preprocessed,phenotyped,spatial_analysis}.h5ad`
- `data/processed/0041534/0041534_{preprocessed,phenotyped,spatial_analysis}.h5ad`
- `data/processed/medulla_zones.csv` (50 rows), `thymic_zone_composition.csv`
  (160 rows), `thymocyte_maturation_gradient.csv` (40 rows),
  `AIRE_neighborhood.csv` (200 rows)
- `data/processed/panel_audit_thymus.csv` (5,062-gene audit)
- `reports/diagnose_0041534_zonation_20260522_151315.md` (final
  diagnostic report; verdict block at the top)
- `figures/{0041435,0041534}/07_thymic_zonation/diagnose_inversion.png`
  (2×2: celltype, recomputed zone, cTEC marker heatmap, mTEC marker
  heatmap) — the bottom-row heatmaps are the definitive evidence
  cells aren't mislabeled.
- Spatial-h5ad backups at `*.bak.preIndetResolve` (rollback path).

**For a future redo of this rebuild:**

```bash
# Prereqs the rebuild needed:
python scripts/compute_panel_detection.py --samples 0041435 0041534 --out-dir /tmp
python scripts/audit_panel_exclusions.py --tissue thymus \
    --samples 0041435 0041534 \
    --output data/processed/panel_audit_thymus.csv --staging-dir /tmp

# Pipeline + post-step chain:
STAGES="01 02 03" TISSUE=thymus SAMPLES="0041435 0041534" \
    nohup bash scripts/run_local_pipeline.sh \
    > logs/pipeline_driver_thymus_$(date +%Y%m%d_%H%M%S).log 2>&1 &

python scripts/refine_rare_celltypes.py --tissue thymus --samples 0041435 0041534
python scripts/resolve_indeterminate.py --samples 0041435 0041534 \
    --target both --write --no-downstream
python scripts/thymic_zonation_analysis.py --samples 0041435 0041534
python scripts/diagnose_0041534_zonation.py
```

**scVI model storage note (clarification from this rebuild session).**
The notebook does NOT save trained scVI models to disk. `vae.train(...)`
runs in-memory; `vae.get_latent_representation()` is copied to
`adata.obsm['X_scvi']` and the model is discarded when the kernel exits.
The only persisted output is the 30-dim latent matrix inside
`_phenotyped.h5ad`. To enable model reuse across sessions, add
`vae.save("data/processed/{sample}/scvi_model/")` after training and a
matching `scvi.model.SCVI.load(...)` branch — neither is in the code today.

### 2026-05-22 (later) — Tissue-dispatch fix for cluster_consensus_annotation plot

**Bug.** `notebooks/02_phenotyping.ipynb` cell 20 hardcoded pancreas-only
`LINEAGE_MAP`, `SPECIFIC_KEEP`, `EXOCRINE`, and `PALETTE` dicts —
shadowing `_markers.LINEAGE_MAP` that cell 1 had already imported. On
thymus runs the third panel was always titled "Specific celltype, no
exocrine" (meaningless — thymus has no exocrine), and `SPECIFIC_KEEP`
didn't include any thymic-specific labels, so thymocytes and TECs were
silently filtered out of panels 2/3 (LINEAGE_MAP `.map()` returned NaN
for every thymic label in panel 1).

**Fix.**
1. Added `SPECIFIC_KEEP`, `PALETTE`, and `THIRD_PANEL = {exclude_types,
   title}` exports to both `scripts/tissue_markers/pancreas.py` and
   `scripts/tissue_markers/thymus.py`. `LINEAGE_MAP` was already exported
   from both. For thymus, the third-panel filter drops `DP_thymocyte +
   DN_thymocyte` (the bulk cortical thymocytes, ~55% of section) so
   TECs, SP4/SP8/Treg, AIRE_mTEC, Hassalls, and rare immune types are
   visible. Panel title is now
   `"Specific celltype, no bulk thymocyte (DP/DN)"`.
2. Refactored `notebooks/02_phenotyping.ipynb` cell 20 to read these
   from `_markers` instead of hardcoding them. Pancreas behavior
   unchanged (same vocab, same panel-3 title `"Specific celltype, no
   exocrine"`).
3. Regenerated `figures/{0041435,0041534}/02_phenotyping/cluster_consensus_annotation.png`
   via standalone `/tmp/regen_consensus_annotation.py` (avoids re-running
   the whole stage 02 notebook).

**Verification.** Re-rendered panel counts:
- 0041435: panel 1 / panel 2 / panel 3 = 542,055 / 542,055 / 375,838
  (panel 3 drops 168K DP + 48K DN thymocytes, ~31% reduction).
- 0041534: 302,056 / 302,056 / 186,629 (panel 3 drops 99K DP + 16K DN,
  ~38% reduction).
- No cells dropped from panel 1 — all thymic types now map to a lineage.

**Comprehensive audit.** An Explore subagent swept all notebooks
(01/02/03), the post-step refinement scripts (`refine_rare_celltypes.py`,
`resolve_indeterminate.py`, `thymic_zonation_analysis.py`), and the
diagnostic `diagnose_0041534_zonation.py` for similar tissue-mismatch
bugs. **The cell-20 issue above was the only one found.** All other
thymus-path code dispatches correctly via `_markers` or is
tissue-agnostic. Pancreas-only stage-07 scripts
(`insulitis_analysis.py`, `figs_07_immune.py`, etc.) are not invoked on
the thymus path so their hardcoding is acceptable.

### 2026-05-22 (later) — Indeterminate purged from all final figures + CSVs

**Bug.** Per project policy #15 (CLAUDE.md), `Indeterminate*` labels
must never appear in final figures or data. The audit found three
output paths that still leaked them:

1. `figures/{sample}/02_phenotyping/predicted_celltypes.png` —
   colored by `predicted_celltype` (1.5–1.8 % of cells in
   `Indeterminate_stromal / _neural / _vascular` buckets).
2. `figures/{sample}/02_phenotyping/scvi_vs_markers.png` — same
   source column.
3. `data/processed/{sample}/{sample}_celltype_assignments.csv` —
   exported the `predicted_celltype` column verbatim (7,939 rows on
   0041435; 5,343 on 0041534).

All three were rendered/written *before* `resolve_indeterminate.py`
ran, so the pre-resolution argmax labels (containing Indeterminate
strings) leaked through.

**Fix.**

- `scripts/rerender_stage02_celltype_figs.py`:
  - Block 3 (predicted_celltypes.png): now colors by **resolved
    `celltype`** + `celltype_confidence` (instead of
    `predicted_celltype` + confidence). Title says
    "Cell type (resolved)".
  - Block 4 (scvi_vs_markers.png): now `leiden_scvi` + **resolved
    `celltype`** instead of leiden + predicted_celltype.
  - Block 1 (cluster_consensus_annotation.png): replaced the
    single-panel `manual_celltype` render with the same 3-panel,
    tissue-dispatched layout that notebook 02 cell 20 produces
    (General lineage | Specific celltype | Specific celltype minus
    dominant bulk).
  - **New** Block 6: writes `{sample}_celltype_assignments.csv` with
    only resolution-clean columns: `celltype`, `celltype_confidence`,
    `celltype_lineage`, `celltype_resolution`,
    `celltype_neighbor_agreement` (when present). An inline assertion
    aborts the run if the resulting CSV still contains any
    `Indeterminate` string.

- `notebooks/02_phenotyping.ipynb`:
  - Cell 17 predicted_celltypes.png block: filters out
    `predicted_celltype.startswith("Indeterminate")` rows before
    `sc.pl.umap` so the inline preview during nbconvert is also
    clean. Title flagged "(Indeterminate hidden)".
  - Cell 18 scvi_vs_markers.png block: same Indeterminate filter
    applied.
  - Cell 24 celltype_assignments.csv export: replaced the
    `[celltype, predicted_celltype, celltype_confidence]` column
    list with `[celltype, celltype_confidence,
    celltype_lineage, celltype_resolution,
    celltype_neighbor_agreement]` (optional columns added when
    present). Cell 24's output is now a *stage-02 snapshot*;
    `rerender_stage02_celltype_figs.py` overwrites it post-resolve.

**Verification.** Final Indeterminate sweep across all thymus
artifacts: **0 Indeterminate strings in any CSV; 0 Indeterminate
classes in any figure legend.** Validated by:

- `data/processed/{0041435,0041534}/_celltype_assignments.csv` →
  zero "Indeterminate" substring matches (assertion inside the
  rerender script also passes).
- All thymus pipeline CSVs (`medulla_zones`,
  `thymic_zone_composition`, `thymocyte_maturation_gradient`,
  `AIRE_neighborhood`, `panel_audit_thymus`, per-sample CSVs from
  stages 02 + 03) → zero Indeterminate rows.
- Both samples' `predicted_celltypes.png` legend now shows 21 clean
  cell types; `scvi_vs_markers.png` shows leiden_scvi vs 21 clean
  cell types; `cluster_consensus_annotation.png` is the 3-panel
  tissue-dispatched layout with zero Indeterminate.
- `celltype` and `celltype_lineage` on both h5ads each contain 21
  distinct labels, none starting with "Indeterminate".

**Audit-only columns retained on h5ad (do NOT plot / export):**
`predicted_celltype` (per-cell argmax pre-resolution),
`celltype_pre_resolution` (snapshot before `resolve_indeterminate`
ran). These remain on `adata.obs` for provenance per CLAUDE.md
decision #15, but no figure renderer or CSV exporter consumes them
anymore.
