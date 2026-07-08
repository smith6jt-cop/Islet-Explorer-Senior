# CLAUDE.md - Pipeline Scripts

Pipeline scripts that build the data files consumed by the Shiny app. The top-level `CLAUDE.md` covers the high-level data lineage and run commands. This file documents the batch-correction pipeline, pseudotime modes, and per-script details.

## Donor Batch Correction & Pseudotime Modes (Phase 17, Apr 2026)

### scVI batch_key bug (root cause)

`single_cell_analysis/SCVI_CODEX_v2.ipynb` cell 16 set `categorical_covariate_keys=["Age", "Gender"]` but **never passed `batch_key="imageid"`**, despite comments claiming otherwise. Result: scVI's batch decoder was never activated and donor-specific technical variance (staining, antibody lot, imaging gain) leaked into `X_scVI_mean`. Symptoms: donor 6356 (ND, 18.8% of cells) formed a standalone Leiden cluster at every resolution (80–92% dominance); Cramer's V donor↔cluster ≈ 0.4 with AMI(donor) > AMI(status).

**Fix:** downstream Harmony correction (avoids scVI retrain on 2.6M cells). Two different θ values for two different goals.

### Step 8: Harmony for Leiden clustering (θ=2)

`scripts/reaggregate_islets.py` Step 8 applies `harmonypy.run_harmony` on `X_scVI_mean` with `batch=imageid`, `theta=2.0` before `sc.pp.neighbors`. Embedding stored as `obsm['X_scVI_harmony']`. Goal: donor-invariant islet cell-types. Post-fix metrics: Cramer's V dropped from 0.38 → 0.135, donor 6356's worst-cluster fraction 0.83 → 0.29, AMI(donor) 0.19 → 0.028. Status-orthogonality is fine/desired (clusters are cell-types, disease lives on pseudotime).

### Step 5: Trajectory pseudotime (raw scVI core, robust medoid root)

Step 5 uses **raw `X_scVI_mean`** (no Harmony) for the default core-only pseudotime. Empirical sweep showed even θ=0.1 compressed ND→T1D span more than it bought in donor noise reduction. Aab+ donor variance accepted as biology (Aab+ donors exhibit a heterogeneous range of islet types).

**Robust root selection** (replaces single-islet `argmax(INS)` over ND): root = medoid of top-50 ND islets by INS in scVI space. Anchored by the centre-of-mass of "healthy ND, INS-high" islets and naturally draws from multiple donors. Current root: index 703, donor 6516, INS 0.663 (99.7th percentile).

### Step 5b: Combined-mode pseudotime (Harmony θ=0.1 + structures α=0.15)

Computes a SECOND pseudotime saved as `dpt_pseudotime_combined`. Input matrix: `[core_scVI(10) | peri_scVI(10) | α_struct · z(struct(16))]` = 36 dims. Structural block: (`prop`, `count`) × (Neural, Blood Vessel, Endothelial, Lymphatic) × (core, peri), z-scored, scaled by `ALPHA_STRUCT = 0.15`. Disease ordering breaks at α ≥ 0.20; α=0.15 is the largest weight preserving ND<Aab+<T1D ordering while strengthening immune/vascular marker correlations (CD3e, CD8a, CD31, CD34, PDPN). Light Harmony (θ=0.1) on the combined matrix; same medoid root as core mode.

### Step 6 validation thresholds

- INS Spearman r < −0.3 (PASS at −0.69 raw / −0.55 combined)
- GCG Spearman r > 0.1 (relaxed from 0.2 because raw donor-amplified ~0.27 dropped post-fix)
- ND<Aab+<T1D mean ordering
- Donor eta² within status: ≤ 0.15 for ND/T1D, **≤ 0.25 for Aab+** (relaxed because Aab+ donors genuinely span heterogeneous islet types)
- 6533 included

### Important: don't claim "scVI is batch-corrected"

Old comments in `reaggregate_islets.py` and `islet_analysis/CLAUDE.md` claimed "scVI latent has donor-level variation corrected out (Age+Gender covariates bijectively map to donors)". This is wrong — Age/Gender as `categorical_covariate_keys` don't activate batch correction; only `batch_key` does. Comments updated Apr 2026.

## Senior trajectory batch-correction + robustness fix (Jul 2026)

The **deployed Senior** trajectory (`data/islets_core_senior.h5ad`, **9,130 islets × 59
markers** — not the stale Junior 5,214/31 the older docs describe; read by the app via the
symlink `data/app_data/adata_ins_root.h5ad → ../islets_core_senior.h5ad`) was found to be a
**donor/slide batch axis, not disease progression**: batch η²=0.988 vs disease η²=0.103;
within-status donor η² ND 0.99 / T1D 0.97; ND<Aab+<T1D ordering FAILED (Aab+ 0.76 > T1D 0.65);
root was a saturated **3-cell** islet. Fixes (land now, execute at the next rebuild):

- **`train_scvi_senior.py`** — the root cause. Now `setup_anndata(batch_key="imageid")` (was
  `categorical_covariate_keys=["Age","Gender"]`, **no batch_key** → batch decoder never fired,
  per-image staining/gain drift leaked into `X_scVI_mean`). Age/Gender dropped (collinear with
  imageid). **Input must be REDSEA+RESTORE-corrected MFI** (`--corrected-mfi-dir data/cells_redsea`
  overlays it by `object_id`, or PHENO's `raw/X` already holds it) so `batch_key` isn't doing
  double duty against spillover REDSEA already removes.
- **`compute_trajectory_senior.py`** (rewritten) — ports the proven Junior recipe:
  (1) **Harmony θ-sweep** `{0,0.1,0.5,1.0}` on `X_scVI_mean` (batch=imageid) → neighbors on
  `X_scVI_harmony_pt` (k=30, cosine); picks the smallest θ passing the gates (harmonypy `Z_corr`
  is `(d,N)` — always `.T`). (2) **Robust root** = medoid of top-50 ND-by-INS islets with
  `total_cells ≥ 20` (kills the 3-cell root). (3) **`low_cell_qc` flag** + a `≥5`-cell
  **sensitivity check** (dataset size unchanged → core↔peri pairing + app keying intact).
  (4) **Hard gates that RAISE before write** (fail → `*.rejected.h5ad` + exit 1, `--force` to
  override): INS ρ<−0.3; GCG ρ>0.1; **T1D>ND** median + one-sided MWU p<0.05 (**Aab+ excluded
  — single donor 6521**); within-status donor η² ≤0.15 (ND,T1D); root sanity. Metrics persist to
  `uns['trajectory_validation']` (surfaced as a validity caption in the Trajectory tab).
  (5) **Combined "Core+Peri" pseudotime** (Step-5b port) → `dpt_pseudotime_combined`
  (`STRUCT_PHENOTYPES` reconciled to the new broad-lineage taxonomy); written only if it passes
  the same gates. (6) Display `X_umap` from raw-marker PCA (not the donor-separated scVI blob).
- **Prototype validation** (on a stand-in composition latent, since the real scVI latent is
  currently missing): θ=1.0 took within-status η² ND 0.54→0.07 / T1D 0.92→0.09, batch η²
  0.99→0.18, ordering ND<Aab+<T1D restored, INS ρ −0.55, root a 41-cell multi-donor medoid.
- **Cohort caveat (surfaced, not fixed):** Aab+ = a single donor (6521, 399 islets) → excluded
  from the hard ordering/η² gates and labelled n=1 in the app.

## 22-donor cohort expansion (Jul 2026) — code prep for the full rerun

Cohort growing **15 → 22 donors** (batch 2: +6450/6505 = Aab+, +6523/6534/6566/6591/6623 = T1D;
new split ND 8 / Aab+ 3 / T1D 11). Plan: `~/.claude/plans/be-sure-the-pseudotime-gentle-brook.md`.
Code landed now; runs when the QuPath 22-donor exports are ready. Confirmed decisions: **58 common
markers** (batch 1 has CD38, batch 2 has b-Catenin1 — a physical panel change → drop both);
**full app rebuild with the 6 broad lineages** (Step-2 fine subclustering deferred); **REDSEA
re-runs all 22**; **archive 15-donor outputs by move + delete the 44 GB `redsea_scratch/`**.

- `scripts/senior/build_phenotyped_h5ad.py` (**new — GAP #1 fix**): no script produced
  `singlecell_protein_phenotyped.h5ad`; this joins `data/phenotype/broad/` into
  `singlecell_protein.h5ad` by `object_id` → `obs['phenotype'] = broad_lineage`. Verified: all
  23.28M cells match (0 Unassigned). **Follow-up:** `compute_islet_markers.py` /
  `build_composition_and_tissue.py` still threshold `X>0.5` (old scimap-rescaled X) — must be
  reconciled to RESTORE `_pos` since X is now raw MFI (RESTORE covers 10 markers; the rest need a
  per-marker threshold).
- `scripts/senior/build_singlecell_anndata.py`: now takes the **marker INTERSECTION** across donor
  partitions (`get_fragments().physical_schema`) → uniform 58-marker panel, no NaN column.
- `scripts/senior/build_viewer_assets.py`: hard-coded 15-donor `DONORS` list (**GAP #2**) replaced
  with `discover_donors()` auto-glob of `parquet/tissue/case_id=*` (falls back to `cells/`).
- `scripts/senior/compute_trajectory_senior.py`: within-status donor η² gate is now **data-driven**
  — gates any status with ≥2 donors, so **Aab+ (now 3 donors) enters the gate** automatically
  (was hard-excluded as a single donor); hard ordering stays ND<T1D, ND<Aab+<T1D is informational.
- `scripts/senior/archive_15donor.sh` (**new**): dry-run by default (`--apply` to execute); moves
  cohort/derived 15-donor outputs → `data/archive_15donor_<date>/`, deletes `redsea_scratch/` +
  the deprecated pre-REDSEA `restore/`/`restore_gated/` family. `cells/`+`cells_redsea/` left in
  place (regenerated over) unless `--archive-cells`.

## Viewer OME-TIFF conversion (Jul 2026)

`scripts/senior/convert_viewer_ometiff.sh` converts the 59-plex `.qptiff` slides
(`/home/smith6jt/IO60panc2nd/Images/`) to pyramidal, tiled OME-TIFF for the app's **Viewer** tab.
bioformats2raw 0.12.1 → raw2ometiff 0.10.0 (zip-app launchers under `scripts/tools/`, **no Docker**,
Java 21); `-s 0 --use-existing-resolutions` (Baseline series only, reuse qptiff pyramid) →
`raw2ometiff --compression LZW` → `/data/islet_ome_tiff/<case_id>.ome.tiff` (+ a
`<case_id>.offsets.json` Viv IFD index via `generate-tiff-offsets`). Idempotent/resumable; intermediate
Zarr on `/data/tmp_b2r` (`/home` is too small). `bash scripts/senior/convert_viewer_ometiff.sh all15`
(skips already-done). The **`Channel_names` sidecar had to be regenerated** to the 59-channel Senior
order (was the 35-plex Junior panel) — full detail incl. per-batch panel differences and the
**HTTP-Range (206) serving caveat** (prod nginx/shiny-server yes, dev `runApp`/httpuv no) in
`app/shiny_app/CLAUDE.md` § "Viewer Tab (Avivator) — producing the OME-TIFFs".

## Pipeline Scripts & Notebooks

- `scripts/fix_tcell_subtypes.py` - Reclassifies CD4/CD8 T cells by marker gating (CD8a/CD4 >= 0.5); `--dry-run` flag for preview
- `scripts/reaggregate_islets.py` - **Primary pipeline**: aggregation (min_cells=0) + trajectory + Leiden in one script
- `scripts/compute_neighborhood_metrics.py` - Computes per-islet peri-islet metrics from single-cell H5AD
- `scripts/extract_per_islet_cells.py` - Extracts per-islet cell CSVs for drill-down viewer
- `scripts/extract_per_donor_tissue.py` - Extracts per-donor tissue CSVs (ALL cells: core+peri+tissue) for Spatial tab scatter
- `scripts/build_h5ad_for_app.py` - Builds enriched H5AD (trajectory + groovy + neighborhood + Leiden)
- `scripts/build_master_excel.py` - Builds master_results.xlsx from groovy TSV exports
- `notebooks/scvi_qc_validation.ipynb` - Validates scVI batch correction (silhouette, LISI, UMAP)
- `notebooks/rebuild_trajectory.ipynb` - Regenerates trajectory from fixed pipeline (superseded by reaggregate_islets.py)
- `islet_analysis/fixed_islet_aggregation.py` - Core aggregation module (islet-level from single-cell, default min_cells=0)

### Spatial Card C — per-cell distance parquet (Jul 2026)

- `scripts/senior/build_cell_distance_parquet.py` — builds `data/parquet/cell_distance/case_id=*`
  (per-cell `donor_status, cell_region, phenotype, dist_islet` + 10 RESTORE `{m}_pos` booleans; 23.3M
  cells, ~244 MB) for the Spatial-tab **Card C right** plot (per-cell signed distance to nearest islet:
  Density by phenotype / RESTORE-marker+ / all, and a marker-Composition small-multiples bar view). It
  joins three sources per donor: **signed `dist_islet` + `object_id`** from `data/cells/donor_id=*`;
  **reliable `phenotype`** from `data/parquet/tissue/case_id=*` by **row order** (its own
  `classification` is ~9–90% NULL — unusable); and the **RESTORE `_pos` calls** for 10 lineage markers
  (INS, GCG, SST, CD3e, CD20, CD163, CD31, SMA, Vimentin, Pan_Cytokeratin) from
  `data/restore_gated_redsea/donor_id=*` by **`object_id`** (calibrated per-image threshold on
  REDSEA-corrected MFI — supersedes the earlier raw-intensity percentile). **Asserts** per-row
  coordinate agreement (`|Δ|<1 µm`) so a future reorder of cells↔tissue fails loudly. Donor→status from
  `islets_core_senior.h5ad` obs. Run: `python scripts/senior/build_cell_distance_parquet.py` (~2 min).
  Registered as the `cell_distance` DuckDB view (`app/shiny_app/R/00_globals.R`). Detail in
  `app/shiny_app/docs/spatial.md`.

## Python Environment

Pipeline scripts use the `scvi-env` conda environment:
```bash
conda activate scvi-env  # scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, harmonypy
```

`harmonypy==0.2.0` is a downstream batch-correction dependency (added Apr 2026). One-time install: `pip install harmonypy` inside `scvi-env`.

## Running Pipeline Scripts

```bash
conda activate scvi-env

# Fix T cell subtypes in canonical H5AD (dry-run first, then apply)
python scripts/fix_tcell_subtypes.py --dry-run
python scripts/fix_tcell_subtypes.py

# Full rebuild from single-cell H5AD (aggregate + trajectory + Leiden, ~15 min)
python scripts/reaggregate_islets.py

# Compute neighborhood metrics (reads 2.6M-cell H5AD, ~3 min)
python scripts/compute_neighborhood_metrics.py

# Extract per-islet cell CSVs (reads 2.6M-cell H5AD, ~5 min)
python scripts/extract_per_islet_cells.py

# Extract per-donor tissue CSVs for Spatial tab scatter (~2 min, independent of islet filtering)
python scripts/extract_per_donor_tissue.py

# Build enriched H5AD (merges trajectory + groovy + neighborhood + Leiden)
python scripts/build_h5ad_for_app.py
```

## RESTORE per-image intensity normalization (Jun 2026)

`scripts/senior/restore_normalize.py` applies **RESTORE** (Chang Lab / OHSU) to the
59-plex protein single-cell data for nine mutually-exclusive markers spanning the broad
pancreas cell types (`DEFAULT_MARKER_PAIRS`, `[target ← reference]`):
- **exocrine** Pan_Cytokeratin ← Vimentin · **mesenchyme** Vimentin ← Pan_Cytokeratin
  (epithelial/stromal Spearman ≈ −0.5, the broadest exclusivity in the tissue)
- **endocrine** INS (β) ← Pan_Cytokeratin · GCG (α) ← Pan_Cytokeratin
- **endothelium** CD31 ← Pan_Cytokeratin · **muscle/pericyte** SMA ← Pan_Cytokeratin
- **immune** CD20 (B) ← Pan_Cytokeratin · CD163 (mac) ← Pan_Cytokeratin · CD3e (T) ← CD163

Pan_Cytokeratin marks ~91% of cells (acinar/ductal) and is cleanly negative in every
non-epithelial population → the abundant universal reference. **CD99 was dropped** (96%+
detectable; not lineage-exclusive in this panel — was a poor islet proxy). Markers vetted
from data: CD56/CD34/EpCAM/SST/IAPP/CD11b/B3TUBB rejected (too broad / weak signal).

**Immune negative-control fix (Jun 2026) — the key correction.** The immune markers were
originally referenced *immune-vs-immune* (CD20←CD3e, CD163←CD3e). That was the bug: CD20/CD163
**background scales with cell brightness/autofluorescence** (CD20 MFI vs cytokeratin/DAPI Pearson
r≈0.5), but the dim, rare T-cell reference never saw the **bright-acinar** background → threshold
too low → **~55%+ of the Immune broad-lineage class were bright acinar FALSE POSITIVES** (donor
6539; `data/immune_false_positive_investigation.png`). Fix: reference CD20/CD163 against
**Pan_Cytokeratin** — the abundant, immune-negative, *bright* acinar population that represents the
brightness-scaling background. Result (`data/restore_reference_fix_validation.png`): bright-acinar
false positives removed (CD20⁺ −87%) AND the disease signal **sharpened** (CD20/CD163 T1D/ND ratio
1.3×→~2.5×); Immune broad composition ND 8.9%→3.9%, T1D 13.9%→10.5% (real infiltration preserved).
**CD3e is the deliberate EXCEPTION — it KEEPS CD163.** CD3e was already clean/specific (never a
bright-acinar false positive) and Pan_Cytokeratin is *positively* correlated with it (autofluorescence
coupling, ref-qc pearson ≈0.39), which drags the SSC threshold too low and **destroys the validated
T-cell infiltration biology** (real SSC cohort: T-gain 4.5×→0.9×; a crude mean+3σ approximation had
wrongly suggested it improved — always validate on the real SSC run, not the approximation). With
CD3e←CD163 the T-gain is preserved (ND 0.9%→T1D 4.1%, 4.5×). **Principle: Pan_Cytokeratin references
the markers whose background is brightness/autofluorescence-driven; a clean, specific marker keeps a
clean, mutually-exclusive (negatively-correlated) reference.** `reference_exclusivity_qc()` writes the
data-driven audit (`data/restore/reference_selection_qc.csv`: svd_ratio σ2/σ1, pearson sign, reference
abundance, peri-islet contamination check). Deliberate deviation from the RESTORE paper's literal
"reject positively-correlated pairs" rule: here a positive CD20↔Pan_CK correlation IS the
autofluorescence we want the negative control to capture, not antibody cross-reactivity.

RESTORE pairs each target with a mutually-exclusive *reference* whose positive cells are,
by biology, target-negative, so the target intensity in that population is background; per
image it splits background vs signal with KMeans/GMM/SSC and sets `threshold = mean + 3·σ`
of the target-negative cluster. Corrects per-image autofluorescence drift; operates on
**RAW MFI** (RESTORE's `idx_select` uses an absolute `>50` floor — do NOT log-transform).

**Vendored, not pip-installed** (its `setup.py` is a broken placeholder template):
```bash
git clone https://github.com/smith6jt-cop/RESTORE.git external/RESTORE   # pinned @ 38df59b
pip install spams-bin   # SSC model needs `spams`; prebuilt wheel matched cp313 in scvi-env
```
`external/` is gitignored. The script injects lightweight stubs for RESTORE's
figure/notebook-only imports (`holoviews`, `selenium`, `tqdm.notebook`) so the vendored
source runs headless **unmodified** with `save_figs=False`, and prepends
`external/RESTORE/python_code` to `sys.path` (resolves both `RESTORE` and the bundled
`ssc` sparse-subspace-clustering namespace package).

**Model choice — SSC (default).** All three models are computed and saved, but on this
panel KMeans/GMM put the abundant immune markers (CD3e/CD20/CD163) at/above the 99.9th
intensity percentile (~0.3–1% positive — implausible). SSC lands at the background→signal
shoulder (plausible per-image fractions). `--model` overrides.

**Robustness guard (default on).** Some image×marker SSC thresholds are degenerate —
abundant markers (Pan_Cytokeratin/Vimentin) overshoot to the tail, INS undershoots
(24/135 in the 9-marker run, e.g. Pan_CK/Vim → 3000–6400 vs ~850 median). `--robust`
(default) imputes any threshold `>factor×` or `<1/factor×` the cohort median (default
factor 3) with the cohort median, logged per imputation; raw per-model thresholds in the
CSV are left untouched. `--no-robust` keeps pure per-image RESTORE thresholds.

**Biological validation (independent).** Per-image RESTORE fractions recover both T1D
hallmarks: β-cell INS⁺ falls ND 3.4% → Aab+ 2.3% → T1D 0.5% (insulin loss), and T-cell
CD3e⁺ rises ND 0.8% → T1D 3.0% (immune infiltration). Note Pan_Cytokeratin⁺ (~20–35%) is
conservative for the exocrine compartment — `mean + 3σ` flags only high-confidence cells.

Mutually-exclusive pairs live in the `DEFAULT_MARKER_PAIRS` constant (`[target,
reference]`); override with `--marker-pairs 'CD3e:CD163,INS:Pan_Cytokeratin'`.

Reads `data/cells/donor_id=*/data_0.parquet` (raw MFI); writes:
- `data/restore/threshs.pkl` — native RESTORE `threshs[batch][scene][marker][model]`
- `data/restore_thresholds.csv` — tidy, all 3 models × 15 images × 9 markers (405 rows)
- `data/restore_gated/donor_id=*/data_0.parquet` — per cell, for each marker: `<m>_pos`,
  `<m>_norm` (= raw / threshold, so threshold → 1.0, harmonized across images), `<m>_log2r`
- `data/restore/qc/*.png` (per-image histograms + threshold + positive-fraction heatmaps)
  and `data/restore/positive_fractions.csv`

```bash
conda activate scvi-env
python scripts/senior/restore_normalize.py                                # full 9-marker run (~5 min)
python scripts/senior/restore_normalize.py --limit-scenes 1 --skip-apply  # single-image dry run
python scripts/senior/restore_normalize.py --reuse-threshs --no-robust    # re-apply, pure thresholds (fast)
python scripts/senior/restore_normalize.py --donors 6539 --no-robust \
    --cells-dir data/cells_redsea --out-dir data/restore_redsea ...       # re-RESTORE on REDSEA-corrected cells
```

## Pixel-level REDSEA spillover correction (Jun 2026)

`scripts/senior/redsea_full.py` (scvi-env; needs `imagecodecs` for LZW qptiff decode) corrects
**lateral signal spillover** across shared cell boundaries — the cause of the "pan-positive
bright" cells that produced the fake "Unassigned" fraction. It reimplements REDSEA (Bai et al.
2021) **scalably** (the reference `redseapy` uses a dense cellNum² matrix + whole-image imread +
pixel loops — impossible at 0.5–2.6M cells/image). Validated on 6539; see the three figures
`data/redsea_{full_validation,residual_investigation,doublet_investigation}.png`.

**Step A — per-cell geometry export (`~/IO60panc2nd/scripts/export_cells_geojson.groovy`).**
Headless QuPath Groovy: `getDetectionObjects()` → `exportObjectsToGeoJson(cells, path,
"FEATURE_COLLECTION", "EXCLUDE_MEASUREMENTS")`. Each feature `id` == `object_id` UUID; polygon
coords are full-res pixels. **Close the QuPath GUI first** (it holds ~400 GB RAM). QuPath `script`
accepts only ONE `--image`, so loop per donor:
```bash
while IFS= read -r img; do
  ~/QuPath/bin/QuPath script --project ~/IO60panc2nd/project.qpproj --image "$img" \
     ~/IO60panc2nd/scripts/export_cells_geojson.groovy
done < donor_images.txt          # -> data/redsea_scratch/geojson/cells__<image>.geojson
```

**Step B — rasterize + compensate (`redsea_full.py`).** Per donor: rasterize GeoJSON → int32
instance mask (skimage.draw, NOT cv2 which can't fill uint32) → read qptiff series0/level0 ONE
channel at a time (`tifffile`) → whole-cell + 1-px-boundary `np.bincount` sums → 8-connected
sparse contact matrix (gap-bridged: cells aren't gapless, ~23% of px are foreground) → compensate
→ `data/cells_redsea/donor_id=*/data_0.parquet` (`object_id` + 59 corrected means). Alignment
verified: full-res means reproduce QuPath's own at r≈0.998 (join on `object_id`, never coords).
```bash
python scripts/senior/redsea_full.py --donor 6539 --keep-mask          # one donor (~8 min, 59 ch)
python scripts/senior/redsea_full.py --donor 6539 --save-intermediates # dump data/edge/contact for fast re-tuning
python scripts/senior/redsea_full.py --donor 6539 --from-intermediates --alpha 0.5  # re-compensate (instant)
python scripts/senior/redsea_full.py --all
```

**Compensation math + the two non-obvious decisions (empirical, validated):**
`corrected_mean = clip(data − α·(F @ edge), 0) / cell_area_px`, where `data`/`edge` are whole-cell
/ 1-px-boundary channel SUMS, `F` = row-normalized contact (`F @ edge` is the *recipient* form —
NOT `F.T`).
- **Subtract-ONLY (`--comp-mode 0`, default), not textbook subtract+reinforce.** Reinforce
  (`+edge`) is a structural **no-op in confluent tissue** at every band thickness — touching cells
  share an ambiguous boundary so `+edge` cancels `−F@edge`. Reinforce only helps background-separated
  MIBI cells.
- **Thin 1-px band (`--edge-radius 0`, default) + `α=1`.** A thick band over-erodes small cells; the
  1-px rim is the natural boundary, and α=1 is full subtraction (no arbitrary tuning).

**Validate (`redsea_full_validate.py`)** with **broad mutually-exclusive compartments**
(Epithelial / Immune / Endocrine / Stromal-vascular) — NOT the 9 individual markers, because
Vimentin/CD31/SMA are not lineage-exclusive (CD31+SMA+Vimentin = a real vessel wall). On 6539:
impossible cross-compartment≥3 0.46→0.05% and multi-immune 0.59→0.06% (~89% removed), single-
compartment cells 90% preserved. Residual co-expression is ~90% real Vimentin⁺ vascular biology;
"Beta+Alpha doublets" are residual islet-core spillover (normal-sized, not merges), 76% resolved by
`argmax(INS,GCG)`. Then re-run `restore_normalize.py --cells-dir data/cells_redsea` (cohort robust)
and proceed to broad-lineage argmax phenotyping.

## Senior phenotyping — broad lineage assignment (Step 1, Jun 2026)

`scripts/senior/assign_broad_lineage.py` (scvi-env) types every one of the 23.3M cells into one of
six **mutually-exclusive** broad lineages with **ZERO "Unassigned"** (the user mandate: fix the
upstream cause, don't invent a junk category). Input: `data/restore_gated_redsea/` (REDSEA-corrected,
re-RESTORE'd `_pos`/`_norm`) + `data/cells/` (region/islet context). Output:
`data/phenotype/broad/donor_id=*/` + `data/phenotype/broad_lineage_composition.png`.

**Hierarchical rule** (not plain argmax — see below): (1) Endocrine if INS|GCG|SST `_pos` (β/α/δ; SST
added 2026-06-29 — δ-cells were else mis-typed Epithelial/Immune, +139k Endocrine cells); (2) Immune
elif CD3e|CD20|CD163 `_pos`; (3) Endothelial elif CD31 `_pos`; (4) else structural background =
argmax of (Pan_Cytokeratin, Vimentin, SMA) `_norm` → Epithelial / **Fibroblast** / Muscle, with
cells below ALL structural thresholds defaulted to Epithelial (`epi_default` flag). Why not a plain
six-way `_norm` argmax: RESTORE's mean+3σ flags only the brightest ~20% of any ABUNDANT marker, so
the dim acinar majority is sub-threshold and a naive argmax scatters them to low-threshold markers
(empirically ~90% of "Immune" calls were such cells). The `epi_default` cells are *validated*
epithelial (median Ker8_18 549 vs 113 in positively-typed cells); acinar/ductal/fibroblast resolution
is deferred to Step-2 clustering. Composition (mean %): ND Epi 70 / Fibroblast 10 / Immune 4 /
Endocrine 4 / Endothelial 9 / Muscle 3; T1D shows the immune gain (Immune 10.5%) + β-loss.

**Taxonomy note — "Fibroblast" (was "Mesenchymal" → "Stromal" → "Fibroblast").** The six classes must
be MUTUALLY-EXCLUSIVE siblings that survive subclassing. "Mesenchymal" names a developmental *state*
(the M in EMT), not a resident compartment; "Stromal" is an *umbrella* that overlaps Endothelial/
Muscle/Immune (stroma = the parenchyma's complement). "Fibroblast" is the exclusive terminal
compartment for the Vimentin⁺ connective-tissue cells (pancreatic stellate cells are its members).
The same rename was propagated to `make_phenotype_rules.py` (now archived under
`archive/phenotyping_legacy/`; the legacy hierarchical rules: the "Stromal" parent dissolved into
Muscle + Fibroblast top-level siblings). The app's
`drilldown_helpers.R` "Stromal" color is provisional/deployed-data and reconciled at the Step-5 app
rebuild (TODO breadcrumb left there). NB watch the numpy fixed-width-string gotcha: build the broad-call
array as `<U16` or assigning a longer label (e.g. "Endothelial") silently truncates.

## False-endocrine fix — hormone over-calling, NOT keratin spillover (Jul 2026)

The Step-1 identity heatmap's **Endocrine row lit for the keratin block** (Pan_CK/Ker8_18) turned out to be
**FALSE endocrine**, not keratin spillover onto real β-cells. RESTORE's per-image hormone threshold lands
**in the noise** for β-loss (T1D) donors — there is no separated bright INS population, so it forces a
threshold into the noise tail and calls thousands of **acinar** cells INS⁺ at `_norm≈1.1` (barely over
threshold), which then show their own real keratin. User-confirmed: donor **6380 has no INS cells**, yet the
un-fixed pipeline called **11,905 INS⁺** there (all `_norm≈1.15`, 100% in tissue, 0 islets with ≥5
co-INS⁺). The keratin-vs-strength relationship is monotonic (median corrected Pan_CK of INS⁺ cells: `_norm`
1–1.5 → **765**, 10+ real β → **126**), so weak calls **are** acinar. **Real β-cells (strong `_norm`) are
already keratin-clean** → the pixel-REDSEA keratin lever solves a non-problem.

**Discriminators (all agree, per-donor):** hormone `_norm` strength (real β 8–25 vs false ~1.1); islet
coherence (real endocrine cluster in islets, false scatter in tissue); per-image separation (p99/p90 of
corrected intensity: ND 6539 INS ~80, β-loss 6380 ~2). Per-donor biology is preserved (hard rule): T1D
**6593 keeps 22,904 strong-`_norm` INS⁺ = real residual β** — never normalized toward a cohort β-null.

**THE FIX = a threshold-relative hormone-strength floor** (per-image-adaptive since `_norm`=raw/threshold):
call Endocrine only when INS/GCG/SST `_norm ≥ K` (rejects the noise-floor false β) **+ bright CD99
`_norm ≥ 3`** for the PP/ε/EC endocrine we lack specific markers for. **Operating point K=5** (user choice),
**committed cohort-wide (22 donors)**:
- `scripts/senior/assign_broad_lineage.py` — gained argparse (`--gated-dir/--out-dir/--extra-gated-dir/
  --hormone-min-norm/--cd99-bright`); the endocrine rule is now `INS/GCG/SST _norm ≥ K | CD99_norm ≥ 3`
  (was `_pos`). Re-typed `data/phenotype/broad` at `--hormone-min-norm 5`.
- `scripts/senior/apply_hormone_floor.py` (**new**) — rewrites `{INS,GCG,SST}_pos = (_norm ≥ K)` in
  `restore_gated_redsea` (leaves `_norm`/`_log2r` and all other columns untouched); re-applyable after any
  RESTORE re-run. Applied at K=5 so `_pos` consumers (Spatial-tab `cell_distance`) see the fix too.
- **Bright CD99** needs `data/restore_gated_redsea_extra` (CD99/B3TUBB/MPO ← Pan_Cytokeratin): produce via
  `restore_normalize.py --marker-pairs 'CD99:Pan_Cytokeratin,B3TUBB:Pan_Cytokeratin,MPO:Pan_Cytokeratin'
  --cells-dir data/cells_redsea --gated-dir data/restore_gated_redsea_extra --no-robust --no-ref-qc`.

**Validation (per-donor, `scripts/senior/redsea_reassess_diag.py` = the yardstick):** 6380 INS-driven endocrine
11,905 → **31**; real β kept (6539 8,084, 6593 22,904); endocrine islet-coherence ND 65→**94%**, T1D 39→
**61%**; GCG(α) coherence 59→94% (removes scattered false α, keeps real); Endocrine-cell **keratin Pan_CK
−62% / Ker8_18 −56%, hormones INS +97% / GCG +71%**; re-plotted heatmap Endocrine row now hormone-defined
(INS/GCG/SST/CD99) with Pan_CK neutral. Figures in `data/redsea_reassess/`
(`false_endocrine_diag`, `endocrine_profile_fix`, `retype_fix_compare`, `redsea_reassess_baseline`).

**Hallmark metric fix:** the CD3e T-gain was hidden by the **mean** (6476 pancreatitis ND has real CD3e⁺
0.37, dominating the ND mean). Use the **median** → both β-loss (INS ND 0.028→T1D 0.005) and CD3e T-gain
(ND 0.007→T1D 0.024) show; 6476 is plotted separately, not hidden. **Pan_CK is NOT a definitive acinar
marker (overly bright)** — do not over-weight "acinar keratin retention" as a guardrail.

**Committed state (deployed data, gitignored):** old typing/gating archived to `data/phenotype/broad.pre_hormonefloor`
and `data/restore_gated_redsea.pre_hormonefloor`; `cell_distance` rebuilt from the floored `_pos`.
`scripts/senior/redsea_full.py` per-channel-α + receptivity-gate (the abandoned keratin lever) is kept
**additive and default bit-identical** (verified mean|diff| 7e-11 vs deployed) — off by default, can revert.
**`marker_taxonomy.py` + `plot_celltype_markers.py` were MOVED `scripts/` → `scripts/senior/`** (their
`Path(__file__).parents[2]/"data"` only resolves to `<repo>/data` from `scripts/senior/`, and they import
`marker_taxonomy` as a sibling). **Follow-up:** the app h5ad still serves the old provisional fine labels; it
picks up the corrected broad lineages at the next `build_phenotyped_h5ad`/Step-2 rebuild. `data/parquet/tissue`
is still 15-donor, so `cell_distance` covers 15 until the 22-donor tissue build.

## Senior phenotyping — QuPath visual inspection (Jun 2026)

To eyeball broad-lineage calls on the images in QuPath 0.7, match by the detection **UUID**
(`PathObject.getID()` == our `object_id`) — exact, no centroid rounding:
- `scripts/senior/export_broad_class_for_qupath.py` → `data/phenotype/qupath_class/pheno_class_<donor>.csv`
  (`object_id, broad_lineage, image`).
- `scripts/groovy/import_broad_lineage.groovy` (QuPath 0.7): reads that folder (hard-coded path,
  `donorsToImport` filter), matches the open image's detections by UUID, sets the PathClass + lineage
  colors directly, self-reports matched/missing. Non-destructive (in-memory unless saved). Replaces the
  fragile CytoMAP centroid-import + separate set-class step.
- `scripts/groovy/view phenotypes.groovy` (rewritten for 0.7): checkbox show/hide per class via
  `OverlayOptions.setPathClassHidden(...)` (the old `hiddenClassesProperty()` was removed in 0.7).
- Legacy `import clusters and cell neighborhoods.groovy` + `set_broad_lineage_class.groovy`
  + `export_phenotype_for_qupath.py` (measurement-based CytoMAP path) were **archived 2026-06-29**
  under `archive/phenotyping_legacy/` — superseded by the UUID path above.
QuPath 0.7 API gotchas (verified via `javap` on the jars): `MeasurementList.put(String,double)` (not
`putMeasurement`), `Dialogs.promptForDirectory` single-arg removed (hard-code the folder),
`PathClass.setColor(int,int,int)`.

## Reproducible pipeline notebook (2026-06-29)

`notebooks/senior_phenotyping_redsea_restore.ipynb` runs the whole pipeline end-to-end as a **thin
orchestrator**: each step shells out to the live `scripts/senior/*.py` (single source of truth — no
logic duplicated), wrapped in an idempotent `run_step()` that **skips any stage whose outputs already
exist** (flip the `FORCE` dict to recompute), bookended by cells that load the parquet outputs and
render QC figures. Kernel: `scvi-env`. The RESTORE step passes the REDSEA dirs explicitly
(`--cells-dir data/cells_redsea --out-dir data/restore_redsea --gated-dir data/restore_gated_redsea`).
The final cell is a **cell-type × marker dotplot** (9 RESTORE markers; dot size = % `_pos`, color =
mean `_norm`) **+ heatmap** (full 59-marker panel, z-scored per marker). Verified end-to-end:
23,325,515 cells, **0 "Unassigned"**, `epi_default` 50.8%.
**Gotcha:** 13/15 `data/cells_redsea/` donors carry a few NaN edge-cells (59–531 each); a plain
`numpy.sum` NaN-poisons the whole marker column, so per-cell aggregation must use `np.nansum` +
finite-count (`restore_gated_redsea` has no NaN). The executed copy (`*_executed.ipynb`) is gitignored.

## Archived / deprecated phenotyping files (2026-06-29)

The predecessors of the unified pipeline were moved to **`archive/phenotyping_legacy/`** (git history
preserved; nothing in `scripts/senior/` imports them). They were the three superseded approaches that
produced the fake "Unassigned" fraction — per-donor **auto-gating** (`auto_gate_per_donor.py`,
`phenotype_scimap.py`), global **Leiden clustering** (`cluster_annotate.py`, `finalize_phenotype.py`),
and the **Astir/consensus** typer (`run_astir.py`, `astir_markers.yaml`, `rules_to_astir_yaml.py`) —
plus their support scripts (`make_phenotype_rules.py`, `marker_reliability_audit.py`,
`redsea_approx.py` table-level proof, `validate_phenotypes.py`, `phenotype_rules_senior.csv`), the
legacy QuPath companions (`export_phenotype_for_qupath.py`, `set_broad_lineage_class.groovy`,
`import clusters and cell neighborhoods.groovy`), and the Junior-era phenotyping notebooks
(`CODEX_Panc_scimap_Analysis.ipynb`, `CODEX_Panc_Analysis.ipynb`, `napari_gating_senior.ipynb`,
`Comb_Xen_Phen.ipynb`, `phenotype_rules1.csv`, `phenotype_rules2.csv`). See
`archive/phenotyping_legacy/README.md` for the full what-replaced-what table. The end-to-end live
pipeline is reproduced in **`notebooks/senior_phenotyping_redsea_restore.ipynb`**.
