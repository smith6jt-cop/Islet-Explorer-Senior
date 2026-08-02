# Scripts guide

## Ownership boundary

`Phenocycler_Analysis/` owns the raw QuPath-to-cell-identity workflow. Root
scripts will consume completed assignments after a run-aware H5AD adapter is
implemented; current app scripts still consume a legacy H5AD. Root scripts
must not invent marker thresholds, positivity frequencies, or alternate
cell-typing rules.

Current source stages:

<<<<<<< HEAD
```text
manifest -> ingest -> geometry -> redsea -> expression -> controls
         -> calibrate -> type -> states -> QuPath export
=======
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
**full app rebuild with the 8 broad lineages** (Neural + Neutrophil added; Step-2 fine subclustering deferred); **REDSEA
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

## TissUUmaps viewer projects (Jul 2026)

`scripts/senior/build_tissuumaps_project.py` turns one donor OME-TIFF into a **TissUUmaps
`.tmap` project** — the alternative Viewer backend (`ISLET_VIEWER_BACKEND=tissuumaps`). Per
donor it extracts each channel to a pyramidal 8-bit TIFF under `channels/<case>/<idx>_<marker>.tif`
(idempotent, skip-if-exists, `--force` to redo), then writes `<case>.tmap` with **layers named
from the slide's own OME-XML** — mandatory because the cohort mixes panels (batch-1 59-ch
SST@35 vs batch-2 58-ch SST@1), so a positional config is wrong for batch 2 — plus
INS/GCG/SST/DAPI default-on in the Avivator colours, `compositeMode: lighter`, `mpp`, and
optional islet-GeoJSON regions + per-cell markers.

Do **not** rely on the fork's `/slide` multi-channel auto-split: it names layers
`<file>_Channel_<i>` (marker identity lost — its OME-XML name regex matches `<Image Name=...>`
first) and converts synchronously inside the HTTP request. Details + the two other bugs found
(region autoLoad blocked by a layer-index prompt unless features carry
`properties.collectionIndex`; greyscale JPEG-in-TIFF unreliable on libvips 8.15.1 → default
`--compression deflate`) in [`docs/tissuumaps_evaluation.md`](../docs/tissuumaps_evaluation.md).

**Coordinate flags matter** — slide is full-res px (`--um-per-px 0.2539`), QuPath GeoJSON is
resolution-#1 px (`--region-scale 2.0`), per-cell CSVs are µm (`coord_factor` derived
automatically). Verify against one donor before batch-converting.

```bash
python scripts/senior/build_tissuumaps_project.py --image /data/islet_ome_tiff/6539.ome.tiff \
    --out-dir /data/tissuumaps --regions data/app_data/json/6539.geojson \
    --markers data/app_data/cells/6539_tissue.csv
python scripts/senior/tissuumaps_smoke_test.py            # 19 end-to-end checks (synthetic)
```

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
>>>>>>> origin/main
```

Run and inspect them through the package:

```bash
conda activate scvi-env
python -m pip install -e Phenocycler_Analysis
python -m phenocycler.pipeline run --config config.ini
python -m phenocycler.pipeline status --config config.ini
```

`senior/run_pipeline.sh` is a convenience wrapper for `run`, `status`, and
`export`; it adds no workflow logic. Set `PHENOCYCLER_CONFIG` only when a
different configuration file is intentional.

See:

- `../Phenocycler_Analysis/README.md`
- `../Phenocycler_Analysis/docs/WORKFLOW.md`
- `../Phenocycler_Analysis/DATA_README.md`

## Downstream scripts and current boundary

The existing aggregation, trajectory, spatial, and application scripts operate
on a legacy single-cell H5AD. They remain useful for the deployed application,
but they are not yet a supported consumer of the content-addressed source
workflow and must not be presented as a current end-to-end rebuild.

Before those scripts are used with new assignments, implement one canonical
run-to-H5AD adapter that reads:

```text
<run>/cells
<run>/expression
<run>/marker_evidence
<run>/assignments
<run>/cell_states
```

The adapter must join on `(donor_id, object_id)`, preserve uncertainty, and
record the run ID and relevant manifest content IDs. Existing H5AD phenotype
columns are not valid calibration or typing inputs. `fix_tcell_subtypes.py` is
legacy H5AD maintenance only and must not modify current assignments.

## Environment separation

- `scvi-env` is the established downstream Scanpy/scVI/Harmony environment.
- `Phenocycler_Analysis/environment.yml` is the focused reproducible source
  workflow environment.
- Do not add embedding, trajectory, Shiny, or app-only dependencies to the
  focused source environment.

## Donor and identity contracts

- Central exclusions for `6457` and `6579` are enforced by
  `phenocycler.cohort`; scripts must not recreate or override the list.
- Join source artifacts by exact `object_id`.
- Do not use row order, rounded centroids, filenames, or provisional
  `classification` values as identity.
- Missing/invalid marker evidence stays unavailable or indeterminate.
- Disease labels and annotations are audit-only.

## Retired phenotyping code

The following approaches are historical and must not be used as live inputs or
fallbacks:

- RESTORE normalization and divisor/floor searches;
- NNMF/GMM/Otsu Gate 2–5 screens;
- ordered residual broad-lineage assignment;
- `data/cells_redsea*`, `data/restore*`, and `data/phenotype/broad*` as current
  source-stage roots;
- Gate 5 dashboards/counterfactuals and anatomy-check release arms;
- hard-coded QuPath import/export paths and the old broad-lineage importer.

The active QuPath handoff is
`Phenocycler_Analysis/scripts/groovy/import_cell_assignments.groovy`, which
uses UUIDs and preserves uncertainty. Historical experiments remain
recoverable from Git history.

After importing, `groovy/view_phenotypes.groovy` can show or hide any exact
broad, specific, ambiguous, or unavailable QuPath class.

## Viewer and application assets

OME-TIFF conversion and browser-specific viewer behavior are downstream
application concerns. Keep channel order explicit, retain HTTP Range support
in production, and follow `../app/shiny_app/CLAUDE.md` for viewer, plotting,
reticulate, and deployment details.

Generated data, figures, run directories, caches, and executed notebooks stay
out of source control.
