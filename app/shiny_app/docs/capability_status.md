# Islet Explorer (Senior) — App Capability Status

**Assessed:** 2026-06-29 · **Method:** the dev server was booted
(`Rscript -e 'shiny::runApp(".", port=7777)'`, `RETICULATE_PYTHON` → `scvi-env`) and every tab was
loaded and screenshotted in a headless browser. Verdicts below are what the **running app actually
did**, not a code read. Screenshots are in `docs/img/capability/`.

## Headline

- The app **boots cleanly** and **all 5 data tabs are functional** on the data currently on disk
  (Viewer needs images — see below).
- **Trajectory works** — an earlier "broken" call was a **premature 8s screenshot**. The tab was
  merely **slow to first render** (~10–30s): `traj_data_clean` ran a 9,130-iteration
  `get_islet_annotations` loop, each scanning the full 9.6k-row spatial lookup (~88M comparisons per
  render). **Fixed 2026-06-29** by vectorizing that lookup into one `match()` (renders instantly).
  The "islet-key join" was a red herring — the keys (`Islet_755` both sides) match perfectly.
- **Viewer is empty** by design (no OME-TIFF images deployed) — it degrades gracefully.
- **Every tab still serves the OLD provisional fine phenotypes** (Beta cell, Alpha cell,
  M2 macrophage, NK cell, Regulatory T cell, Stromal, …). The new REDSEA/RESTORE/broad-lineage
  outputs are **intentionally not wired into the app** — they live only in the phenotyping notebook
  (`notebooks/senior_phenotyping_redsea_restore.ipynb`). No app data rebuild was requested.

Startup log confirmed: DuckDB views `islets` + `tissue` registered; 4 mirai workers; 9,677 islet
spatial records; 35 channel names; `www/local_images` absent.

## Per-tab matrix

| Tab | Verdict | Evidence (screenshot) | Notes |
|-----|---------|-----------------------|-------|
| **Plot** | ✅ Works | `plot.png` | Sidebar filters, "Insulin⁺ fraction vs Islet Size" scatter (clean ND↘ / Aab+↗ / T1D-low biology), donor-group violin, AI panel — all render. Drill-down controls present (9,706 per-islet CSVs + 15 GeoJSONs on disk). |
| **Trajectory** | ✅ Works (was slow) | `trajectory.png` (8s, mid-load) → after 30s the pseudotime scatter (INS vs pseudotime, ND/Aab+/T1D + trend lines) + both UMAPs render fully. Slow first render was vectorized away 2026-06-29. |
| **Viewer** | ⚠️ Empty (graceful) | `viewer.png` | "No images found. Place OME-TIFF files in `www/local_images/` or set `LOCAL_IMAGE_ROOT`." UI renders; no images deployed. Not a bug — missing optional data. |
| **Statistics** | ✅ Works | `statistics.png` | Full narrative + Configure Analysis controls + Run Statistics / Download CSV. Result panels (Hypothesis Testing, Effect-Size Forest, AUC) populate **on demand** after Run Statistics; the donor-aware (N=15) mixed-effects workflow is wired. |
| **Spatial** | ✅ Works | `spatial.png` | Tissue scatter renders **145,651 foreground cells** for donor 115 via the **DuckDB tissue backend** (no per-donor CSV fallback needed); Leiden UMAP (9,130 islets, clusters 0–14), Donor-Status UMAP, and "Mean Phenotype Composition by Cluster" all render. |
| **Drill-Down** *(embedded in Plot/Trajectory)* | ✅ Works via Plot | — | Per-islet segmentation + single-cell overlay is driven by a click on a scatter point. Functional from **Plot** (cell CSVs + GeoJSONs present); unusable from **Trajectory** only because that scatter is empty. |
| **AI Assistant** *(Plot side-panel)* | ✅ Renders | `plot.png` | "Panc Floyd" panel + model selector (Navigator `gpt-oss-20b`) shown and ready. Live chat needs UF Navigator API creds in `.Renviron` (not validated here). |

## The Trajectory "failure" — what it actually was (and the fix)

Not a bug in the data or the join — a **performance** problem plus a premature screenshot. Diagnosis
(2026-06-29):

- The trajectory data is fully valid: all **9,130** islets pass the filters for INS (`dpt_pseudotime`
  finite 0–1, INS expression finite, `islet_key`/`case_id` present). The server log even shows
  `[traj_data_clean] Successfully processed 9130 observations for INS`.
- The islet keys match perfectly on both sides (`Islet_755`), so there is **no join defect** —
  `prep_data()` already calls `add_islet_key()` on `comp`, and `traj()` synthesizes
  `combined_islet_id`. The "0 islet metrics" line in the first screenshot was the **initial frame**
  before `input$traj_feature` propagated.
- The real cost was `traj_data_clean` calling `get_islet_annotations()` **once per islet** (9,130×),
  each scanning the full 9.6k-row `islet_spatial_lookup` — ~88M comparisons **on every render**
  (feature change, filter change). First render took ~10–30s, so the 8s screenshot caught an empty
  tab.

**Fix (`mod_trajectory_server.R`):** replaced that per-row loop with a single vectorized `match()`
against the lookup (identical result, O(n)). The full render (pseudotime scatter + both UMAPs) was
verified working via screenshot; the `match()` returns the same `centroid_x/y_um`/`has_spatial`
columns as the loop (same first-match semantics, same lookup columns), so it only changes speed.

## Other observations

- **No Excel fallback.** `data/app_data/master_results.xlsx` is absent, so `load_master_auto()` has a
  single point of failure: if the H5AD / reticulate-`anndata` path ever fails, the app won't load.
  (It loaded fine here with `RETICULATE_PYTHON` pointed at `scvi-env`.)
- **DuckDB backend is doing real work** — the Spatial tissue scatter (145k+ cells) comes from
  `data/parquet/tissue` via partition pruning, confirming `ISLET_USE_DUCKDB=TRUE` is engaged.
- **OLD-vs-new phenotypes.** The Spatial phenotype list still shows the full provisional taxonomy
  (Acinar, Ductal, CD4/CD8 T cell, M2 macrophage, NK cell, Regulatory T cell, Stromal, …). The
  notebook's broad-lineage output (6 classes) is a separate, QC'd artifact and is not surfaced in the
  app by design.

## Reproduce

```bash
cd app/shiny_app
RETICULATE_PYTHON=/home/smith6jt/miniconda3/envs/scvi-env/bin/python \
  Rscript -e 'shiny::runApp(".", port = 7777, host = "127.0.0.1")'
# then load http://127.0.0.1:7777 and click through Plot / Trajectory / Viewer / Statistics / Spatial
```
