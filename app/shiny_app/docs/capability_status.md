# Islet Explorer (Senior) — App Capability Status

**Assessed:** 2026-06-29 · **Method:** the dev server was booted
(`Rscript -e 'shiny::runApp(".", port=7777)'`, `RETICULATE_PYTHON` → `scvi-env`) and every tab was
loaded and screenshotted in a headless browser. Verdicts below are what the **running app actually
did**, not a code read. Screenshots are in `docs/img/capability/`.

## Headline

- The app **boots cleanly** and **4 of 5 tabs are functional** on the data currently on disk.
- **Trajectory is broken** (data-join failure — see below). Static code review had called it
  "functional"; the running app says otherwise. This is the one regression to fix.
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
| **Trajectory** | ❌ **Broken** | `trajectory.png` | Controls + legend render, but the panel reads **"Loaded AnnData (9130 obs). Joined 0 rows (0 unique islets) of 0 islet metrics. Source: adata_ins_root.h5ad"** and the scatter shows **"No trajectory data available"**; both UMAPs are empty. The islet metrics join returns 0 rows. |
| **Viewer** | ⚠️ Empty (graceful) | `viewer.png` | "No images found. Place OME-TIFF files in `www/local_images/` or set `LOCAL_IMAGE_ROOT`." UI renders; no images deployed. Not a bug — missing optional data. |
| **Statistics** | ✅ Works | `statistics.png` | Full narrative + Configure Analysis controls + Run Statistics / Download CSV. Result panels (Hypothesis Testing, Effect-Size Forest, AUC) populate **on demand** after Run Statistics; the donor-aware (N=15) mixed-effects workflow is wired. |
| **Spatial** | ✅ Works | `spatial.png` | Tissue scatter renders **145,651 foreground cells** for donor 115 via the **DuckDB tissue backend** (no per-donor CSV fallback needed); Leiden UMAP (9,130 islets, clusters 0–14), Donor-Status UMAP, and "Mean Phenotype Composition by Cluster" all render. |
| **Drill-Down** *(embedded in Plot/Trajectory)* | ✅ Works via Plot | — | Per-islet segmentation + single-cell overlay is driven by a click on a scatter point. Functional from **Plot** (cell CSVs + GeoJSONs present); unusable from **Trajectory** only because that scatter is empty. |
| **AI Assistant** *(Plot side-panel)* | ✅ Renders | `plot.png` | "Panc Floyd" panel + model selector (Navigator `gpt-oss-20b`) shown and ready. Live chat needs UF Navigator API creds in `.Renviron` (not validated here). |

## The Trajectory failure — actionable detail

The trajectory module loads `adata_ins_root.h5ad` (→ `islets_core_senior.h5ad`, **9,130 obs**) but the
join to the per-islet metrics yields **0 rows / 0 unique islets**, so there is nothing to plot. The
same 9,130-islet data renders fine on the **Spatial** tab's Leiden UMAP, so the AnnData itself is
intact — the defect is in the **islet-key join** between the trajectory source and the metrics table
(an islet-ID/key mismatch in the Senior data, or the joined metric columns are absent). Fixing the
join key (or rebuilding the trajectory h5ad so its islet index matches the metrics) should restore
the tab. This is the first thing to check before any app data refresh.

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
