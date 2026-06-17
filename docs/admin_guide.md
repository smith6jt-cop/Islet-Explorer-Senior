# Islet Explorer — Administration & Operations Guide

Developer- and maintainer-facing reference for rebuilding the app data, deploying the Shiny app, and troubleshooting. End-user feature documentation lives in [user_guide.md](user_guide.md).

**Production URL**: `http://<server-ip>:8080/islet-explorer/`

---

## Data Pipeline

### Rebuilding the App Data

The canonical pipeline rebuilds `data/islet_explorer.h5ad` from source:

```bash
# Activate the Python environment
conda activate scvi-env

# Step 1: Reaggregate islets + trajectory + Leiden (~15 min)
python scripts/reaggregate_islets.py
# → islet_analysis/islets_core_fixed.h5ad (5,214 islets)
# → data/adata_ins_root.h5ad (+ pseudotime + UMAP)
# → islet_analysis/islets_core_clustered.h5ad (+ Leiden at 4 resolutions)

# Step 2: Compute neighborhood metrics (from 2.6M-cell single-cell H5AD)
python scripts/compute_neighborhood_metrics.py
# → data/neighborhood_metrics.csv (5,214 rows × ~100 cols, per-phenotype enrich_z_* and min_dist_*)

# Step 3: Extract per-islet cell CSVs (for drill-down viewer)
python scripts/extract_per_islet_cells.py
# → data/cells/*.csv (5,214 files, ~203 MB)

# Step 4: Extract per-donor tissue CSVs (for Spatial tab scatter, independent of islet filtering)
python scripts/extract_per_donor_tissue.py
# → data/donors/*.csv (15 files, ~78 MB)

# Step 5: Build enriched H5AD (trajectory + groovy + neighborhood + Leiden + donor metadata)
python scripts/build_h5ad_for_app.py
# → data/islet_explorer.h5ad (~70 MB)
```

### Full Lineage

```
CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, canonical single-cell)
  ↓ scripts/reaggregate_islets.py (min_cells=0, require_paired=True)
islets_core_fixed.h5ad  (5,214 islets, proteins + scVI embeddings + trajectory + Leiden)
  ↓ scripts/build_h5ad_for_app.py (+ groovy + neighborhood + donor metadata)
data/islet_explorer.h5ad  (complete app data, incl. Leiden clustering)
```

Branch pipelines:
- `scripts/compute_neighborhood_metrics.py` — peri-islet composition, immune metrics, enrichment z-scores, distances
- `scripts/extract_per_islet_cells.py` — individual cell CSVs for drill-down
- `scripts/extract_per_donor_tissue.py` — per-donor tissue CSVs for Spatial tab scatter
- Leiden clustering from `islet_analysis/islets_core_clustered.h5ad` (4 resolutions) — merged by `build_h5ad_for_app.py`

### Upstream Data Sources

| Source | Location | Description |
|--------|----------|-------------|
| Single-cell H5AD | `single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad` | 2.6M cells, scVI-embedded (note: scVI training was not batch_key-aware; downstream Harmony in reaggregate_islets.py corrects donor variance for clustering and combined-mode pseudotime), phenotyped |
| Groovy TSV exports | `~/panc_CODEX/results/groovy_exports/` *(Junior legacy — Senior uses no groovy)* | QuPath islet measurements (15 donors × 4 types); Senior derives composition/markers directly from the phenotyped single-cell H5AD (`~/IO60panc2nd/senior_data/singlecell_protein_phenotyped.h5ad`) instead |
| GeoJSON boundaries | `data/json/*.geojson` | Islet segmentation polygons |
| Spatial lookup | `data/islet_spatial_lookup.csv` | Islet centroid coordinates |

See `data/DATA_PROVENANCE.md` for the complete lineage documentation including validation results.

---

## Deployment & Administration

### Architecture

```
nginx (port 8080)
  ↓ reverse-proxy /islet-explorer/
shiny-server (port 3838)
  ↓ symlink /srv/shiny-server/islet-explorer → app/shiny_app/
Shiny app (R worker process)
```

### Updating the App

Code changes take effect when shiny-server spawns a new R worker. No restart is needed for code-only changes.

If a stale R worker persists with old code:
```bash
# Find and kill stale Shiny R workers
ps aux | grep R | grep shiny
kill <PID>
```

Note: `sudo systemctl restart shiny-server` is not available — use the worker kill approach above.

### Local Development

```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
# Then open http://localhost:7777 in browser
```

Port 7777 is for local testing only and is not accessible to end users.

### R Dependencies

```
shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite,
RColorBrewer, scales, anndata, reticulate, lmerTest, lme4, emmeans
```

### Python Dependencies (pipeline only)

```bash
conda activate scvi-env
# scanpy, anndata, scvi-tools, sklearn, scib-metrics, scipy, pandas, numpy
```

---

## Troubleshooting

### App fails to load

1. Check that `data/islet_explorer.h5ad` exists and is ~70 MB
2. Check that R packages are installed: `library(anndata); library(reticulate)`
3. Check that Python is accessible via reticulate: `reticulate::py_available()`
4. Fall back to Excel: ensure `data/master_results.xlsx` exists (reduced functionality)

### No phenotype or demographic filters

The app is running from Excel fallback. Rebuild the H5AD:
```bash
conda activate scvi-env
python scripts/build_h5ad_for_app.py
```

### No neighborhood metrics in Plot selector

1. Verify `data/neighborhood_metrics.csv` exists (5,214 rows)
2. Rebuild the H5AD: `python scripts/build_h5ad_for_app.py`
3. Kill stale R workers and refresh the browser

### Single-cell drill-down not available

1. Verify `data/cells/` directory contains ~5,023 CSV files
2. If missing, regenerate: `python scripts/extract_per_islet_cells.py`

### Segmentation viewer shows "No GeoJSON found"

1. Check that `data/json/` or `data/gson/` has `.geojson` / `.geojson.gz` files
2. Case ID zero-padding: GeoJSON files use 4-digit padded IDs (`0112.geojson`), data uses unpadded (`112`). The app handles this automatically via `sprintf("%04d", ...)` fallback.

### Statistics tab shows no results

The Statistics tab uses data from the Plot sidebar. Ensure:
1. A feature is selected in the Plot sidebar
2. At least 2 disease groups are checked
3. The diameter range includes some islets

### Spatial tab tissue scatter is empty

1. Verify `data/donors/` directory contains 15 CSV files
2. If missing, regenerate: `python scripts/extract_per_donor_tissue.py`
3. If Leiden panel says "not available", rebuild the H5AD with Leiden data: `python scripts/build_h5ad_for_app.py` (requires `islet_analysis/islets_core_clustered.h5ad`)

### AI assistant shows error or no response

1. **"LLM key not found"**: Set `KEY=your-api-key` and `BASE=https://api.ai.it.ufl.edu` in `~/.Renviron` or `app/shiny_app/.Renviron`
2. **"key not allowed to access model"**: The selected model isn't available. Verify available models with `GET /v1/models`. Current models: `gpt-oss-20b` and `gpt-oss-120b`.
3. **Authentication error (401)**: Confirm the API key is valid and has no extra whitespace
4. **Timeout or empty response**: The Large model (120b) needs time for reasoning. Try the Fast model (20b) for quicker responses, or increase the timeout.
5. **"httr2 package required"**: Install with `install.packages('httr2')` and restart the app
6. **Debug mode**: Set `DEBUG_CREDENTIALS=1` in `.Renviron` and check the R console output for credential loading details

### Leiden UMAP shows a blob

The Leiden UMAP visualization uses raw marker PCA coordinates (same as the Trajectory tab). If it shows a blob, the clustered H5AD may have outdated scVI-based UMAP coordinates. Fix:
1. Copy visualization UMAP from `data/adata_ins_root.h5ad` into `islet_analysis/islets_core_clustered.h5ad`
2. Rebuild: `python scripts/build_h5ad_for_app.py`
3. Kill stale R workers and refresh

### Drill-down does nothing or Viewer flickers on macOS (fixed Jun 2026)

Earlier builds had two macOS-only problems: clicking a plot point did nothing (macOS suppresses plotly's `plotly_click` event), and the Viewer tab flickered/reloaded continuously (the WebGL `<iframe>` was destroyed and recreated on every reactive update — heavy flicker under Firefox's WebRender compositor on macOS, Mozilla bug #1555544). Both are fixed:

- **Drill-down** uses a native-click bridge (`islet_click_bridge()` in `app/shiny_app/R/00_globals.R`) that fires on a native mouse click even when `plotly_click` is suppressed; both the Plot and Trajectory scatters route through the shared `select_islet_from_key()` resolver.
- **Viewer** mounts the iframe once and updates its `src` in place via `session$sendCustomMessage`, so the WebGL canvas is never recreated.

If these symptoms reappear on macOS, confirm the deployment is running the latest `app/shiny_app/R/` code (kill stale R workers, hard-refresh the browser). On the Plot tab, enable **Show individual points** so there are per-islet points to click. If the Viewer still flickers specifically on a Retina display, the residual cause is deck.gl's device-pixel backing buffer — rebuild Avivator (`www/avivator`) with `useDevicePixels: false` (the bundle only reads `image_url` + `channel_config` query params, so it cannot be toggled at runtime).
