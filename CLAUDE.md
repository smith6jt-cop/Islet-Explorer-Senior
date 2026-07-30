# CLAUDE.md - Islet Explorer

This file provides guidance to Claude Code when working with this repository.

## Project Overview

Interactive Shiny app for exploring pancreatic islet CODEX multiplexed imaging data from nPOD donors. Tracks insulin loss and immune infiltration across ND -> Aab+ -> T1D progression. Features islet-level aggregated views, pseudotime trajectory, spatial neighborhood analysis, and single-cell drill-down.

## ⚠️ Working rule — no frequency/uniformity assumptions

**Never assume marker or cell-type frequencies, and never treat per-donor variation as an error to
normalize toward a cohort norm — donor biology is genuinely variable.** e.g. donor **6476** (lone
**pancreatitis** donor) has real extreme CD3e; **CD3e is sparser than CD163** (not the reverse);
positivity is **RESTORE-threshold-defined**, not a raw-MFI floor. Validate a per-donor threshold/call
against **that donor's own data** (its histograms, negative-control population, spatial coherence),
not against other donors. Justify method changes on methodological grounds (stability/representativeness),
not on "this donor should look like the others."

**In-flight work:** the 15→22 donor rebuild is mid-pipeline — see [`docs/HANDOFF_22donor_rebuild.md`](docs/HANDOFF_22donor_rebuild.md) for current state, next steps, and the note that `data/` + raw inputs are NOT in git.

## Sub-directory CLAUDE.md files

- `app/shiny_app/CLAUDE.md` — per-tab feature details (Drill-Down, Trajectory, Plot, AI Assistant) + conventions including the May 2026 font-standardisation system
  - `app/shiny_app/docs/spatial.md` — Spatial Neighborhood Analysis deep reference (split out May 2026)
  - `app/shiny_app/docs/statistics.md` — Statistics Tab deep reference (split out May 2026)
- `scripts/CLAUDE.md` — pipeline scripts, batch correction, pseudotime modes, **the Senior phenotyping pipeline (REDSEA spillover → RESTORE normalization → broad-lineage assignment → QuPath inspection)**, and the **Senior trajectory batch-correction + robustness fix (Jul 2026)** — scVI `batch_key`, Harmony θ-sweep, robust root, hard gates, combined mode
  - (`islet_analysis/` was removed; its content folded into `scripts/senior/` + `scripts/CLAUDE.md`)

## App Architecture (Modular - Feb 2026)

The Shiny app lives in `app/shiny_app/` and uses a modular architecture:

```
app/shiny_app/
  app.R                        # Entrypoint (~510 lines): UI layout, CSS, server wiring, root-level outputs
  R/
    00_globals.R               # Package loads, constants, feature flags, paths (h5ad_path, master_path)
    data_loading.R             # load_master_auto(), prep_data(), bin_islet_sizes(), H5AD loading
    utils_stats.R              # summary_stats(), per_bin_anova(), per_bin_kendall(), cohens_d(), eta_squared(), pairwise_wilcox()
    utils_safe_join.R          # safe_left_join(), add_islet_key(), compute_diameter_um()
    segmentation_helpers.R     # GeoJSON cache, build_segmentation_base_plot(), render_islet_segmentation_plot()
    drilldown_helpers.R        # PHENOTYPE_COLORS, load_islet_cells(), render_islet_drilldown_plot(), render_drilldown_summary()
    viewer_helpers.R           # Viewer backend (Avivator | TissUUmaps) URL builders, channel config, environment detection
    ai_helpers.R               # OpenAI credential loading, call_openai_chat()
    mod_plot_ui.R / mod_plot_server.R           # Plot tab (scatter, distribution, outliers, segmentation+drilldown panel)
    mod_trajectory_ui.R / mod_trajectory_server.R  # Trajectory tab (UMAP, heatmap, pseudotime, segmentation+drilldown panel)
    mod_viewer_ui.R / mod_viewer_server.R       # Viewer tab (OME-TIFF iframe; ISLET_VIEWER_BACKEND=avivator|tissuumaps)
    mod_statistics_ui.R / mod_statistics_server.R  # Statistics tab (5-section narrative)
    spatial_helpers.R           # Per-donor tissue CSV loader with new.env() caching
    mod_spatial_ui.R / mod_spatial_server.R     # Spatial tab (3-panel: sidebar controls, tissue scatter, Leiden panel)
    mod_ai_assistant_ui.R / mod_ai_assistant_server.R  # AI chat panel
  www/                         # Static assets (images, logos)
  app_original.R               # Backup of pre-modularization monolithic app (5,397 lines)
```

### Data Directory Layout (Senior, relocated Jun 2026)

`data/` is now a **real directory at the repo root** (was `~/IO60panc2nd/senior_data`, one level up; moved in and gitignored — never committed; ~21 GB). It holds the raw/intermediate pipeline files at top level (`singlecell_protein.h5ad`, `cells/` parquet, `islet_explorer_senior.h5ad`, `cluster_*`, `marker_audit_*`, …) **plus** a curated **`data/app_data/`** subdirectory that is the app's view (per-islet drill-down CSVs, GeoJSONs, `islet_spatial_lookup.csv`, and relative symlinks `islet_explorer.h5ad`→`../islet_explorer_senior.h5ad`, `adata_ins_root.h5ad`→`../islets_core_senior.h5ad`, `parquet`, `phenotype_rules.csv`). **The app reads from `data/app_data/`** (R paths use `file.path("..","..","data","app_data", …)`); the **Senior pipeline scripts** read/write `data/` top-level via `Path(__file__).resolve().parents[2] / "data"` (move-proof, no absolute paths).

### Data Loading (H5AD with Excel Fallback)

The app uses `load_master_auto()` which tries H5AD first, falls back to Excel:
- `h5ad_path` -> `data/app_data/islet_explorer.h5ad` (Senior; -> `../islet_explorer_senior.h5ad`)
- `master_path` -> `data/app_data/master_results.xlsx` (Excel fallback; absent in Senior)

Both produce the same `list(markers, targets, comp, lgals3, phenotypes, donor_demographics, neighborhood)` structure consumed by `prep_data()`. The `phenotypes`, `donor_demographics`, and `neighborhood` elements are extracted from H5AD `.obs` and are `NULL` when loading from Excel (graceful fallback).

### Performance: Reticulate Bridge Optimization (Mar 2026)

**Critical**: `ad$uns[[key]]` crosses the Python-R bridge on EVERY call (~221ms each). With ~62 groovy columns × 3 sheets = ~186 calls, this took ~16s. Fix: cache `uns <- ad$uns` once (0.2s) then index the R list. `reconstruct_groovy_df_from_list(uns, sheet)` replaces per-key bridge crossings.

**Deferred loading**: `annotations.tsv` (72 MB, 576K rows) was loaded eagerly at startup but only serves as a fallback when `islet_spatial_lookup` doesn't have the islet (never happens). Now lazy-loaded via `.seg_lazy` environment on first fallback access.

Total startup: **~22s → ~3.5s**.

### Shared Reactive State (wired in app.R server)

- `prepared()` - core data reactive from H5AD or Excel, consumed by Plot, Trajectory, Statistics, Spatial
- `selected_islet` - reactiveVal with `case_id`, `islet_key`, `centroid_x/y`, `donor_status`, `donor_age`, `donor_gender`
- `forced_image` - reactiveVal, written by Trajectory, read by Viewer
- `active_tab` - `reactive(input$tabs)`, passed to Plot and Trajectory to prevent duplicate output IDs
- `active_donor_colors` - reactive from `donor_palette_name()`, synced across 3 non-namespaced inputs: `sidebar_donor_palette` (Plot sidebar), `traj_donor_palette` (Trajectory), `spatial_donor_palette` (Spatial)

### Critical: Root-Level Outputs + Active Tab Guard

The following render outputs are defined in **app.R** at the root level, NOT inside any module:
- `output$islet_segmentation_view` - GeoJSON boundary plot (shared by Plot + Trajectory)
- `output$islet_drilldown_view` - Single-cell spatial overlay plot
- `output$islet_drilldown_summary` - Cell composition bar chart
- `output$islet_drilldown_table` - Core/peri cell counts

Both Plot and Trajectory tabs embed these using non-namespaced `plotOutput("islet_segmentation_view")` etc. **IMPORTANT**: Each module's `segmentation_viewer_panel` renderUI must guard with `if (active_tab() != "Plot") return(NULL)` (or `"Trajectory"` respectively). Without this guard, both renderUIs fire when `selected_islet()` changes, creating duplicate DOM IDs — Shiny can only bind one, causing the other to grey out. Shiny `renderUI` fires based on reactive dependencies, NOT tab visibility.

### Non-Namespaced Inputs for Root-Level Outputs

The segmentation panel `renderUI` in both modules generates non-namespaced inputs:
- `drilldown_view_mode` - "Boundaries" or "Single Cells" toggle
- `drilldown_color_by` - "phenotype" or a marker column name
- `drilldown_show_peri` - Whether to show peri-islet cells

These are read by the root-level `renderPlot` outputs in `app.R`.

### Plotly Click Events in Modules

- Plot module: `source = ns("plot_scatter")` with `event_register("plotly_click")`
- Trajectory module: `source = ns("traj_scatter")` with `event_register("plotly_click")`
- Both use namespaced sources; `event_data()` calls must use matching `ns()` prefix
- **macOS native-click bridge (Jun 2026)**: `plotly_click` is suppressed on macOS, so both scatters also pipe `islet_click_bridge(ns("<plot|traj>_native_click"))` (native mousedown/mouseup → `Shiny.setInputValue`) after `event_register`. The new server observers and the legacy plotly_click handler both call the shared `select_islet_from_key()` (`segmentation_helpers.R`). See `app/shiny_app/CLAUDE.md` § "macOS Native-Click Bridge + Shared Resolver".

## Key Data Files

- `data/islet_explorer.h5ad` - **Primary app data** (~70 MB, 5,214 islets, groovy + trajectory + donor + neighborhood + Leiden)
- `data/master_results.xlsx` - Aggregated islet-level data (composition, markers, targets) -- Excel fallback
- `data/adata_ins_root.h5ad` - Trajectory h5ad (5,214 islets, 31 vars, DPT pseudotime, UMAP)
- `data/neighborhood_metrics.csv` - Per-islet peri-islet metrics (5,214 rows, 62 columns)
- `data/cells/*.csv` - Per-islet single-cell CSVs for drill-down (~5,214 files, ~203 MB total)
- `data/donors/*.csv` - Per-donor tissue-wide cell CSVs for Spatial tab scatter (15 files, ~78 MB total)
- `data/islet_spatial_lookup.csv` - Centroid coordinates for segmentation viewer
- `data/json/*.geojson` / `data/gson/*.geojson.gz` - QuPath segmentation boundaries
- `data/DATA_PROVENANCE.md` - Documents canonical H5AD lineage and data sources

## Data Pipeline

### Canonical H5AD Lineage

```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, 31 markers)
  |
  |-- scripts/reaggregate_islets.py  (one-stop: aggregate + trajectory + Leiden)
  |     -> islet_analysis/islets_core_fixed.h5ad  (5,214 islets, min_cells=0, require_paired=True)
  |     -> islet_analysis/islets_peri_fixed.h5ad   (5,214 peri-islet)
  |     -> islet_analysis/islets_merged_fixed.h5ad (5,214 paired core+peri)
  |     -> data/adata_ins_root.h5ad  (+ core + combined pseudotimes, UMAP, X_scVI_harmony_pt, X_scVI_combined_harmony_pt)
  |     -> islet_analysis/islets_core_clustered.h5ad  (+ Leiden at 4 resolutions, Harmony θ=2 on imageid)
  |
  |-- scripts/compute_neighborhood_metrics.py
  |     -> data/neighborhood_metrics.csv  (5,214 rows, ~100 cols — per-phenotype enrich_z and min_dist)
  |
  |-- scripts/extract_per_islet_cells.py
  |     -> data/cells/*.csv  (5,214 files, 37 cols: coords, phenotype, region, 31 markers)
  |
  |-- scripts/extract_per_donor_tissue.py
  |     -> data/donors/*.csv  (15 files, 5 cols: X/Y coords, phenotype, cell_region, islet_name)
  |
  +-- scripts/build_h5ad_for_app.py (trajectory + groovy + donor + neighborhood + Leiden)
        -> data/islet_explorer.h5ad  (all app data in one file, incl. leiden_* + leiden_umap_*)
```

`scripts/reaggregate_islets.py` replaces the manual notebook workflow (aggregation + trajectory + Leiden in one script). `islet_analysis/islets_core_clustered.h5ad` provides Leiden clustering (4 resolutions) + UMAP coords, merged by `build_h5ad_for_app.py` step 4.7.

See `scripts/CLAUDE.md` for per-script details, batch correction (Phase 17), and the two-pseudotime architecture (`dpt_pseudotime` core-only vs `dpt_pseudotime_combined` core+peri+structures).

### Senior phenotyping pipeline (current, Jun 2026) — supersedes the provisional auto-labels

The Senior 59-plex data is being **re-phenotyped from scratch** because global clustering produced a large fake "Unassigned" fraction caused by **upstream artifacts**, fixed in order **REDSEA → RESTORE → broad lineage** (full detail in `scripts/CLAUDE.md`):

```
data/cells/donor_id=*/  (raw 59-marker MFI + object_id UUID, 23.3M cells)
  |-- QuPath per-cell GeoJSON export (~/IO60panc2nd/scripts/export_cells_geojson.groovy)
  |-- scripts/senior/redsea_full.py        pixel-level lateral-spillover correction (subtract-only, 1px band, α=1)
  |     -> data/cells_redsea/donor_id=*/    (59 spillover-corrected means; alignment r≈0.998)
  |-- scripts/senior/restore_normalize.py   per-image + autofluorescence normalization (mutually-exclusive refs)
  |     -> data/restore_redsea/, data/restore_gated_redsea/donor_id=*/  ({m}_pos/_norm/_log2r)
  |     -> Pan_Cytokeratin is the universal negative control; CD20/CD163←Pan_CK, CD3e←CD163 (immune fix)
  +-- scripts/senior/assign_broad_lineage.py  every cell typed, 0 Unassigned (hierarchical _pos + structural argmax)
        -> data/phenotype/broad/donor_id=*/  (8 classes: Epithelial/Fibroblast/Immune/Endocrine/Endothelial/Muscle/Neural/Neutrophil)
        -> QuPath inspection: export_broad_class_for_qupath.py + scripts/groovy/import_broad_lineage.groovy
```

**Status (updated 2026-07-08):** REDSEA + RESTORE done & validated (**22 donors**); broad lineage (Step 1) done and **re-typed to 8 classes** (adds Neural←B3TUBB, Neutrophil←MPO, + a bright-CD99 Endocrine gate) via the **`Phenocycler_Analysis` submodule** — now the canonical steps-1–5 engine (`python -m phenocycler.pipeline`, driven by repo-root `config.ini` / `scripts/senior/run_pipeline.sh`). The parity gate confirmed the submodule is byte-identical to the committed `assign_broad_lineage.py`. **Next: Step 2 per-lineage subclustering.** Validation figures live under `data/` (`redsea_*`, `broad_lineage_composition.png`, `phenotype/celltype_marker_{dotplot,heatmap}.png`, `redsea_reassess/`). The app still serves the **provisional** un-QC'd labels until the Step-5 app rebuild.

**Reproducible notebook (Jun 2026):** `notebooks/senior_phenotyping_redsea_restore.ipynb` runs the full pipeline end-to-end (idempotent thin orchestrator over the live `scripts/senior/*.py`) and ends with a **cell-type × marker dotplot + heatmap**. The superseded gating/clustering/Astir scripts were archived to **`archive/phenotyping_legacy/`** (see `scripts/CLAUDE.md`). App tab-by-tab capability is tracked in `app/shiny_app/docs/capability_status.md`.

## Viewer backend (Avivator | TissUUmaps, Jul 2026)

The Viewer tab can render through **Avivator** (default) or **TissUUmaps** — flip with
`ISLET_VIEWER_BACKEND=tissuumaps`. TissUUmaps adds region drawing, GeoJSON import/export and
per-cell marker overlays, and tiles images itself (no HTTP-Range requirement, so the Viewer
also works under `shiny::runApp`). It consumes **`.tmap` projects** built per donor by
`scripts/senior/build_tissuumaps_project.py`, which names each channel layer from the slide's
own OME-XML — required because the cohort mixes panels (batch-1 SST@35 vs batch-2 SST@1).
Findings, costs and the deployment recipe: [`docs/tissuumaps_evaluation.md`](docs/tissuumaps_evaluation.md);
app wiring in `app/shiny_app/CLAUDE.md` § "Viewer Tab — pluggable backend".

## Spatial / Statistics / Trajectory / AI tabs

Per-tab implementation details — Spatial Neighborhood Analysis (Phase 7), Single-Cell Drill-Down (Phase 8), Statistics Tab pseudoreplication fix (Phase 16), Cell-Count-Weighted Trajectory (Phase 12), Pseudotime Heatmap, AI Assistant (Phase 10) — live in `app/shiny_app/CLAUDE.md`.

## Important Conventions

- `ggplot2::coord_sf()` and `ggplot2::geom_sf()` - these are ggplot2 functions, NOT sf functions
- `PIXEL_SIZE_UM = 0.3774` - micrometers per pixel conversion constant
- UMAP "Selected Feature" uses continuous viridis (inferno) colormap scaled to data min/max
- Donor status colors: ND = steel blue (#4477AA), Aab+ = burnt umber (#CC6633), T1D = forest green (#228833) — Paul Tol bright variant, colorblind-safe, centralized as `DONOR_COLORS` in `00_globals.R`
- Phenotyping uses Rules1 (`data/phenotype_rules.csv`) - 21 phenotypes, 18 markers. CSV headers use adata.var_names (CD8a, HLADR, PDPN — not CD8, HLA-DR, Podoplanin)
- **Donor metadata**: Comes from `islets_core_fixed.h5ad` obs, NOT from `CODEX_Pancreas_Donors.xlsx` (different cohort)
- **User-facing terminology**: Use "Sex" not "Gender" in all UI labels, plot titles, and methods text. The underlying data column remains `gender` for backwards compatibility.
- **H5AD obs index**: `islets_core_fixed.h5ad` index name is `islet_id` -- same as column, use `reset_index(drop=True)`
- **Case ID zero-padding**: GeoJSON files use `0112.geojson` (4-digit padded), data uses `112` (unpadded). `load_case_geojson()` and click handlers use `sprintf("%04d", ...)` fallback
- **Log-scale with zeros**: Use `scales::pseudo_log_trans(base=10)` instead of `scale_y_log10()` -- zeros map to 0 (visible) instead of -Infinity (dropped)
- **Plot defaults**: Point size = 3.0, transparency = 0.6
- **Plot fonts & downloads**: Every plot uses `theme_islet()` (ggplot) and `plotly_islet_fonts()` + `plotly_islet_config()` (plotly) — single source of truth in `app/shiny_app/R/00_globals.R`. Base 16pt / titles 18pt / axes 14-16pt / legends 14-16pt / heatmap rows 12pt min. Plotly modebar PNG downloads ship 2800 × 1800 px (scale=2). ggplot `renderPlot` outputs that should be downloadable are refactored into a `reactive()` consumed by both `renderPlot` and a `downloadHandler` using `ggsave(..., width=10, height=7, dpi=300)` (heatmaps 12 × 8). See `app/shiny_app/CLAUDE.md` § "Plot Fonts & Publication-Quality Downloads".
- **Phenotype marker rules in UI**: `data/phenotype_rules.csv` is parsed once at startup into `PHENOTYPE_RULES_CSV` in `00_globals.R`, with hierarchical composition (sub-phenotypes inherit parent constraints) and an alias map for the four PHENOTYPE_COLORS-vs-CSV naming mismatches (`CD8a Tcell` ↔ `CD8 T cell`, etc.). Surfaced via (1) inline blue hint under every phenotype `selectInput` showing the active rule, (2) global "Phenotype rules — marker definitions" button at the top of every tab opening a reference modal listing all 21 phenotypes + markers, (3) expandable rules list in the drill-down composition panel naming each phenotype present in the clicked islet. See `app/shiny_app/CLAUDE.md` § "Phenotype Rules — Marker Definitions in the UI".
- **Peri-islet data guard**: Always check `total_cells_peri > 0` before using peri metrics
- **Column name sanitization**: Phenotype names use `_` for spaces, `plus` for `+` in peri-islet columns
- **Single-cell Parent column**: `Islet_N` = core cells, `Islet_N_exp20um` = peri-islet cells
- **macOS plotly clicks & Viewer iframe (Jun 2026)**: `plotly_click` is suppressed on macOS — drill-down uses a native-click bridge (`islet_click_bridge()` in `00_globals.R`). The Viewer iframe is mounted once and its `src` updated in place via `session$sendCustomMessage` to avoid recreating the WebGL canvas (macOS WebRender flicker). Details + the `reactiveVal`-dedup gotcha in `app/shiny_app/CLAUDE.md`.

App-specific UI conventions (CSS overflow, font minimums, large-scatter rendering, palette syncing, reticulate bridge, synthetic-union rows, macOS click bridge & Viewer iframe persistence) live in `app/shiny_app/CLAUDE.md`.

## Deployment Architecture

- **Production URL**: `http://<server-ip>:8080/islet-explorer/`
- nginx (port 8080) reverse-proxies `/islet-explorer/` -> shiny-server (port 3838)
- shiny-server serves from symlink: `/srv/shiny-server/islet-explorer` -> `app/shiny_app/`
- Changes take effect when shiny-server spawns a fresh R worker (no restart needed; kill stale workers if needed)
- Dev server (`Rscript -e 'shiny::runApp(".", port=7777)'`) is for local testing only

## Running the App

Production (always use this):
```
http://<server-ip>:8080/islet-explorer/
```

Development only:
```bash
cd app/shiny_app
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Running Pipeline Scripts

See `scripts/CLAUDE.md` for the full per-script command list. Quick rebuild:

```bash
conda activate scvi-env
python scripts/reaggregate_islets.py          # ~15 min: aggregation + trajectory + Leiden
python scripts/compute_neighborhood_metrics.py # ~3 min
python scripts/build_h5ad_for_app.py          # ~1 min: merge into app H5AD
```

## Analysis Pipeline

See `scripts/CLAUDE.md` for the Python pipeline; phenotyping steps 1–5 (REDSEA → RESTORE → hormone floor → 8-class broad lineage) live in the `Phenocycler_Analysis` submodule. (The old `islet_analysis/` dir was removed; its content folded into `scripts/senior/` + the submodule.)
See `data/DATA_PROVENANCE.md` for full data lineage documentation.
See `docs/user_guide.md` for end-user documentation.

## Cross-Session Knowledge

Cross-session development knowledge for this project lives in Claude's auto-memory (see `MEMORY.md`). The earlier `Skills_Registry/` submodule + `/advise` + `/retrospective` workflow was retired Apr 2026 in favour of the auto-memory mechanism.

## R Dependencies

shiny, shinyjs, plotly, ggplot2, dplyr, tidyr, readxl, sf, jsonlite, RColorBrewer, scales, anndata, reticulate, httr2, stringr, lmerTest, lme4, emmeans
