# Islet Explorer Senior

The big-data, multimodal successor to Islet Explorer — an interactive R Shiny application for exploring nPOD pancreas data at scale (PhenoCycler/CODEX protein now; Xenium 5K spatial transcriptomics and post-Xenium IF to follow). Backed by a DuckDB + Parquet data layer (with mirai async and optional rdeck WebGL) so it scales from thousands of islets toward millions of cells. Tracks insulin loss, immune infiltration, and cellular microenvironment changes across the ND → Aab+ → T1D disease progression.

<div align="center">

[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue.svg)](https://github.com/smith6jt-cop/Islet-Explorer-Senior)
[![R Shiny](https://img.shields.io/badge/R-Shiny-blue.svg)](https://shiny.rstudio.com/)

</div>

## Overview

Islet Explorer is a modular R Shiny application serving **two deliberately separate datasets**, one
at a time, selected by the `ISLET_DOMAIN_CONFIG` environment variable. They share five donors and no
vocabulary, so they are **never merged**: a cross-donor comparison spanning both would compare
methods, not biology.

| domain | config | data root | cohort |
|---|---|---|---|
| legacy CODEX | `apps/islet/config/islet_default.yml` (the default) | `data/app_data/` | 5,214 islets, 15 nPOD donors, 31 protein markers |
| PhenoCycler run cohort | `apps/islet/config/run20_default.yml` | `data/development/typing/app/` | 20 nPOD donors, subpanel classification |

Every data path resolves through `ieCore::domain_data_root()`; a path that escapes the active
domain's root serves the wrong cells for a right-looking islet, and the testthat suite asserts
against exactly that. The sections below describe the **legacy CODEX domain** unless they say
otherwise.

### Key Capabilities

- **Tabs**: `app.R` registers nine hidden nav panels — Plot, Trajectory, Statistics, Spatial,
  Neighborhoods, CrossModal, LymphNode, Viewer, plus Inflammation, which is hidden unless
  `SHOW_INFLAMMATION_TAB=1`. The AI panel lives inside Plot and is off unless `ISLET_ENABLE_AI=TRUE`
  (hard-blocked in production regardless)
- **Single-cell drill-down**: Click any islet to view segmentation boundaries and cell composition
- **Pseudotime trajectory**: Diffusion pseudotime (DPT) with INS root, cell-count-weighted visualization
- **Spatial neighborhood analysis**: Immune infiltration, enrichment z-scores, distance metrics
- **Leiden clustering**: 4 resolution levels with UMAP visualization
- **AI-powered analysis**: Chat assistant via UF Navigator API
- **Configurable palettes**: Paul Tol (default), Bright, Okabe-Ito, Diverging — synced across all tabs
- **Demographic filters**: Age, gender, autoantibody status

### Disease Groups

| Group | Color (default) | Description |
|-------|-----------------|-------------|
| **ND** | Steel Blue (#4477AA) | Non-diabetic control donors |
| **Aab+** | Burnt Umber (#CC6633) | Autoantibody-positive (pre-T1D) |
| **T1D** | Forest Green (#228833) | Type 1 diabetic donors |

## Quick Start

```bash
# Install R dependencies
Rscript scripts/install_shiny_deps.R

# Run the application — from apps/islet, which must stay exactly two directories deep
cd apps/islet && ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  Rscript -e 'shiny::runApp(".", port = 7777)'
```

Omit `ISLET_DOMAIN_CONFIG` for the legacy CODEX domain. Use **R 4.5.1** with the `~/R/…/4.5`
library; `/usr/bin/Rscript` is 4.6.1 and silently lacks the dependencies.

**Production URL**: `http://<server-ip>:8080/islet-explorer/`

The development server is **not** a deployment test: Avivator needs HTTP Range (206), which
`shiny::runApp` does not serve. See `deploy/DEPLOY.md`.

## Installation

### Prerequisites

**R Shiny App (R 4.5.1, library `~/R/…/4.5`):**

| Package | Purpose |
|---------|---------|
| `shiny`, `shinyjs` | App framework |
| `readxl` | Excel fallback data loading |
| `dplyr`, `tidyr`, `stringr` | Data wrangling |
| `ggplot2`, `plotly` | Visualization |
| `broom`, `jsonlite` | Stats and JSON |
| `yaml` | Validated islet/lymph DOMAIN configuration |
| `scales` | Axis formatting |
| `sf` | GeoJSON segmentation viewer (optional) |
| `anndata`, `reticulate` | H5AD loading, trajectory (optional) |
| `httr2` | AI assistant API calls (optional) |
| `vitessceR` | Embedded OME-TIFF viewer (optional) |

**Python Pipeline (conda `scvi-env`):**
- `scanpy`, `anndata`, `scvi-tools`, `sklearn`, `scib-metrics`, `scipy`

### Install

```bash
# R dependencies (automated)
Rscript scripts/install_shiny_deps.R

# Optional: Install Bioconductor packages for trajectory + viewer
INSTALL_VITESSCER=1 Rscript scripts/install_shiny_deps.R

# Optional: Avivator image viewer static bundle
./scripts/install_avivator.sh
```

## Application Tabs

### Plot Tab

Primary exploration interface with sidebar controls shared with Statistics.

- **Scatter mode**: Feature expression vs islet diameter (31 markers, targets)
- **Composition mode**: Cell type proportions with 4 option groups:
  - Hormone Fractions (INS, GCG, SST)
  - Cell Type Proportions (21 phenotypes)
  - Peri-Islet Proportions (20 μm expansion zone)
  - Immune Metrics (infiltration ratios, density)
- **Filters**: Donor status, autoantibody (GADA, IA2A, ZnT8A, IAA, mIAA), age range, gender, diameter range
- **Normalization**: None, global z-score, robust per-donor (median/MAD×1.4826)
- **Click interaction**: Click any point to open islet segmentation + single-cell drill-down

### Trajectory Tab

Pseudotime analysis computed via diffusion pseudotime (DPT) with INS as root anchor.

- **UMAP scatter**: Color by donor status, selected feature (viridis inferno), or pseudotime
- **Cell-count weighting**: Point size by `sqrt(total_cells)` — honestly represents measurement quality
- **Weighted LOESS trends**: `log1p(total_cells)` weights prevent noisy small islets from distorting fits
- **Multi-feature heatmap**: Z-scored expression along pseudotime, clamped [-2.5, 2.5]
- **Click to drill-down**: Opens segmentation viewer and can jump to OME-TIFF Viewer tab

### Statistics Tab

Five-section narrative layout with shared sidebar from Plot tab:

1. **Configure Analysis** — Test type, alpha, outlier handling, bin width, diameter range
2. **Primary Results** — Hypothesis table + forest plot (Cohen's d with CI)
3. **Size-Dependent Patterns** — BH-corrected stratified test heatmap + Kendall τ trend analysis
4. **Confounders & Deeper Analysis** — Demographics + AUC analysis
5. **Methods Reference** — Dynamic statistical methods description

### Spatial Tab

Three-panel layout with neighborhood analysis cards:

- **Tissue scatter** (ggplot2, ~177K cells/donor): Foreground/background layering, brush-to-zoom, phenotype filtering
- **Leiden panel**: UMAP of 5,214 islets colored by Leiden cluster (4 resolutions: 0.3/0.5/0.8/1.0) + stacked bar chart
- **Neighborhood analysis cards**:
  - **Immune Infiltration**: Violin plots + peri vs core scatter
  - **Immune Cell Enrichment**: Grouped bar chart (7 immune types × 3 stages) + heatmap
  - **Immune Proximity**: Distance box plots + enrichment vs distance scatter

### Viewer Tab

Embedded OME-TIFF viewer on either of two backends, chosen by `ISLET_VIEWER_BACKEND`
(`avivator` default, `tissuumaps`). The cohort mixes panels, so default-on channels resolve **by
marker name, never by index**, against each image's own OME-XML.

### Neighborhoods, CrossModal and LymphNode Tabs

Neighborhoods holds the spatial neighbourhood cards. CrossModal and LymphNode are **thin card
registries** — `IE_CROSSMODAL_CARDS` / `IE_LYMPHNODE_CARDS` in `apps/islet/R/mod_*_ui.R`. Read the
current card list from those files; see `apps/islet/CLAUDE.md`.

### AI Assistant

Chat interface powered by UF Navigator API (`gpt-oss-20b` / `gpt-oss-120b` reasoning model). Supports streaming responses with token budget tracking.

### Single-Cell Drill-Down

Available from both Plot and Trajectory tabs when clicking an islet:

- **Segmentation view**: GeoJSON boundary polygons from QuPath
- **Single-cell overlay**: 5,214 per-islet CSV files with 37 columns (coordinates, phenotype, 31 markers)
- **Composition summary**: Horizontal bar chart + cell count table (core/peri breakdown)
- **Color by**: Phenotype (21 categories) or any protein marker

## Data

### Primary Data File (legacy CODEX domain)

`data/app_data/islet_explorer.h5ad` — single enriched file containing all app data for that domain:

- 5,214 islets from 15 nPOD donors
- 31 protein markers (mean expression per islet)
- Pseudotime trajectory (DPT with INS root, scVI-corrected)
- 21 phenotype proportions (from single-cell phenotyping)
- Donor demographics (age, gender, autoantibody status)
- 62 neighborhood metrics (peri-islet composition, immune infiltration, enrichment z-scores, distances)
- Leiden clustering at 4 resolutions + UMAP coordinates
- QuPath groovy data (targets, markers, composition, LGALS3)

**Fallback**: `data/app_data/master_results.xlsx` (4 sheets: Islet_Markers, Islet_Targets, Islet_Composition, LGALS3) — no phenotype proportions, demographics, or neighborhood metrics.

### Supporting Data Files

All of these live under the domain root — `data/app_data/` for legacy CODEX,
`data/development/typing/app/` for the run cohort. Never hardcode either; resolve through
`ieCore::domain_data_file()`.

| File | Description |
|------|-------------|
| `adata_ins_root.h5ad` | Trajectory source (legacy: 5,214 islets, pseudotime, UMAP) |
| `cells/*.csv` | Per-islet single-cell CSVs (legacy: 5,214 files) |
| `donors/*.csv` | Per-donor tissue CSVs for Spatial tab (legacy: 15 files) |
| `neighborhood_metrics.csv` | Per-islet peri-islet metrics (legacy: 5,214 rows, 62 columns) |
| `islet_spatial_lookup.csv` | Islet centroid coordinates |
| `json/*.geojson` | QuPath segmentation boundaries |
| `gson/*.geojson.gz` | Compressed segmentation boundaries |
| `phenotype_rules.csv` | Phenotype rules parsed once at startup (legacy: 19 phenotypes, 15 markers) |

### Data Conventions

**Islet keying**: All joins use (`Case ID`, `Donor Status`, `islet_key`). `islet_key` is derived from region names (e.g., `Islet_200` from `Islet_200_core`).

**Region types**: `islet_core`, `islet_band`, `islet_union` — no other spellings.

**Islet diameter**: `islet_diam_um = 2 * sqrt(core_region_um2 / pi)`

**Region synthesis**: If union rows are missing — Targets: synthesize `count` only. Markers: synthesize `n_cells` and `pos_count` as sums, compute `pos_frac`.

**Factor ordering**: ND < Aab+ < T1D (always).

**Peri-islet expansion**: every expansion is 20 µm — derive it from the signed distance, not from
the annotation name. Say *peri-islet expansion*, not ring or band, outside the legacy region labels
above.

## Data Pipeline

The upstream source workflow is the manifest-controlled
QuPath → ingest → geometry QC → REDSEA → expression pipeline in
[core/phenocycler/README.md](core/phenocycler/README.md). It stops after
expression: marker calling and typing are blocked until a verified typing
workflow exists (root `CLAUDE.md`). Run it from `core/phenocycler/`, which owns
the only `config.ini`:

```bash
python -m pip install -e core/phenocycler
cd core/phenocycler
python -m phenocycler.pipeline run --config config.ini --pipelined --through expression
python -m phenocycler.pipeline status --config config.ini
```

What runs today, lane by lane, is mapped in
[docs/WORKFLOW_2026-09-02_MAP.md](docs/WORKFLOW_2026-09-02_MAP.md).

That workflow writes to stable roots under `core/phenocycler/data/{cells,
geometry_qc,redsea,expression}/`. It does not currently produce verified
marker calls or cell identities. **The content-addressed run framework that
formerly wrote these outputs under
`core/phenocycler/data/runs/<run_id>/` was removed on 2026-09-16**; there is
no run ID, stage manifest, or fingerprint check any more (see
[core/phenocycler/README.md](core/phenocycler/README.md#outputs)). Runs made
before 2026-09-16 remain, untouched, in that historical `data/runs/<run_id>/`
tree.

**Cell typing is a development lane, and it is not part of that workflow.** The subpanel
classification lives in `core/phenocycler/scripts/typing/` and writes **outside** the package, to
the repo-root `data/development/typing/`: `typing/cells/` was rewritten on 2026-09-19 with the
islet chain applied to all 20 donors — 25,230,109 cells, eight classes plus `Excluded` (measured
2026-09-24). Nothing there is a release, the
`NO_VERIFIED_TYPING` latch is unmoved, and `typing/subpanel/subpanel_calls.parquet` is a SUPERSEDED
earlier generation. How it works is `core/phenocycler/scripts/typing/RECONSTRUCTION.md`; read it
rather than re-deriving it.

The application aggregation below is a **separate, legacy lane**. It consumes a
legacy single-cell H5AD and is not a supported bridge from the current
PhenoCycler workflow's output; the two are not connected today. It does not define
marker thresholds or cell identities either.

```
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad  (2.6M cells, 31 markers)
  │
  ├── scripts/reaggregate_islets.py          → islet-level aggregation + trajectory + Leiden
  │     → islet_analysis/islets_core_fixed.h5ad      (5,214 islets)
  │     → data/adata_ins_root.h5ad                    (+ pseudotime + UMAP)
  │     → islet_analysis/islets_core_clustered.h5ad   (+ Leiden at 4 resolutions)
  │
  ├── scripts/compute_neighborhood_metrics.py → data/neighborhood_metrics.csv
  ├── scripts/extract_per_islet_cells.py      → data/cells/*.csv
  ├── scripts/extract_per_donor_tissue.py     → data/donors/*.csv
  │
  └── scripts/build_h5ad_for_app.py           → data/islet_explorer.h5ad (all app data)
```

### Running Pipeline Scripts

```bash
conda activate scvi-env

# Full rebuild from single-cell H5AD (~15 min)
python scripts/reaggregate_islets.py

# Compute peri-islet neighborhood metrics (~3 min)
python scripts/compute_neighborhood_metrics.py

# Extract per-islet single-cell CSVs (~5 min)
python scripts/extract_per_islet_cells.py

# Extract per-donor tissue CSVs for Spatial tab (~2 min)
python scripts/extract_per_donor_tissue.py

# Build enriched app H5AD (merges all sources)
python scripts/build_h5ad_for_app.py

# Build Excel fallback (from groovy TSV exports)
python scripts/build_master_excel.py
```

## Project Structure

```
Islet-Explorer-Senior/
├── apps/R/ieCore/              # Shared domain, theme, stats, and join package
├── apps/islet/                 # Modular R Shiny application
│   ├── app.R                      # Entrypoint: UI layout, CSS, server wiring
│   ├── config/islet_default.yml   # Downstream display/domain contract
│   ├── R/
│   │   ├── 00_00_iecore.R         # Loads the one shared helper implementation
│   │   ├── 00_globals.R           # Package loads and app runtime paths
│   │   ├── data_loading.R         # load_master_auto(), prep_data(), H5AD loading
│   │   ├── segmentation_helpers.R # GeoJSON cache, segmentation plot builders
│   │   ├── drilldown_helpers.R    # Single-cell overlay, phenotype colors
│   │   ├── spatial_helpers.R      # Per-donor tissue CSV loader with caching
│   │   ├── viewer_helpers.R       # Avivator URL builders, channel config
│   │   ├── ai_helpers.R           # OpenAI/UF Navigator API, streaming
│   │   ├── mod_plot_*.R           # Plot tab module (UI + server)
│   │   ├── mod_trajectory_*.R     # Trajectory tab module
│   │   ├── mod_statistics_*.R     # Statistics tab module
│   │   ├── mod_spatial_*.R        # Spatial tab module
│   │   ├── mod_viewer_*.R         # Viewer tab module
│   │   └── mod_ai_assistant_*.R   # AI chat panel module
│   └── www/                       # Static assets (logo, avivator/)
├── apps/lymph/                 # Lymph config + synthetic schema fixture; no app fork
├── data/app_data/                 # Legacy CODEX domain root (gitignored)
├── data/development/              # Mutable review artifacts, incl. typing/ and typing/app (gitignored)
├── scripts/                       # App-data and viewer-asset scripts; see scripts/CLAUDE.md
│   └── senior/                    # The nine layer builders + build_layers.py, which owns their order
├── islet_analysis/                # Generated islet-level H5ADs (gitignored)
├── core/phenocycler/              # Subtree: QuPath→expression workflow, development lanes, integration
├── core/panels/                   # Subtree: XeniumPanelExplorer — the one panel-taxonomy authority
├── analysis/xenium/               # Xenium 5K lane (absorbed 2026-08-28)
├── analysis/follicle/             # Spleen follicle lane (absorbed 2026-09-03)
├── deploy/                        # nginx, shiny-server, TissUUmaps container recipes
├── docs/                          # Workflow map, handoffs, dated measurement records
└── CLAUDE.md                      # Developer reference (architecture, conventions)
```

## Development

### Running Locally

```bash
cd apps/islet && ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  Rscript -e 'shiny::runApp(".", port = 7777)'
```

The app must be launched from `apps/islet` — the path resolution assumes exactly two directories
deep. Add `NOT_CRAN=true ... Rscript tests/tab_smoke.R /tmp/tabshots` for the shinytest2 tab smoke
driver (it lives in `tests/`, not `tests/testthat/`).

### Deployment

- **nginx** (port 8080) reverse-proxies `/islet-explorer/` → **shiny-server** (port 3838)
- Symlink: `/srv/shiny-server/islet-explorer` → `apps/islet/`
- Code changes take effect when shiny-server spawns a fresh R worker (no restart needed)

### Testing & Diagnostics

```bash
R CMD INSTALL apps/R/ieCore
Rscript apps/islet/tests/parse_check.R
Rscript apps/islet/tests/run_tests.R
Rscript scripts/diagnose_islets.R            # Islet data quality
Rscript scripts/diagnose_marker_regions.R    # Marker region checks
Rscript scripts/na_audit.R                   # Missing value audit
./scripts/clear_shiny_cache.sh               # Clear Shiny cache
```

### Environment Variables

| Variable | Purpose |
|----------|---------|
| `LOCAL_IMAGE_ROOT` | Directory containing OME-TIFF files for Viewer tab |
| `KEY` | UF Navigator API key (for AI assistant) |
| `BASE` | API base URL (defaults to `https://api.openai.com/v1`) |
| `RETICULATE_PYTHON` | Python binary for reticulate/anndata |
| `DEBUG_CREDENTIALS` | Set to `1` for verbose credential loading |
| `VIEWER_DEBUG` | Set to `1` for Viewer tab debug output |

### Working with Images

```bash
./scripts/convert_to_ometiff.sh input.tiff output.ome.tiff
./scripts/convert_to_omezarr.sh input.tiff output.zarr
python scripts/annotate_ome_tiff_channels.py input.ome.tiff channel_names.txt
python scripts/auto_color_ome_tiff_channels.py input.ome.tiff
```

## Contributing

### Conventions

- **Data keys**: Preserve (`Case ID`, `Donor Status`, `islet_key`) three-key system
- **Region labels**: Only `islet_core`, `islet_band`, `islet_union`
- **Factor ordering**: ND < Aab+ < T1D
- **Normalization**: `none`, `global z-score`, `robust per-donor` (median/MAD×1.4826)
- **Statistical models**: Global `value ~ donor_status + islet_diam_um`; pairwise on residuals; BH correction
- **Large scatter plots (>50K points)**: Use `ggplot2::renderPlot()` not plotly
- **Coordinate convention**: `coord_fixed() + scale_y_reverse()` for tissue scatter (microscopy convention)

### Architecture Notes

- The Plot sidebar is shared with Statistics via `conditionalPanel`
- Root-level outputs (`islet_segmentation_view`, `islet_drilldown_view`, etc.) are defined in `app.R` and shared by Plot + Trajectory modules with active-tab guards to prevent duplicate DOM IDs
- Donor palette is synced across 3 non-namespaced `selectInput`s via `observeEvent` + `updateSelectInput`

See [CLAUDE.md](CLAUDE.md) for full architectural documentation.

## Troubleshooting

**"Package not found" errors:**
```bash
Rscript scripts/install_shiny_deps.R
```

**"Cannot find islet_explorer.h5ad" / "Cannot find master_results.xlsx":**
- Run the data pipeline scripts (see [Data Pipeline](#data-pipeline))
- Or ensure they exist under the active domain root — `data/app_data/` by default

**Trajectory/demographics not available:**
- Install `anndata`: `BiocManager::install('anndata')`
- Ensure `islet_explorer.h5ad` exists under the active domain root (the Excel fallback lacks these features)

**Image viewer not loading:**
```bash
./scripts/install_avivator.sh
./scripts/verify_static_images.sh
```

**App performance issues:**
```bash
./scripts/clear_shiny_cache.sh
```

## License

This project is licensed under the [MIT License](LICENSE). Copyright (c) 2026 smith6jt-cop.

## Citation

If you use Islet Explorer in your research, please cite this repository:

```
Islet Explorer: Interactive Shiny application for exploring human pancreatic islet
CODEX multiplexed imaging data
https://github.com/smith6jt-cop/Islet-Explorer
```

## Acknowledgments

- Developed for analysis of human pancreatic islet CODEX multiplexed imaging data from nPOD donors
- Built with R Shiny, ggplot2, plotly, and sf
- Avivator/Viv viewer for OME-TIFF visualization
- AI assistant powered by UF Navigator API
- Statistical analysis with R stats and broom

---

**Version**: 2.0 (modular architecture)
**Last Updated**: 2026-09-18
