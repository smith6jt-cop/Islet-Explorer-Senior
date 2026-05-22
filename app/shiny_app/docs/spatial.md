# Spatial Neighborhood Analysis — Detailed Reference

Detailed reference for the Spatial tab. Linked from `app/shiny_app/CLAUDE.md`.

## Spatial Neighborhood Analysis (Phase 7, Feb 2026 — extended Apr 2026)

### Peri-Islet Metrics (~100 columns, all phenotypes systematic)

| Category | Columns | Description |
|----------|---------|-------------|
| Peri-islet composition | `peri_prop_*`, `peri_count_*`, `total_cells_peri` | Proportion & count for each of 21 phenotypes in 20µm expansion zone |
| Aggregate immune metrics (legacy) | `immune_frac_peri/core`, `immune_ratio`, `cd8_to_macro_ratio`, `tcell_density_peri` | Immune fractions, ratios, density. Columns retained for compatibility but no longer surfaced in any UI. |
| Enrichment z-scores (per phenotype) | `enrich_z_<phenotype>` | Poisson z (observed vs tissue-baseline) for **all 20 phenotypes** (excludes "Unknown"). Extended Apr 2026 from a hand-curated 7-immune subset. |
| Distance metrics (per phenotype) | `min_dist_<phenotype>`, `min_dist_immune_mean` | Min distance from islet core centroid to nearest peri cell of each phenotype, plus the aggregate "any immune" mean across the 7 IMMUNE_TYPES. |

### Data Flow

`compute_neighborhood_metrics.py` reads single-cell H5AD -> `neighborhood_metrics.csv` (5,214 rows) -> merged into `islet_explorer.h5ad` `.obs` by `build_h5ad_for_app.py` (step 4.5) -> extracted in `data_loading.R` `load_master_h5ad()` -> merged into `comp` by `prep_data()` -> available in Plot composition selector (4 option groups) + Statistics tab + Spatial tab.

### Key Details

- All 5,214 islets have peri-islet data (100% coverage with min_cells=0 + require_paired=True)
- Column naming sanitizes phenotype names: spaces -> `_`, `+` -> `plus` (e.g., `peri_prop_ECADplus`)
- Immune signal validated: T1D immune_frac_peri (0.155) > Aab+ (0.106) > ND (0.069)
- **Plot composition selector** has 3 option groups: Hormone Positivity, Phenotype Proportions, Peri-Islet Proportions. The 5 aggregate ratio metrics were removed (Apr 2026) in favour of phenotype-driven Spatial UI selectors.
- **`IMMUNE_TYPES`** in `compute_neighborhood_metrics.py` (CD8a Tcell, CD4 Tcell, T cell, B cell, Macrophage, APCs, Immune) is intentionally narrow — drives only aggregates where "what counts as immune" needs definition. **`per_type_phenotypes`** (all `obs.phenotype.unique()` minus `PER_TYPE_EXCLUDE = {"Unknown"}`) drives `enrich_z_*` and `min_dist_*` per-phenotype loops.

### Spatial Tab Layout (Phase 9+13+15, mod_spatial_ui/server.R)

Row 1: Controls sidebar (col-2) + Tissue Scatter (col-6) + Leiden Panel (col-4). Row 2-4: Neighborhood Analysis Cards (A/B/C).

1. **Controls Sidebar** (col-2) -- donor selector, color-by (phenotype/Leiden), Leiden resolution, region filter (All/Core+Peri/Core), "Color background cells" checkbox, phenotype filter (checkboxes with All/None links), donor status checkboxes, phenotype + donor palette selectors, download button.
2. **Tissue Scatter** (col-6) -- ggplot2 `renderPlot` (NOT plotly) showing ~177K cells/donor. Background tissue cells in light grey or colored by phenotype (toggle). Foreground (core/peri) colored by phenotype or Leiden cluster. `coord_fixed() + scale_y_reverse()`. Height: 800px. Supports brush-to-zoom.
3. **Leiden Panel** (col-4) -- plotly UMAP of 5,214 islets colored by selected Leiden resolution (0.3/0.5/0.8/1.0) + stacked bar chart of mean phenotype composition per cluster. UMAP uses raw marker PCA visualization coords (same as trajectory).
4. **Neighborhood Analysis Cards** (3 sections, conditionally rendered via `has_neighborhood()` guard). Apr 2026: all four card selectors are now phenotype-driven and dynamically populated from `prepared()$comp` columns.
   - **Global controls bar**: Min cells/islet filter, point size, opacity, live islet count.
   - **Card A: Peri-Islet Phenotype Enrichment** -- Bar chart (donor-level mean/median enrichment z-score per status for selected phenotype) + peri vs core proportion scatter. Inputs: `infiltration_phenotype`, `scatter_phenotype`.
   - **Card B: Phenotype Composition & Enrichment** -- Grouped bar chart + heatmap across all 20 phenotypes × 3 stages. Region toggle: peri enrichment z-scores (diverging blue→red), core proportions (sequential white→red), peri proportions (sequential white→red).
   - **Card C: Phenotype Proximity to Islet** -- Distance box plots (per phenotype + "Immune (all)" aggregate) + signed-distance KDE. Inputs: `distance_metric` (built from `min_dist_*` columns) and `kde_immune_type` (built from `immune_distance_kde.csv`).
   - Section headings reuse `section_heading()` pattern (gradient pill badge + title + subtitle).
   - Distance NA values: caused by zero immune cells in peri-islet zone (20.1% NA for any immune, 88.6% for CD8+ T-cells). Biological, not a bug. Min-cells filter at ≥50 reduces NA from 20% to 6%.

Wired as `spatial_server("spatial", prepared, active_palette, active_donor_colors)`.

#### Neighborhood Cards Server Architecture
- `nbr_comp()` shared reactive -- filters `prepared()$comp` by `input$groups` + `input$nbr_min_cells`. Reused by all 6 outputs.
- `phenotype_options()` reactive (Apr 2026) -- builds a list of phenotype records from `prepared()$comp` columns. Each entry: `label`, `safe`, `enrich`, `prop_core`, `peri_prop`, `peri_count`, `min_dist`. Recovers original phenotype label from sanitised column suffix by trying common transforms (`_` → space, `plus` → `+`).
- `phenotype_choices()` reactive -- named vector for `selectInput` choices: names = labels, values = safe-names. All four card selectors populate via `updateSelectInput` from this.
- `enrich_summary()` intermediate reactive -- iterates over `phenotype_options()` per cell type × donor status. Respects `input$enrich_region`, `input$enrich_clip`, `input$enrich_stat`.
- All 6 outputs use `donor_colors_reactive()` for palette syncing.
- **Plotly categorical axis ordering**: Must set `categoryorder = "array", categoryarray = c("ND", "Aab+", "T1D")` on x-axis. Plotly defaults to alphabetical, which puts "Aab+" before "ND".
- **KDE selector race fix (May 2026)**: The `kde_immune_type` `updateSelectInput` MUST live inside the same observer that updates `distance_metric` (which depends on `phenotype_choices()`). When it lived in its own `observe({})` it sometimes fired before the Card C `renderUI` had mounted the selectInput on the client, leaving the dropdown empty. Sharing the reactive dependency guarantees the update fires after the card's UI is in the DOM. See `mod_spatial_server.R` ~line 898.

#### Supporting Files
- `spatial_helpers.R` -- `donor_tissue_available()`, `get_available_donors()`, `load_donor_tissue(imageid)` with `new.env()` caching
- `data/donors/{imageid}.csv` -- 15 per-donor CSVs (5 cols: X_centroid, Y_centroid, phenotype, cell_region, islet_name), ~78 MB total
- `islets_core_clustered.h5ad` -- Leiden clustering source (4 resolutions: leiden_0.3/0.5/0.8/1.0)

#### Tissue Scatter Design
- Uses `ggplot2::renderPlot` NOT plotly -- 177K points would freeze plotly
- Foreground/background layering: tissue cells at `size=0.15, alpha=0.3` in grey; core/peri at `size=0.4, alpha=0.6` in color
- Leiden coloring maps `islet_name -> cluster` via islet-level lookup from `prepared()$comp`
- **Brush-to-zoom**: `plotOutput(brush=brushOpts(...), dblclick=...)` + `reactiveValues(xmin,xmax,ymin,ymax)`. When zoomed, switches from `coord_fixed()` to `coord_cartesian(xlim, ylim)` with `sort()` on ylim for reversed y-axis.
- **Phenotype filter**: `checkboxGroupInput` with All/None action links. Dynamically populated from donor's unique phenotypes.
- Donor 6533 has 191 islets (9,815 core + 13,840 peri cells) — fully integrated after patching Parent annotations.

