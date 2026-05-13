# CLAUDE.md - Shiny App (Islet Explorer)

Deep details for the modular Shiny app under `app/shiny_app/`. The top-level `CLAUDE.md` covers the entrypoint architecture, root-level outputs, shared reactive state, and data loading. This file documents per-tab features and conventions.

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

## Single-Cell Drill-Down (Phase 8, Feb 2026)

### Per-Islet Cell Data

`extract_per_islet_cells.py` reads single-cell H5AD -> outputs `data/cells/{imageid}_Islet_{N}.csv`:
- 5,214 files, ~203 MB total
- 37 columns: `X_centroid`, `Y_centroid`, `phenotype`, `cell_region` (core/peri), `Cell Area`, `Nucleus Area`, + 31 protein markers
- File naming matches `combined_islet_id` from `islet_spatial_lookup.csv`

### Islet Viewer Panel (Segmentation + Drill-Down)

When an islet is clicked, the viewer panel shows with controls on the left sidebar:
- **Title**: `Islet_N (case_id, status, age, gender)` — donor info from `prepared()$comp`
- **Left sidebar (col-2)**: Layer toggles (Single Cells, Peri Boundary, Structures), Color-by dropdown, Phenotype Palette, Show peri-islet cells.
- **Center plot (col-6)**: Segmentation map with GeoJSON boundaries + optional single-cell overlay
- **Right sidebar (col-4)**: Cell composition bar chart + cell count table (when Single Cells enabled)
- `selected_islet()` carries `case_id`, `islet_key`, `centroid_x/y`, `donor_status`, `donor_age`, `donor_gender`

### Coordinate Alignment

Cell centroids (`X_centroid`, `Y_centroid`) are in micrometers. GeoJSON polygons are in pixels. Conversion: `x_px = X_centroid / PIXEL_SIZE_UM`. Matches `render_islet_segmentation_plot()` coordinate system.

### drilldown_helpers.R

- `PHENOTYPE_COLORS` -- 21 phenotypes with distinctive hex palette (endocrine=warm, immune=cool, structural=neutral)
- `drilldown_available()` -- checks if `data/cells/` exists with CSVs
- `load_islet_cells(imageid, islet_key)` -- CSV loader with `new.env()` caching
- `render_islet_drilldown_plot(info, cells, color_by, show_peri)` -- cells over `build_segmentation_base_plot()`. Uses `fill` aesthetic (shape 21) for phenotype/marker coloring, `alpha` for core/peri (0.9/0.4). `colour` aesthetic reserved for structure boundaries.
- `render_drilldown_summary(cells)` -- horizontal bar chart of phenotype composition

### segmentation_helpers.R Refactor

`build_segmentation_base_plot(info)` extracted from `render_islet_segmentation_plot()`. Returns ggplot with GeoJSON polygon layers + coord_sf (WITHOUT crosshairs or title). Reused by both `render_islet_segmentation_plot()` and `render_islet_drilldown_plot()`.

**Structure Legend**: combines all visible polygons into one sf with a `structure_type` column, uses `aes(colour = structure_type)` + `scale_colour_manual(name = "Structures")` with `key_glyph = "path"`. Display names: Islet, Peri Boundary, Nerve, Capillary, Lymphatic. The `colour` aesthetic is reserved for structures; drilldown cell scatter uses `fill` (shape 21) to avoid scale conflict.

## Statistics Tab (Phase 6+14+16, Feb-Mar 2026)

### Pseudoreplication Fix (Phase 16)

The app has 5,214 islets from only **15 donors** (5 per group). Islets within a donor are correlated. Phase 16 fixes this textbook pseudoreplication by:

1. **Donor-level aggregation**: All tests use `aggregate_to_donor()` to compute per-donor means (N=15), not islet-level (N=5,214). Produces realistic p-values (0.01-0.05 range, not 1e-50).
2. **Mixed-effects sensitivity**: `lmer_test_donor()` runs `lmerTest::lmer(value ~ donor_status + (1 | case_id))` on islet-level data with donor as random intercept. ICC quantifies donor-level clustering.
3. **Min-cells filter**: Available on all 4 analysis tabs — Plot, Trajectory, Statistics, Spatial. All use `total_cells_core + total_cells_peri >= threshold` (Trajectory uses `total_cells` directly). Defaults to 1.
4. **Normality testing**: Shapiro-Wilk on donor-level means per group, with auto-suggestion for test type.
5. **No-binning option**: `stats_no_binning` checkbox hides Section 3 (size-dependent analysis).
6. **Bin slider fix**: Max increased from 75 to 150, step=5. Heatmap x-axis labels angled (-45°).
7. **Per-bin donor-level**: `per_bin_donor_anova()` and `per_bin_donor_kendall()` aggregate to donor means within each bin. Bins with <2 donors per group are skipped.
8. **Demographics redesign**: Age scatter shows all islets (not donor means) with Pearson r. Sex analysis uses box plot faceted by sex. Autoantibody analysis (Aab+ only): per-donor AAb profile table + feature by AAb count box plot. Use 'Sex' not 'Gender' in user-facing labels.

### Key Utilities (`utils_stats.R`)
- `aggregate_to_donor(rdf, agg_fn)` -- Groups by `Case ID` + `donor_status`, computes mean/median of `value`, preserves age/gender
- `lmer_test_donor(rdf)` -- Mixed-effects model, returns `list(p_value, icc, fit, error_msg)`. Uses `lmerTest::anova()` (Satterthwaite F-test). Returns diagnostic `error_msg` on failure (singular fit, convergence, missing packages).
- `per_bin_donor_anova(df, ...)` -- Donor-level per-bin ANOVA. Skips bins where any group has <2 donors.
- `per_bin_donor_kendall(df, ...)` -- Donor-level per-bin Kendall tau. Skips bins with <3 total donors.
- `normality_tests(ddf)` -- Shapiro-Wilk per group, returns data.frame with `W`, `p_value`, `is_normal`
- `cohens_d(x, y)` -- returns `list(d, ci_lo, ci_hi)` (NOT `ci_lower`/`ci_upper`)
- `eta_squared(fit)` -- from `anova()` output
- `pairwise_wilcox(df, group_col, value_col)` -- returns data.frame with columns `group1`, `group2`, `p_value` (NOT `p.adj`/`statistic`)
- `flag_outliers(df, value_col, group_col, method, threshold, transform, exclude_zeros, enabled)` -- Single source of truth for outlier flagging. Appends `is_outlier` (logical) and `outlier_score` (z or fence-distance) columns. Attributes carry method/threshold/enabled for downstream summaries.
- `summarize_outlier_filter(df)` -- Compact summary list for badges: `n_total`, `n_excluded`, `n_kept`, `method`, `threshold`, `enabled`, `label`.
- `outlier_title_suffix(summary)` -- Formats the `" — N outliers excluded (label)"` string for plot titles. Returns `""` when no exclusion.
- `GLOBAL_OUTLIER_TOOLTIP` -- String constant used as `title=` attribute on every outlier control across tabs.

### Unified Outlier Handling (May 2026)
- One non-namespaced checkbox per tab (`sidebar_remove_outliers`, `traj_remove_outliers`, `spatial_remove_outliers`) — synced bidirectionally via three `observeEvent`s in `app.R`. The synced reactive `remove_outliers_global()` is injected into every module server (`plot_server`, `trajectory_server`, `statistics_server`, `spatial_server`) as the `remove_outliers` argument.
- Default: ON. Toggling propagates instantly to all four tabs.
- Plot, Trajectory, Statistics, Spatial Card C: `flag_outliers(method = "zscore", threshold = 3)` per donor status. Card C transforms via `log1p` to handle the long-tailed distance distribution.
- Spatial Card A: `flag_outliers(method = "iqr", threshold = 1.5, exclude_zeros = TRUE)` per donor status. Keeps the previous Tukey logic because phenotype enrichment is zero-inflated.
- Visual surfacing: every affected plot shows a **count badge** in its title or subtitle and renders excluded points as **shape=4 grey × ghost markers** (Plot, Trajectory, Spatial Card C). Stats forest plot title shows the badge. Statistics Methods Reference paragraph names the rule and cites the exact `n_excluded` for the current Run. Download CSV header records `Outlier rule: …` and the count.
- Spatial Card B's "Clip z > 5" checkbox is a **display saturation only** (caps the colour scale) — left alone, NOT an outlier control. Documented as such in `docs/user_guide.md`.
- Trajectory bug fix: previously discarded outliers silently before plotting and logged only to console; now drawn as ghosts and counted in a badge above the scatter.

### R Package Dependencies for Phase 16
- `lmerTest` (3.2-1) -- Mixed-effects model with Satterthwaite df
- `lme4` (1.1.37) -- `VarCorr()` for ICC computation

### Shared Sidebar Architecture
The Plot sidebar (mode, feature, region, donor status, AAb, age, gender filters) is visible on both the Plot and Statistics tabs. Achieved via `conditionalPanel` condition `"input.tabs == 'Plot' || input.tabs == 'Statistics'"` and matching JS `adjustLayout()` logic in `app.R`. The Statistics module consumes `plot_returns$raw_df` and `plot_returns$summary_df` directly -- no data duplication.

### 5-Section Narrative Layout

1. **Configure Analysis** -- Overview banner, Run/CSV buttons, test type, alpha, outliers, min-cells filter, normality test results, no-binning checkbox, bin width, diameter range.
2. **Primary Results** -- Pseudoreplication info banner (ICC, mixed-effects p) + hypothesis table (col-4) + forest plot (col-4) + AUC (col-4). Cards equal-height via `display: flex` row.
3. **Size-Dependent Patterns** -- Wrapped in `conditionalPanel` for no-binning toggle. Donor-level stratified tests heatmap + trend analysis.
4. **Confounders & Deeper Analysis** -- Demographics: donor summary table, age scatter (islet-level, all points), sex box plot (faceted by sex), autoantibody profile (Aab+ only), covariate model (donor-level)
5. **Methods Reference** -- Pseudoreplication, donor-level aggregation, ICC, mixed-effects, normality testing

### Multiple Testing Correction
Per-bin ANOVA and Kendall τ p-values are BH-corrected across bins via `p.adjust(method = "BH")`. Heatmap shows `-log10(q)` (corrected) with star annotations using corrected values.

### Plotly Legend Positioning
Forest plot uses `theme(legend.position = "none")` in ggplot (prevents ggplotly auto-placement), then `layout(margin = list(b = 60), legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))` in plotly. The `margin(b=60)` creates room below the x-axis title to prevent overlap. The trend plot uses `margin(b=100)` with legend at `y = -0.3`.

## Trajectory Tab Features

### Cell-Count-Weighted Trajectory (Phase 12, Mar 2026)

#### Problem
60.8% of islets have ≤10 cells — their aggregated scVI embeddings are noisy. DPT operates on 10-dim scVI latent space but the scatter shows single-marker expression from `.X`. An islet can be "far" in latent space (unusual across ALL 31 markers → high DPT) while having normal INS expression. Real biology, but shouldn't get equal visual weight.

#### Solution: Cell-Count-Aware Visualization
- **Point sizing**: `sqrt(total_cells)` — area proportional to cell count. Default in trajectory scatter.
- **Weighted LOESS**: `weight = log1p(total_cells)` — well-measured islets drive trends. Raw counts too skewed (median=9, max=1,902) → `log1p` compresses 1,902:1 ratio to ~11:1.
- **Hover tooltip**: Shows "Cells: N".
- **UI controls**: Point size selector: Cell Count (default) / Islet Diameter / Uniform.

#### Key Details
- `total_cells` is in `adata_ins_root.h5ad` `.obs` (all 5,214 islets, min=1, max=1,902)
- Added to `traj_data_clean()` result_df with NA guard (fills missing with 1)
- `log1p(total_cells)` for LOESS weights — prevents extreme islets from distorting trend lines
- All three trend line variants weighted: overall, by-donor (donor_id coloring), by-donor (donor_status coloring)
- **Critical**: Raw `total_cells` as LOESS weight causes wild curves — a 1,902-cell islet dominates the fit. Always use `log1p()`.

### Pseudotime Mode Selector (Phase 17, Apr 2026)

`mod_trajectory_ui.R` has a "Pseudotime mode" dropdown with two options: "Core only" (default, `dpt_pseudotime`) and "Core + Peri" (`dpt_pseudotime_combined`). `mod_trajectory_server.R::traj_with_pseudotime()` switches `obs$pseudotime` based on `input$pt_mode`. The selector controls X-axis (pseudotime) only — Y-axis features remain the 31-marker panel. See `scripts/CLAUDE.md` for the upstream pipeline that produces these two pseudotimes.

### Pseudotime Heatmap: Dominant Status × Feature Intensity (May 2026)

The strip heatmap below the trajectory scatter uses a 2-D encoding per pseudotime bin:

- **Hue** = dominant donor_status, computed as the argmax of per-status NORMALIZED density (`count_in_bin / total_for_status`). Normalization corrects for cohort imbalance (ND ~2300, Aab+ ~1668, T1D ~1246 islets) so dominance reflects enrichment relative to each status's pool, not raw count.
- **Intensity** = bin-mean of selected feature, normalized via 5–95th percentile clip, blended in `[INTENSITY_FLOOR, 1]` toward the donor color so even the dimmest bin is a clear tint (never white). `INTENSITY_FLOOR = 0.50`.
- **Renderer**: ggplot + `geom_tile` + `scale_fill_identity` (NOT plotly heatmap). Guarantees bars span [0, 1]. Plotly previously misaligned bin widths with the scatter axis margins.
- **Hover**: Shows pseudotime, dominant status, raw counts per status (n_ND, n_Aab, n_T1D), normalized density per status (%), and feature mean.

Replaces the previous "average donor_status code (0/1/2)" heatmap, which conflated cohort-level imbalance with per-bin enrichment and produced a single colorscale gradient that hid expression intensity entirely.

### Multi-Feature Trajectory Heatmap
Z-scored expression along pseudotime. Clamped [-2.5, 2.5], blue-white-red, dynamic height.

## Plot Tab Features (Phase 5, Feb 2026)

### Phenotype Composition Explorer
Plot tab Composition mode exposes 21 cell-type proportions (`prop_*` from H5AD `.obs`) + 21 peri-islet proportions (`peri_prop_*`) + 5 immune metrics alongside 3 hormone fractions. Grouped `selectInput` with 4 option groups.

### Age & Sex Demographic Filters
Age slider and sex checkboxes in Plot sidebar (H5AD only). `renderUI` returns `NULL` when absent.

### Segmentation + Drill-Down Pattern
Click a point -> `selected_islet()` updates -> embedded panel renders inline with Boundaries/Single Cells toggle. Close button sets `selected_islet(NULL)`.

## AI Assistant (Phase 10, Feb 2026)

### Architecture

The AI chat panel (`mod_ai_assistant_ui/server.R`) uses the University of Florida Navigator AI Toolkit via an OpenAI-compatible API. Wired as `ai_assistant_server("ai_chat")` in `app.R`.

**Key files:**
- `ai_helpers.R` -- Credential loading, API calls (streaming + non-streaming), model selection, error handling
- `mod_ai_assistant_ui.R` -- Chat panel UI with model picker, textarea, send/reset buttons
- `mod_ai_assistant_server.R` -- Chat history management, streaming callback, token budget tracking, `strip_markdown()` for plain-text display

### UF Navigator API

- **Endpoint**: `https://api.ai.it.ufl.edu/v1/chat/completions`
- **Available models** (as of Feb 2026):
  - `gpt-oss-20b` -- Fast model (default)
  - `gpt-oss-120b` -- Large reasoning model (chain-of-thought in `reasoning_content`, answer in `content`)
- **Authentication**: Bearer token from `KEY` env var (loaded from `.Renviron`)
- **Base URL**: `BASE` env var in `.Renviron` (defaults to `https://api.openai.com/v1` if unset)

### Reasoning Model Behavior (gpt-oss-120b)

The 120b model is a reasoning model (similar to DeepSeek R1):

- **Non-streaming**: Response has `message.content` (final answer) + `message.reasoning_content` (chain-of-thought). Content may be `null` if `max_tokens` is too low for reasoning + answer.
- **Streaming**: Chunks first arrive with `delta.reasoning_content` only (no `delta.content`). After reasoning completes, `delta.content` chunks appear with the actual answer. The streaming parser shows "Thinking..." during the reasoning phase.

### Credential Loading

`ai_helpers.R` searches for `.Renviron` in priority order:
1. `app/shiny_app/.Renviron` (deployment)
2. Script directory `.Renviron`
3. `R_ENVIRON_USER` env var
4. `~/.Renviron`
5. `~/.Renviron.local`
6. `$HOME/.Renviron`

Required env vars: `KEY` (API key), `BASE` (API base URL, optional).

### Model Fallback Logic

`select_openai_models()` returns a candidate list. For each candidate, `call_openai_chat()` tries streaming first, then non-streaming. If a model returns HTTP 400/401/404 with "model" in the error, it falls back to the next candidate. The fallback model defaults to `gpt-oss-120b` (env var `OPENAI_DEFAULT_MODEL`).

**Critical**: UF Navigator returns HTTP **401** (not 404) for model access denied. The fallback logic must include 401 in its status code checks.

### Important AI Conventions

- **Model names must match API**: Query `GET /v1/models` to verify available models.
- **Reasoning models need sufficient tokens**: `max_output_tokens` must be large enough for reasoning + answer.
- **UF Navigator uses `/chat/completions` only**: The code detects `api.ai.it.ufl.edu` in the base URL and skips the `/responses` endpoint.
- **Streaming parser must handle `reasoning_content`**: For reasoning models, `delta.reasoning_content` chunks arrive first — use a "Thinking..." indicator.
- **Token budget**: Controlled by `OPENAI_TOKEN_BUDGET` env var. Chat panel shows cumulative usage and blocks requests when budget exhausted.
- **Debug mode**: Set `DEBUG_CREDENTIALS=1` env var for verbose credential loading logs.
- **Plain-text display**: `strip_markdown()` in `mod_ai_assistant_server.R` removes markdown formatting from AI responses before display.

## App-Specific Conventions

- **CSS overflow for cards with dropdowns**: `selectInput` menus extend below their container. Add `overflow: visible;` to card styles or dropdowns get clipped behind cards.
- **Font size minimums**: Use `h5` (not `h6`) with explicit `font-size: 15px` for panel headings. Legend items minimum 14-15px.
- **Large scatter plots (>50K points)**: Use `ggplot2::renderPlot()` NOT `plotlyOutput()` -- plotly cannot handle 100K+ points interactively. Use small `size` (0.3-0.5) and low `alpha` (0.3-0.6).
- **Tissue scatter coordinate convention**: `coord_fixed() + scale_y_reverse()` matches microscopy convention (y increases downward, spatial proportions preserved)
- **Leiden cluster mapping for cells**: Islet-level Leiden assignments map to single cells via `islet_name` column lookup. Tissue background cells without an islet get `cluster = "tissue"` (grey).
- **Per-donor tissue CSVs**: `data/donors/{imageid}.csv` with 5 columns (X_centroid, Y_centroid, phenotype, cell_region, islet_name). `cell_region` = core/peri/tissue.
- **Leiden in H5AD**: 4 resolution columns (`leiden_0.3/0.5/0.8/1.0`) + 2 UMAP coords (`leiden_umap_1/2`). Extracted via `^leiden_` regex in `data_loading.R`.
- **Leiden UMAP uses visualization coords**: `islets_core_clustered.h5ad` X_umap is copied from trajectory's raw marker PCA UMAP, NOT computed separately from scVI (which produces a blob). Leiden clustering itself still uses scVI neighbors.
- **Donor palette syncing**: 3 non-namespaced `selectInput`s (`sidebar_donor_palette`, `traj_donor_palette`, `spatial_donor_palette`) synced bidirectionally via `observeEvent` + `updateSelectInput` in `app.R`. All feed `donor_palette_name()` -> `active_donor_colors()`.
- **Trajectory legend uses reactive colors**: `donor_colors_reactive()` NOT hardcoded hex. Changes when palette selector updates.
- **Spatial background cell coloring**: `color_background` checkbox. When checked + phenotype mode, background cells get phenotype colors at alpha=0.25.
- **Reticulate bridge crossing**: `ad$uns[[key]]` costs ~221ms per call. ALWAYS cache `uns <- ad$uns` as an R list first, then index. `reconstruct_groovy_df_from_list()` uses cached list; legacy `reconstruct_groovy_df()` is a thin wrapper.
- **Deferred heavy file loading**: `annotations.tsv` (72 MB) lazy-loaded via `.seg_lazy` environment in `segmentation_helpers.R`.
- **Synthetic union rows need density**: `prep_data()` synthesizes `islet_union` rows by summing core+band. Must compute `area_um2`, `region_um2`, and `area_density = area_um2/region_um2` — setting these to `NA` causes Microenvironment Core+Peri to drop 12/15 donors (default metric is Density, NA filtered by `!is.na(value)`).
- **`updateSelectInput` after `renderUI`**: Server-side updates targeting an input that lives inside a `renderUI`-mounted block must share a reactive dependency with the renderUI, or fire only after the UI has mounted on the client. Putting the update inside its own naked `observe({})` can race the UI mount and silently leave the dropdown empty. Pattern: nest the `updateSelectInput` inside the same observer that builds the choices reactive.
