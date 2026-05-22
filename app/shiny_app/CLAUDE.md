# CLAUDE.md - Shiny App (Islet Explorer)

Deep details for the modular Shiny app under `app/shiny_app/`. The top-level `CLAUDE.md` covers the entrypoint architecture, root-level outputs, shared reactive state, and data loading. This file documents per-tab features and conventions.

## Spatial Neighborhood Analysis

See [`docs/spatial.md`](docs/spatial.md) for full details — peri-islet metrics, tab layout (Phase 9+13+15), neighborhood-card server architecture, KDE selector race fix, tissue scatter design.

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

## Statistics Tab

See [`docs/statistics.md`](docs/statistics.md) for full details — pseudoreplication fix (Phase 16), `utils_stats.R` reference, unified outlier handling, zero-value exclusion, sex-stratified BH correction, 5-section narrative layout, multiple-testing correction.

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
Plot tab Cell Populations mode exposes 21 cell-type proportions (`prop_*` from H5AD `.obs`) + 21 peri-islet proportions (`peri_prop_*`) + 21 core+peri combined proportions (`combined_prop_*`, computed on the fly in `raw_df_base()`) alongside 3 hormone fractions. Grouped `selectInput` with up to 4 option groups (Hormone Positivity, Core, Peri-Islet, Core+Peri). The Combined group's choices are the intersection of `prop_*` and `peri_prop_*` stems — one entry per phenotype that has both region values.

### Age & Sex Demographic Filters
Age slider and sex checkboxes in Plot sidebar (H5AD only). `renderUI` returns `NULL` when absent.

### Segmentation + Drill-Down Pattern
Click a point -> `selected_islet()` updates -> embedded panel renders inline with Boundaries/Single Cells toggle. Close button sets `selected_islet(NULL)`.

### Metric Label Vocabulary (May 2026)
Sidebar radio labels and axis titles are unit-explicit — never use the bare word "density" without naming the denominator. Canonical phrasing:
- **Cell Populations Counts** → `"Cells per islet (core)"` / `"Cells per islet (peri)"` / `"Cells per islet (core+peri)"`.
- **Cell Populations Density** → selector `"Cells per mm² of core area"` / `"… peri area"` / `"… core+peri area"`; y-axis `"Cells per mm² (core area)"` / `"(peri area)"` / `"(core+peri area)"`. Denominator is the QuPath-measured annotation area (`core_region_um2`, `peri_region_um2`, or their sum), joined into `comp` by `prep_data()` from `targets_all` rows. Peri-islet area comes from the `_exp20um` annotation in `annotations.tsv`.
- **Microenvironment** → metric radio hidden (only one option). Y-axis names both the class and the region: `"Nerve area (% of core area)"` / `"Capillary area (% of peri-islet area)"` / `"Lymphatic area (% of core+peri area)"`. **Value computed from `100 * area_um2 / region_um2`** (% of the selected region's QuPath annotation area occupied by the structure). Counts are not exposed — these targets (Nerve / Capillary / Lymphatic) are area-valued (a single nerve may segment into several pieces), so object counts depend on segmentation morphology, not biology.
- **Core+Peri combined proportions** → computed on the fly: `value = (prop_X · total_cells_core + peri_prop_X · total_cells_peri)` divided by either `(total_cells_core + total_cells_peri)` (Percentage) or `(core_region_um2 + peri_region_um2)` (Density). NA totals/areas are treated as 0 so islets with only one zone still contribute.
- **Spatial Card C distance KDE** y-axis is `"Cells per µm of distance (KDE × N)"`, not a per-area density.

A `helpText` block under the Plot metric selector restates the active mode's denominator. Keep this phrasing if you add new metrics. The y-axis class name is parsed from `input$which` (encoded as `className|region`).

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
- **Donor palette syncing**: 3 non-namespaced `selectInput`s (`sidebar_donor_palette`, `traj_donor_palette`, `spatial_donor_palette`) synced bidirectionally via `observeEvent` + `updateSelectInput` in `app.R`. All feed `donor_palette_name()` -> `active_donor_colors()`. Default selection: **Bright** (`DONOR_COLORS_BRIGHT` = `#3366CC` / `#FF6600` / `#109618`); the legacy "Paul Tol" palette remains available as a choice. The fallback `get_donor_color_palette()` default arg in `00_globals.R` is also `"Bright"` so headless/test contexts match.
- **Trajectory legend uses reactive colors**: `donor_colors_reactive()` NOT hardcoded hex. Changes when palette selector updates.
- **Spatial background cell coloring**: `color_background` checkbox. When checked + phenotype mode, background cells get phenotype colors at alpha=0.25.
- **Reticulate bridge crossing**: `ad$uns[[key]]` costs ~221ms per call. ALWAYS cache `uns <- ad$uns` as an R list first, then index. `reconstruct_groovy_df_from_list()` uses cached list; legacy `reconstruct_groovy_df()` is a thin wrapper.
- **Deferred heavy file loading**: `annotations.tsv` (72 MB) lazy-loaded via `.seg_lazy` environment in `segmentation_helpers.R`.
- **Synthetic union rows need density**: `prep_data()` synthesizes `islet_union` rows by summing core+band. Must compute `area_um2`, `region_um2`, and `area_density = area_um2/region_um2` — setting these to `NA` causes Microenvironment Core+Peri to drop 12/15 donors (default metric is Density, NA filtered by `!is.na(value)`).
- **`updateSelectInput` after `renderUI`**: Server-side updates targeting an input that lives inside a `renderUI`-mounted block must share a reactive dependency with the renderUI, or fire only after the UI has mounted on the client. Putting the update inside its own naked `observe({})` can race the UI mount and silently leave the dropdown empty. Pattern: nest the `updateSelectInput` inside the same observer that builds the choices reactive.

## Plot Fonts & Publication-Quality Downloads (May 2026)

Single source of truth in `R/00_globals.R`:

- `ISLET_FONT_SIZES` — named list with `base=16`, `title=18`, `subtitle=14`, `axis_title=16`, `axis_text=14`, `legend_title=16`, `legend_text=14`, `strip=14`, `dense_min=12` (heatmap rows).
- `theme_islet(base_size = 16)` — ggplot theme extending `theme_minimal()`. Use on every `renderPlot` / inline `ggplot()` that is also handed to `ggplotly()`.
- `plotly_islet_fonts()` — returns a named list of `list(size = N)` font specs (keys: `global`, `title`, `axis_title`, `axis_tick`, `legend_title`, `legend_text`, `colorbar_title`). Splice into `layout(font=..., title=list(font=...), xaxis=list(title=list(font=...), tickfont=...), yaxis=..., legend=list(font=..., title=list(font=...)))`. `ggplotly()` does **not** fully respect ggplot themes — the explicit `layout(...)` overlay is required after every `ggplotly()` call.
- `plotly_islet_config(p, filename_root, width=1400, height=900, scale=2)` — pipe with `%>%` after `layout()`. Sets `toImageButtonOptions` so the modebar camera icon exports a 2800 × 1800 px PNG and disables the plotly logo. Apply to **every** plotly output except deliberate UI strips like the trajectory pseudotime heatmap (`mod_trajectory_server.R:1284` uses `config(displayModeBar = FALSE)` and is left as-is).
- `ISLET_PNG_DIMS` (10 × 7 in @ 300 dpi → 3000 × 2100 px) and `ISLET_HEATMAP_PNG_DIMS` (12 × 8 in) for `ggsave()` in `downloadHandler()` `content` functions.
- `islet_png_filename(stem)` — `<stem>_YYYYMMDD_HHMMSS.png`.

**Pattern for ggplot downloads** — when a `renderPlot` output should be downloadable, refactor it into a `reactive()` and a thin `renderPlot({ plot_reactive() })` so both the renderer and the `downloadHandler` consume the same ggplot object. Then `ggsave(file, plot = plot_reactive(), width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height, dpi = 300, units = "in")`. Examples: `mod_trajectory_server.R` UMAP donor/feature + multi-feature heatmap (`traj_umap_donor_plot()`, `traj_umap_feature_plot()`, `traj_multi_heatmap_plot()`); `mod_spatial_server.R` tissue scatter and donor UMAP (`spatial_tissue_scatter_plot()`, `spatial_umap_donor_plot()`); `app.R` segmentation + drilldown (`islet_segmentation_view_plot()`, `islet_drilldown_summary_plot()`).

**Download button placement** — non-namespaced `downloadButton("dl_islet_segmentation_view", ...)` for root-level outputs (segmentation, drilldown) lives inside each module's `segmentation_viewer_panel` renderUI alongside the plotOutput; namespaced `downloadButton(ns("dl_..."), ...)` for module-local plots lives in the module UI immediately under the plotOutput. Render handlers for the non-namespaced download buttons sit at root in `app.R`; namespaced ones live in the module server.

**Don't** mix the helper with manual `font = list(size = X)` per-key calls — the helper is the entire convention. When you need a smaller font for dense rows (e.g., heatmap text annotations) reach for `ISLET_FONT_SIZES$dense_min`, not a fresh literal.
