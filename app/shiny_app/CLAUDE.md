# CLAUDE.md - Shiny App (Islet Explorer)

Deep details for the modular Shiny app under `app/shiny_app/`. The top-level `CLAUDE.md` covers the entrypoint architecture, root-level outputs, shared reactive state, and data loading. This file documents per-tab features and conventions.

**Per-tab capability status** (what actually works on the current data, verified by launching the dev server + screenshotting every tab, 2026-06-29) lives in [`docs/capability_status.md`](docs/capability_status.md). Headline: all 5 data tabs work; Viewer is empty-but-graceful (no OME-TIFFs); **Trajectory works** but was slow to first render — vectorized 2026-06-29 (see below). All tabs still serve the OLD provisional fine-phenotype labels by design — the new broad-lineage outputs are notebook-only.

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

### macOS Native-Click Bridge + Shared Resolver (Jun 2026)

**Problem**: On macOS, plotly's click/drag discrimination suppresses `plotly_click` entirely, so the server `observeEvent(event_data("plotly_click", ...))` never fires — clicking a point did nothing (no notification). Both Plot and Trajectory drill-down were affected (Trajectory's key→customdata→coord fallbacks all live *inside* the plotly_click observer, so a suppressed event bypassed them too). Works on Linux/Windows; macOS-only — no platform-conditional code, so the cause is the macOS browser event layer.

**Fix** — `islet_click_bridge(p, input_id)` in `00_globals.R`:
- `htmlwidgets::onRender` that listens for native DOM `mousedown`/`mouseup` (which always fire, immune to plotly's drag suppression). On a click-like gesture (movement < 6px) it reads the hovered islet's key and pushes `{key, nonce}` to a Shiny input via `Shiny.setInputValue(..., {priority:'event'})`.
- Key source: the ggplot `key` aes (`"case_id|islet_key"`), read as `pt.data.key[pt.pointNumber]` from the cached `plotly_hover` event AND from plotly's live `gd._hoverdata` at mousedown (belt-and-suspenders).
- Piped after each scatter's `event_register("plotly_click")`: Plot uses `ns("plot_native_click")`, Trajectory `ns("traj_native_click")`. The legacy `plotly_click` handler is kept for platforms where it already works.
- **`gd.__isletBridge` attach-once guard**: `onRender` re-runs on every plot re-render and would stack listeners (one click → N fires). Guard at the top: `if (gd.__isletBridge) return; gd.__isletBridge = true;`.

**Gotcha**: ggplotly surfaces the `key` aes to the browser ONLY on the **individual-point marker traces** (`gd.data[i].key`), NOT on summary-point/line traces — so the Plot tab needs "Show individual points" ON for a click to resolve a key. Verified headlessly (chromote): trace keys like `112|Islet_1`; a real native gesture fires exactly one `[PLOT CLICK] native bridge` / `[TRAJ CLICK] native bridge` and renders the panel.

**Shared resolver** — `select_islet_from_key(case_id, islet_key, prepared, selected_islet)` in `segmentation_helpers.R`: the common tail of both drill-down paths (SF/lookup guards, centroid lookup with zero-padded case-id fallback, donor-metadata join, `selected_islet(list(...))`). Called by both the legacy plotly_click observers and the new native-click observers in Plot + Trajectory.

## Statistics Tab

See [`docs/statistics.md`](docs/statistics.md) for full details — pseudoreplication fix (Phase 16), `utils_stats.R` reference, unified outlier handling, zero-value exclusion, sex-stratified BH correction, 5-section narrative layout, multiple-testing correction.

## Trajectory Tab Features

> **PERF FIX (2026-06-29, Senior data):** the Trajectory tab was **slow to first render** (~10–30s),
> not broken — `traj_data_clean` called `get_islet_annotations()` once per islet (9,130×), each
> scanning the full 9.6k-row `islet_spatial_lookup` (~88M comparisons **per render**). Replaced with a
> single vectorized `match()` (identical result, O(n)) → renders instantly. NB the data + islet-key
> join were always fine (keys are `Islet_755` on both sides; `prep_data()` adds `islet_key`, `traj()`
> synthesizes `combined_islet_id`); the "0 islet metrics" in an early screenshot was just the
> pre-`input$traj_feature` frame. See [`docs/capability_status.md`](docs/capability_status.md).

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

## Viewer Tab (Avivator) — iframe Persistence (Jun 2026)

The Viewer embeds a local Avivator build (`www/avivator`, Vite bundle) in an `<iframe>` to render OME-TIFFs. Files: `mod_viewer_server.R`, `mod_viewer_ui.R`, `viewer_helpers.R`.

**Problem (macOS-only flashing)**: a single `renderUI` emitted the dropdown AND the iframe together, so any reactive invalidation destroyed + recreated the `<iframe>` (a WebGL canvas). Repeated canvas add/remove flickers heavily under Firefox's WebRender compositor on macOS (Mozilla bug #1555544); Retina also doubles deck.gl's drawing buffer (`useDevicePixels`).

**Architecture (the fix)**:
- **Split outputs**: `output$image_selector` (dropdown; depends on the stable available-images list) and `output$viewer_frame` (iframe; mounts once).
- **Mount-once iframe**: `viewer_frame` depends on `has_image_stable()` (flip-only boolean, below) + `current_tab()`, and reads the URL via `isolate(viewer_info()$iframe_src)` so it does NOT re-render when the image URL changes.
- **In-place src update**: `observeEvent(viewer_info()$iframe_src, ...)` → `session$sendCustomMessage("islet_update_viewer_src", list(id, src))`. A JS handler registered once in `mod_viewer_ui.R` sets `iframe.src` only when it differs from a `data-relsrc` attribute (dedup). Switching images updates the existing node — never recreates it (verified via chromote: same DOM node, src updated).
- **`probe_viewer_asset` gated** behind `VIEWER_DEBUG_ENABLED` (it's a blocking `curl` HEAD); `session$clientData` URL reads `isolate()`d (origin is fixed per session).

**CRITICAL gotcha — flip-only boolean**: a plain `reactive()` boolean invalidates its dependents whenever its INPUTS change (Shiny propagates on invalidation, NOT value-equality). So a `has_image` reactive reading `input$selected_image` made `viewer_frame` re-render (recreate the iframe) on every image switch. Fix: mirror it into a `reactiveVal` via `observe({ has_image_stable(isTRUE(has_image_now())) })` — `reactiveVal` dedups, so it only propagates when the boolean actually FLIPS (no-image → image), keeping the iframe alive across image switches.

**Secondary mitigation** (if macOS flicker persists after the above): rebuild Avivator with `useDevicePixels: false`. The bundle only reads `image_url` + `channel_config` query params, so this can't be toggled at runtime.

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
- **`reactiveVal` to stabilise a `renderUI` boolean**: a `reactive()` boolean re-fires its dependents whenever its INPUTS change, even if the value is unchanged — so a `renderUI` keyed on it rebuilds (destroying child DOM, including iframes/WebGL canvases). Mirror the boolean into a `reactiveVal` via an `observe`; `reactiveVal` dedups, so the `renderUI` only re-renders when the boolean truly flips. See § "Viewer Tab (Avivator) — iframe Persistence".
- **Robust plotly clicks**: `plotly_click` can be suppressed on macOS (clicking does nothing, no notification). Use `islet_click_bridge()` (native mousedown/mouseup → `Shiny.setInputValue`) for click-to-select, and guard `onRender` with an attach-once flag so listeners don't stack. See § "macOS Native-Click Bridge".

## Phenotype Rules — Marker Definitions in the UI (May 2026)

Source of truth: `data/phenotype_rules.csv`. Parsed once at startup in
`R/00_globals.R` into `PHENOTYPE_RULES_CSV` (19 phenotypes, marker columns
`CD99, INS, GCG, SST, SMA, CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR, PDPN,
CD34, CD31, B3TUBB, CD56`). The composer walks each phenotype's parent chain
(e.g. `CD8 T cell` → `T cell` → `Immune`) and unions positive/negative
constraints, then strips redundant `anypos` groups that are supersets of a
stricter constraint downstream. Four residual phenotypes (`ECAD+`, `APCs`,
`Structural`, `Unknown`) have no CSV row — they live in `PHENOTYPE_RESIDUALS`
with hand-written descriptions.

Naming aliases: PHENOTYPE_COLORS (drilldown_helpers.R) uses
`CD8a Tcell` / `CD4 Tcell` / `Blood Vessel` / `SMA+`; the CSV uses
`CD8 T cell` / `CD4 T cell` / `Blood vessel` / `Smooth Muscle`. The
`PHENOTYPE_ALIASES` map normalises PHENOTYPE_COLORS → CSV before lookup.

Helpers:

- `get_phenotype_rule(name)` — returns `list(parent, pos, neg, anypos_groups,
  description)` for either CSV form or PHENOTYPE_COLORS form. Handles aliases
  and residuals.
- `format_phenotype_rule(name, style = "inline" | "html" | "plain")` — formats
  the rule as `"CD3e+, CD8a+  —  CD4−, CD20−"` (inline), styled HTML
  (`<strong>` positives, muted negatives), or plain text for plot subtitles.
- `extract_phenotype_from_column(col)` — resolves a sanitised column name
  (e.g. `prop_CD8a_Tcell`, `peri_prop_Beta_cell`, `min_dist_Macrophage`,
  `enrich_z_B_cell`) back to a phenotype label by stripping known prefixes
  and reversing `_`→space / `plus`→`+`.
- `phenotype_hint_ui(name, show_modal_link = "<id>")` — small blue inline
  hint block embedded under each phenotype selector. The trailing
  "see all rules" actionLink triggers the same modal as the tab-top button.
- `phenotype_rules_button(input_id)` — pill-styled actionLink placed at the
  top of every tab. Each tab uses a unique non-namespaced id
  (`show_phenotype_rules_<tab>`); a single shared loop in `app.R`
  registers an `observeEvent` per id that fires `showModal(phenotype_rules_modal_ui())`.
- `phenotype_rules_modal_ui()` — returns the modal: a sortable table of all
  23 phenotype names (19 CSV + 4 residuals) with category, positive markers,
  negative markers, and notes.

Per-tab integration:

- **Plot tab** (`mod_plot_server.R::output$which_ui`) — when in `"Cell
  Populations"` mode, parses `input$which` via `extract_phenotype_from_column()`
  and emits `phenotype_hint_ui()` directly under the metric `selectInput`.
- **Spatial Card A** — `output$infiltration_pheno_hint` and
  `output$scatter_pheno_hint` look up the selected `safe` name via
  `phenotype_options()` and emit the hint UI.
- **Spatial Card C** — `output$distance_pheno_hint` and `output$kde_pheno_hint`
  follow the same pattern. The aggregate `"Immune (all)"` selection emits a
  custom description rather than a per-phenotype rule.
- **Drill-down composition panel** — root-level `output$islet_drilldown_rules`
  in `app.R` renders an expandable `<details>` element listing colour-coded
  swatches + rules for every phenotype present in the clicked islet. Embedded
  via `uiOutput("islet_drilldown_rules")` inside both the Plot and Trajectory
  `segmentation_viewer_panel` renderUIs.

When adding a new phenotype-selecting input anywhere in the app, follow the
existing pattern: resolve the selected value to a label via
`extract_phenotype_from_column()` (or your existing `phenotype_options()`
record), then emit `phenotype_hint_ui(label, show_modal_link = ns("…"))` and
register an `observeEvent(input$…, showModal(phenotype_rules_modal_ui()))` to
wire the inline link.

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
