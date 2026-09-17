# Islet Explorer — coding agent instructions

The authoritative guides are root `CLAUDE.md` for repository scope and provenance,
and `apps/islet/CLAUDE.md` for app behavior, UI, plotting, and deployment details.
Follow those guides when changing the systems summarized here. Paths are repo-relative.

## Repository map

- `core/phenocycler/`: upstream QuPath measurements, preserving UUIDs as `object_id`.
  The current workflow ends at expression: ingest → geometry QC → REDSEA → expression.
  Run from this directory with its sole `config.ini`; never add a root `config.ini`.
- `analysis/xenium/`: separate Xenium 5K spatial-transcriptomics lane; use `scvi-env`.
  Preserve its numbered decisions #1–#17. `core/panels` is the sole panel authority;
  generated panel snapshots must not be hand-edited or committed as duplicate sources.
- `apps/islet/`: modular R Shiny app for aggregated islet, spatial, trajectory,
  single-cell, and image views.
- `apps/R/ieCore`: shared R package and sole home of domain loading, plot typography,
  safe joins, and statistical helpers. Implement shared helpers here, once.
- `apps/lymph/`: second DOMAIN config (`config/lymphoid_default.yml`) plus a
  synthetic fixture. It is deliberately not an app: keep it free of `R/` and `app.R`,
  as asserted by `apps/islet/tests/testthat/test-lymph-data-contract.R`.

## Modular app layout

`apps/islet/app.R` wires the app; approximately 25 files under `apps/islet/R/`
hold the modules and helpers. Put tab changes in the corresponding module.

- `00_00_iecore.R`: shared-package and DOMAIN bootstrap; must remain the first
  auto-sourced app file.
- `00_globals.R`: shared constants, state, and centralized donor colors;
  `00_audit.R`: audit helpers.
- `data_loading.R`: application data loading and preparation.
- `mod_*_ui.R` / `mod_*_server.R`: paired modules for `plot`, `statistics`,
  `inflammation`, `spatial`, `trajectory`, `viewer`, and `ai_assistant`.
- `viewer_helpers.R`, `segmentation_helpers.R`, `drilldown_helpers.R`,
  `spatial_helpers.R`, `ai_helpers.R`, `rdeck_helpers.R`: specialized app helpers.

## DOMAIN, data, and provenance

- The default DOMAIN contract is `apps/islet/config/islet_default.yml`, selected
  through `ISLET_DOMAIN_CONFIG` and validated fail-closed by `ie_domain_validate()`.
  Do not bypass validation or supply silent defaults for an incomplete contract.
- DOMAIN pixel size, palettes, and phenotype vocabulary are presentation metadata;
  they must never define upstream marker positivity or cell identity.
- The app reads `data/app_data/`, preferring DuckDB/parquet, then H5AD, then Excel.
  Preserve this loading priority when changing `data_loading.R`.
- The legacy application H5AD is not connected to current content-addressed source
  runs. Do not present its rebuild as a supported bridge from current Phenocycler.
  A future adapter must join by donor and `object_id` and preserve run/manifest IDs.
- Donor `6457` is centrally excluded in `core/phenocycler/phenocycler/cohort.py`.
  Derive from `DONOR_EXCLUSIONS`, never from a literal, and no lane may restate the list.
  **`6579`'s exclusion was lifted 2026-09-17** and it is recorded in `REINSTATED` with why;
  the pancreatitis it was excluded for is a disease context, not a reason to drop it.
- There are no verified production target/reference pairs. Pair-dependent controls,
  calibration, typing, state calling, and assignment export remain blocked upstream.
  Never infer positivity from expected frequencies, uniformity, or disease ordering.
  Existing classifications and reviews are audit context, not training or release gates.
- A missing measurement is unavailable, not zero or negative. Do not feed provisional
  or stale labels into application scripts or revive retired upstream workflows.
  Legacy roots still consumed by app Card C remain distinct from current source runs.

## App coding contracts

- Preserve donor status order `ND < Aab+ < T1D`. The default palette is **Bright**:
  ND `#3366CC`, Aab+ `#FF6600`, T1D `#109618`, centralized in `00_globals.R`.
  Use `active_donor_colors()` and existing reactive palette helpers instead of
  hardcoding colors in plots; keep the sidebar, trajectory, and spatial selectors synced.
- Use **Sex** in user-facing labels; the stored column remains `gender`.
- Core parents are `Islet_N`; peri-islet parents are `Islet_N_exp20um`.
  Check `total_cells_peri > 0` before using peri-islet metrics.
- Synthetic `islet_union` rows in `prep_data()` sum core + band and must include
  `area_um2`, `region_um2`, and `area_density = area_um2 / region_um2`.
- Name metric units and denominators explicitly. Microenvironment values are
  `100 * area_um2 / region_um2`; cell density uses cells per mm² of the selected region.
- Use `scales::pseudo_log_trans(base = 10)` when zeros must remain visible.
  Tissue scatter uses `coord_fixed() + scale_y_reverse()`; convert cell coordinates
  from micrometers to polygon pixels using DOMAIN `PIXEL_SIZE_UM`.
- Guard root-level outputs by the active tab and namespace module-local inputs.
  Coordinate `updateSelectInput()` with `renderUI` mounting to avoid empty dropdowns.
- Reuse ieCore typography through `theme_islet()`, `plotly_islet_fonts()`, and
  `plotly_islet_config()`; apply explicit font layout after `ggplotly()`.
  Build downloadable ggplots once in a reactive shared by rendering and `downloadHandler`.
- Preserve donor-aware statistics and the pseudoreplication correction; reuse ieCore
  statistics helpers and the app guide's statistics conventions.
- Cache reticulate `ad$uns` as an R list before repeated indexing. Preserve deferred
  loading of heavy data and vectorized lookups.
- Keep `islet_click_bridge()` and its attach-once guard for reliable plot selection.
  Preserve the Viewer iframe across image switches; update its `src` in place and
  stabilize presence booleans with `reactiveVal` to prevent DOM/WebGL recreation.

## Viewer and development

- `ISLET_VIEWER_BACKEND` selects default `avivator` or `tissuumaps`.
  Avivator streams local OME-TIFFs and needs HTTP Range (206); ordinary `runApp`
  serving is insufficient for its image rendering. TissUUmaps serves `.tmap` projects
  through a separate service and works with the development Shiny server.
- Resolve default channels by marker name from each image's own OME-XML; panels vary.
  Preserve INS/red, GCG/blue, SST/yellow, DAPI/grey defaults; avoid positional indices.
- The AI panel is opt-in via `ISLET_ENABLE_AI=TRUE`; public deployments omit it.
- Install R manually before running the app or its tests. Development launch:

  ```bash
  cd apps/islet
  Rscript -e 'shiny::runApp(".", port = 7777)'
  ```

- Production uses nginx and shiny-server; follow `deploy/DEPLOY.md` and the app guide.
- Keep generated data, run manifests, output-bearing notebooks, caches, credentials,
  and local model/image paths out of Git. The repo uses subtrees, not submodules.
