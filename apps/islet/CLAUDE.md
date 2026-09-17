# CLAUDE.md — Shiny app (`apps/islet`)

Per-tab conventions. The root `CLAUDE.md` covers the entrypoint, the two domains, and the data
boundary. Longer write-ups live in [`docs/spatial.md`](docs/spatial.md),
[`docs/statistics.md`](docs/statistics.md) and
[`docs/capability_status.md`](docs/capability_status.md).

## Domain boundary

`apps/R/ieCore` is the sole implementation of domain loading, plot typography, safe joins and the
statistical helpers. `R/00_00_iecore.R` must stay the first auto-sourced file; it loads the DOMAIN
and fails closed on an incomplete contract.

`PIXEL_SIZE_UM`, the grouping palettes and the phenotype display vocabulary come from the DOMAIN.
They are **presentation metadata** and must never define upstream positivity or cell identity.

**Every data path resolves through `domain_data_root()` / `domain_data_file()`.** A path that
escapes the active domain's root serves the wrong cells for a right-looking islet, because the two
cohorts share five donors and an islet key exists on both sides. `tests/testthat/test-domain-data-paths.R`
asserts against that leak and pins the CODEX literals.

**A cohort fact written as a literal is a domain leak too, and the suite cannot see it.** The app
asserted "Aab+ is a single donor (6521)", "donors (N=15 total)", "N=5 per group", "all fifteen",
`choices = c("ND", "Aab+", "T1D")` and "(scVI)" — every one a fact about the CODEX cohort, printed
as a measurement of whichever dataset is loaded. Every path resolved and every column existed, so
238 assertions passed while the run cohort's pages named a donor it does not contain. Derive these
from the loaded data or the DOMAIN; never write a cohort's size, its donor ids or its typing method
into a string.

`tests/testthat/tab_smoke.R` is what catches them: it drives Chrome through all five tabs and
reports per tab the R output errors, the JS console errors, and how many outputs exist but drew
nothing. Run it whenever a tab's data changes — it needs `NOT_CRAN=true`, or `AppDriver$new()`
skips with "Reason: On CRAN", which reads like a failure and is really a refusal to run.

```bash
NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml Rscript tests/tab_smoke.R /tmp/tabshots
```

## Tabs

Plot, Trajectory, Viewer, Statistics, Spatial — plus Inflammation, hidden unless
`SHOW_INFLAMMATION_TAB=1`, and the AI panel inside Plot, off unless `ISLET_ENABLE_AI=TRUE` and
hard-blocked in production regardless.

**Plot.** Cell-population proportions (`prop_*`, `peri_prop_*`, and a combined form computed on the
fly), hormone fractions, and Microenvironment. Labels are unit-explicit — never the bare word
"density" without naming its denominator. Counts read "Cells per islet (core)"; densities name the
region whose QuPath-measured area is the denominator. Microenvironment is **area-valued**
(`100 × area_um2 / region_um2`): a nerve may segment into several pieces, so object counts measure
morphology, not biology, and are not exposed.

**Trajectory.** Loads its own H5AD separately from the app's main one, resolved through the domain
root with an `ADATA_PATH` fallback. Point size defaults to `sqrt(total_cells)`; LOESS weights are
`log1p(total_cells)` — raw counts let one 1,900-cell islet dominate the fit. The validity caption
reads `uns['trajectory_validation']` and is NULL-guarded throughout.

**Spatial.** Tissue scatter and the Leiden UMAP come from the DuckDB views; the neighbourhood cards
need the per-region columns. Card C's per-cell distance view reads the `cell_distance` view on a
symlog axis with bin edges anchored at 0, the islet border.

**Statistics.** Entirely downstream of the Plot sidebar's current selection and shares its DOM. It
keeps its own outlier control rather than the app-wide one, deliberately.

**Drill-down** is embedded in Plot and Trajectory: a click sets `selected_islet()` and the panel
renders inline. Cell centroids are micrometres and outlines are pixels — convert with
`PIXEL_SIZE_UM`. The `colour` aesthetic is reserved for structure outlines; cell scatter uses `fill`
(shape 21) to avoid a scale conflict.

## Viewer

Two backends, chosen by `ISLET_VIEWER_BACKEND` (`avivator` default, `tissuumaps`). All the iframe
machinery is shared; only the URL builder and the dropdown source differ.

**The load-bearing constraint for both:** the cohort mixes panels — batch 1 is 59-channel with
SST@35, batch 2 puts SST@1 and some slides are 58-channel — so default-on markers are emitted **by
name**, never by index, and resolved against each image's own OME-XML. Do not regress this, and do
not use the TissUUmaps fork's `/slide` auto-split, which names layers `<file>_Channel_<i>` and loses
marker identity.

**Avivator requires HTTP Range (206).** Production nginx/shiny-server serves it; `shiny::runApp`
returns 200 with the whole body, so images render only under the deployed server. TissUUmaps tiles
images itself and works under either.

Two Shiny traps this tab taught, both general:

- A `reactive()` boolean invalidates its dependents whenever its *inputs* change, so a `renderUI`
  keyed on one rebuilds and destroys child DOM — including a WebGL canvas. Mirror it into a
  `reactiveVal`, which dedups and only propagates on a real flip.
- `plotly_click` can be suppressed on macOS. `islet_click_bridge()` listens for native
  mousedown/mouseup instead; guard `onRender` with an attach-once flag or listeners stack.

## Conventions

- `selectInput` menus extend past their container — give card styles `overflow: visible`.
- Scatter plots above ~50k points use `renderPlot`, not plotly.
- Tissue scatter uses `coord_fixed() + scale_y_reverse()`, the microscopy convention.
- `updateSelectInput` targeting a `renderUI`-mounted input must share a reactive dependency with
  that `renderUI`, or it races the mount and silently leaves the dropdown empty.
- Reticulate: `ad$uns[[key]]` costs ~220 ms per crossing. Cache `uns` as an R list once, then index.
- Typography is `ISLET_FONT_SIZES` / `theme_islet()` / `plotly_islet_fonts()` / `plotly_islet_config()`.
  `ggplotly()` does not fully respect ggplot themes, so the explicit `layout()` overlay after every
  `ggplotly()` call is required, not optional.
- For a downloadable plot, build the ggplot in a `reactive()` and have both the renderer and the
  `downloadHandler` consume it — never build it twice.

## Phenotype rules in the UI

`phenotype_rules.csv` in the active domain's root is the source of truth, parsed once at startup.
The composer walks each phenotype's parent chain and unions its constraints. Every phenotype-
selecting input emits `phenotype_hint_ui()` under it and links to the shared modal; register the
`observeEvent` for any new tab's button, since the shared loop is keyed by id.
