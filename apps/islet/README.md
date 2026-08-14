# Shiny App: Islet Area Distributions

Interactive Shiny app to explore islet area distributions across donor status (ND → Aab+ → T1D) with:

- Band-region Islet_Markers: % positive per islet vs. size (mean ± SE)
- Band-region Islet_Targets: density/count vs. size (mean ± SE)
- Islet composition vs. size (INS/GCG/SST fractions)
- One-way ANOVA per size bin across donor groups (BH-adjusted p-values)

Data source: `../data/master_results.xlsx` (generated from TSVs in `../data/results` with donor key `CODEX_Pancreas_Donors.xlsx`).

## Run

Prereqs: R >= 4.2 with packages: `shiny`, `readxl`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `plotly`, `broom`.

```r
shiny::runApp("shiny")
```

Or from the shell:

```bash
R -q -e 'shiny::runApp("shiny", launch.browser=TRUE)'
```

## Notes

- Islet diameter is standardized to CORE area: derived from `Islet_Targets` rows where `type == "islet_core"` (core `region_um2`). Diameter = `2 * sqrt(core_area/pi)`.
- Band-region rows are used for markers (`region_type == "islet_band"`) and targets (`type == "islet_band"`).
- Binning is by diameter (µm); default bin width = 50 µm. Adjust in the sidebar.
- Markers use `% positive` (100 × `pos_frac`). Targets use `area_density` by class.
- Composition uses `*_any / cells_total` for INS/GCG/SST.
- Stats tab runs one-way ANOVA per bin across groups (ND, Aab+, T1D) and reports BH-adjusted p-values.

