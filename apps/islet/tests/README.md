# Shiny App Tests & CI

Automated checks for the R Shiny app (`apps/islet/`). Run in CI by
`.github/workflows/islet-r.yml` on every push/PR that touches `apps/islet/**`,
`apps/R/ieCore/**`, or `apps/lymph/**`.

## What runs

1. **Parse check** — `tests/parse_check.R` syntax-parses the ieCore package,
   every app `R/*.R` file, and `app.R`.
2. **Unit and contract tests** — `tests/run_tests.R` covers the shared ieCore
   statistics/domain helpers, exact islet/lymph configs, the committed synthetic
   lymph fixture, and UI/server wiring for every tab. These need no scientific
   app data and no Python bridge.

3. **Tab smoke driver** — `tests/tab_smoke.R`. **Manual/local only** (not in CI, and not
   part of the testthat suite — it lives in `tests/`, not `tests/testthat/`, so `test_dir()`
   never picks it up). It drives Chrome through the tabs named in its own `TABS` vector and
   reports, per tab, the R output errors, the JS console errors, and how many outputs exist
   but drew nothing. It needs `NOT_CRAN=true`, or `AppDriver$new()` skips with
   "Reason: On CRAN", which reads like a failure and is really a refusal to run.

   ```bash
   NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml \
     Rscript tests/tab_smoke.R /tmp/tabshots
   ```

4. **Statistics audit harness** — `tests/manual/stats_audit.R`. **Manual/local only**
   (not in CI): it drives the *real* `plot_server` → `statistics_server` chain
   headlessly with `shiny::testServer` against the live H5AD, and pins each of the
   six defects from the Jul 2026 correctness pass. This is the layer that catches
   silent no-ops and size-dependent method switches — the class of bug the unit
   tests structurally cannot see, because those live in module reactives rather
   than in `utils_stats.R`. Skips cleanly (exit 0) when data or packages are absent.

## Run locally

```bash
cd apps/islet
Rscript tests/parse_check.R      # syntax check
Rscript tests/run_tests.R        # unit tests

# Statistics audit (needs data/ + reticulate + a working H5AD path)
export ISLET_AUDIT_CACHE=/tmp/islet_audit_cache.rds   # optional: persist the ~40s load
Rscript tests/manual/stats_audit.R
```

Expected on healthy data (2026-08-14): **26 checks, 0 failures**, reproducing the
audit's own numbers — min-cells sweeping 5382 → 4246 → 689 islets, pairwise
`df ≈ 12.1` (Satterthwaite, not `Inf`), and the donor-level AAb test at
`p = 0.6703` against the islet-level `1.6e-15` it replaced.

### Two traps if you write more testServer code here

- **`testServer` starts every unset input at `NULL`, and `NULL` is not always
  "unfiltered".** `raw_df_base()` drops *every* `Aab+` islet when
  `input$aab_flags` is empty. The live UI ships all five boxes ticked, so this
  never bites a user — but a harness that omits it silently runs on a two-group
  cohort while appearing to request three. `stats_audit.R` sets `aab_flags`
  explicitly and asserts all three donor groups arrived before testing anything.
- **The two modules are chained**, so one `testServer` call cannot drive both:
  `statistics_server` consumes `plot_returns$raw_df_prezero`. Stage 1 drives
  `plot_server` and grabs `session$returned$raw_df_prezero()` as a plain data
  frame; stage 2 passes it in as a constant reactive. `stats_run()` is an
  `eventReactive(ignoreInit = TRUE)` — fire it with `setInputs(run_tests = 1)`.

Dependencies for the tests: `testthat`, `dplyr`, `tidyr`, `tibble`, `lme4`,
`lmerTest`, `yaml`, and `readxl`. Install with:

```r
install.packages(c("testthat", "dplyr", "tidyr", "tibble", "lme4", "lmerTest",
                   "yaml", "readxl"))
```

## Why CI does not launch the full app

The complete app needs (a) heavy runtime deps — `sf` (GDAL/GEOS/PROJ),
`reticulate` + `anndata` (a Python bridge), `shiny`, `plotly` — and (b) the data
artifacts in `data/` (H5AD ~70 MB, per-islet/per-donor CSVs ~280 MB, GeoJSON,
`phenotype_rules.csv`). `data/` is **gitignored**, so it is absent in CI and in a
fresh clone, and `R/00_globals.R` reads some of it at startup. CI therefore
parse-checks the app files and unit-tests the data-independent utilities.

To launch the full app locally (with the data present):

```bash
cd apps/islet
# ALWAYS name the domain. With ISLET_DOMAIN_CONFIG unset the app falls back to the
# legacy CODEX contract and serves the wrong dataset under right-looking islet keys --
# the two domains share five donors, so nothing downstream catches it.
ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  Rscript -e 'shiny::runApp(".", port = 7777)'
```

R **4.5.1** with the `~/R/…/4.5` library. `/usr/bin/Rscript` is 4.6.1 and silently
lacks the dependencies.

## Adding tests

Drop `test-*.R` files in `tests/testthat/`. Dependencies and the code under test
are loaded once by `tests/testthat/helper-setup.R` (sourced automatically before
the test files). Keep new tests data-independent so they run in CI.

## Installing R

R must be installed manually if it is not already available, including in Claude
Code web sessions. This repository does not provide a `SessionStart` hook that
installs R. After installing R, install the test packages listed above before
running the checks. Launching the full app also requires its runtime dependencies
and the gitignored `data/` artifacts described above.
