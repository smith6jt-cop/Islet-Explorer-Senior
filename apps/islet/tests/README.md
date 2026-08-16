# Shiny App Tests & CI

Automated checks for the R Shiny app (`apps/islet/`). Run in CI by
`.github/workflows/shiny-app-ci.yml` on every push/PR that touches
`apps/islet/**`.

## What runs

1. **Parse check** — `tests/parse_check.R` syntax-parses every `R/*.R` file plus
   `app.R`. Catches syntax errors without loading any package or data.
2. **Unit tests** — `tests/run_tests.R` runs the `testthat` suite in
   `tests/testthat/`. Currently covers `R/utils_stats.R` — the pure statistical
   helpers behind the Statistics tab (donor-level aggregation, mixed-effects test,
   per-bin ANOVA / Kendall τ / pairwise, Cohen's d, η², normality, outlier
   flagging). These need no app data and no Python bridge.

3. **Statistics audit harness** — `tests/manual/stats_audit.R`. **Manual/local only**
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
`lmerTest`. Install with:

```r
install.packages(c("testthat", "dplyr", "tidyr", "tibble", "lme4", "lmerTest"))
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
Rscript -e 'shiny::runApp(".", port = 7777)'
```

## Adding tests

Drop `test-*.R` files in `tests/testthat/`. Dependencies and the code under test
are loaded once by `tests/testthat/helper-setup.R` (sourced automatically before
the test files). Keep new tests data-independent so they run in CI.

## Auto-install R in Claude Code web sessions (SessionStart hook)

A fresh Claude Code on the web container has **no R**, so the checks above can't run
until R is installed. This repo ships a `SessionStart` hook that installs R + the
test packages automatically at the start of every web session:

- `.claude/hooks/session-start.sh` — installs `r-base-core` plus the precompiled
  `r-cran-*` Debian packages (testthat, dplyr, tidyr, tibble, lme4, lmerTest). It is
  idempotent (skips what's already present), remote-only (`$CLAUDE_CODE_REMOTE`
  guard, so it's a no-op on local machines), and always exits 0 so it never blocks
  session startup.
- `.claude/settings.json` — registers the hook under `hooks.SessionStart`.

The hook runs **synchronously**, so the session waits until R is ready (no race
where the agent runs tests before the packages exist). The `r-cran-*` debs install
fast on Ubuntu and bundle their own system dependencies (BLAS/LAPACK, gfortran
runtime) — no compilation. Verified by running the full suite (73 assertions pass as of
2026-08-14; the count drifts as tests are added — re-run rather than trusting this number,
and note the suite was red from `d403563` to `2026-08-14` because that commit widened
`cohens_d()`'s return value without updating its test) after
a clean install.

The hook only enables the automated **checks**. The **full app** still cannot be
launched in a web session — it needs the gitignored `data/` artifacts (see above).
To disable the hook, remove the `SessionStart` entry from `.claude/settings.json`.
