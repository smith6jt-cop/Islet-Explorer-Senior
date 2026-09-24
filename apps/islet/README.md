# Islet Explorer application

This is the monorepo's modular R Shiny consumer for islet, spatial,
trajectory, viewer, statistics, inflammation, and optional AI views.

## Two domains, selected by `ISLET_DOMAIN_CONFIG`

The app serves **two deliberately separate datasets**. They share five donors and
no vocabulary, so they are never merged — a cross-donor comparison spanning both
would compare methods, not biology.

| domain | config | data root |
|---|---|---|
| legacy CODEX, 15 donors | `config/islet_default.yml` | `../../data/app_data` |
| PhenoCycler run cohort, 20 donors (17 clean + 3 marginal) | `config/run20_default.yml` | `../../data/development/typing/app` |

`ISLET_DOMAIN_CONFIG` picks one; with it unset the app falls back to the CODEX
contract. Every data path resolves through `ieCore::domain_data_root()` /
`domain_data_file()`, and `tests/testthat/test-domain-data-paths.R` asserts that a
path cannot escape the active domain's root.

The run cohort's display vocabulary has **nine entries and eight cell types**: the
eight classes the subpanel classification decides, plus `Unresolved`, which is an
*absence* rather than a class. **No cell carries `Unresolved` today** — the
refinement pass that writes it was ruled in on 2026-09-23 and is wired as the
`build_cells.py --refine` join (into `refined_class`, beside `cell_class`), and that
join has not run; `typing/app` is built from `cell_class`, eight classes plus
`Excluded`. It is declared ahead of its data so a card meets a named
absence instead of an unmapped value. See the header of `config/run20_default.yml`.

These configs are downstream display contracts. They must not be used to define
upstream marker calls, positivity or cell identities.

## Running it

Shared domain, typography, safe-join and statistical code lives once in
`../R/ieCore`; `R/00_00_iecore.R` loads it before the tab modules. Run from this
directory so the relative data roots above stay exact. **R 4.5.1 with the
`~/R/…/4.5` library** — `/usr/bin/Rscript` is 4.6.1 and silently lacks the
dependencies.

```bash
cd apps/islet
Rscript tests/parse_check.R
Rscript tests/run_tests.R

# run cohort
ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  Rscript -e 'shiny::runApp(".", port = 7777)'

# legacy CODEX domain, with the heavy backends off
ISLET_USE_DUCKDB=FALSE ISLET_USE_MIRAI=FALSE ISLET_ENABLE_AI=FALSE \
  Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

The development server is not a deployment test: Avivator needs HTTP Range, which
`runApp` does not serve.

Per-tab conventions are in `CLAUDE.md`; longer write-ups in `docs/`. Deployment
and the deliberately separate publication build are documented in
`../../deploy/DEPLOY.md`. The lymphoid configuration/fixture is under `../lymph`;
it is not a second app fork.
