# Islet Explorer application

This is the monorepo's modular R Shiny consumer for islet, spatial,
trajectory, viewer, statistics, inflammation, and optional AI views.

Application-wide display configuration lives in
`config/islet_default.yml`. Shared domain, typography, safe-join, and
statistical code lives once in `../R/ieCore`; `R/00_00_iecore.R` loads it before
the tab modules. These are downstream display contracts and must not be used to
define upstream marker calls or cell identities.

The app reads external/ignored data through `../../data/app_data`. Run it from
this directory so that path remains exact:

```bash
cd apps/islet
Rscript tests/parse_check.R
Rscript tests/run_tests.R
ISLET_USE_DUCKDB=FALSE ISLET_USE_MIRAI=FALSE ISLET_ENABLE_AI=FALSE \
  Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
```

Deployment and the deliberately separate publication build are documented in
`../../deploy/DEPLOY.md`. The lymphoid configuration/fixture is under
`../lymph`; it is not a second app fork.
