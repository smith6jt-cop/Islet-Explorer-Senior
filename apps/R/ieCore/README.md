# ieCore

`ieCore` is the single shared R implementation for application-level domain
configuration, plot typography, safe joins, and statistical helpers.

It is deliberately downstream-only. A `DOMAIN` file controls display labels,
palette vocabulary, and the pixel-to-micrometre conversion used by the viewer;
it cannot define marker positivity, cell identity, donor inclusion, or upstream
scientific thresholds.

The islet app loads these sources from `apps/islet/R/00_00_iecore.R` before
Shiny auto-sources its modules. The directory is also an installable R package:

```bash
R CMD INSTALL apps/R/ieCore
Rscript apps/islet/tests/run_tests.R
```

Domain contracts:

- `apps/islet/config/islet_default.yml` — the legacy CODEX cohort, 15 donors,
  `data/app_data`
- `apps/islet/config/run20_default.yml` — the PhenoCycler run cohort, 20 donors
  (17 clean + 3 marginal), `data/development/typing/app`. Its display vocabulary
  is the eight classes the subpanel classification decides plus one declared
  *absence*, `Unresolved`, which no producer writes today
- `apps/lymph/config/lymphoid_default.yml` — the byte-preserved LymphVizkit
  contract and its synthetic fixture

`ISLET_DOMAIN_CONFIG` selects which islet contract the app loads. The two islet
domains share five donors and no vocabulary and are never merged.

The loader validates required scalars, unique group/phenotype vocabularies,
palette coverage, a positive pixel size, and exactly declared region keys. It
fails closed on an incomplete contract.
