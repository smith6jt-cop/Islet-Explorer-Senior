# R-lane ownership disposition

This is the execution ledger for consolidation Phase 5 / plan Track D.2. The
source table is `docs/CONSOLIDATION.md` § R lane. No source repository was
modified while resolving it.

| Family | Final owner / disposition | Evidence |
|---|---|---|
| domain configuration, pixel size, grouping palettes, phenotype vocabulary | `apps/R/ieCore` engine; exact configs under `apps/islet/config` and `apps/lymph/config` | derived from `LymphVizkit@820bb38`; exact YAML hash in `apps/R/ieCore/IMPORT_PROVENANCE.md` |
| `utils_stats` | `apps/R/ieCore/R/utils_stats.R` | Junior correctness harvest already landed at `259a1d5`; existing 73 assertions retained, DOMAIN ordinal coding added |
| `utils_safe_join` | `apps/R/ieCore/R/utils_safe_join.R` | one implementation; unit-key extraction now reads the active DOMAIN |
| plot typography/theme | `apps/R/ieCore/R/theme.R` | moved without changing the established Islet font values; app source smoke passes |
| `ai_helpers`, AI modules | `apps/islet/R` | current Senior implementation retained; no lymph fork copied |
| `segmentation_helpers` | `apps/islet/R` | Senior geometry/data paths retained; `PIXEL_SIZE_UM` now comes from DOMAIN |
| `viewer_helpers`, `mod_viewer_server`, `rdeck_helpers` | `apps/islet/R` | Senior-only TissUUmaps/MinIO/OME-Zarr/rdeck capabilities retained |
| `spatial_helpers` | `apps/islet/R` | Junior phenotype-order fix already landed at `621e0a4`; Senior DuckDB backend retained |
| `mod_statistics_server` | `apps/islet/R` | Junior implementation already harvested at `60b470a`; tab contract and stats suite pass |
| `mod_trajectory_server` | `apps/islet/R` | reconciled Senior implementation retained; tab contract passes |
| `mod_spatial_server` | `apps/islet/R` | Senior rdeck path and current download behavior retained; tab contract passes |
| `drilldown_helpers` | `apps/islet/R` | UI palette now reads DOMAIN; legacy and Senior labels remain explicit display compatibility values |
| `mod_inflammation_*` | `apps/islet/R` | Modernist CLR module already harvested at `661fa27`; remains feature-gated |
| Modernist shell (`cohort_filter`, `mod_cohort`, rail, bslib shell) | deliberately superseded | Junior merged it at `f827106`, explicitly reverted it at `9cf3d50`, then released the old layout at `f2a5c54`; D.2 does not reverse that decision silently |
| `app/shiny_redesign` | absent / superseded | no runnable duplicate imported |
| Lymph app fork | absent / replaced by config + data | `apps/lymph` contains no `app.R` or `R/`; contract test enforces this |
| Junior review-hardening `e589961` | portable controls reconciled in D.4 | Senior now has fixed-schema validation audit events, browser-input/path allowlists, production-safe AI/error handling, sanitized launch flags, and accessibility fixes. Junior's OCI evidence pipeline is deliberately not copied because it depends on Junior's release backend; see `SOURCE_RECONCILIATION_2026-09-04.md` |

Exit evidence:

- `R CMD build apps/R/ieCore` followed by `R CMD check` on the source tarball
  finishes `Status: OK`.
- `Rscript apps/islet/tests/parse_check.R` parses the package and app sources.
- `Rscript apps/islet/tests/run_tests.R` passes 188 assertions, including
  exact islet/lymph domain checks, the synthetic lymph schema, one-definition
  checks, UI/server wiring for Plot, Trajectory, Viewer, Statistics, Spatial,
  Inflammation, and AI, and the D.4 production security boundaries.
- A production-flagged app-source smoke with DuckDB, mirai, and AI disabled loads the islet
  config from `apps/islet/config/islet_default.yml` and returns a server
  function.
