# Lymph domain lane

This directory is intentionally **not another Shiny application**. It contains:

- the byte-preserved LymphVizkit domain contract in
  `config/lymphoid_default.yml`; and
- a 588 KB synthetic, non-scientific fixture under `data/app_data/` for schema
  and tab-contract testing.

Shared behavior lives once in `apps/R/ieCore`. To inspect this domain contract:

```bash
Rscript -e 'source("apps/R/ieCore/R/domain_config.R"); d <- ie_load_domain("apps/lymph/config/lymphoid_default.yml"); print(domain_group_levels(d))'
```

The fixture demonstrates file shape only. It is not evidence for biological
thresholds, identities, disease effects, or a production lymphoid workflow.
Source hashes are recorded in `apps/R/ieCore/IMPORT_PROVENANCE.md`.
