# Every data path the app reads must resolve inside the ACTIVE domain's root.
#
# Before these existed, only master/h5ad/parquet were domain-aware and nine other
# constants were hard-coded to ../../data/app_data. A run-cohort session would then read
# the CODEX cohort's cells, outlines and rules -- and the two cohorts share five donors, so
# an islet key exists on both sides and the app would serve the wrong cells for a
# right-looking islet. That is a correctness bug, not a cosmetic one, which is why the
# leak is asserted against rather than the happy path.

CODEX_ROOT <- file.path("..", "..", "data", "app_data")

.app_data_paths <- function() {
  list(
    cells       = file.path(domain_data_root(), "cells"),
    donors      = file.path(domain_data_root(), "donors"),
    json        = file.path(domain_data_root(), "json"),
    gson        = file.path(domain_data_root(), "gson"),
    zarr        = file.path(domain_data_root(), "zarr_images"),
    annotations = file.path(domain_data_root(), "annotations.tsv"),
    lookup      = file.path(domain_data_root(), "islet_spatial_lookup.csv"),
    rules       = file.path(domain_data_root(), "phenotype_rules.csv"),
    adata       = file.path(domain_data_root(), "adata_ins_root.h5ad"),
    parquet     = file.path(domain_data_root(), "parquet"),
    master      = domain_data_file("master", "master_results.xlsx"),
    h5ad        = domain_data_file("h5ad", "islet_explorer.h5ad")
  )
}

test_that("a domain with no data: block resolves to the historical CODEX literals", {
  ie_set_domain(ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml")))
  paths <- .app_data_paths()

  # These are the exact strings the constants carried before they were routed through
  # domain_data_root(); the CODEX build must be byte-identical in behaviour.
  expect_identical(domain_data_root(), CODEX_ROOT)
  expect_identical(paths$cells,       file.path(CODEX_ROOT, "cells"))
  expect_identical(paths$donors,      file.path(CODEX_ROOT, "donors"))
  expect_identical(paths$json,        file.path(CODEX_ROOT, "json"))
  expect_identical(paths$gson,        file.path(CODEX_ROOT, "gson"))
  expect_identical(paths$annotations, file.path(CODEX_ROOT, "annotations.tsv"))
  expect_identical(paths$lookup,      file.path(CODEX_ROOT, "islet_spatial_lookup.csv"))
  expect_identical(paths$rules,       file.path(CODEX_ROOT, "phenotype_rules.csv"))
  expect_identical(paths$adata,       file.path(CODEX_ROOT, "adata_ins_root.h5ad"))
  expect_identical(paths$master,      file.path(CODEX_ROOT, "master_results.xlsx"))
  expect_identical(paths$h5ad,        file.path(CODEX_ROOT, "islet_explorer.h5ad"))
})

test_that("no app data path escapes the active domain's root", {
  for (config in c("apps/islet/config/islet_default.yml",
                   "apps/islet/config/run20_default.yml")) {
    if (!file.exists(suppressWarnings(tryCatch(.find_repo_file(config),
                                               error = function(e) "")))) next
    ie_set_domain(ie_load_domain(.find_repo_file(config)))
    root <- domain_data_root()
    for (name in names(.app_data_paths())) {
      path <- .app_data_paths()[[name]]
      if (is.null(path)) next
      expect_true(startsWith(path, root),
                  info = sprintf("%s: %s escaped domain root %s", config, name, root))
    }
  }
})

test_that("the run cohort declares its own root, distinct from the CODEX one", {
  config <- tryCatch(.find_repo_file("apps/islet/config/run20_default.yml"),
                     error = function(e) NULL)
  skip_if(is.null(config), "run cohort config not present")

  ie_set_domain(ie_load_domain(config))
  expect_false(identical(domain_data_root(), CODEX_ROOT))
  expect_true(grepl("typing", domain_data_root(), fixed = TRUE))
  # The two cohorts share five donors; serving one from the other's root is the failure
  # this separation exists to prevent.
  expect_false(startsWith(domain_data_root(), CODEX_ROOT))
})
