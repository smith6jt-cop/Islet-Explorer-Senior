# testthat helper — sourced automatically before any test-*.R file. It loads
# the shared ieCore sources exactly as the Shiny bootstrap does, without
# loading the heavy app or scientific data.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# base R >= 4.4 ships `%||%`; utils_stats.R's summarize_outlier_filter() uses it.
# Define a fallback on the global environment (where source() places utils_stats.R
# functions) so the suite runs under older R too.
if (!exists("%||%", mode = "function")) {
  assign("%||%", function(a, b) if (is.null(a)) b else a, envir = globalenv())
}

# Locate a repository file regardless of the working directory test_dir() uses.
.find_repo_file <- function(rel) {
  cands <- c(rel,
             file.path("apps", "islet", rel),
             file.path("..", rel),
             file.path("../..", rel),
             file.path("../../..", rel),
             file.path("../../../..", rel))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) {
    stop(sprintf("Could not locate %s from wd=%s", rel, getwd()))
  }
  normalizePath(hit[1])
}

.find_source <- function(rel) .find_repo_file(file.path("apps", "islet", rel))

.iecore_dir <- dirname(.find_repo_file("apps/R/ieCore/R/domain_config.R"))
for (.file in sort(list.files(.iecore_dir, pattern = "[.]R$", full.names = TRUE))) {
  sys.source(.file, envir = globalenv())
}

.islet_domain <- ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml"))
ie_set_domain(.islet_domain)

sys.source(.find_source("R/00_audit.R"), envir = globalenv())
