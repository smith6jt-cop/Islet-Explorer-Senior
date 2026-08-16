#!/usr/bin/env Rscript
# Run the Shiny app's unit test suite (testthat).
#
#   cd apps/islet && Rscript tests/run_tests.R      # from the app dir
#   Rscript apps/islet/tests/run_tests.R            # from the repo root
#
# The suite covers R/utils_stats.R only — pure statistical helpers that need no
# app data or Python bridge. Dependencies loaded + code sourced by
# tests/testthat/helper-setup.R.

suppressPackageStartupMessages(library(testthat))

cands <- c("tests/testthat", "apps/islet/tests/testthat")
test_path <- cands[dir.exists(cands)]
if (length(test_path) == 0) {
  stop("Cannot find tests/testthat from wd = ", getwd())
}

testthat::test_dir(test_path[1], stop_on_failure = TRUE, reporter = "progress")
