# Parity: the R filter and the Python authority must agree on the full donor list, and both
# must exclude 6457 and 6579. `ie_cohort.R` only ever reads the roster built by
# `scripts/senior/build_cohort_roster.py`; this test checks that roster against a *live* call
# into `phenocycler.cohort` so a stale JSON cannot pass silently.

.cohort_authority_path <- .find_repo_file("core/phenocycler/phenocycler/cohort.py")
.cohort_repo_root <- dirname(dirname(dirname(dirname(.cohort_authority_path))))
.cohort_json_path <- .find_repo_file("apps/islet/config/cohort.json")

#' Call the Python authority directly and return its eligible-donor answer for `donors`.
.python_eligible_donors <- function(donors, repo_root) {
  python_bin <- Sys.which("python3")
  if (!nzchar(python_bin)) skip("python3 not available to run the parity check")

  phenocycler_dir <- normalizePath(file.path(repo_root, "core", "phenocycler"))
  script <- sprintf(
    paste(
      "import sys, json",
      "sys.path.insert(0, %s)",
      "from phenocycler.cohort import filter_eligible_donors",
      "donors = json.loads(sys.argv[1])",
      "print(json.dumps(filter_eligible_donors(donors)))",
      sep = "\n"
    ),
    jsonlite::toJSON(phenocycler_dir, auto_unbox = TRUE)
  )
  script_file <- tempfile(fileext = ".py")
  on.exit(unlink(script_file), add = TRUE)
  writeLines(script, script_file)

  donors_json <- as.character(jsonlite::toJSON(donors))
  out <- system2(
    python_bin, c(shQuote(script_file), shQuote(donors_json)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("python authority call failed: ", paste(out, collapse = "\n"))
  }
  jsonlite::fromJSON(out[length(out)])
}

test_that("cohort.json names 6457 and 6579 as excluded, with a measured reason each", {
  cohort <- ie_cohort_load(.cohort_json_path)

  expect_true(all(c("6457", "6579") %in% ie_cohort_excluded_ids(cohort)))
  expect_true(all(c("6457", "6579") %in% cohort$all_donors))

  reasons <- setNames(cohort$excluded_donors$reason, cohort$excluded_donors$donor_id)
  expect_true(nzchar(reasons[["6457"]]))
  expect_true(nzchar(reasons[["6579"]]))
})

test_that("ie_cohort_filter() removes 6457 and 6579 from a donor vector", {
  cohort <- ie_cohort_load(.cohort_json_path)

  filtered <- ie_cohort_filter(c("6457", "6579", "6539", "115"), cohort = cohort)
  expect_setequal(filtered, c("6539", "115"))
  expect_false(any(c("6457", "6579") %in% filtered))
})

test_that("ie_cohort_filter() normalises donor ids before comparing", {
  cohort <- ie_cohort_load(.cohort_json_path)

  # '0112' and '112' are one donor (docs/ISLET_LAYERS_SPEC.md ss 0.3); leading zeros on an
  # eligible donor must not make it look excluded, and vice versa.
  expect_identical(ie_cohort_norm_donor("06457"), "6457")
  expect_true(ie_cohort_is_excluded("06457", cohort = cohort))
  expect_identical(ie_cohort_filter("06457", cohort = cohort), character(0))
})

test_that("the R filter and the live Python authority agree on the full donor list", {
  cohort <- ie_cohort_load(.cohort_json_path)

  r_eligible <- ie_cohort_filter(cohort$all_donors, cohort = cohort)
  py_eligible <- .python_eligible_donors(cohort$all_donors, .cohort_repo_root)

  expect_setequal(r_eligible, py_eligible)
  expect_setequal(r_eligible, cohort$eligible_donors)

  expect_false(any(c("6457", "6579") %in% r_eligible))
  expect_false(any(c("6457", "6579") %in% py_eligible))
})
