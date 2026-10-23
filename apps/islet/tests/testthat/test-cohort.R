# Parity: the R filter and the Python authority must agree on the full donor list, and both
# must exclude 6457 and 6579. `ie_cohort.R` only ever reads the roster built by
# `scripts/senior/build_cohort_roster.py`; this test checks that roster against a *live* call
# into `phenocycler.cohort` so a stale JSON cannot pass silently.

.cohort_authority_path <- .find_repo_file("core/multiplex5000/multiplex5000/cohort.py")
.cohort_repo_root <- dirname(dirname(dirname(dirname(.cohort_authority_path))))
.cohort_json_path <- .find_repo_file("apps/islet/config/cohort.json")

#' Call the Python authority directly and return its eligible-donor answer for `donors`.
# The exclusion list itself, asked of the same authority. Naming donors here is the restatement
# CLAUDE.md forbids -- and it went stale the day 6579 was reinstated, failing this file rather
# than anything real.
.python_excluded_donors <- function(repo_root) {
  python_bin <- Sys.which("python3")
  if (!nzchar(python_bin)) skip("python3 not available to run the parity check")
  phenocycler_dir <- normalizePath(file.path(repo_root, "core", "multiplex5000"))
  script <- sprintf(
    paste("import sys, json",
          "sys.path.insert(0, %s)",
          "from phenocycler.cohort import DONOR_EXCLUSIONS",
          "print(json.dumps(sorted(DONOR_EXCLUSIONS)))",
          sep = "\n"),
    shQuote(phenocycler_dir, type = "cmd"))
  out <- system2(python_bin, c("-c", shQuote(script)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    skip(paste("python authority unavailable:", paste(out, collapse = " ")))
  }
  jsonlite::fromJSON(tail(out, 1))
}

.python_eligible_donors <- function(donors, repo_root) {
  python_bin <- Sys.which("python3")
  if (!nzchar(python_bin)) skip("python3 not available to run the parity check")

  phenocycler_dir <- normalizePath(file.path(repo_root, "core", "multiplex5000"))
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

test_that("cohort.json names every excluded donor, with a measured reason each", {
  cohort <- ie_cohort_load(.cohort_json_path)
  excluded <- .python_excluded_donors(.cohort_repo_root)

  # Derived from the authority, never named here. This file used to assert the pair
  # {6457, 6579} and so could not survive 6579 being reinstated on 2026-09-17 -- it failed on
  # a correct change, which is the failure mode "no lane may restate the list" exists to stop.
  expect_setequal(ie_cohort_excluded_ids(cohort), excluded)
  expect_true(all(excluded %in% cohort$all_donors))

  reasons <- setNames(cohort$excluded_donors$reason, cohort$excluded_donors$donor_id)
  for (donor in excluded) expect_true(nzchar(reasons[[donor]]))
})

test_that("ie_cohort_filter() removes exactly the excluded donors from a donor vector", {
  cohort <- ie_cohort_load(.cohort_json_path)
  excluded <- .python_excluded_donors(.cohort_repo_root)
  keep <- c("6539", "115")

  filtered <- ie_cohort_filter(c(excluded, keep), cohort = cohort)
  expect_setequal(filtered, keep)
  expect_false(any(excluded %in% filtered))
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

  # Derived, never named. Whoever is in DONOR_EXCLUSIONS is out of both lists and nobody else is,
  # so a reinstatement flows through instead of breaking this file.
  excluded <- .python_excluded_donors(.cohort_repo_root)
  expect_true(length(excluded) > 0)
  expect_false(any(excluded %in% r_eligible))
  expect_false(any(excluded %in% py_eligible))
  expect_setequal(setdiff(cohort$all_donors, r_eligible), intersect(cohort$all_donors, excluded))
})

# The platforms each exclusion covers, asked of the same authority.
.python_excluded_platforms <- function(repo_root) {
  python_bin <- Sys.which("python3")
  if (!nzchar(python_bin)) skip("python3 not available to run the parity check")
  phenocycler_dir <- normalizePath(file.path(repo_root, "core", "multiplex5000"))
  script <- sprintf(
    paste("import sys, json",
          "sys.path.insert(0, %s)",
          "from phenocycler.cohort import DONOR_EXCLUSIONS, PLATFORMS, excluded_platforms",
          "print(json.dumps({'platforms': list(PLATFORMS),",
          "                  'covers': {d: list(excluded_platforms(d)) for d in DONOR_EXCLUSIONS}}))",
          sep = "\n"),
    jsonlite::toJSON(phenocycler_dir, auto_unbox = TRUE))
  script_file <- tempfile(fileext = ".py")
  on.exit(unlink(script_file), add = TRUE)
  writeLines(script, script_file)
  out <- system2(python_bin, shQuote(script_file), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
    skip(paste("python authority unavailable:", paste(out, collapse = " ")))
  }
  jsonlite::fromJSON(tail(out, 1), simplifyVector = FALSE)
}

test_that("cohort.json records the platforms each exclusion covers, as the authority does", {
  cohort <- ie_cohort_load(.cohort_json_path)
  authority <- .python_excluded_platforms(.cohort_repo_root)

  expect_identical(cohort$platforms, unlist(authority$platforms))
  for (donor in names(authority$covers)) {
    expect_setequal(cohort$excluded_platforms[[donor]], unlist(authority$covers[[donor]]))
  }
})

test_that("a platform-scoped filter keeps 6457's Xenium and refuses its PhenoCycler", {
  # The maintainer, 2026-10-22: "the xenium for 6457 is not excluded".
  cohort <- ie_cohort_load(.cohort_json_path)
  expect_identical(cohort$excluded_platforms[["6457"]], "phenocycler")

  expect_identical(ie_cohort_filter(c("06457", "115"), cohort = cohort, platform = "xenium"),
                   c("6457", "115"))
  expect_identical(ie_cohort_filter(c("06457", "115"), cohort = cohort, platform = "if"),
                   c("6457", "115"))
  expect_identical(ie_cohort_filter(c("06457", "115"), cohort = cohort, platform = "phenocycler"),
                   "115")
  expect_identical(ie_cohort_filter(c("06457", "115"), cohort = cohort), "115")
  expect_false(ie_cohort_is_excluded("6457", cohort = cohort, platform = "xenium"))
  expect_true(ie_cohort_is_excluded("6457", cohort = cohort, platform = "phenocycler"))
  expect_true(ie_cohort_is_excluded("6457", cohort = cohort))
})

test_that("a misspelt platform is refused, and an unscoped entry covers every platform", {
  cohort <- ie_cohort_load(.cohort_json_path)
  expect_error(ie_cohort_filter("6457", cohort = cohort, platform = "Xenium"), "unknown platform")

  roster <- jsonlite::fromJSON(.cohort_json_path, simplifyVector = FALSE)
  roster$excluded_donors[[1]]$platforms <- NULL
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  jsonlite::write_json(roster, path, auto_unbox = TRUE)
  unscoped <- ie_cohort_load(path)
  donor <- unscoped$excluded_donors$donor_id[[1]]
  for (platform in unscoped$platforms) {
    expect_true(ie_cohort_is_excluded(donor, cohort = unscoped, platform = platform))
  }
})
