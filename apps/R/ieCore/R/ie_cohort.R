# One cohort definition, read by both languages.
#
# `core/phenocycler/phenocycler/cohort.py` is the authority for donor exclusions (pinned by
# `core/phenocycler/tests/current/test_cohort.py`). `scripts/senior/build_cohort_roster.py`
# reads that authority and `core/phenocycler/cohort-spec.json` to emit
# `apps/islet/config/cohort.json`. This file reads that JSON and filters a donor vector -- it
# never retypes an exclusion or its reason.
#
# `donor_id` is a string and is NOT normalised across trees: `0112` and `112` are one donor
# (docs/ISLET_LAYERS_SPEC.md ss 0.3). Every donor id this file writes or compares is normalised
# with `ie_cohort_norm_donor()` first.
#
# WHAT THIS ROSTER IS, AND WHAT IT IS NOT. `eligible_donors` is the acquisition roster minus the
# central exclusions -- nothing more. It is NOT the typed cohort, and the two differ. Measured
# 2026-09-18: the built `cohort.json` carries 23 candidates and 22 eligible (only 6457 excluded;
# 6579's exclusion was lifted 2026-09-17), while `data/development/typing/cells/` holds 20 donor
# partitions -- 6533 and 6579 are eligible and have no typed partition at all. So a card that
# needs "donors this build can show" must ask the loaded data or its own artifact, and use this
# only to REMOVE an excluded donor from a list it already has.

#' Default location of the built cohort roster, relative to a repo root.
#' @export
IE_COHORT_JSON_REL <- file.path("apps", "islet", "config", "cohort.json")

#' Normalise a donor id the way `scripts/senior/build_cohort_roster.py` does.
#'
#' `"0112" -> "112"`. Never `as.integer()`: `6533` must not become `6533.0` through a
#' numeric column, and a donor id is never a number to this function, only a string.
#' @export
ie_cohort_norm_donor <- function(donor_id) {
  v <- trimws(as.character(donor_id))
  stripped <- sub("^0+", "", v)
  ifelse(!nzchar(stripped), "0", stripped)
}

#' Load the built cohort roster.
#'
#' @param path path to `cohort.json`; defaults to `IE_COHORT_JSON_REL` under `repo_root`
#' @param repo_root repo root to resolve `path` against when `path` is relative and does not
#'   already exist relative to the working directory
#' @return a list with `all_donors`, `excluded_donors` (a data.frame of `donor_id`/`reason`),
#'   and `eligible_donors`, all normalised
#' @export
ie_cohort_load <- function(path = IE_COHORT_JSON_REL, repo_root = NULL) {
  resolved <- path
  if (!file.exists(resolved) && !is.null(repo_root)) {
    resolved <- file.path(repo_root, path)
  }
  if (!file.exists(resolved)) {
    stop("ie_cohort_load(): no cohort roster at '", resolved, "'. Build it with ",
         "'python scripts/senior/build_cohort_roster.py'.", call. = FALSE)
  }
  raw <- jsonlite::fromJSON(resolved, simplifyVector = TRUE)

  excluded_donors <- data.frame(
    donor_id = ie_cohort_norm_donor(raw$excluded_donors$donor_id),
    reason = as.character(raw$excluded_donors$reason),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      all_donors = ie_cohort_norm_donor(raw$all_donors),
      excluded_donors = excluded_donors,
      eligible_donors = ie_cohort_norm_donor(raw$eligible_donors),
      sources = raw$sources,
      generated_at = raw$generated_at
    ),
    class = "ie_cohort"
  )
}

#' The excluded donor ids in a loaded cohort roster.
#' @export
ie_cohort_excluded_ids <- function(cohort) {
  if (!inherits(cohort, "ie_cohort")) {
    stop("ie_cohort_excluded_ids(): not an ie_cohort (load one with ie_cohort_load())",
         call. = FALSE)
  }
  cohort$excluded_donors$donor_id
}

#' Filter a donor vector down to the eligible cohort.
#'
#' Donor exclusion is a decision made once, at the build side (`build_cohort_roster.py`
#' reading `phenocycler.cohort.DONOR_EXCLUSIONS`) -- this only ever applies that decision to a
#' vector the caller already has. It never retypes an exclusion or invents a new one.
#'
#' @param donors character vector of donor ids, any padding
#' @param cohort an `ie_cohort` from [ie_cohort_load()]; loaded with defaults if omitted
#' @return `donors`, normalised, with excluded donor ids removed
#' @export
ie_cohort_filter <- function(donors, cohort = ie_cohort_load()) {
  normalised <- ie_cohort_norm_donor(donors)
  normalised[!normalised %in% ie_cohort_excluded_ids(cohort)]
}

#' Is this a normalised donor id excluded from the cohort?
#' @export
ie_cohort_is_excluded <- function(donor_id, cohort = ie_cohort_load()) {
  ie_cohort_norm_donor(donor_id) %in% ie_cohort_excluded_ids(cohort)
}
