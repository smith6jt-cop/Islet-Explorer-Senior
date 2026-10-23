# One cohort definition, read by both languages.
#
# `core/multiplex5000/multiplex5000/cohort.py` is the authority for donor exclusions (pinned by
# `core/multiplex5000/tests/current/test_cohort.py`). `scripts/senior/build_cohort_roster.py`
# reads that authority and `core/multiplex5000/cohort-spec.json` to emit
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
#
# AN EXCLUSION CAN COVER ONE PLATFORM ONLY. Each excluded entry records the platforms its
# exclusion covers (`platforms`): 6457's covers its PhenoCycler data (maintainer, 2026-10-22:
# "the xenium for 6457 is not excluded"). A caller filtering a list of one platform's donors
# passes `platform =`; with none, every excluded donor is removed, as before. An entry that does
# not say which platforms it covers is read as covering all of them, and a roster too old to
# name its platforms refuses a per-platform question rather than answering it.

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
#' @return a list with `platforms`, `all_donors`, `excluded_donors` (a data.frame of
#'   `donor_id`/`reason`), `excluded_platforms` (donor id -> the platforms its exclusion
#'   covers), and `eligible_donors`, all normalised
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
  platforms <- as.character(raw$platforms %||% character(0))
  # Fail closed: an entry that does not say which platforms it covers covers all of them.
  covered <- raw$excluded_donors$platforms
  excluded_platforms <- lapply(seq_len(nrow(excluded_donors)), function(i) {
    entry <- if (is.null(covered)) NULL else covered[[i]]
    if (length(entry) == 0) platforms else as.character(entry)
  })
  names(excluded_platforms) <- excluded_donors$donor_id

  structure(
    list(
      platforms = platforms,
      all_donors = ie_cohort_norm_donor(raw$all_donors),
      excluded_donors = excluded_donors,
      excluded_platforms = excluded_platforms,
      eligible_donors = ie_cohort_norm_donor(raw$eligible_donors),
      sources = raw$sources,
      generated_at = raw$generated_at
    ),
    class = "ie_cohort"
  )
}

#' The excluded donor ids in a loaded cohort roster.
#'
#' @param platform `NULL` for every excluded donor, or one of `cohort$platforms` for the donors
#'   whose exclusion covers that platform
#' @export
ie_cohort_excluded_ids <- function(cohort, platform = NULL) {
  if (!inherits(cohort, "ie_cohort")) {
    stop("ie_cohort_excluded_ids(): not an ie_cohort (load one with ie_cohort_load())",
         call. = FALSE)
  }
  ids <- cohort$excluded_donors$donor_id
  if (is.null(platform)) return(ids)
  # A misspelt platform matches no scoped exclusion and would readmit the donor: refuse it.
  if (length(platform) != 1 || !platform %in% cohort$platforms) {
    stop("ie_cohort_excluded_ids(): unknown platform '", paste(platform, collapse = ", "),
         "'; one of ", paste(cohort$platforms, collapse = ", "), call. = FALSE)
  }
  ids[vapply(ids, function(id) platform %in% cohort$excluded_platforms[[id]], logical(1))]
}

#' Filter a donor vector down to the eligible cohort.
#'
#' Donor exclusion is a decision made once, at the build side (`build_cohort_roster.py`
#' reading `phenocycler.cohort.DONOR_EXCLUSIONS`) -- this only ever applies that decision to a
#' vector the caller already has. It never retypes an exclusion or invents a new one.
#'
#' @param donors character vector of donor ids, any padding
#' @param cohort an `ie_cohort` from [ie_cohort_load()]; loaded with defaults if omitted
#' @param platform `NULL` to remove every excluded donor, or the platform the donors' data
#'   comes from, to remove only those whose exclusion covers it
#' @return `donors`, normalised, with excluded donor ids removed
#' @export
ie_cohort_filter <- function(donors, cohort = ie_cohort_load(), platform = NULL) {
  normalised <- ie_cohort_norm_donor(donors)
  normalised[!normalised %in% ie_cohort_excluded_ids(cohort, platform)]
}

#' Is this a normalised donor id excluded from the cohort (on `platform`, when given)?
#' @export
ie_cohort_is_excluded <- function(donor_id, cohort = ie_cohort_load(), platform = NULL) {
  ie_cohort_norm_donor(donor_id) %in% ie_cohort_excluded_ids(cohort, platform)
}
