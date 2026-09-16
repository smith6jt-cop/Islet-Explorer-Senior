# ---------- Shared cohort filter: pure implementation ----------
# One definition of "which islets are in scope", used by every tab. The module
# in R/mod_cohort.R is a thin reactive wrapper around these.
#
# WHY THIS FILE IS PURE: it is the single place the row subset is defined, so it
# must be drivable from shiny::testServer AND from a plain Rscript with no
# session. That is the harness that caught the silent no-ops in the first place.
#
# PORTED FROM JUNIOR (f827106) with ONE substantive change: every `GROUP_LEVELS`
# became `domain_group_levels()`. Junior has one cohort and can afford the
# literal; this app carries two with no shared vocabulary, so a literal here
# would filter the run cohort against the CODEX cohort's disease labels -- and
# because the filter is "keep rows whose status is in this set", a mismatch
# empties every panel rather than erroring.
#
# WHY EVERY FILTER RETURNS AN AUDIT ROW: every lookup in this app is by exact
# column name, so a missing or renamed column DISABLES a filter rather than
# erroring. That has already cost this project twice -- the AAb `.x`/`.y` join
# collision, and the min-cells filter that held Capillary at 5,380 islets and
# p = 0.3476 from min-cells 1 through 200 while the CSV header still reported the
# value the user had set. So cohort_apply() records, per filter: was it
# requested, which column did it bind to, did it apply, how many rows went. The
# root CLAUDE.md rule -- "any new filter guarded on a column's presence must
# surface the skip, never swallow it" -- is mechanised here rather than
# re-implemented per tab.

COHORT_COL_CANDIDATES <- list(
  status     = c("Donor Status", "donor_status", "status"),
  donor      = c("Case ID", "case_id", "donor_id", "imageid"),
  age        = c("age", "donor_age"),
  sex        = c("gender", "sex", "donor_sex"),
  diameter   = c("islet_diam_um"),
  cells_core = c("total_cells_core"),
  cells_peri = c("total_cells_peri")
)

COHORT_AAB_COLS <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")

cohort_resolve_cols <- function(df) {
  cn <- colnames(df %||% data.frame())
  lapply(COHORT_COL_CANDIDATES, function(cands) {
    hit <- intersect(cands, cn)
    if (length(hit)) hit[[1]] else NA_character_
  })
}

cohort_spec_default <- function(age_range = c(0, 100)) {
  list(groups = domain_group_levels(), age = age_range, sex = NULL,
       aab_flags = COHORT_AAB_COLS, min_cells = 1, diam = c(0, 500),
       outlier_enabled = FALSE, outlier_z = 3.0)
}

.cohort_audit_row <- function(filter, requested, applied, bound, n_before, n_after, note = "") {
  data.frame(filter = filter, requested = requested, applied = applied,
             bound_col = bound %||% NA_character_,
             n_before = n_before, n_after = n_after,
             removed = n_before - n_after, note = note, stringsAsFactors = FALSE)
}

#' Apply the cohort row subset to any islet-level frame.
#'
#' ORDER IS FIXED AND LOAD-BEARING: groups -> age -> sex -> AAb -> min-cells ->
#' diameter. Outlier flagging is NOT done here; it runs in each module AFTER
#' this, which is what makes "N islets excluded" count against the tested set
#' rather than a superset of it.
#'
#' A NOTE IS ABOUT THE FRAME, NOT THE REQUEST. The audit's `note` used to read
#' "no sex column in this frame" whenever the sex filter was simply set to
#' "all" -- true of the request, false of the frame, and the frame is what the
#' reader is being told about. It is empty now unless the column really is
#' absent.
#'
#' @param cells_expr "core_peri" (default) or "core". Trajectory used to bind
#'   min-cells to `total_cells`, which is core-only -- verified 5214/5214 islets
#'   match core and 0/5214 match core+peri, so the same slider gave 130 islets
#'   there and 309 everywhere else.
cohort_apply <- function(df, spec, cells_expr = c("core_peri", "core")) {
  cells_expr <- match.arg(cells_expr)
  if (is.null(df) || !is.data.frame(df)) return(df)
  cols  <- cohort_resolve_cols(df)
  audit <- list()
  n0 <- nrow(df)

  # ---- 1. donor status ----
  req_groups <- !setequal(spec$groups %||% domain_group_levels(), domain_group_levels())
  if (!is.na(cols$status)) {
    n <- nrow(df)
    df <- df[as.character(df[[cols$status]]) %in% (spec$groups %||% domain_group_levels()), , drop = FALSE]
    audit[[length(audit) + 1]] <- .cohort_audit_row("groups", TRUE, TRUE, cols$status, n, nrow(df))
  } else {
    audit[[length(audit) + 1]] <- .cohort_audit_row("groups", req_groups, FALSE, NA, n0, n0,
      "no donor-status column in this frame")
  }

  # ---- 2. donor age ----
  if (!is.null(spec$age) && length(spec$age) == 2 && !is.na(cols$age)) {
    n <- nrow(df); a <- suppressWarnings(as.numeric(df[[cols$age]]))
    df <- df[is.finite(a) & a >= spec$age[1] & a <= spec$age[2], , drop = FALSE]
    audit[[length(audit) + 1]] <- .cohort_audit_row("age", TRUE, TRUE, cols$age, n, nrow(df))
  } else {
    audit[[length(audit) + 1]] <- .cohort_audit_row("age", !is.null(spec$age), FALSE, NA,
      nrow(df), nrow(df),
      if (is.na(cols$age)) "no age column in this frame" else "")
  }

  # ---- 3. sex ----
  # NA-SAFE. `df[[col]] == value` yields NA for a missing sex, and NA-indexing a
  # data frame returns a row of NAs rather than dropping it. `%in%` returns FALSE.
  if (!is.null(spec$sex) && length(spec$sex) && !is.na(cols$sex)) {
    n <- nrow(df)
    df <- df[as.character(df[[cols$sex]]) %in% spec$sex, , drop = FALSE]
    audit[[length(audit) + 1]] <- .cohort_audit_row("sex", TRUE, TRUE, cols$sex, n, nrow(df))
  } else {
    audit[[length(audit) + 1]] <- .cohort_audit_row("sex", !is.null(spec$sex), FALSE, NA,
      nrow(df), nrow(df),
      if (is.na(cols$sex)) "no sex column in this frame" else "")
  }

  # ---- 4. autoantibody flags (Aab+ donors only) ----
  # Semantics ported verbatim from the Plot sidebar: an UNCHECKED flag EXCLUDES
  # donors carrying it, and no flags at all excludes every Aab+ donor.
  aab_avail <- intersect(COHORT_AAB_COLS, colnames(df))
  aab_req   <- !setequal(spec$aab_flags %||% character(0), COHORT_AAB_COLS)
  if (length(aab_avail) && !is.na(cols$status)) {
    n <- nrow(df)
    is_aabp <- as.character(df[[cols$status]]) == "Aab+"
    if (any(is_aabp)) {
      keep_aabp <- rep(TRUE, sum(is_aabp))
      if (!length(spec$aab_flags)) {
        keep_aabp[] <- FALSE
      } else {
        unchecked <- intersect(setdiff(COHORT_AAB_COLS, spec$aab_flags), aab_avail)
        if (length(unchecked)) {
          m <- as.data.frame(df[is_aabp, unchecked, drop = FALSE])
          for (cc in names(m)) m[[cc]] <- as.logical(m[[cc]])
          keep_aabp <- !(rowSums(m, na.rm = TRUE) > 0)
        }
      }
      keep <- rep(TRUE, nrow(df)); keep[is_aabp] <- keep_aabp
      df <- df[keep, , drop = FALSE]
    }
    audit[[length(audit) + 1]] <- .cohort_audit_row("aab", TRUE, TRUE,
      paste(aab_avail, collapse = "+"), n, nrow(df))
  } else {
    audit[[length(audit) + 1]] <- .cohort_audit_row("aab", aab_req, FALSE, NA,
      nrow(df), nrow(df),
      if (!length(aab_avail)) "no AAb_* columns in this frame"
      else "no donor-status column, so Aab+ donors cannot be identified")
  }

  # ---- 5. min cells / islet ----
  mc <- suppressWarnings(as.numeric(spec$min_cells %||% 1))
  # A cleared numericInput sends NA, which used to crash the Statistics run.
  if (!(length(mc) == 1 && is.finite(mc) && mc >= 1)) mc <- 1
  mc_req <- mc > 1
  need <- if (cells_expr == "core") "cells_core" else c("cells_core", "cells_peri")
  have <- all(!is.na(unlist(cols[need])))
  if (mc_req && have) {
    n <- nrow(df)
    tc <- suppressWarnings(as.numeric(df[[cols$cells_core]]))
    tc <- ifelse(is.finite(tc), tc, 0)
    if (cells_expr == "core_peri") {
      tp <- suppressWarnings(as.numeric(df[[cols$cells_peri]]))
      tc <- tc + ifelse(is.finite(tp), tp, 0)
    }
    df <- df[tc >= mc, , drop = FALSE]
    audit[[length(audit) + 1]] <- .cohort_audit_row("min_cells", TRUE, TRUE,
      paste(unlist(cols[need]), collapse = " + "), n, nrow(df),
      if (cells_expr == "core") "core only" else "")
  } else {
    audit[[length(audit) + 1]] <- .cohort_audit_row("min_cells", mc_req, !mc_req, NA,
      nrow(df), nrow(df),
      if (mc_req) "per-islet cell counts unavailable in this frame" else "")
  }

  # ---- 6. islet diameter ----
  if (!is.null(spec$diam) && length(spec$diam) == 2 && !is.na(cols$diameter)) {
    n <- nrow(df); d <- suppressWarnings(as.numeric(df[[cols$diameter]]))
    df <- df[is.finite(d) & d >= spec$diam[1] & d <= spec$diam[2], , drop = FALSE]
    audit[[length(audit) + 1]] <- .cohort_audit_row("diameter", TRUE, TRUE, cols$diameter, n, nrow(df))
  } else {
    audit[[length(audit) + 1]] <- .cohort_audit_row("diameter", !is.null(spec$diam), FALSE, NA,
      nrow(df), nrow(df),
      if (is.na(cols$diameter)) "no islet_diam_um column in this frame" else "")
  }

  attr(df, "cohort_audit") <- do.call(rbind, audit)
  attr(df, "cohort_spec")  <- spec
  df
}

#' Donor keep-list for CELL-LEVEL frames (Spatial tissue CSVs, cell_distance).
#' Those carry case_id and donor_status but no age/sex/AAb/diameter, so the
#' donor-scope filters resolve against the islet table and land as an id
#' membership test. min-cells and diameter are islet-level and never apply.
cohort_donor_ids <- function(donor_tbl, spec) {
  if (is.null(donor_tbl) || !nrow(donor_tbl)) return(character(0))
  keep <- cohort_apply(donor_tbl, utils::modifyList(spec, list(min_cells = 1, diam = NULL)))
  cols <- cohort_resolve_cols(keep)
  if (is.na(cols$donor)) return(character(0))
  unique(as.character(keep[[cols$donor]]))
}

#' Red banner for any filter the user asked for that this frame could not apply.
cohort_audit_badge <- function(df, context = "") {
  a <- attr(df, "cohort_audit")
  if (is.null(a)) return(NULL)
  bad <- a[a$requested & !a$applied, , drop = FALSE]
  if (!nrow(bad)) return(NULL)
  shiny::tags$div(
    class = "ie-callout ie-callout-warn",
    shiny::tags$strong("Cohort filter not applied here. "),
    if (nzchar(context)) paste0("(", context, ") "),
    paste(sprintf("%s — %s", bad$filter, bad$note), collapse = "; "),
    " The value you set is recorded in the CSV header with a NOT APPLIED marker."
  )
}

#' THE single source of CSV provenance for cohort scope. Generated from the spec
#' plus the audit, so an export cannot omit a filter the panel owns and cannot
#' silently report one that was skipped.
cohort_provenance <- function(spec, audit = NULL) {
  sx <- c(F = "Female", M = "Male")
  mark <- function(name) {
    if (is.null(audit)) return("")
    r <- audit[audit$filter == name, , drop = FALSE]
    if (nrow(r) && r$requested[1] && !r$applied[1])
      paste0("  *** NOT APPLIED: ", r$note[1], " ***") else ""
  }
  ex <- setdiff(COHORT_AAB_COLS, spec$aab_flags %||% character(0))
  c("[COHORT PANEL] Filters below apply to every analysis tab.",
    paste0("Donor groups: ", paste(spec$groups, collapse = ", "), mark("groups")),
    paste0("Donor age range: ", spec$age[1], " - ", spec$age[2], " years", mark("age")),
    paste0("Sex: ", if (is.null(spec$sex)) "all"
           else paste(ifelse(spec$sex %in% names(sx), sx[spec$sex], spec$sex), collapse = ", "),
           mark("sex")),
    paste0("AAb filter: ",
           if (!length(spec$aab_flags)) "all Aab+ donors EXCLUDED"
           else if (length(ex)) paste("excluding donors positive for", paste(ex, collapse = ", "))
           else "all Aab+ donors included", mark("aab")),
    paste0("Min cells per islet: ", spec$min_cells, mark("min_cells")),
    paste0("Islet diameter range: ", spec$diam[1], " - ", spec$diam[2], " um", mark("diameter")),
    paste0("Outlier rule: ",
           if (isTRUE(spec$outlier_enabled))
             paste0("ON, robust median/MAD |z| > ", spec$outlier_z, " per donor status")
           else "OFF"))
}
