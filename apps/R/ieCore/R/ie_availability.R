# Five states, because "missing" is four different things and collapsing them lies.
#
# The repo rule is "a missing measurement is unavailable, not zero or negative". Honouring it
# needs more than a NA check, because the distinctions live in the DATA and are easy to erase:
#
#   measured     the value exists
#   blank        never typed. In cells_xen this is lineage_native == "" -- all 960,284 cells
#                of donor 115 and 17 other donors, because the typer has not been run on them
#   unavailable  measured, then dropped at QC. lineage_native == "Unavailable" -- 53,626 cells
#                on 6436 and 52,533 on 6442. Something WAS measured here; it did not survive
#   withheld     measurable and deliberately not released. The post-Xenium IF hormone calls
#                sit below the 0.75 AUROC bar they were preregistered at (6436 INS 0.5143 /
#                GCG 0.5220; 6442 0.6571 / 0.7276)
#   failed       attempted and failed. Donor 115's restain gate at 7.198 um
#
# THE TRAP, measured: in `lineage_common` the QC-dropped cells and the never-typed donors are
# BYTE-IDENTICAL empty strings. Deriving identity state from that column silently converts a
# measured QC drop into an absence of typing. `lineage_native` is the only column that keeps
# them apart, so it is the only one this file reads.

#' @export
IE_STATES <- c("measured", "blank", "unavailable", "withheld", "failed")

#' An availability value: a state, a reason, and optionally a count
#' @export
ie_avail <- function(state, reason = NA_character_, n = NA_integer_, source = NA_character_) {
  if (length(state) != 1L || !state %in% IE_STATES) {
    stop("ie_avail(): state must be one of ", paste(IE_STATES, collapse = ", "),
         "; got ", paste(format(state), collapse = ", "), call. = FALSE)
  }
  if (state != "measured" && (is.na(reason) || !nzchar(reason))) {
    stop("ie_avail(): state '", state, "' requires a reason. An absence without a recorded ",
         "reason is indistinguishable from an oversight.", call. = FALSE)
  }
  structure(list(state = state, reason = reason, n = n, source = source),
            class = "ie_avail")
}

#' @export
ie_is_measured <- function(x) inherits(x, "ie_avail") && identical(x$state, "measured")

#' Identity state for a Xenium/PhenoCycler cell table
#'
#' Reads `lineage_native` ONLY. See the note at the top of this file: `lineage_common` cannot
#' distinguish a QC drop from an untyped donor, and reading it is the specific error this
#' function exists to prevent.
#' @export
ie_identity_state <- function(lineage_native) {
  if (is.null(lineage_native)) return(factor(character(0), levels = IE_STATES))
  v <- as.character(lineage_native)
  out <- ifelse(is.na(v) | !nzchar(v), "blank",
         ifelse(v == "Unavailable", "unavailable", "measured"))
  factor(out, levels = IE_STATES)
}

#' Human label for a state
#'
#' Never returns "Unknown". `Unknown` is a real phenotype in the DOMAIN display vocabulary with
#' its own colour (#CCCCCC); using it for an absence would put a measurement and a hole in the
#' same bucket on the same legend.
#' @export
ie_state_label <- function(state, n = NA_integer_) {
  base <- switch(as.character(state),
    measured    = "Measured",
    blank       = "Not typed",
    unavailable = "Measured, dropped at QC",
    withheld    = "Withheld",
    failed      = "Failed",
    stop("ie_state_label(): unknown state '", state, "'", call. = FALSE))
  if (identical(as.character(state), "measured") && !is.na(n)) {
    return(paste0(base, " (n = ", format(n, big.mark = ","), ")"))
  }
  base
}

#' Render an availability value as one line of text
#'
#' `ie_avail("measured", n = 0)` and `ie_avail("unavailable", ...)` must never produce the same
#' string: a measured zero is a finding, an unavailable is a hole, and a reader who cannot tell
#' them apart will treat the hole as a zero.
#' @export
ie_state_text <- function(x) {
  if (!inherits(x, "ie_avail")) stop("ie_state_text(): not an ie_avail", call. = FALSE)
  if (identical(x$state, "measured")) {
    n <- if (is.na(x$n)) "" else paste0(" — n = ", format(x$n, big.mark = ","))
    return(paste0(ie_state_label("measured"), n))
  }
  paste0(ie_state_label(x$state), " — ", x$reason)
}

#' Colours for the five states
#'
#' Deliberately outside the DOMAIN phenotype palette: an absence must not be able to collide
#' with a cell type on a shared legend.
#' @export
ie_state_palette <- function() {
  c(measured    = "#2B7A0B",
    blank       = "#9AA0A6",
    unavailable = "#495057",
    withheld    = "#7A6A00",
    failed      = "#A61B1B")
}

#' Guard a card on a measured value
#'
#' Returns a NAMED NOTE when the value is not measured -- never NULL. Returning NULL inside a
#' renderUI makes the card vanish with no trace, which is how a reader concludes a card has
#' nothing to say rather than that its input is missing. The note names the guard so the
#' message can be traced to the code that emitted it.
#' @export
ie_require_measured <- function(x, guard, artifact = NULL, builder = NULL) {
  if (ie_is_measured(x)) return(NULL)
  if (missing(guard) || !nzchar(guard)) {
    stop("ie_require_measured(): every guard needs a name", call. = FALSE)
  }
  reason <- if (inherits(x, "ie_avail")) ie_state_text(x) else "no availability recorded"
  parts <- c(reason)
  if (!is.null(artifact)) parts <- c(parts, paste0("Artifact: ", artifact))
  if (!is.null(builder)) parts <- c(parts, paste0("Rebuild with: ", builder))
  structure(list(guard = guard, lines = parts), class = "ie_guard_note")
}

#' @export
format.ie_guard_note <- function(x, ...) {
  paste0(paste(x$lines, collapse = " "), "  [", x$guard, "]")
}

#' @export
print.ie_guard_note <- function(x, ...) {
  cat(format(x), "\n"); invisible(x)
}
