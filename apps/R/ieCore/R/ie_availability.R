# Six states, because "missing" is five different things and collapsing them lies.
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
#   uncalled     typed, and the typing resolved no identity. NOT an absence of measurement and
#                NOT an absence of typing -- the producer looked, and recorded that it found
#                nothing. Added 2026-09-18 for the refinement pass's
#                `Unresolved (no lineage evidence)` and the Xenium releases' `Uncalled_*`
#                family, both of which this function previously read as `measured`
#
# THE TRAP, measured: in `lineage_common` the QC-dropped cells and the never-typed donors are
# BYTE-IDENTICAL empty strings. Deriving identity state from that column silently converts a
# measured QC drop into an absence of typing. `lineage_native` is the only column that keeps
# them apart, so it is the only one this file reads.

#' The availability states, in the order a legend walks them.
#'
#' `uncalled` is APPENDED, not inserted: these are factor levels and a card legend enumerates
#' them in order, so no existing level may move.
#' @export
IE_STATES <- c("measured", "blank", "unavailable", "withheld", "failed", "uncalled")

#' Protein-side labels that are the producer's own recorded non-call.
#'
#' The Xenium releases name their non-calls by a prefix convention that
#' `scripts/senior/xenium_release.py::classify_label` already reads (`Unavailable*` ->
#' unavailable, `Uncalled*` / `*_uncalled` -> uncalled). The protein refinement pass does NOT
#' follow that convention -- its label is prose -- so the one protein spelling is enumerated
#' here rather than matched by shape. A label that is not on this list and does not carry one
#' of the two Xenium prefixes is a CALL.
#' @export
IE_PROTEIN_NONCALL_LABELS <- c("Unresolved (no lineage evidence)")

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
#'
#' THE RULE IS NOT INVENTED HERE. `scripts/senior/xenium_release.py::classify_label` is this
#' repository's authority on how a release label is read, and it is copied rather than
#' paraphrased: `Unavailable*` is the release saying the answer is unavailable; `Uncalled*` and
#' `*_uncalled` are "looked at, nothing resolved"; everything else is a call, INCLUDING a
#' measured negative such as `Non_immune`, which is an answer and must never be collapsed into
#' an absence.
#'
#' WHAT THIS FIXES. The equality test it replaces matched the exact string `"Unavailable"` and
#' nothing else, so `Unavailable_no_nucleus` -- the label the pancreas and pLN releases actually
#' write, and the one cell this repository excludes from every calculation -- read as
#' `measured`. So did all four of `Uncalled_no_naming_evidence`, `Uncalled_no_enrichment`,
#' `Uncalled_below_cluster_qc` and `Uncalled_ambiguous`, and so would the protein refinement's
#' `Unresolved (no lineage evidence)`. Each of those is a producer's honest absence being
#' counted as an identity.
#' @export
ie_identity_state <- function(lineage_native) {
  if (is.null(lineage_native)) return(factor(character(0), levels = IE_STATES))
  v <- as.character(lineage_native)
  blank <- is.na(v) | !nzchar(v)
  # `FALSE & NA` is FALSE in R, so a missing value cannot leak out of these guards.
  unavail <- !blank & startsWith(v, "Unavailable")
  uncalled <- !blank & !unavail &
    (startsWith(v, "Uncalled") | endsWith(v, "_uncalled") | v %in% IE_PROTEIN_NONCALL_LABELS)
  out <- ifelse(blank, "blank",
         ifelse(unavail, "unavailable",
         ifelse(uncalled, "uncalled", "measured")))
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
    # Not "Unresolved" and not "Unknown": this names the STATE OF THE CALL, and the phrase is
    # `xenium_release.classify_label`'s own ("looked at, no call resolved").
    uncalled    = "Typed, no identity resolved",
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

#' Colours for the six states
#'
#' Deliberately outside the DOMAIN phenotype palette: an absence must not be able to collide
#' with a cell type on a shared legend.
#'
#' `uncalled` is the one exception, and it is an exception on purpose -- see
#' `ie_absence_class_color()`. A DOMAIN whose vocabulary contains a class that IS an absence
#' takes its colour from here, so the absence looks the same wherever it is drawn.
#' @export
ie_state_palette <- function() {
  c(measured    = "#2B7A0B",
    blank       = "#9AA0A6",
    unavailable = "#495057",
    withheld    = "#7A6A00",
    failed      = "#A61B1B",
    # ACHROMATIC on purpose: R == G == B, no hue at all. Every cell-type colour in every DOMAIN
    # palette carries a hue, so "no hue" is the visual statement that there is no identity here,
    # and it survives being read next to a saturated legend. Mid-grey rather than near-white so
    # the cells are still visible: a class drawn in a colour nobody can see is a worse lie than
    # one drawn in the wrong colour. 5.1:1 on white, 4.1:1 on black. Distinct from `blank`
    # (#9AA0A6) and `unavailable` (#495057), which are both lighter/darker AND cool-cast.
    uncalled    = "#6E6E6E")
}

#' The colour a DOMAIN phenotype takes when the phenotype IS an absence
#'
#' One authority, so a config literal and this palette cannot drift. The YAML cannot call R, so
#' the domain file carries the hex and the testthat suite asserts it against this.
#' @export
ie_absence_class_color <- function() unname(ie_state_palette()[["uncalled"]])

#' Which of a DOMAIN's display classes record an absence rather than a cell type
#'
#' Takes the domain's declared `phenotypes.residuals` list (see `domain_phenotype_residuals()`)
#' and returns the names of the entries that declare `absence: true`. Keeping the question here
#' rather than in a module means a tab asks "is this an absence" instead of matching a label by
#' spelling, and a second domain that grows one only has to declare it.
#'
#' A class is an absence only where the DOMAIN says so. Nothing is inferred from the name: a
#' phenotype called `Unknown` in the legacy CODEX vocabulary is a real upstream residual class,
#' not a recorded non-call, and guessing from the string would put a measurement and a hole in
#' the same bucket.
#' @export
ie_absence_phenotypes <- function(residuals) {
  if (is.null(residuals) || !length(residuals)) return(character(0))
  out <- character(0)
  for (entry in residuals) {
    if (!is.list(entry) || !isTRUE(entry[["absence"]])) next
    name <- entry[["name"]]
    if (is.null(name) || length(name) != 1L || is.na(name) || !nzchar(as.character(name))) next
    out <- c(out, as.character(name))
  }
  unique(out)
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
