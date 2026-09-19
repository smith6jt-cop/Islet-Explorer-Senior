# A selector offers only what it can actually show.
#
# The maintainer's rule, 2026-09-18: "Only donors that have data ready should be available to choose
# in the app ... The same goes for similar situations elsewhere."
#
# The post-Xenium IF overlay was offering all 25 donors of its availability layer and then rendering
# a named note for 23 of them, because the choices came from every row of the artifact rather than
# from the rows that are measured. An unavailable state is the right thing to RECORD; it is the
# wrong thing to put in a dropdown. A reader should not have to click through 23 donors to find the
# two with data, and a selector that offers a donor is a claim that the donor has something to show.
#
# This is not the same as hiding the absence. The count of what was withheld, and why, still belongs
# on the card -- see `ie_selectable_note()`. What changes is that it is stated once, as a fact about
# the cohort, instead of 23 times as a fact about whatever the reader happened to click.

#' The donors a card can actually draw, from its own availability table.
#'
#' @param data a data.frame with a donor column and a state column.
#' @param state_col the column carrying the availability value, one of [IE_STATES]. Defaults to
#'   `"state"`. (Five states until 2026-09-18, six since: `uncalled` was appended.)
#' @param donor_col the donor identifier column. Defaults to `"donor_id"`.
#' @param states which states count as showable. Defaults to `"measured"` alone -- `blank` means the
#'   lane ran and found nothing, which is a result worth recording and not something to plot.
#' @param cohort optional cohort object; when it exposes `donor_ids()`, the result is intersected
#'   with the active domain's roster, so a donor present in a shared artifact but outside this
#'   domain is not offered. The two domains share five donors and no vocabulary.
#' @return a character vector of donor ids, sorted, possibly empty.
#' @export
ie_selectable_donors <- function(data, state_col = "state", donor_col = "donor_id",
                                 states = "measured", cohort = NULL) {
  if (is.null(data) || !nrow(data)) return(character())
  if (!all(c(state_col, donor_col) %in% names(data))) {
    stop(sprintf("ie_selectable_donors(): need columns %s and %s", donor_col, state_col))
  }
  keep <- as.character(data[[state_col]]) %in% states
  ids <- unique(as.character(data[[donor_col]][keep]))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!is.null(cohort) && is.function(cohort$donor_ids)) {
    ids <- intersect(ids, as.character(cohort$donor_ids()))
  }
  sort(ids)
}


#' One line saying what was left out of a selector, and why.
#'
#' Rendered once beside the control. Without it, narrowing the choices would hide the absence
#' rather than state it -- which is the opposite of what the [IE_STATES] vocabulary is for.
#'
#' @param data the same table passed to `ie_selectable_donors()`.
#' @param offered the ids that survived.
#' @inheritParams ie_selectable_donors
#' @return a `shiny::tags$p`, or `NULL` when nothing was withheld.
#' @export
ie_selectable_note <- function(data, offered, state_col = "state", donor_col = "donor_id") {
  if (is.null(data) || !nrow(data)) return(NULL)
  ids <- unique(as.character(data[[donor_col]]))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  missing <- setdiff(ids, offered)
  if (!length(missing)) return(NULL)
  states <- as.character(data[[state_col]][match(missing, as.character(data[[donor_col]]))])
  tally <- sort(table(states[!is.na(states) & nzchar(states)]), decreasing = TRUE)
  shiny::tags$p(
    class = "ie-card-note", style = "font-size: 12px",
    sprintf("%d of %d donors are not offered: %s. They are recorded, not missing.",
            length(missing), length(ids),
            paste(sprintf("%d %s", as.integer(tally), names(tally)), collapse = ", "))
  )
}
