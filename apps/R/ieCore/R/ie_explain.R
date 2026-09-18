# Explanatory text, collapsed by default behind a descriptive toggle.
#
# The app had 4,665 words of standing explanation across 53 blocks in 12 files, printed above the
# figures on every render -- Statistics alone carried 2,329. None of it was wrong and none of it is
# dropped here. What was wrong is that a caveat which is always true was given the same weight, and
# the same vertical space, as a number that changes.
#
# The rule this encodes, from the maintainer 2026-09-18: if a large block of explanatory text is
# needed, it is HIDDEN BY DEFAULT behind a descriptive icon that shows and hides it.
#
# `<details>`/`<summary>` rather than a modal or a `bslib` accordion, deliberately:
#   * no server round-trip and no `observeEvent`, so a card gains no state and no namespaced id.
#     The modal this replaces needed one observer per card and was already a coordination cost.
#   * it is a native disclosure widget, so it is keyboard reachable and screen readers announce the
#     expanded/collapsed state without an `aria-expanded` we would have to keep in sync ourselves.
#   * it survives a `renderUI` rebuild, which a modal does not.
#
# MEASURED STATE DOES NOT GO IN HERE. A count, a donor list, a withheld reason, an unavailable
# state -- anything that is a fact about what the reader is looking at right now -- stays visible.
# This is for text that would read identically on every donor and every render.

#' A collapsed block of explanatory text with a descriptive toggle.
#'
#' @param ... the explanation. Any shiny tags; passed through untouched.
#' @param label short, specific summary text shown beside the icon. "How this is measured" beats
#'   "More info": the reader decides whether to open it from this string alone.
#' @param icon_char the disclosure glyph. Defaults to an information mark. The caret is supplied by
#'   the browser's own `<summary>` marker, so this is the subject, not the direction.
#' @param open `TRUE` to start expanded. Use only where a reader who skips the text would
#'   misread the figure -- which is rare, and is not the same as the text being important.
#' @return a `shiny::tags$details` element.
#' @export
ie_explain <- function(..., label = "How this is measured", icon_char = "ⓘ", open = FALSE) {
  if (!nzchar(label)) stop("ie_explain(): label must say what is inside, not 'More info'")
  summary_tag <- shiny::tags$summary(
    class = "ie-explain-summary",
    style = paste("cursor: pointer; font-size: 12px; color: #5a6270;",
                  "padding: 2px 0; list-style-position: outside"),
    title = paste0(label, " (click to show or hide)"),
    shiny::tags$span(`aria-hidden` = "true", style = "margin-right: 4px", icon_char),
    label
  )
  shiny::tags$details(
    class = "ie-explain",
    style = "margin: 4px 0 8px 0",
    if (isTRUE(open)) shiny::tagAppendAttributes(summary_tag, open = NA) else summary_tag,
    shiny::tags$div(
      class = "ie-explain-body",
      style = paste("font-size: 12px; line-height: 1.45; color: #4a5160;",
                    "padding: 6px 10px; margin-top: 4px;",
                    "border-left: 2px solid #d8dce3"),
      ...
    )
  ) -> el
  if (isTRUE(open)) shiny::tagAppendAttributes(el, open = NA) else el
}


#' How many words of ALWAYS-VISIBLE prose a UI object carries.
#'
#' Counts text outside any `ie_explain()` block, which is what a reader sees before opening
#' anything. Used by the guard test so the app cannot quietly regrow its walls of text.
#'
#' @param tag a shiny tag or tagList.
#' @return integer word count of the visible text.
#' @export
ie_visible_words <- function(tag) {
  drop_explained <- function(x) {
    if (inherits(x, "shiny.tag")) {
      if (identical(x$name, "details") &&
          isTRUE(grepl("ie-explain", paste(x$attribs$class, collapse = " ")))) {
        return(NULL)
      }
      x$children <- lapply(x$children, drop_explained)
      return(x)
    }
    if (is.list(x)) return(lapply(x, drop_explained))
    x
  }
  visible <- drop_explained(tag)
  txt <- paste(as.character(visible), collapse = " ")
  txt <- gsub("<[^>]*>", " ", txt)          # tags out; their text stays
  txt <- gsub("&[a-z]+;|&#[0-9]+;", " ", txt)
  length(strsplit(trimws(gsub("\\s+", " ", txt)), " ")[[1]])
}
