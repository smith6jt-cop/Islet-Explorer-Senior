# ---------- Shared UI blocks (Modernist) ----------
# These replace four copies of `section_heading()` and two of the "About this
# tab" banner, each of which hardcoded the same blue gradient inline. Inline
# `style=` beats any stylesheet without `!important`, so every one of those was
# a hole the design system could not reach through.
#
# The fourth section_heading copy is the one worth naming: it lived INSIDE a
# renderUI in mod_spatial_server.R, with the honest comment "inline -- cannot
# share with UI function directly". It is the copy most likely to be missed,
# and it kept painting blue gradients on the Spatial tab after the other three
# were fixed.

#' Numbered section heading.
#'
#' @param step   short marker, usually a number. Sequence is real here: these
#'   number the order an analysis is read in, not decoration.
#' @param tip    optional native-tooltip text on a "?" affordance.
ie_section_heading <- function(step, title, subtitle = NULL, tip = NULL) {
  div(
    class = "ie-section",
    div(
      class = "ie-section-row",
      span(class = "ie-section-num", step),
      span(class = "ie-section-title", title),
      if (!is.null(tip) && nzchar(tip))
        span(class = "ie-section-tip", "?", title = tip)
    ),
    if (!is.null(subtitle)) p(class = "ie-section-sub", subtitle)
  )
}

#' Callout block. `tone` picks the left rule colour, nothing else.
#'
#' info  -- neutral orientation ("About this tab")
#' note  -- a methodological aside
#' warn  -- something the reader must not miss
ie_callout <- function(..., tone = c("info", "note", "warn")) {
  div(class = paste0("ie-callout ie-callout-", match.arg(tone)), ...)
}

#' Small titled panel, used for the Trajectory legend box and similar.
ie_panel_note <- function(title, ...) {
  div(class = "ie-legend-box",
      if (!is.null(title)) div(class = "ie-legend-title", title), ...)
}
