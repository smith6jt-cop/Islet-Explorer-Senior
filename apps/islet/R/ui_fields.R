# ---------------------------------------------------------------------------
# Modernist form primitives
#
# Extracted from design/delivered/shiny/R/mod_cohort.R. Only the two UI
# primitives are taken here; `cohortFilterUI` / `cohortFilterServer` /
# `normalize_features` are NOT ported yet -- the shared cohort bar rewires
# raw_df_base(), stats_run() and nbr_comp() plus four sync-observer triads, and
# bundling that into the shell port would make the shell diff unreviewable. It
# lands as its own phased migration. The file is named for what it contains so
# nothing implies a cohort module already exists.
#
# GROUP_LEVELS lives in R/theme_modernist.R (it is a data-level constant, used
# by the scales), not here.
# ---------------------------------------------------------------------------

#' Label + control cell. Used by the cohort bar and per-panel control rows.
ie_field <- function(label, ..., width = NULL) {
  div(
    class = "ie-field",
    style = if (!is.null(width)) paste0("width:", width),
    tags$span(class = "ie-field-label", label),
    ...
  )
}

#' Hand-rolled segmented control.
#'
#' Shiny's RadioInputBinding binds any `.shiny-input-radiogroup` container and
#' reads `input:radio[name=<container id>]:checked`, so this gets full markup
#' control and a pure-CSS active state with no extra package.
#'
#' *** SAFE-USAGE RULE -- verified against the shiny 1.13.0 source, not assumed ***
#'   updateRadioButtons(session, id, selected = <value>)   SAFE
#'   updateRadioButtons(session, id, choices  = <...>)     BROKEN
#'
#' `RadioInputBinding.receiveMessage` replaces options by removing
#' `div.shiny-options-group` and `label.radio` and appending the new set. This
#' markup contains NEITHER, so nothing is removed and a second, Bootstrap-styled
#' option block is APPENDED. Worse, the appended `<input checked>` shares the
#' group `name`, so the browser unchecks the hand-rolled one -- the control then
#' renders visually EMPTY while `input$<id>` still returns a sensible value, so
#' it reads as a CSS bug and gets misdiagnosed. `label = ` is a silent no-op for
#' the same reason (getLabelNode looks for `label[for=<id>]`, which is absent).
#'
#' To change a segmented control's OPTIONS, re-render it via uiOutput/renderUI.
#' Real radioButtons() and selectInput() are unaffected.
ie_segmented <- function(inputId, choices, selected = NULL) {
  values <- if (is.null(names(choices))) unlist(choices, use.names = FALSE)
            else unname(unlist(choices))
  labels <- if (is.null(names(choices))) values else names(choices)
  if (is.null(selected)) selected <- values[[1]]
  div(
    id = inputId, class = "shiny-input-radiogroup ie-seg",
    role = "radiogroup", `aria-label` = inputId,
    Map(function(v, l) {
      tags$label(
        class = "ie-seg-opt",
        tags$input(
          type = "radio", name = inputId, value = v,
          checked = if (identical(as.character(v), as.character(selected))) "checked"
        ),
        # <input> then <span> so `input:checked + span` styles the active
        # option -- see the :has() note in R/ui_rail.R.
        tags$span(l)
      )
    }, values, labels)
  )
}
