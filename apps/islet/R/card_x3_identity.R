# Hormone identity is separate from measured IF intensity. Never turn withheld
# calls into negative cells, derive thresholds, or label a hormone as a cell type.

card_x3_identity_ui <- function(ns) {
  shiny::tags$div(class = "ie-card", style = "overflow: visible;",
    shiny::tags$h3(class = "ie-card-title", "IF hormone identity"),
    shiny::selectInput(ns("x3_identity_donor"), "Donor", choices = character()),
    shiny::uiOutput(ns("x3_identity_content")))
}

ie_restain_identity <- function(row) {
  do.call(rbind, lapply(c("INS", "GCG"), function(marker) {
    prefix <- paste0("hormone_positivity_", marker, "_")
    fields <- paste0(prefix, c("state", "reason", "auroc", "bar"))
    if (!all(fields %in% names(row))) stop("missing hormone availability fields")
    availability <- ie_avail(row[[fields[[1]]]][[1]], reason = row[[fields[[2]]]][[1]])
    # Missing validation stays an explicit unavailable state, never a numeric zero.
    number <- function(x) if (is.na(x)) "Unavailable" else formatC(x, digits = 4, format = "f")
    data.frame(Hormone = marker, State = ie_state_text(availability),
      AUROC = number(row[[fields[[3]]]][[1]]),
      `Prerecorded AUROC threshold` = number(row[[fields[[4]]]][[1]]), check.names = FALSE)
  }))
}

card_x3_identity_server <- function(input, output, session, prepared = NULL, cohort = NULL,
                                    active_tab = NULL, tab_key = NULL, ...) {
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                             identical(active_tab(), tab_key))
  selection <- ie_restain_reactives(input, output, session, "x3_identity", cohort, active)
  identity <- shiny::reactive({
    selected <- selection()
    if (!ie_is_measured(selected$availability)) return(selected)
    tryCatch(list(availability = selected$availability, data = ie_restain_identity(selected$data)),
      error = function(e) ie_restain_failure(paste0("invalid_hormone_artifact: ",
                                                   conditionMessage(e)), "failed"))
  })
  output$x3_identity_content <- shiny::renderUI({
    d <- identity()
    if (!ie_is_measured(d$availability)) return(ie_restain_note(d$availability, "x3_identity"))
    shiny::tagList(shiny::tableOutput(session$ns("x3_identity_table")),
      shiny::tags$p(class = "ie-card-note",
        "Measured IF intensity does not establish hormone positivity. Withheld calls stay withheld; unavailable calls stay unavailable. No positivity thresholds are fitted by this card."))
  })
  output$x3_identity_table <- shiny::renderTable({
    d <- identity()
    shiny::validate(shiny::need(ie_is_measured(d$availability),
                               "IF hormone identity unavailable [x3_identity]"))
    d$data
  }, striped = TRUE, rownames = FALSE)
  invisible(list(active = active, selection = selection, identity = identity))
}
