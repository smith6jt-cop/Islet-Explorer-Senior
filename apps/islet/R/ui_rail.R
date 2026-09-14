# ---------------------------------------------------------------------------
# Left rail + hidden tabset
#
# Ported from Junior's f827106 (which ported it from design/delivered/shiny)
# with three Senior changes, all noted in place:
#
#   1. The `%||%` redefinition is dropped, as it was in Junior: Shiny sources
#      every R/*.R into the GLOBAL environment, which sits ahead of
#      package:base, so a copy here would shadow base's (>= 4.4) for the whole
#      app -- and mod_statistics_server.R depends on `%||%` NOT catching NA.
#   2. rail_row() wraps its visible content in a span, so the active state is a
#      plain sibling combinator instead of `:has()`. `:has()` needs Firefox
#      121+; on ESR 115 the rail would render with NO active-tab indicator at
#      all while still working, and nothing would log.
#   3. SENIOR: the rail's VALUES are the existing tab labels ("Plot",
#      "Trajectory", ...), not new lowercase keys. Junior renamed them and then
#      had to chase every `active_tab() != "Trajectory"`, every
#      `condition = "input.tabs == 'Plot'"` and viewer_server's own comparison.
#      Keeping the values means the navset keeps the id `tabs` and NOTHING
#      downstream of it changes. The rail is a new way to reach the same tabs.
#
# The rail is a hand-rolled Shiny radio group: one <input type="radio"> per
# destination inside a .shiny-input-radiogroup container. Shiny's
# RadioInputBinding finds any such container and reads
# `input:radio[name=<container id>]:checked` -- a document-wide name lookup --
# so this binds with no extra package and no custom input binding. That buys:
#   - arbitrary markup per row (numeral, label, subtitle, section headings)
#   - active state in pure CSS, so switching tabs costs no re-render, no custom
#     message handler and no actionButton counter bookkeeping
#
#   railUI("rail", IE_PIPELINE, IE_REFERENCE)
#   railServer("rail", "tabs")   # keeps navset_hidden in sync
# ---------------------------------------------------------------------------

rail_row <- function(group, key, selected = FALSE, num = NULL,
                     label = NULL, sub = NULL, marker = NULL) {
  tags$label(
    class = "ie-nav",
    tags$input(
      type = "radio", name = group, value = key,
      checked = if (selected) "checked"
    ),
    # WRAPPER SPAN, not in the delivered file -- see note 2 in the header.
    tags$span(
      class = "ie-nav-body",
      tags$span(class = "ie-nav-num", num %||% marker),
      tags$span(
        class = "ie-nav-text",
        tags$span(label),
        if (!is.null(sub)) tags$span(class = "ie-nav-sub", sub)
      )
    )
  )
}

rail_heading <- function(text) div(class = "ie-rail-heading", text)

#' @param pipeline,reference lists of `list(key, num/marker, label, sub)`.
#'   `key` is the value the navset receives, so it must match a
#'   `nav_panel_hidden()` value exactly.
railUI <- function(id, pipeline, reference = list(),
                   reference_heading = "The images themselves",
                   brand = "Islet Explorer", logo = "logo.png",
                   meta = NULL, help = NULL, footer = NULL) {
  ns <- NS(id)
  first <- pipeline[[1]]$key

  tags$nav(
    class = "ie-rail",
    div(
      class = "ie-rail-brand",
      # Logo only. The mark already carries the wordmark, so a text line under
      # it read as a duplicate; `brand` survives as the alt text.
      div(
        class = "ie-brandmark",
        if (!is.null(logo))
          tags$img(class = "ie-brandlogo", src = logo, alt = paste(brand, "logo"))
        else div(class = "ie-brandtext", brand)
      ),
      if (!is.null(meta)) div(class = "ie-rail-meta", HTML(meta))
    ),

    # A single radio group spans BOTH sections, so exactly one row is ever lit
    # -- the Viewer is a destination like any other, it is only grouped apart.
    div(
      id = ns("tab"), class = "shiny-input-radiogroup",
      role = "radiogroup", `aria-label` = "Analysis pipeline",

      rail_heading("Analysis pipeline"),
      lapply(pipeline, function(it) {
        rail_row(ns("tab"), it$key, identical(it$key, first),
                 num = it$num, label = it$label, sub = it$sub)
      }),

      if (length(reference)) tagList(
        div(class = "ie-rail-divider"),
        rail_heading(reference_heading),
        lapply(reference, function(it) {
          rail_row(ns("tab"), it$key, FALSE,
                   marker = "◇", label = it$label, sub = it$sub)
        })
      )
    ),

    # Help lives here rather than at the top of each tab. Senior mounts the
    # phenotype-rules pill on five tabs; one place, always reachable, costs no
    # vertical space in any panel.
    if (!is.null(help)) tagList(
      div(class = "ie-rail-divider"),
      rail_heading("Help"),
      div(class = "ie-rail-help", help)
    ),

    if (!is.null(footer)) div(class = "ie-rail-footer", footer)
  )
}

#' A help entry in the rail. An actionLink, so the server opens a modal.
rail_help_link <- function(inputId, label, sub = NULL) {
  actionLink(inputId, class = "ie-rail-helplink",
    label = tagList(
      tags$span(class = "ie-nav-num", "?"),
      tags$span(class = "ie-nav-text",
        tags$span(label),
        if (!is.null(sub)) tags$span(class = "ie-nav-sub", sub))
    ))
}

#' A labelled control in the rail footer.
rail_palette_select <- function(inputId, label, choices, selected) {
  div(
    class = "ie-rail-footer-field",
    tags$label(class = "ie-rail-footer-label", `for` = inputId, label),
    selectInput(inputId, NULL, choices = choices, selected = selected,
                width = "100%")
  )
}

#' Static label/value pair for the rail footer (kept for non-interactive notes).
rail_palette_note <- function(label, value) {
  div(
    div(class = "ie-rail-footer-label", label),
    div(class = "ie-rail-footer-value", value)
  )
}

#' @param tabset_id outer id of the bslib::navset_hidden holding the panels
#' @return reactive giving the active key
railServer <- function(id, tabset_id) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$tab, {
      # rootScope() is REQUIRED: without it nav_select resolves `tabset_id`
      # inside this module's namespace and silently does nothing.
      bslib::nav_select(tabset_id, selected = input$tab,
                        session = session$rootScope())
    }, ignoreInit = TRUE)
    reactive(input$tab)
  })
}
