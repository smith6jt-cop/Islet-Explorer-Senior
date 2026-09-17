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
    # No numeral and no marker glyph. The rail is a list of destinations, not an
    # ordered procedure, and the leading column of digits was the loudest thing
    # in it. `num` and `marker` stay in the signature so a caller passing them
    # is not an error; they are simply not drawn.
    tags$span(
      class = "ie-nav-body",
      tags$span(
        class = "ie-nav-text",
        tags$span(label),
        if (!is.null(sub)) tags$span(class = "ie-nav-sub", sub)
      )
    )
  )
}

# rail_heading() removed with the two section headings it drew.

#' @param pipeline,reference lists of `list(key, label, sub)`.
#'   `key` is the value the navset receives, so it must match a
#'   `nav_panel_hidden()` value exactly.
railUI <- function(id, pipeline, reference = list(),
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

      # No section headings. "Analysis pipeline" and "The images themselves"
      # labelled two groups that the divider already separates.
      lapply(pipeline, function(it) {
        rail_row(ns("tab"), it$key, identical(it$key, first),
                 label = it$label, sub = it$sub)
      }),

      if (length(reference)) tagList(
        div(class = "ie-rail-divider"),
        lapply(reference, function(it) {
          rail_row(ns("tab"), it$key, FALSE, label = it$label, sub = it$sub)
        })
      )
    ),

    # Help behind a single "?" -- a native <details>, so it opens and closes with
    # no JS and no state to keep. It was a labelled section taking four lines of
    # rail for something read once.
    if (!is.null(help)) tags$details(
      class = "ie-rail-helpmenu",
      tags$summary(class = "ie-rail-helpbtn", title = "Help", "?"),
      div(class = "ie-rail-help", help)
    ),

    # Pinned to the BOTTOM of the rail by `margin-top: auto` on .ie-rail-footer,
    # so it sits in the lower left whatever is above it.
    if (!is.null(footer)) div(class = "ie-rail-footer", footer)
  )
}

#' A help entry in the rail. An actionLink, so the server opens a modal.
rail_help_link <- function(inputId, label, sub = NULL) {
  actionLink(inputId, class = "ie-rail-helplink",
    label = tagList(
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
