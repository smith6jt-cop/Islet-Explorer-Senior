#' Overview module — landing page summarising the panel audit, with a
#' tissue picker in the sidebar.
#'
#' Renders the active tissue's `subpanel_summary.csv` as a `DT::datatable`
#' and a `plotly` bar chart of `n_genes`. Selecting a row sets
#' `app_state$selected_subpanel` and `app_state$nav_target = "panel_browser"`,
#' which the app-level server watches to switch tabs and pre-select the panel.
#'
#' Switching the tissue dropdown updates `app_state$selected_tissue`, which
#' the top-level `panels` reactive depends on.

overview_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::h5("Tissue"),
      shiny::selectInput(ns("tissue"), label = NULL,
                         choices = NULL, selectize = FALSE),
      shiny::p(shiny::em("Switching tissue reloads subpanels, audit, and ",
                         "the optional custom panel. The currently-loaded ",
                         "Xenium dataset is preserved.")),
      shiny::hr(),
      shiny::h5("Audit overview"),
      shiny::uiOutput(ns("summary_source")),
      shiny::p("Click a row to open the panel browser pre-filtered to ",
               "that subpanel."),
      shiny::uiOutput(ns("counts"))
    ),
    shiny::uiOutput(ns("env_card")),
    bslib::card(
      bslib::card_header("Genes per subpanel"),
      plotly::plotlyOutput(ns("nbar"), height = "320px")
    ),
    bslib::card(
      bslib::card_header("Subpanels"),
      DT::DTOutput(ns("summary"))
    )
  )
}

overview_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    # Populate the tissue selector once at startup.
    shiny::observe({
      tids <- available_tissues()
      if (!length(tids)) return()
      current <- shiny::isolate(app_state$selected_tissue)
      if (is.null(current) || !(current %in% tids)) current <- tids[1]
      shiny::updateSelectInput(session, "tissue",
                               choices  = tids,
                               selected = current)
    }, priority = 100)

    # Push user picks into the app-wide reactive (drives the panels reactive).
    shiny::observeEvent(input$tissue, {
      sel <- input$tissue
      if (is.null(sel) || !nzchar(sel)) return()
      if (!identical(sel, app_state$selected_tissue)) {
        app_state$selected_tissue <- sel
      }
    }, ignoreInit = TRUE)

    output$env_card <- shiny::renderUI({
      errs <- app_state$env_errors   %||% character()
      warn <- app_state$env_warnings %||% character()
      if (!length(errs) && !length(warn)) return(NULL)
      bullets <- function(items) shiny::tags$ul(
        lapply(items, function(s) shiny::tags$li(s)))
      cls <- if (length(errs)) "alert alert-danger" else "alert alert-warning"
      shiny::div(class = cls,
        if (length(errs)) shiny::tagList(
          shiny::strong("Environment errors"), bullets(errs)),
        if (length(warn)) shiny::tagList(
          shiny::strong("Environment warnings"), bullets(warn)))
    })

    output$summary_source <- shiny::renderUI({
      p <- panels()
      shiny::p("Source: ",
               shiny::code(file.path("data", "tissues", p$tissue$id,
                                     "subpanel_summary.csv")))
    })

    output$counts <- shiny::renderUI({
      p <- panels()
      custom_lbl <- custom_panel_label(app_state$custom_panel_status, p)
      tagList <- shiny::tagList(
        shiny::p(shiny::strong("Tissue: "), p$tissue$display_name),
        shiny::p(shiny::strong("Subpanels indexed: "),
                 nrow(p$summary)),
        shiny::p(shiny::strong("Subpanel CSVs loaded: "),
                 length(p$subpanels)),
        shiny::p(shiny::strong("5K reference genes: "),
                 nrow(p$reference_5k))
      )
      if (!is.null(p$custom) && nrow(p$custom) > 0L) {
        tagList <- shiny::tagAppendChild(tagList,
          shiny::p(shiny::strong(sprintf("%s genes: ", custom_lbl)),
                   nrow(p$custom)))
      } else {
        tagList <- shiny::tagAppendChild(tagList,
          shiny::p(shiny::em("No custom panel for this tissue.")))
      }
      tagList
    })

    output$summary <- DT::renderDT({
      df <- panels()$summary
      DT::datatable(
        df,
        rownames  = FALSE,
        selection = list(mode = "single", selected = NULL),
        options   = list(pageLength = 25,
                         order = if (nrow(df) > 0) list(list(0, "asc"))
                                 else list()),
        class     = "stripe hover compact"
      )
    })

    output$nbar <- plotly::renderPlotly({
      df <- panels()$summary
      if (!nrow(df)) {
        return(plotly::plotly_empty(type = "bar") |>
          plotly::layout(annotations = list(list(
            text = "No subpanels indexed for this tissue.",
            showarrow = FALSE, x = 0.5, y = 0.5,
            xref = "paper", yref = "paper"))))
      }
      df <- df[order(df$n_genes, decreasing = TRUE), , drop = FALSE]
      df$subpanel <- factor(df$subpanel, levels = df$subpanel)
      if (!"description" %in% names(df)) df$description <- ""
      df$description[is.na(df$description)] <- ""
      plotly::plot_ly(
        df,
        x      = ~n_genes,
        y      = ~subpanel,
        type   = "bar",
        orientation = "h",
        hoverinfo   = "text",
        text   = ~paste0(subpanel, " — ", n_genes, " genes",
                         ifelse(nzchar(description),
                                paste0("<br>", description), ""))
      ) |>
        plotly::layout(
          margin = list(l = 220),
          xaxis  = list(title = "n_genes"),
          yaxis  = list(title = "", autorange = "reversed"),
          showlegend = FALSE
        )
    })

    shiny::observeEvent(input$summary_rows_selected, {
      idx <- input$summary_rows_selected
      if (length(idx) != 1L) return()
      sp <- panels()$summary$subpanel[idx]
      key <- resolve_subpanel_key(sp, panels())
      app_state$selected_subpanel <- if (is.na(key)) sp else key
      app_state$nav_target        <- "Panel Browser"
    }, ignoreInit = TRUE)
  })
}
