#' Panel-vs-Data Compare module — coverage and detection across subpanels.
#'
#' For every subpanel (plus custom-100 and the full xenium5k_in_audit)
#' compute fraction present, fraction passing the user-set detection cutoff,
#' and detection-percentile / expression summaries on the loaded dataset.
#' Renders a coverage `DT`, a `plotly` bar of fraction passing, and a
#' `geom_tile` heatmap of the top-N genes per subpanel.

panel_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::sliderInput(ns("min_det"),
                         "Min detection % (loaded data)",
                         min = 0, max = 100, value = 0, step = 0.5),
      shiny::sliderInput(ns("topn"),
                         "Top-N genes per subpanel (heatmap)",
                         min = 3, max = 30, value = 10, step = 1),
      shiny::checkboxInput(ns("hide_5k"),
                           "Hide xenium5k_in_audit row in summary",
                           value = TRUE),
      shiny::hr(),
      shiny::p(shiny::em("Filter state is preserved in app_state ",
                         "so other modules can read it."))
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        bslib::card_header("Coverage by subpanel"),
        shiny::uiOutput(ns("nodata_msg")),
        DT::DTOutput(ns("coverage"))
      ),
      bslib::card(
        bslib::card_header("Fraction of subpanel genes passing the cutoff"),
        plotly::plotlyOutput(ns("frac_bar"), height = "420px")
      )
    ),
    bslib::card(
      bslib::card_header("Top-N genes per subpanel — detection %"),
      shiny::p("Rows are subpanels; columns are the top-N genes by ",
               "detection % in the loaded dataset. White cells = gene ",
               "below cutoff or absent from data."),
      shiny::plotOutput(ns("heatmap"), height = "1100px")
    )
  )
}

panel_compare_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    # Persist filter state to app_state so M5/M8 can pick it up.
    shiny::observe({
      app_state$compare_min_det <- input$min_det
      app_state$compare_topn    <- input$topn
    })

    have_xen <- shiny::reactive(!is.null(app_state$xen))

    output$nodata_msg <- shiny::renderUI({
      if (have_xen()) return(NULL)
      shiny::div(class = "alert alert-info",
                 "Load a Xenium dataset on the ",
                 shiny::strong("Load Xenium"),
                 " tab to populate this comparison.")
    })

    coverage_df <- shiny::reactive({
      shiny::req(have_xen())
      compute_subpanel_coverage(
        panels(), app_state$xen,
        min_detection_pct = input$min_det,
        custom_label      = custom_panel_label(app_state$custom_panel_status,
                                               panels())
      )
    })

    output$coverage <- DT::renderDT({
      df <- coverage_df()
      if (isTRUE(input$hide_5k)) {
        df <- df[df$subpanel != "xenium5k_in_audit", , drop = FALSE]
      }
      # round numeric cols for display
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(x) signif(x, 4))
      DT::datatable(
        df, rownames = FALSE, filter = "top",
        options = list(pageLength = 25, scrollX = TRUE,
                       order = list(list(5, "desc"))),
        class = "stripe hover compact nowrap"
      ) |>
        DT::formatPercentage(c("frac_present", "frac_passing"), 1)
    })

    output$frac_bar <- plotly::renderPlotly({
      df <- coverage_df()
      df <- df[df$subpanel != "xenium5k_in_audit", , drop = FALSE]
      df <- df[order(df$frac_passing, decreasing = TRUE), , drop = FALSE]
      df$subpanel <- factor(df$subpanel, levels = df$subpanel)
      plotly::plot_ly(
        df,
        x    = ~frac_passing,
        y    = ~subpanel,
        type = "bar",
        orientation = "h",
        text = ~sprintf("%d / %d (%.1f%%) passing",
                        n_passing, n_genes, 100 * frac_passing),
        hoverinfo = "text"
      ) |>
        plotly::layout(
          margin = list(l = 220),
          xaxis = list(title = "fraction passing", tickformat = ".0%"),
          yaxis = list(title = "", autorange = "reversed"),
          showlegend = FALSE
        )
    })

    output$heatmap <- shiny::renderPlot({
      shiny::req(have_xen())
      df <- top_n_by_subpanel(panels(), app_state$xen,
                              n = input$topn,
                              min_detection_pct = input$min_det)
      shiny::validate(shiny::need(
        nrow(df) > 0,
        "No genes pass the cutoff in any subpanel."
      ))

      # Stable column ordering: rank within subpanel (1..N).
      df <- df[order(df$subpanel, -df$detection_pct), ]
      df$col <- stats::ave(df$detection_pct, df$subpanel,
                           FUN = function(v) seq_along(v))

      ggplot2::ggplot(df,
                      ggplot2::aes(x = col, y = subpanel,
                                   fill = detection_pct)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
        ggplot2::geom_text(ggplot2::aes(label = gene),
                           size = 4, colour = "black") +
        ggplot2::scale_fill_viridis_c(option = "C",
                                      name = "detection %",
                                      limits = c(0, 100)) +
        ggplot2::scale_x_continuous(breaks = NULL) +
        ggplot2::scale_y_discrete(limits = rev) +
        ggplot2::labs(x = sprintf("rank within subpanel (top-%d)",
                                  input$topn),
                      y = NULL) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
    })
  })
}
