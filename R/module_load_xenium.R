#' Load Xenium module — directory / .rds / .h5ad chooser + cache via qs2.
#'
#' Three input paths (CLAUDE.md §2): a Xenium output bundle directory,
#' a saved Seurat `.rds` / `.qs2`, or an AnnData `.h5ad` (TODO).
#' After load, renders cell/gene counts, the panel-validation report,
#' and a rasterized spatial preview.

load_xenium_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::h5("Load a Xenium dataset"),
      shinyFiles::shinyDirButton(ns("bundle_dir"),
                                 "Pick Xenium bundle directory",
                                 title  = "Select Xenium output bundle"),
      shiny::br(), shiny::br(),
      shinyFiles::shinyFilesButton(ns("rds_file"),
                                   "Pick saved Seurat (.rds / .qs2)",
                                   title  = "Select saved Seurat object",
                                   multiple = FALSE),
      shiny::br(), shiny::br(),
      shiny::textInput(ns("manual_path"),
                       "Or paste a path",
                       placeholder = "/path/to/xenium_outs"),
      shiny::actionButton(ns("load_manual"), "Load from path",
                          class = "btn-primary"),
      shiny::hr(),
      shiny::h5("Demo mode"),
      shiny::p(shiny::em("No data on hand? Load a small synthetic ",
                         "dataset (1.5k cells × 800 genes sampled from ",
                         "the panel) to try every tab end-to-end.")),
      shiny::actionButton(ns("load_demo"), "Load demo dataset",
                          class = "btn-secondary",
                          icon  = shiny::icon("flask")),
      shiny::hr(),
      shiny::checkboxInput(ns("force_refresh"),
                           "Bypass cache (force re-ingest)",
                           value = FALSE),
      shiny::p(shiny::em("Caches are written to "),
               shiny::code(app_paths$cache),
               shiny::em(" keyed by absolute path + mtime."))
    ),
    bslib::card(
      bslib::card_header("Load status"),
      shiny::uiOutput(ns("status"))
    ),
    bslib::card(
      bslib::card_header("Panel validation"),
      shiny::uiOutput(ns("validation"))
    ),
    bslib::card(
      bslib::card_header("Custom panel coverage in loaded data"),
      shiny::uiOutput(ns("custom_coverage_summary")),
      shiny::div(style = "margin: 6px 0 10px;",
                 shiny::downloadButton(
                   ns("dl_missing"),
                   "Download missing genes (CSV)",
                   class = "btn-sm btn-outline-secondary")),
      DT::DTOutput(ns("missing_table"))
    ),
    bslib::card(
      bslib::card_header("Spatial preview"),
      shiny::p("Rasterized scatter of cell centroids ",
               "(", shiny::code("x_centroid"), " vs ",
               shiny::code("y_centroid"), ")."),
      plotly::plotlyOutput(ns("spatial"), height = "520px")
    )
  )
}

load_xenium_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {
    roots <- c(home = fs::path_home(),
               root = "/",
               cwd  = normalizePath("."))
    shinyFiles::shinyDirChoose(input, "bundle_dir", roots = roots,
                               allowDirCreate = FALSE)
    shinyFiles::shinyFileChoose(input, "rds_file", roots = roots,
                                filetypes = c("rds", "qs", "qs2"))

    chosen_path <- shiny::reactiveVal(NULL)
    # Wrap path in a list with a per-click nonce so re-clicks with the same
    # path still invalidate downstream reactives (plain reactiveVal compares
    # by identity and would no-op on the second click).
    set_chosen <- function(p) {
      chosen_path(list(path = p, nonce = stats::runif(1)))
    }

    shiny::observeEvent(input$bundle_dir, {
      sel <- shinyFiles::parseDirPath(roots, input$bundle_dir)
      if (length(sel) == 1L && nzchar(sel)) set_chosen(as.character(sel))
    })
    shiny::observeEvent(input$rds_file, {
      sel <- shinyFiles::parseFilePaths(roots, input$rds_file)
      if (nrow(sel) == 1L) set_chosen(as.character(sel$datapath[1]))
    })
    shiny::observeEvent(input$load_manual, {
      p <- trimws(input$manual_path %||% "")
      if (nzchar(p)) set_chosen(path.expand(p))
    })
    shiny::observeEvent(input$load_demo, {
      set_chosen(demo_path_sentinel())
    })

    waiter::useWaiter()

    xen_load <- shiny::reactive({
      cp <- chosen_path()
      shiny::req(cp, cp$path)
      p <- cp$path
      is_demo <- identical(p, demo_path_sentinel())
      w <- waiter::Waiter$new(
        html = shiny::tagList(
          waiter::spin_dots(),
          shiny::h4(if (is_demo) "Building demo dataset…"
                    else "Loading Xenium dataset…")
        ),
        color = "rgba(0,0,0,0.55)"
      )
      w$show()
      on.exit(w$hide(), add = TRUE)
      tryCatch(
        if (is_demo) load_demo_xenium(panels())
        else load_xenium(p, refresh = isTRUE(input$force_refresh)),
        error = function(e) structure(list(error = conditionMessage(e)),
                                       class = "load_xenium_error")
      )
    })

    shiny::observeEvent(xen_load(), {
      out <- xen_load()
      if (inherits(out, "load_xenium_error")) {
        app_state$xen <- NULL
        shiny::showNotification(
          paste0("Load failed: ", out$error),
          type = "error", duration = 10)
        return()
      }
      app_state$xen      <- out
      app_state$xen_path <- chosen_path()$path
      shiny::showNotification(
        sprintf("Loaded %d cells × %d genes.", ncol(out), nrow(out)),
        type = "message", duration = 5)
    }, ignoreNULL = FALSE)

    output$status <- shiny::renderUI({
      out <- xen_load()
      if (inherits(out, "load_xenium_error")) {
        return(shiny::div(class = "alert alert-danger",
                          shiny::strong("Load failed: "), out$error))
      }
      shiny::req(out)
      cache_info <- attr(out, "load_xenium_cache")
      cp <- chosen_path()
      is_demo <- identical(cp$path, demo_path_sentinel())
      path_line <- if (is_demo) {
        shiny::p(shiny::strong("Source: "),
                 shiny::span(class = "badge bg-warning text-dark",
                             "Demo (synthetic)"))
      } else {
        shiny::p(shiny::strong("Path: "), shiny::code(cp$path))
      }
      shiny::tagList(
        path_line,
        shiny::p(shiny::strong("Class: "), class(out)[1]),
        shiny::p(shiny::strong("Cells: "),  ncol(out)),
        shiny::p(shiny::strong("Genes: "),  nrow(out)),
        shiny::p(shiny::strong("Cache: "),
                 if (isTRUE(cache_info$hit)) "hit (read from qs2)"
                 else paste0("miss (wrote ", basename(cache_info$file %||% ""), ")"))
      )
    })

    output$validation <- shiny::renderUI({
      out <- xen_load()
      if (inherits(out, "load_xenium_error") || is.null(out)) {
        return(shiny::p(shiny::em("Load a dataset to validate.")))
      }
      p <- panels()
      rep <- panel_validate(rownames(out), p)
      cls <- if (rep$pct_reference_covered >= 0.95) "alert-success"
             else if (rep$pct_reference_covered >= 0.5) "alert-warning"
             else "alert-info"
      shiny::div(class = paste("alert", cls),
                 panel_validate_summary(
                   rep,
                   custom_label = custom_panel_label(
                     app_state$custom_panel_status, p)))
    })

    # Custom-panel rows whose `gene` is not in the loaded dataset.
    # `panels()$custom` is the (possibly overridden) uploaded panel; NULL
    # for tissues without a default and no upload.
    custom_coverage <- shiny::reactive({
      shiny::req(app_state$xen)
      p <- panels()
      cust <- p$custom
      if (is.null(cust) || nrow(cust) == 0L) return(NULL)
      data_genes <- rownames(app_state$xen)
      is_present <- cust$gene %in% data_genes
      list(
        label      = custom_panel_label(app_state$custom_panel_status, p),
        total      = nrow(cust),
        n_present  = sum(is_present),
        n_missing  = sum(!is_present),
        missing_df = cust[!is_present, , drop = FALSE]
      )
    })

    output$custom_coverage_summary <- shiny::renderUI({
      if (is.null(app_state$xen)) {
        return(shiny::p(shiny::em(
          "Load a Xenium dataset to see custom-panel coverage.")))
      }
      cc <- custom_coverage()
      if (is.null(cc)) {
        return(shiny::p(shiny::em("No custom-panel data available.")))
      }
      cls <- if (cc$n_missing == 0L)               "alert-success"
             else if (cc$n_missing / cc$total <= 0.1) "alert-info"
             else                                     "alert-warning"
      shiny::div(
        class = paste("alert", cls),
        sprintf(
          "%d of %d genes from %s are present in the loaded dataset (%d missing).",
          cc$n_present, cc$total, cc$label, cc$n_missing))
    })

    output$missing_table <- DT::renderDT({
      shiny::req(app_state$xen)
      cc <- custom_coverage()
      shiny::req(cc)
      df <- cc$missing_df
      if (nrow(df) == 0L) {
        return(DT::datatable(
          data.frame(message = "All custom-panel genes are present in the loaded data."),
          rownames = FALSE,
          options  = list(dom = "t"),
          class    = "stripe compact"))
      }
      DT::datatable(
        df, rownames = FALSE, filter = "top",
        options = list(pageLength = 25, scrollX = TRUE),
        class   = "stripe hover compact nowrap")
    })

    output$dl_missing <- shiny::downloadHandler(
      filename = function() {
        lbl <- if (is.null(app_state$xen)) "custom"
               else custom_panel_label(app_state$custom_panel_status,
                                       panels())
        sprintf("missing_in_data_%s_%s.csv",
                lbl, format(Sys.time(), "%Y%m%dT%H%M%S"))
      },
      content = function(file) {
        cc <- custom_coverage()
        df <- if (is.null(cc)) data.frame() else cc$missing_df
        utils::write.csv(df, file, row.names = FALSE)
      })

    output$spatial <- plotly::renderPlotly({
      out <- xen_load()
      shiny::req(out)
      if (inherits(out, "load_xenium_error")) return(NULL)
      md <- out@meta.data
      if (!all(c("x_centroid", "y_centroid") %in% names(md))) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
                 plotly::layout(annotations = list(list(
                   text = "x_centroid / y_centroid not in meta.data",
                   showarrow = FALSE, x = 0.5, y = 0.5,
                   xref = "paper", yref = "paper"))))
      }
      n  <- nrow(md)
      sz <- if (n > 50000L) 2 else if (n > 10000L) 3 else 4
      xr <- range(md$x_centroid, na.rm = TRUE)
      yr <- range(md$y_centroid, na.rm = TRUE)
      plotly::plot_ly(
        md, x = ~x_centroid, y = ~y_centroid,
        type = "scattergl", mode = "markers",
        marker = list(size = sz, opacity = 0.7,
                      color = "rgba(31,119,180,0.6)"),
        hoverinfo = "skip"
      ) |>
        plotly::layout(
          xaxis = list(title = "x_centroid",
                       range = xr,
                       scaleanchor = "y", scaleratio = 1,
                       constrain = "domain"),
          yaxis = list(title = "y_centroid",
                       range = yr,
                       constrain = "domain"),
          showlegend = FALSE
        )
    })
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
