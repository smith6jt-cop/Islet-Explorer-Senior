#' Cluster module — the core analytical tab (CLAUDE.md §5 tab 5, §6).
#'
#' Sidebar binds 1:1 to [run_cluster_pipeline()] arguments. A "Run" button
#' triggers the pipeline (with a `waiter` overlay); after the run lands in
#' `app_state$xen_clustered`, the resolution slider reveals UMAP/spatial/
#' cluster-size/sample-stack plots plus a collapsible run-log panel.

cluster_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::h5("Feature set"),
      shiny::selectizeInput(ns("subpanels"), "Subpanels (union → variable features)",
                            choices = NULL, multiple = TRUE,
                            options = list(placeholder = "(none)")),
      shiny::selectizeInput(ns("extra_genes"), "Add individual genes",
                            choices = NULL, multiple = TRUE,
                            options = list(placeholder = "type to search")),
      shiny::checkboxInput(ns("use_all_features"),
                           "Use all genes (override subpanel selection)",
                           value = FALSE),

      shiny::hr(),
      shiny::h5("Cell filters"),
      shiny::sliderInput(ns("nCount"),  "nCount range",
                         min = 0, max = 1, value = c(0, 1)),
      shiny::sliderInput(ns("nFeature"),"nFeature range",
                         min = 0, max = 1, value = c(0, 1)),
      shiny::selectizeInput(ns("meta_filter_col"), "Optional metadata filter (column)",
                            choices = NULL,
                            options = list(placeholder = "(none)",
                                           allowEmptyOption = TRUE)),
      shiny::selectizeInput(ns("meta_filter_vals"), "...keep these values",
                            choices = NULL, multiple = TRUE),

      shiny::hr(),
      shiny::h5("Normalization"),
      shiny::radioButtons(ns("norm_method"), NULL, inline = TRUE,
                          choices = c("LogNormalize" = "LogNormalize",
                                      "SCT"          = "SCT",
                                      "skip"         = "skip"),
                          selected = "LogNormalize"),
      shiny::checkboxInput(ns("do_scale"),
                           shiny::span(
                             "Scale data ",
                             shiny::tags$abbr(
                               title = paste("Disabling scaling is sometimes",
                                             "preferred when low-count genes",
                                             "produce projection artifacts."),
                               "(?)")),
                           value = TRUE),

      shiny::hr(),
      shiny::h5("PCA"),
      shiny::sliderInput(ns("npcs"), "Number of PCs",
                         min = 5, max = 50, value = 30, step = 1),

      shiny::hr(),
      shiny::h5("Batch correction"),
      shiny::selectInput(ns("batch"), "Method",
                         choices = c("None", "Harmony"), selected = "None"),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'Harmony'", ns("batch")),
        shiny::selectizeInput(ns("batch_var"), "Group by",
                              choices = NULL),
        shiny::sliderInput(ns("harmony_theta"), "theta",
                           min = 0.5, max = 4, value = 2, step = 0.1)
      ),

      shiny::hr(),
      shiny::h5("Neighbors / UMAP"),
      shiny::numericInput(ns("k_param"),     "k.param",     30, min = 5,  max = 200),
      shiny::numericInput(ns("n_neighbors"), "n.neighbors", 30, min = 5,  max = 200),
      shiny::numericInput(ns("min_dist"),    "min.dist",    0.3, min = 0,
                          max = 1, step = 0.05),
      shiny::selectInput(ns("metric"), "metric",
                         choices = c("cosine", "euclidean"),
                         selected = "cosine"),

      shiny::hr(),
      shiny::h5("Clustering"),
      shiny::selectInput(ns("cluster_algorithm"), "Algorithm",
                         choices = c("Louvain", "Leiden"),
                         selected = "Louvain"),
      shiny::sliderInput(ns("res_range"), "Resolution range",
                         min = 0.1, max = 2.0, value = c(0.2, 1.0),
                         step = 0.1),
      shiny::sliderInput(ns("res_step"), "Step",
                         min = 0.05, max = 0.5, value = 0.1, step = 0.05),
      shiny::numericInput(ns("seed"), "Random seed", 1, min = 0, max = .Machine$integer.max),

      shiny::hr(),
      shiny::actionButton(ns("run"), "Run pipeline",
                          class = "btn-primary", icon = shiny::icon("play"))
    ),
    bslib::card(
      bslib::card_header("Run status"),
      shiny::uiOutput(ns("status"))
    ),
    shiny::conditionalPanel(
      sprintf("output['%s']", ns("clustered_ready")),
      bslib::layout_columns(
        col_widths = c(8, 4),
        bslib::card(
          bslib::card_header("UMAP — coloured by selected resolution"),
          shiny::sliderInput(ns("active_res"), "Resolution",
                             min = 0.1, max = 2.0, value = 0.5, step = 0.1),
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("umap"), height = "520px"),
            type = 8, color = "#3498db")
        ),
        bslib::card(
          bslib::card_header("Cluster size"),
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("cluster_bar"), height = "520px"),
            type = 8, color = "#3498db")
        )
      ),
      bslib::layout_columns(
        col_widths = c(7, 5),
        bslib::card(
          bslib::card_header("Spatial scatter (rasterized for n > 50k)"),
          shinycssloaders::withSpinner(
            shiny::plotOutput(ns("spatial"), height = "520px"),
            type = 8, color = "#3498db")
        ),
        bslib::card(
          bslib::card_header("Cluster × sample composition"),
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("sample_stack"), height = "520px"),
            type = 8, color = "#3498db")
        )
      ),
      bslib::accordion(
        id = ns("acc"), open = FALSE,
        bslib::accordion_panel("Run log", shiny::verbatimTextOutput(ns("run_log")))
      )
    )
  )
}

cluster_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    # --- Populate UI choices reactively from the loaded data ----------
    shiny::observe({
      shiny::updateSelectizeInput(session, "subpanels",
                                  choices = sort(panels()$meta$subpanel_keys),
                                  server  = TRUE)
    })

    shiny::observe({
      x <- app_state$xen
      shiny::req(x)
      shiny::updateSelectizeInput(session, "extra_genes",
                                  choices = rownames(x), server = TRUE)
      assay <- SeuratObject::DefaultAssay(x)
      ncol_var  <- paste0("nCount_",   assay)
      nfeat_var <- paste0("nFeature_", assay)
      md <- x@meta.data
      if (ncol_var %in% names(md)) {
        rng <- range(md[[ncol_var]],   na.rm = TRUE)
        shiny::updateSliderInput(session, "nCount",
                                 min = rng[1], max = rng[2], value = rng)
      }
      if (nfeat_var %in% names(md)) {
        rng <- range(md[[nfeat_var]], na.rm = TRUE)
        shiny::updateSliderInput(session, "nFeature",
                                 min = rng[1], max = rng[2], value = rng)
      }
      cats <- names(md)[vapply(md, function(v) {
        is.character(v) || is.factor(v)
      }, logical(1))]
      shiny::updateSelectizeInput(session, "meta_filter_col",
                                  choices = c("(none)", cats), server = TRUE)
      shiny::updateSelectizeInput(session, "batch_var",
                                  choices = cats, server = TRUE,
                                  selected = if ("orig.ident" %in% cats)
                                    "orig.ident" else cats[1])
    })

    shiny::observe({
      x <- app_state$xen
      col <- input$meta_filter_col
      shiny::req(x, col)
      if (identical(col, "(none)") || !nzchar(col)) {
        shiny::updateSelectizeInput(session, "meta_filter_vals",
                                    choices = character())
        return()
      }
      vals <- sort(unique(as.character(x@meta.data[[col]])))
      shiny::updateSelectizeInput(session, "meta_filter_vals",
                                  choices = vals, server = TRUE)
    })

    # --- Run pipeline on click ---------------------------------------
    shiny::observeEvent(input$run, {
      if (is.null(app_state$xen)) {
        shinyWidgets::sendSweetAlert(
          session = session,
          title   = "No dataset loaded",
          text    = paste("Load a Xenium run on the 'Load Xenium' tab",
                          "before running the cluster pipeline."),
          type    = "info",
          btn_labels = "OK"
        )
        return()
      }
      x <- app_state$xen
      opts <- list(
        subpanels         = input$subpanels         %||% character(),
        extra_genes       = input$extra_genes       %||% character(),
        use_all_features  = isTRUE(input$use_all_features),
        nCount_min        = input$nCount[1],
        nCount_max        = input$nCount[2],
        nFeature_min      = input$nFeature[1],
        nFeature_max      = input$nFeature[2],
        meta_filter_col   = if (identical(input$meta_filter_col, "(none)")) NULL
                            else input$meta_filter_col,
        meta_filter_vals  = input$meta_filter_vals,
        norm_method       = input$norm_method,
        do_scale          = isTRUE(input$do_scale),
        npcs              = input$npcs,
        batch             = input$batch,
        batch_var         = input$batch_var,
        harmony_theta     = input$harmony_theta,
        k_param           = input$k_param,
        n_neighbors       = input$n_neighbors,
        min_dist          = input$min_dist,
        metric            = input$metric,
        cluster_algorithm = input$cluster_algorithm,
        resolutions       = seq(input$res_range[1], input$res_range[2],
                                by = input$res_step),
        seed              = as.integer(input$seed)
      )
      out <- shiny::withProgress(
        message = "Cluster pipeline",
        detail  = "starting…",
        value   = 0,
        expr = {
          tryCatch(
            run_cluster_pipeline(
              x, panels(), opts,
              progress = function(amount, message, detail = NULL) {
                shiny::incProgress(amount = amount,
                                   message = message,
                                   detail = detail)
              }
            ),
            error = function(e) e
          )
        }
      )
      if (inherits(out, "error")) {
        app_state$cluster_error <- conditionMessage(out)
        return()
      }
      app_state$cluster_error <- NULL
      app_state$xen_clustered <- out
      # snap the active-resolution slider to a value that exists
      avail <- cluster_resolution_values(out)
      if (length(avail)) {
        shiny::updateSliderInput(session, "active_res",
                                 min = min(avail), max = max(avail),
                                 step = stats::median(diff(sort(avail))) %|na|% 0.1,
                                 value = avail[ceiling(length(avail) / 2)])
      }
    })

    # Clustree tab can request the active resolution be set here.
    shiny::observeEvent(app_state$cluster_jump_res, {
      x <- app_state$xen_clustered
      shiny::req(x)
      avail <- cluster_resolution_values(x)
      if (!length(avail)) return()
      r <- avail[which.min(abs(avail - app_state$cluster_jump_res))]
      shiny::updateSliderInput(session, "active_res", value = r)
      app_state$cluster_jump_res <- NULL
    }, ignoreNULL = TRUE)

    output$status <- shiny::renderUI({
      if (!is.null(app_state$cluster_error)) {
        return(shiny::div(class = "alert alert-danger",
                          shiny::strong("Pipeline error: "),
                          app_state$cluster_error))
      }
      x <- app_state$xen_clustered
      if (is.null(x)) {
        msg <- if (is.null(app_state$xen))
          "Load a Xenium dataset, then configure and click Run."
        else
          "Configure controls in the sidebar, then click Run."
        return(shiny::div(class = "alert alert-info", msg))
      }
      hist <- x@misc$pipeline_history[[x@misc$last_run_id]]
      shiny::div(
        class = "alert alert-success",
        shiny::strong("Pipeline complete. "),
        sprintf("%d cells × %d features, %d resolutions clustered (%s).",
                hist$n_cells_out, hist$n_features,
                length(hist$resolutions), hist$algorithm))
    })

    output$clustered_ready <- shiny::reactive({
      !is.null(app_state$xen_clustered)
    })
    shiny::outputOptions(output, "clustered_ready", suspendWhenHidden = FALSE)

    # --- Plotting -----------------------------------------------------
    active_col <- shiny::reactive({
      x <- app_state$xen_clustered
      shiny::req(x)
      avail <- cluster_resolution_values(x)
      if (!length(avail)) return(NULL)
      r <- avail[which.min(abs(avail - input$active_res))]
      sprintf("seurat_clusters_res_%g", r)
    })

    output$umap <- plotly::renderPlotly({
      x <- app_state$xen_clustered
      col <- active_col()
      shiny::req(x, col)
      umap <- as.data.frame(SeuratObject::Embeddings(x, "umap"))
      names(umap) <- c("UMAP_1", "UMAP_2")
      umap$cluster <- factor(x@meta.data[[col]])
      umap$cell    <- rownames(umap)
      plotly::plot_ly(
        umap, x = ~UMAP_1, y = ~UMAP_2,
        color = ~cluster, type = "scattergl", mode = "markers",
        marker = list(size = 4, opacity = 0.8),
        text = ~paste0("cell: ", cell, "<br>cluster: ", cluster),
        hoverinfo = "text"
      ) |>
        plotly::layout(legend = list(title = list(text = "cluster")))
    })

    output$cluster_bar <- plotly::renderPlotly({
      x <- app_state$xen_clustered
      col <- active_col()
      shiny::req(x, col)
      tab <- as.data.frame(table(cluster = x@meta.data[[col]]))
      tab$cluster <- factor(tab$cluster, levels = tab$cluster)
      plotly::plot_ly(tab, x = ~cluster, y = ~Freq, type = "bar",
                      color = ~cluster) |>
        plotly::layout(showlegend = FALSE,
                       xaxis = list(title = "cluster"),
                       yaxis = list(title = "n cells"))
    })

    output$spatial <- shiny::renderPlot({
      x <- app_state$xen_clustered
      col <- active_col()
      shiny::req(x, col)
      md <- x@meta.data
      shiny::validate(shiny::need(
        all(c("x_centroid", "y_centroid") %in% names(md)),
        "x_centroid / y_centroid not in meta.data — load a Xenium bundle to populate them."
      ))
      df <- data.frame(x = md$x_centroid, y = md$y_centroid,
                       cluster = factor(md[[col]]))
      g <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = cluster))
      pt <- if (nrow(df) > 50000L && requireNamespace("ggrastr", quietly = TRUE)) {
        ggrastr::geom_point_rast(size = 0.5, alpha = 0.7,
                                 raster.dpi = 150)
      } else {
        ggplot2::geom_point(size = 0.5, alpha = 0.7)
      }
      g + pt + ggplot2::coord_fixed() +
        ggplot2::scale_colour_manual(values = scales::hue_pal()(
          length(levels(df$cluster)))) +
        ggplot2::theme_minimal()
    })

    output$sample_stack <- plotly::renderPlotly({
      x <- app_state$xen_clustered
      col <- active_col()
      shiny::req(x, col)
      md <- x@meta.data
      sample_col <- if ("orig.ident" %in% names(md)) "orig.ident" else NULL
      if (is.null(sample_col) ||
          length(unique(md[[sample_col]])) < 2L) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") |>
                 plotly::layout(annotations = list(list(
                   text = "Need ≥ 2 distinct orig.ident values for stacked composition.",
                   showarrow = FALSE, x = 0.5, y = 0.5,
                   xref = "paper", yref = "paper"))))
      }
      tab <- as.data.frame(table(cluster = md[[col]],
                                 sample  = md[[sample_col]]))
      plotly::plot_ly(tab, x = ~cluster, y = ~Freq, color = ~sample,
                      type = "bar") |>
        plotly::layout(barmode = "stack",
                       xaxis = list(title = "cluster"),
                       yaxis = list(title = "n cells"))
    })

    output$run_log <- shiny::renderPrint({
      x <- app_state$xen_clustered
      shiny::req(x)
      h <- x@misc$pipeline_history[[x@misc$last_run_id]]
      h$opts$resolutions <- paste(round(h$opts$resolutions, 3), collapse = ", ")
      str(h, max.level = 2L, give.attr = FALSE)
    })
  })
}

# small helper: NA-safe %||%
`%|na|%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
