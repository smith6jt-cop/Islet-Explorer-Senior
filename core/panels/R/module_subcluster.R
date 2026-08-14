#' Subcluster module — drill into a chosen cluster, rerun the pipeline,
#' and stack the results so the user can pop back (CLAUDE.md §5 tab 6).
#'
#' Inherits the pipeline of §6 via [run_cluster_pipeline()]. State lives
#' in `app_state$cluster_stack`, a list whose head ("top of stack") is
#' the currently-active view. Each entry:
#'   list(obj, label, parent_res, parent_cluster, opts)
#'
#' The Cluster tab's `xen_clustered` is the immutable root run; this
#' module never mutates it. When `xen_clustered` changes, the app server
#' resets the stack with that root as its only entry.

subcluster_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::h5("Drill into"),
      shiny::selectInput(ns("parent_res"), "Parent resolution",
                         choices = NULL),
      shiny::selectizeInput(ns("parent_cluster"),
                            "Cluster(s) to keep",
                            choices = NULL, multiple = TRUE,
                            options = list(placeholder = "(pick one or more)")),
      shiny::hr(),
      shiny::h5("Feature set (optional)"),
      shiny::selectizeInput(ns("subpanels"),
                            "Subpanels (union → variable features)",
                            choices = NULL, multiple = TRUE,
                            options = list(placeholder = "(inherit parent)")),
      shiny::checkboxInput(ns("use_all_features"),
                           "Use all genes",
                           value = TRUE),
      shiny::hr(),
      shiny::h5("Pipeline options (slim)"),
      shiny::numericInput(ns("npcs"), "PCs", 20, min = 5, max = 50),
      shiny::sliderInput(ns("res_range"), "Resolution range",
                         min = 0.1, max = 2.0, value = c(0.2, 1.0),
                         step = 0.1),
      shiny::sliderInput(ns("res_step"), "Step",
                         min = 0.05, max = 0.5, value = 0.2, step = 0.05),
      shiny::numericInput(ns("seed"), "Seed", 1, min = 0, max = .Machine$integer.max),
      shiny::hr(),
      shiny::actionButton(ns("drill"), "Drill in",
                          class = "btn-primary",
                          icon = shiny::icon("magnifying-glass-plus")),
      shiny::actionButton(ns("back"), "Back",
                          class = "btn-outline-secondary",
                          icon = shiny::icon("arrow-left"))
    ),
    bslib::card(
      bslib::card_header("Stack"),
      shiny::uiOutput(ns("stack_breadcrumb"))
    ),
    bslib::card(
      bslib::card_header("Status"),
      shiny::uiOutput(ns("status"))
    ),
    shiny::conditionalPanel(
      sprintf("output['%s']", ns("active_ready")),
      bslib::layout_columns(
        col_widths = c(7, 5),
        bslib::card(
          bslib::card_header("UMAP — coloured by selected resolution"),
          shiny::sliderInput(ns("active_res"), "Resolution",
                             min = 0.1, max = 2.0, value = 0.5, step = 0.1),
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("umap"), height = "520px"),
            type = 8, color = "#3498db")
        ),
        bslib::card(
          bslib::card_header("Spatial scatter"),
          shinycssloaders::withSpinner(
            shiny::plotOutput(ns("spatial"), height = "520px"),
            type = 8, color = "#3498db")
        )
      )
    )
  )
}

subcluster_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    active_obj <- shiny::reactive({
      stk <- app_state$cluster_stack
      if (!length(stk)) return(NULL)
      stk[[length(stk)]]$obj
    })

    have_active <- shiny::reactive(!is.null(active_obj()))

    output$active_ready <- shiny::reactive(have_active())
    shiny::outputOptions(output, "active_ready", suspendWhenHidden = FALSE)

    # Populate the parent-resolution + parent-cluster choices off the
    # current top of stack.
    shiny::observe({
      x <- active_obj()
      shiny::req(x)
      cols <- cluster_resolution_columns(x)
      vals <- cluster_resolution_values(x)
      ord <- order(vals); cols <- cols[ord]; vals <- vals[ord]
      shiny::updateSelectInput(session, "parent_res",
                               choices = setNames(cols, sprintf("res = %g", vals)))
      shiny::updateSelectizeInput(session, "subpanels",
                                  choices = sort(panels()$meta$subpanel_keys),
                                  server  = TRUE)
    })

    shiny::observe({
      x   <- active_obj()
      col <- input$parent_res
      shiny::req(x, col, col %in% names(x@meta.data))
      cl <- sort(unique(as.character(x@meta.data[[col]])))
      shiny::updateSelectizeInput(session, "parent_cluster",
                                  choices = cl, server = TRUE,
                                  selected = cl[1])
    })

    # Drill in: build a subset Seurat from the parent, rerun pipeline,
    # push onto the stack.
    waiter::useWaiter()

    shiny::observeEvent(input$drill, {
      x   <- active_obj()
      col <- input$parent_res
      pcl <- input$parent_cluster
      shiny::req(x, col, length(pcl))

      keep <- which(as.character(x@meta.data[[col]]) %in% pcl)
      if (length(keep) < 30L) {
        app_state$subcluster_error <-
          sprintf("Only %d cells in cluster(s) %s — too few to recluster.",
                  length(keep), paste(pcl, collapse = ","))
        return()
      }
      sub <- x[, keep]

      # Drop parent's res cols so the run starts clean.
      for (c0 in cluster_resolution_columns(sub)) {
        sub@meta.data[[c0]] <- NULL
      }

      opts <- list(
        subpanels         = input$subpanels         %||% character(),
        use_all_features  = isTRUE(input$use_all_features) ||
                            !length(input$subpanels),
        npcs              = input$npcs,
        resolutions       = seq(input$res_range[1], input$res_range[2],
                                by = input$res_step),
        seed              = as.integer(input$seed)
      )

      w <- waiter::Waiter$new(
        html  = shiny::tagList(waiter::spin_dots(),
                               shiny::h4("Reclustering subset…")),
        color = "rgba(0,0,0,0.55)")
      w$show(); on.exit(w$hide(), add = TRUE)

      out <- tryCatch(run_cluster_pipeline(sub, panels(), opts),
                      error = function(e) e)
      if (inherits(out, "error")) {
        app_state$subcluster_error <- conditionMessage(out)
        return()
      }
      app_state$subcluster_error <- NULL

      label <- sprintf("%s = {%s}",
                       sub("^seurat_clusters_", "", col),
                       paste(pcl, collapse = ","))
      stk <- app_state$cluster_stack %||% list()
      stk[[length(stk) + 1L]] <- list(
        obj            = out,
        label          = label,
        parent_res     = col,
        parent_cluster = pcl,
        opts           = opts
      )
      app_state$cluster_stack <- stk

      avail <- cluster_resolution_values(out)
      if (length(avail)) {
        shiny::updateSliderInput(session, "active_res",
                                 min = min(avail), max = max(avail),
                                 step = stats::median(diff(sort(avail))) %|na|% 0.1,
                                 value = avail[ceiling(length(avail) / 2)])
      }
    })

    shiny::observeEvent(input$back, {
      stk <- app_state$cluster_stack %||% list()
      if (length(stk) > 1L) {
        app_state$cluster_stack <- stk[-length(stk)]
        app_state$subcluster_error <- NULL
      }
    })

    output$stack_breadcrumb <- shiny::renderUI({
      stk <- app_state$cluster_stack %||% list()
      if (!length(stk)) {
        return(shiny::p(shiny::em(
          "No clustered run yet. Visit the Cluster tab.")))
      }
      labels <- vapply(stk, function(e) e$label %||% "(unnamed)", character(1))
      shiny::tags$ol(class = "breadcrumb",
        lapply(seq_along(labels), function(i) {
          cls <- if (i == length(labels)) "breadcrumb-item active" else "breadcrumb-item"
          shiny::tags$li(class = cls,
                         sprintf("%d. %s (%d cells)",
                                 i, labels[i], ncol(stk[[i]]$obj)))
        }))
    })

    output$status <- shiny::renderUI({
      err <- app_state$subcluster_error
      if (!is.null(err)) {
        return(shiny::div(class = "alert alert-danger", err))
      }
      stk <- app_state$cluster_stack %||% list()
      if (!length(stk)) {
        return(shiny::div(class = "alert alert-info",
                          "Run the Cluster tab to seed the stack."))
      }
      top <- stk[[length(stk)]]
      shiny::div(class = "alert alert-success",
                 sprintf("Active view: %s (%d cells, %d res cols).",
                         top$label %||% "root",
                         ncol(top$obj),
                         length(cluster_resolution_columns(top$obj))))
    })

    active_col <- shiny::reactive({
      x <- active_obj()
      shiny::req(x)
      avail <- cluster_resolution_values(x)
      if (!length(avail)) return(NULL)
      r <- avail[which.min(abs(avail - input$active_res))]
      sprintf("seurat_clusters_res_%g", r)
    })

    output$umap <- plotly::renderPlotly({
      x <- active_obj()
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

    output$spatial <- shiny::renderPlot({
      x <- active_obj()
      col <- active_col()
      shiny::req(x, col)
      md <- x@meta.data
      shiny::validate(shiny::need(
        all(c("x_centroid", "y_centroid") %in% names(md)),
        "x_centroid / y_centroid not in meta.data."))
      df <- data.frame(x = md$x_centroid, y = md$y_centroid,
                       cluster = factor(md[[col]]))
      g <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = cluster))
      pt <- if (nrow(df) > 50000L &&
                requireNamespace("ggrastr", quietly = TRUE)) {
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
  })
}
