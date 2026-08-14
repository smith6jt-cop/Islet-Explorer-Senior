#' Clustree module — tree, overlay, stability summary (CLAUDE.md §7).
#'
#' Read-only with respect to clustering: consumes `seurat_clusters_res_*`
#' columns produced by [run_cluster_pipeline()]. Provides a "Use this
#' resolution" button that signals the Cluster tab to set its active
#' resolution, via `app_state$cluster_jump_res`.

clustree_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      shiny::h5("Resolutions"),
      shiny::checkboxGroupInput(ns("res_cols"),
                                "Include resolutions",
                                choices  = NULL),

      shiny::hr(),
      shiny::h5("Node colour"),
      shiny::radioButtons(ns("node_colour_kind"), NULL,
                          choices = c("cluster id" = "cluster",
                                      "sc3 stability" = "sc3_stability",
                                      "gene expression" = "gene")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'gene'", ns("node_colour_kind")),
        shiny::selectizeInput(ns("gene"), "Gene",
                              choices = NULL,
                              options = list(placeholder = "type a gene"))
      ),

      shiny::hr(),
      shiny::h5("Edges"),
      shiny::radioButtons(ns("edge_colour"), "Edge colour",
                          choices = c("count" = "count",
                                      "in_prop" = "in_prop"),
                          inline = TRUE),
      shiny::sliderInput(ns("prop_filter"),
                         "Edge prop_filter",
                         min = 0, max = 1, value = 0.1, step = 0.05),

      shiny::hr(),
      shiny::h5("Node size"),
      shiny::radioButtons(ns("node_size"), NULL,
                          choices = c("cell count" = "size",
                                      "sc3 stability" = "sc3_stability")),

      shiny::hr(),
      shiny::h5("Use selection"),
      shiny::numericInput(ns("jump_res"),
                          "Jump to resolution (Cluster tab)",
                          value = 0.5, min = 0, max = 5, step = 0.1),
      shiny::actionButton(ns("jump"), "Use this resolution",
                          class = "btn-secondary",
                          icon = shiny::icon("share")),

      shiny::hr(),
      shiny::h5("Downloads"),
      shiny::downloadButton(ns("dl_pdf"), "Tree PDF"),
      shiny::downloadButton(ns("dl_png"), "Tree PNG"),
      shiny::downloadButton(ns("dl_csv"), "Edge table CSV")
    ),
    bslib::card(
      bslib::card_header("Status"),
      shiny::uiOutput(ns("status"))
    ),
    bslib::card(
      bslib::card_header("Resolution-tree"),
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("tree"), height = "720px"),
        type = 8, color = "#3498db")
    ),
    bslib::card(
      bslib::card_header("Tree on UMAP"),
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("overlay"), height = "640px"),
        type = 8, color = "#3498db")
    ),
    bslib::card(
      bslib::card_header("Stability summary"),
      DT::DTOutput(ns("stability"))
    )
  )
}

clustree_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    have_xen <- shiny::reactive(!is.null(app_state$xen_clustered))

    output$status <- shiny::renderUI({
      if (!have_xen()) {
        return(shiny::div(class = "alert alert-info",
                          "Run the Cluster pipeline to populate the tree."))
      }
      x <- app_state$xen_clustered
      cols <- cluster_resolution_columns(x)
      vals <- cluster_resolution_values(x)
      shiny::div(class = "alert alert-success",
                 sprintf("%d resolutions available: %s.",
                         length(cols),
                         paste(round(vals, 3), collapse = ", ")))
    })

    # --- Populate sidebar choices when a clustered object lands -------
    shiny::observe({
      shiny::req(have_xen())
      x <- app_state$xen_clustered
      cols <- cluster_resolution_columns(x)
      vals <- cluster_resolution_values(x)
      ord <- order(vals)
      cols <- cols[ord]; vals <- vals[ord]
      labels <- setNames(cols, sprintf("res = %g", vals))
      shiny::updateCheckboxGroupInput(session, "res_cols",
                                       choices  = labels,
                                       selected = cols)
      shiny::updateSelectizeInput(session, "gene",
                                  choices = rownames(x),
                                  server  = TRUE)
      mid <- vals[ceiling(length(vals) / 2)]
      shiny::updateNumericInput(session, "jump_res", value = mid,
                                min = min(vals), max = max(vals),
                                step = 0.1)
    })

    # --- Build a temporary Seurat with only the chosen res cols -------
    sub_obj <- shiny::reactive({
      shiny::req(have_xen())
      x <- app_state$xen_clustered
      cols <- input$res_cols
      shiny::req(length(cols) >= 2L)
      keep_meta <- setdiff(cluster_resolution_columns(x), cols)
      if (length(keep_meta)) {
        x@meta.data[keep_meta] <- NULL
      }
      x
    })

    tree_plot <- shiny::reactive({
      x <- sub_obj()
      shiny::req(x)
      args <- list(
        x          = x,
        prefix     = "seurat_clusters_res_",
        edge_arrow = FALSE,
        edge_width = 1
      )
      args$node_colour <- switch(input$node_colour_kind,
        cluster        = "seurat_clusters_res_",
        sc3_stability  = "sc3_stability",
        gene           = if (nzchar(input$gene %||% "")) input$gene else NULL)
      if (is.null(args$node_colour)) args$node_colour <- "seurat_clusters_res_"
      args$node_size  <- if (identical(input$node_size, "sc3_stability"))
        "sc3_stability" else "size"
      args$edge_colour <- input$edge_colour
      args$prop_filter <- input$prop_filter
      tryCatch(do.call(clustree::clustree, args),
               error = function(e) {
                 ggplot2::ggplot() +
                   ggplot2::annotate("text", x = 0.5, y = 0.5,
                     label = paste("clustree error:", conditionMessage(e))) +
                   ggplot2::theme_void()
               })
    })

    overlay_plot <- shiny::reactive({
      x <- sub_obj()
      shiny::req(x)
      shiny::validate(shiny::need(
        "umap" %in% SeuratObject::Reductions(x),
        "UMAP reduction missing — re-run the Cluster pipeline."))
      um <- SeuratObject::Embeddings(x, "umap")
      df <- x@meta.data
      df$UMAP_1 <- um[, 1]; df$UMAP_2 <- um[, 2]
      tryCatch(
        clustree::clustree_overlay(
          df,
          prefix      = "seurat_clusters_res_",
          x_value     = "UMAP_1",
          y_value     = "UMAP_2",
          edge_arrow  = FALSE,
          prop_filter = input$prop_filter
        ),
        error = function(e) {
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
              label = paste("clustree_overlay error:", conditionMessage(e))) +
            ggplot2::theme_void()
        })
    })

    output$tree    <- shiny::renderPlot(tree_plot())
    output$overlay <- shiny::renderPlot(overlay_plot())

    output$stability <- DT::renderDT({
      shiny::req(have_xen())
      df <- compute_stability_summary(app_state$xen_clustered)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(v) signif(v, 4))
      DT::datatable(df, rownames = FALSE,
                    options = list(pageLength = 25,
                                   order = list(list(0, "asc"))),
                    class = "stripe hover compact")
    })

    # Cross-tab handoff
    shiny::observeEvent(input$jump, {
      app_state$cluster_jump_res <- input$jump_res
      app_state$nav_target       <- "Cluster"
    })

    # Downloads
    output$dl_pdf <- shiny::downloadHandler(
      filename = function() sprintf("clustree_%s.pdf",
                                    format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        ggplot2::ggsave(file, tree_plot(), device = "pdf",
                        width = 10, height = 9)
      })
    output$dl_png <- shiny::downloadHandler(
      filename = function() sprintf("clustree_%s.png",
                                    format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        ggplot2::ggsave(file, tree_plot(), device = "png",
                        width = 10, height = 9, dpi = 200)
      })
    output$dl_csv <- shiny::downloadHandler(
      filename = function() sprintf("clustree_edges_%s.csv",
                                    format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        shiny::req(have_xen())
        utils::write.csv(clustree_edge_table(app_state$xen_clustered),
                         file, row.names = FALSE)
      })
  })
}
