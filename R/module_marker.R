#' Markers module — presto::wilcoxauc tables, dotplot, heatmap, features.
#'
#' Operates on whichever Seurat object the user picks: the immutable
#' Cluster-tab root or the current top of the Subcluster stack. Caches
#' the marker table per-(object, group_col) so successive UI tweaks don't
#' recompute. Top-N table feeds the ranked bar, the heatmap, and the
#' DotPlot. The FeaturePlot panel is independent — pick any gene from the
#' assay to colour the UMAP and the spatial scatter.

marker_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      shiny::h5("Source"),
      shiny::radioButtons(ns("source"), NULL, inline = TRUE,
                          choices = c("Root (Cluster tab)" = "root",
                                      "Stack top (Subcluster)" = "top")),
      shiny::selectInput(ns("group_col"), "Grouping (resolution column)",
                         choices = NULL),
      shiny::hr(),
      shiny::h5("Filters"),
      shiny::sliderInput(ns("min_pct_in"), "Minimum pct_in",
                         min = 0, max = 100, value = 10, step = 1),
      shiny::sliderInput(ns("max_padj"), "Maximum padj",
                         min = 0, max = 0.5, value = 0.05, step = 0.01),
      shiny::sliderInput(ns("topn"), "Top-N per cluster",
                         min = 3, max = 30, value = 10, step = 1),
      shiny::hr(),
      shiny::h5("FeaturePlot gene"),
      shiny::selectizeInput(ns("feature_gene"), NULL,
                            choices = NULL,
                            options = list(placeholder = "type a gene")),
      shiny::hr(),
      shiny::actionButton(ns("compute"), "Compute markers",
                          class = "btn-primary",
                          icon = shiny::icon("flask")),
      shiny::downloadButton(ns("dl_csv"), "Download markers CSV")
    ),
    bslib::card(
      bslib::card_header("Status"),
      shiny::uiOutput(ns("status"))
    ),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        bslib::card_header("Top-N markers (filtered)"),
        DT::DTOutput(ns("top_table"))
      ),
      bslib::card(
        bslib::card_header("Ranked AUC — top-N per cluster"),
        shinycssloaders::withSpinner(
          plotly::plotlyOutput(ns("rank_bar"), height = "520px"),
          type = 8, color = "#3498db")
      )
    ),
    bslib::card(
      bslib::card_header("DotPlot — average expression × pct_in"),
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("dotplot"), height = "560px"),
        type = 8, color = "#3498db")
    ),
    bslib::card(
      bslib::card_header("Heatmap — mean expression of top-N markers per cluster"),
      shinycssloaders::withSpinner(
        shiny::plotOutput(ns("heatmap"), height = "640px"),
        type = 8, color = "#3498db")
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        bslib::card_header(shiny::textOutput(ns("feature_caption_umap"),
                                              inline = TRUE)),
        plotly::plotlyOutput(ns("feature_umap"), height = "440px")
      ),
      bslib::card(
        bslib::card_header(shiny::textOutput(ns("feature_caption_spatial"),
                                              inline = TRUE)),
        shiny::plotOutput(ns("feature_spatial"), height = "440px")
      )
    )
  )
}

marker_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    active_obj <- shiny::reactive({
      if (identical(input$source, "top")) {
        stk <- app_state$cluster_stack %||% list()
        if (length(stk)) stk[[length(stk)]]$obj else NULL
      } else {
        app_state$xen_clustered
      }
    })

    have_active <- shiny::reactive(!is.null(active_obj()))

    # Populate group-col + feature-gene off the active object
    shiny::observe({
      x <- active_obj()
      shiny::req(x)
      cols <- cluster_resolution_columns(x)
      vals <- cluster_resolution_values(x)
      ord <- order(vals); cols <- cols[ord]; vals <- vals[ord]
      shiny::updateSelectInput(session, "group_col",
        choices  = setNames(cols, sprintf("res = %g", vals)),
        selected = cols[ceiling(length(cols) / 2)])
      shiny::updateSelectizeInput(session, "feature_gene",
        choices = rownames(x), server = TRUE)
    })

    # Cache per (digest of active object's last_run_id + col). Lifted
    # to app_state so the Export module can include the most-recent
    # markers in the session report.
    markers_df <- shiny::eventReactive(input$compute, {
      x <- active_obj()
      shiny::req(x, input$group_col)
      key   <- paste(x@misc$last_run_id %||% "noid", input$group_col,
                     sep = "::")
      cache <- app_state$markers_cache %||% list()
      if (!is.null(cache[[key]])) {
        app_state$markers_last_key <- key
        return(cache[[key]])
      }
      m <- tryCatch(compute_markers(x, input$group_col),
                    error = function(e) e)
      if (inherits(m, "error")) {
        app_state$markers_error <- conditionMessage(m)
        return(NULL)
      }
      app_state$markers_error <- NULL
      cache[[key]] <- m
      app_state$markers_cache    <- cache
      app_state$markers_last_key <- key
      m
    }, ignoreInit = TRUE)

    top_df <- shiny::reactive({
      m <- markers_df()
      shiny::req(m)
      top_markers(m, n = input$topn,
                  min_pct_in = input$min_pct_in / 100,
                  max_padj   = input$max_padj)
    })

    output$status <- shiny::renderUI({
      if (!have_active()) {
        return(shiny::div(class = "alert alert-info",
                          "Run the Cluster pipeline (or the Subcluster tab) to populate markers."))
      }
      err <- app_state$markers_error
      if (!is.null(err)) {
        return(shiny::div(class = "alert alert-danger", err))
      }
      m <- markers_df()
      if (is.null(m)) {
        return(shiny::div(class = "alert alert-secondary",
                          "Click ", shiny::strong("Compute markers"),
                          " to run presto::wilcoxauc on the chosen group."))
      }
      tdf <- top_df()
      shiny::div(class = "alert alert-success",
        sprintf("%d markers across %d clusters; %d in the filtered top-N.",
                nrow(m), length(unique(m$group)), nrow(tdf)))
    })

    output$top_table <- DT::renderDT({
      df <- top_df()
      shiny::req(nrow(df) > 0)
      num <- vapply(df, is.numeric, logical(1))
      df[num] <- lapply(df[num], function(v) signif(v, 4))
      DT::datatable(df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 25, scrollX = TRUE),
                    class = "stripe hover compact nowrap")
    })

    output$rank_bar <- plotly::renderPlotly({
      df <- top_df()
      shiny::req(nrow(df) > 0)
      df$label <- paste0(df$group, " · ", df$feature)
      df <- df[order(df$group, -df$auc), ]
      df$label <- factor(df$label, levels = rev(df$label))
      plotly::plot_ly(df, x = ~auc, y = ~label, color = ~factor(group),
                      type = "bar", orientation = "h",
                      hoverinfo = "text",
                      text = ~sprintf("%s · %s — AUC %.3f, logFC %.3f, padj %.3g",
                                       group, feature, auc, logFC, padj)) |>
        plotly::layout(showlegend = TRUE,
                       legend = list(title = list(text = "cluster")),
                       xaxis = list(title = "AUC"),
                       yaxis = list(title = ""))
    })

    output$dotplot <- shiny::renderPlot({
      x   <- active_obj()
      df  <- top_df()
      shiny::req(x, nrow(df) > 0, input$group_col)
      genes <- unique(df$feature)
      Idents_old <- SeuratObject::Idents(x)
      SeuratObject::Idents(x) <- as.factor(x@meta.data[[input$group_col]])
      on.exit(SeuratObject::Idents(x) <- Idents_old, add = TRUE)
      Seurat::DotPlot(x, features = genes, dot.scale = 4) +
        ggplot2::theme(axis.text.x =
                       ggplot2::element_text(angle = 45, hjust = 1))
    })

    output$heatmap <- shiny::renderPlot({
      x   <- active_obj()
      df  <- top_df()
      shiny::req(x, nrow(df) > 0, input$group_col)
      genes <- unique(df$feature)
      grp   <- as.character(x@meta.data[[input$group_col]])
      data_mat <- SeuratObject::GetAssayData(x,
                                             assay = SeuratObject::DefaultAssay(x),
                                             layer = "data")[genes, , drop = FALSE]
      avg <- vapply(split(seq_len(ncol(data_mat)), grp), function(idx) {
        Matrix::rowMeans(data_mat[, idx, drop = FALSE])
      }, numeric(length(genes)))
      long <- as.data.frame.table(as.matrix(avg),
                                  responseName = "avg_expr",
                                  stringsAsFactors = FALSE)
      names(long)[1:2] <- c("gene", "cluster")
      long$gene <- factor(long$gene, levels = genes)
      ggplot2::ggplot(long,
                      ggplot2::aes(x = cluster, y = gene, fill = avg_expr)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
        ggplot2::scale_fill_viridis_c(option = "C",
                                      name = "avg expr") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(size = 8))
    })

    # FeaturePlot — UMAP + spatial coloured by gene expression
    feature_data <- shiny::reactive({
      x <- active_obj()
      g <- input$feature_gene
      shiny::req(x, g, g %in% rownames(x))
      assay <- SeuratObject::DefaultAssay(x)
      vals <- SeuratObject::GetAssayData(x, assay = assay, layer = "data")[g, ]
      list(gene = g, vals = as.numeric(vals))
    })

    output$feature_caption_umap <- shiny::renderText({
      fd <- feature_data(); shiny::req(fd)
      sprintf("UMAP — %s expression", fd$gene)
    })
    output$feature_caption_spatial <- shiny::renderText({
      fd <- feature_data(); shiny::req(fd)
      sprintf("Spatial — %s expression", fd$gene)
    })

    output$feature_umap <- plotly::renderPlotly({
      x  <- active_obj()
      fd <- feature_data()
      shiny::req(x, fd, "umap" %in% SeuratObject::Reductions(x))
      um <- as.data.frame(SeuratObject::Embeddings(x, "umap"))
      names(um) <- c("UMAP_1", "UMAP_2")
      um$expr <- fd$vals
      plotly::plot_ly(um, x = ~UMAP_1, y = ~UMAP_2,
                      color = ~expr, colors = viridisLite::viridis(50),
                      type = "scattergl", mode = "markers",
                      marker = list(size = 4, opacity = 0.8),
                      hoverinfo = "text",
                      text = ~sprintf("%s = %.3g", fd$gene, expr))
    })

    output$feature_spatial <- shiny::renderPlot({
      x  <- active_obj()
      fd <- feature_data()
      shiny::req(x, fd)
      md <- x@meta.data
      shiny::validate(shiny::need(
        all(c("x_centroid", "y_centroid") %in% names(md)),
        "x_centroid / y_centroid not in meta.data."))
      df <- data.frame(x = md$x_centroid, y = md$y_centroid,
                       expr = fd$vals)
      g <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = expr))
      pt <- if (nrow(df) > 50000L &&
                requireNamespace("ggrastr", quietly = TRUE)) {
        ggrastr::geom_point_rast(size = 0.5, alpha = 0.85,
                                  raster.dpi = 150)
      } else {
        ggplot2::geom_point(size = 0.5, alpha = 0.85)
      }
      g + pt + ggplot2::coord_fixed() +
        ggplot2::scale_colour_viridis_c(option = "C",
                                         name = sprintf("%s expr", fd$gene)) +
        ggplot2::theme_minimal()
    })

    output$dl_csv <- shiny::downloadHandler(
      filename = function() sprintf("markers_%s.csv",
                                     format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        m <- markers_df(); shiny::req(m)
        utils::write.csv(m, file, row.names = FALSE)
      })
  })
}
