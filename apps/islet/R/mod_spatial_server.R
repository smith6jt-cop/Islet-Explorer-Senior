# ---------- Spatial Tab module server ----------
# Exports: spatial_server(id, prepared)
#
# Dependencies:
#   prepared()$comp — must contain peri_prop_*, immune_*, enrich_z_*, leiden_* columns
#   prepared()$neighborhood — raw neighborhood data with leiden_umap_1/2
#   PHENOTYPE_COLORS — from drilldown_helpers.R
#   donor_tissue_available(), get_available_donors(), load_donor_tissue() — from spatial_helpers.R

#' @param cohort the shared cohort filter from cohortFilterServer(), or NULL.
#'   NULL keeps this module's own inputs, which is what the testthat suite and
#'   any standalone use get; app.R always passes one.
spatial_server <- function(id, prepared, palette = reactive(PHENOTYPE_COLORS),
                           donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                           remove_outliers = reactive(FALSE),
                           outlier_threshold = reactive(3.0),
                           cohort = NULL) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Donor status, from the cohort panel when it is mounted and from this
    # tab's own checkbox group otherwise. One reader, so the five sites below
    # cannot drift apart from each other the way the three outlier mirrors did.
    active_groups <- reactive({
      if (!is.null(cohort)) cohort$groups() else input$groups
    })

    # ---- Outlier-table reactiveVals (parity with Plot/Trajectory) ----
    # Populated inside the Card A and Card C renderPlotly bodies, consumed by
    # the namespaced "Show excluded-outlier table" toggles below each card.
    infiltration_outliers <- reactiveVal(NULL)
    distance_outliers <- reactiveVal(NULL)

    # 20-color qualitative palette for Leiden clusters
    leiden_palette <- c(
      "#1f77b4", "#ff7f0e", "#66c2a5", "#d62728", "#8da0cb",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
      "#9467bd", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
    )

    # ---- Helpers: column identification ----
    get_nbr_columns <- function(comp_cols) {
      list(
        peri_prop  = grep("^peri_prop_", comp_cols, value = TRUE),
        peri_count = grep("^peri_count_", comp_cols, value = TRUE),
        immune     = intersect(c("immune_frac_peri", "immune_frac_core", "immune_ratio",
                                 "cd8_to_macro_ratio", "tcell_density_peri"), comp_cols),
        enrich     = grep("^enrich_z_", comp_cols, value = TRUE),
        distance   = grep("^min_dist_", comp_cols, value = TRUE),
        leiden     = grep("^leiden_[0-9]", comp_cols, value = TRUE),
        leiden_umap = intersect(c("leiden_umap_1", "leiden_umap_2"), comp_cols)
      )
    }

    # Offer only the disease statuses this dataset actually has. Both status filters were
    # hard-coded c("ND", "Aab+", "T1D"); under the run cohort, which has no Aab+ donor, that
    # invites the reader to filter on a category whose box can never contain a point. These
    # update statically-mounted inputs, so they cannot race a renderUI mount.
    statuses_present <- reactive({
      pd <- prepared()
      if (is.null(pd$comp) || !"Donor Status" %in% colnames(pd$comp)) return(domain_group_levels())
      present <- unique(as.character(pd$comp$`Donor Status`))
      intersect(domain_group_levels(), present)
    })

    observe({
      present <- statuses_present()
      if (!length(present)) return()
      for (id in c("groups", "dist_status")) {
        current <- isolate(input[[id]])
        updateCheckboxGroupInput(session, id, choices = present,
                                 selected = if (is.null(current)) present
                                            else intersect(current, present),
                                 inline = identical(id, "dist_status"))
      }
    })

    has_neighborhood <- reactive({
      pd <- prepared()
      if (is.null(pd$comp)) return(FALSE)
      nbr <- get_nbr_columns(colnames(pd$comp))
      length(nbr$peri_prop) > 0
    })

    has_leiden <- reactive({
      pd <- prepared()
      if (is.null(pd$comp)) return(FALSE)
      nbr <- get_nbr_columns(colnames(pd$comp))
      length(nbr$leiden) > 0
    })

    has_tissue <- reactive({
      donor_tissue_available()
    })

    # ---- Donor selector ----
    output$donor_selector <- renderUI({
      donors <- get_available_donors()
      if (length(donors) == 0) {
        return(div(style = "font-style: italic; padding-top: 8px",
                   "No per-donor tissue data available."))
      }
      pd <- prepared()
      # Build friendly labels: imageid (Donor Status)
      donor_labels <- sapply(donors, function(d) {
        if (!is.null(pd$comp) && "Donor Status" %in% colnames(pd$comp) &&
            "Case ID" %in% colnames(pd$comp)) {
          ds_match <- pd$comp$`Donor Status`[pd$comp$`Case ID` == as.integer(d)]
          ds <- if (length(ds_match) > 0) ds_match[1] else "?"
          paste0(d, " (", ds, ")")
        } else {
          d
        }
      })
      choices <- setNames(donors, donor_labels)
      selectInput(ns("donor"), "Donor", choices = choices, selected = donors[1])
    })

    # Accept browser input only when it names a server-discovered donor
    # partition, before it reaches a path, query, label, or filename.
    selected_donor <- reactive({
      donor_id <- input$donor
      if (is.null(donor_id)) return(NULL)
      if (!islet_valid_scalar_choice(
        donor_id, get_available_donors(), "spatial_donor_invalid"
      )) return(NULL)
      donor_id
    })

    # ---- Leiden resolution selector ----
    output$leiden_res_selector <- renderUI({
      if (!has_leiden()) return(NULL)
      nbr <- get_nbr_columns(colnames(prepared()$comp))
      # Available resolutions (e.g., leiden_0.3 -> "0.3")
      res_labels <- gsub("^leiden_", "", nbr$leiden)
      choices <- setNames(nbr$leiden, res_labels)
      selectInput(ns("leiden_res"), "Leiden Resolution",
                  choices = choices, selected = nbr$leiden[2] %||% nbr$leiden[1])
    })

    # ---- Donor cells reactive (load tissue CSV for selected donor) ----
    donor_cells <- reactive({
      donor_id <- selected_donor()
      req(donor_id)
      load_donor_tissue(donor_id)
    })

    # ---- Donor status for selected donor ----
    donor_status <- reactive({
      donor_id <- selected_donor()
      req(donor_id)
      pd <- prepared()
      if (!is.null(pd$comp) && "Donor Status" %in% colnames(pd$comp) &&
          "Case ID" %in% colnames(pd$comp)) {
        ds <- pd$comp$`Donor Status`[pd$comp$`Case ID` == as.integer(donor_id)]
        if (length(ds) > 0) return(ds[1])
      }
      "?"
    })

    # Leiden mapping for tissue scatter: map islet_name -> cluster for core/peri cells
    islet_leiden_map <- reactive({
      req(has_leiden())
      pd <- prepared()
      leiden_col <- input$leiden_res
      req(leiden_col, leiden_col %in% colnames(pd$comp))

      # Build islet_key -> cluster mapping for the selected donor
      donor_id <- selected_donor()
      req(donor_id)

      comp <- pd$comp
      donor_mask <- comp$`Case ID` == as.integer(donor_id)
      if (!any(donor_mask)) return(NULL)

      sub <- comp[donor_mask, c("islet_key", leiden_col), drop = FALSE]
      setNames(as.character(sub[[leiden_col]]), as.character(sub$islet_key))
    })

    # ---- Phenotype filter ----
    # Choices come from the CANONICAL vocabulary (`canonical_phenotype_order()`,
    # backed by PHENOTYPE_COLORS), never from `unique(donor_cells()$phenotype)`.
    # Deriving them per-donor is why "CD4 Tcell" could vanish: the default donor
    # (112) has zero CD4 cells, as do 6533 and 6534, so the checkbox simply did
    # not exist on first load. Phenotypes with no cells for the current donor are
    # still listed, so the list never changes shape when the donor changes.
    output$phenotype_filter <- renderUI({
      cells <- tryCatch(donor_cells(), error = function(e) NULL)
      if (is.null(cells)) return(NULL)
      phenos <- canonical_phenotype_order(cells$phenotype)
      if (length(phenos) == 0) return(NULL)
      # Preserve previous selections if available (donor change shouldn't reset)
      prev <- isolate(input$pheno_filter)
      if (!is.null(prev)) {
        selected_phenos <- intersect(prev, phenos)
      } else {
        hide_default <- c("Acinar", "Ductal", "Endocrine", "T cell", "Unknown")
        selected_phenos <- setdiff(phenos, hide_default)
      }
      tagList(
        div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 4px",
          tags$label("Phenotypes", style = "font-weight: 600; font-size: 13px; margin: 0"),
          actionLink(ns("pheno_all"), "All", style = "font-size: 11px"),
          span("|"),
          actionLink(ns("pheno_none"), "None", style = "font-size: 11px")
        ),
        div(style = "column-count: 2; column-gap: 8px",
          checkboxGroupInput(ns("pheno_filter"), NULL,
                             choices = phenos, selected = selected_phenos,
                             inline = FALSE)
        )
      )
    })

    observeEvent(input$pheno_all, {
      cells <- tryCatch(donor_cells(), error = function(e) NULL)
      if (!is.null(cells)) {
        # Same canonical list the checkboxes were built from, so "All" selects
        # every box that exists rather than only the ones this donor populates.
        phenos <- canonical_phenotype_order(cells$phenotype)
        updateCheckboxGroupInput(session, "pheno_filter", selected = phenos)
      }
    })

    observeEvent(input$pheno_none, {
      updateCheckboxGroupInput(session, "pheno_filter", selected = character(0))
    })

    # ---- Zoom state for tissue scatter ----
    scatter_zoom <- reactiveValues(xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL)

    observeEvent(input$scatter_brush, {
      brush <- input$scatter_brush
      if (!is.null(brush)) {
        scatter_zoom$xmin <- brush$xmin
        scatter_zoom$xmax <- brush$xmax
        # y is reversed in the plot, so brush coords are already in data space
        scatter_zoom$ymin <- brush$ymin
        scatter_zoom$ymax <- brush$ymax
      }
    })

    observeEvent(input$scatter_dblclick, {
      scatter_zoom$xmin <- NULL
      scatter_zoom$xmax <- NULL
      scatter_zoom$ymin <- NULL
      scatter_zoom$ymax <- NULL
    })

    observeEvent(input$scatter_reset_zoom, {
      scatter_zoom$xmin <- NULL
      scatter_zoom$xmax <- NULL
      scatter_zoom$ymin <- NULL
      scatter_zoom$ymax <- NULL
    })

    # Reset zoom when donor changes
    observeEvent(input$donor, {
      scatter_zoom$xmin <- NULL
      scatter_zoom$xmax <- NULL
      scatter_zoom$ymin <- NULL
      scatter_zoom$ymax <- NULL
    })

    # Phase 3: pick the rendering backend for the tissue scatter.
    # ggplot2 is the default (handles ~200K cells comfortably and supports
    # the R-side brush zoom). When `ISLET_USE_RDECK=TRUE` and the rdeck
    # package is installed, the WebGL backend is used instead. The switch
    # happens in renderUI so existing outputs do not change.
    output$tissue_scatter_container <- renderUI({
      if (is_rdeck_available()) {
        rdeck::rdeckOutput(ns("tissue_scatter_deck"), height = "700px")
      } else {
        plotOutput(ns("tissue_scatter"), height = "700px",
                   dblclick = ns("scatter_dblclick"),
                   brush = brushOpts(id = ns("scatter_brush"), resetOnNew = TRUE))
      }
    })

    # ==== Card 2 (rdeck variant): WebGL tissue scatter for million-cell donors ====
    output$tissue_scatter_deck <- if (requireNamespace("rdeck", quietly = TRUE)) {
      rdeck::renderRdeck({
        cells <- donor_cells()
        req(cells)

        # Apply region + phenotype filters (same semantics as ggplot path)
        region_mode <- input$region_filter %||% "all"
        selected_phenos <- input$pheno_filter
        plot_df <- cells
        if (!is.null(selected_phenos) && length(selected_phenos) > 0) {
          plot_df <- plot_df[plot_df$phenotype %in% selected_phenos, , drop = FALSE]
        }
        if (region_mode == "core") {
          plot_df <- plot_df[plot_df$cell_region == "core", , drop = FALSE]
        } else if (region_mode == "core_peri") {
          plot_df <- plot_df[plot_df$cell_region %in% c("core", "peri"), , drop = FALSE]
        }
        if (nrow(plot_df) == 0) return(NULL)
        render_tissue_deck(plot_df, palette = palette())
      })
    } else {
      NULL
    }

    # ==== Card 2: Tissue Scatter Plot (ggplot2, rasterized) ====
    spatial_tissue_scatter_plot <- reactive({
      cells <- donor_cells()
      req(cells)

      # Copy so we can filter
      plot_df <- cells

      # Region filter: highlight selected region, dim others
      region_mode <- input$region_filter %||% "all"

      # Determine coloring mode
      color_mode <- input$color_by %||% "phenotype"
      use_leiden <- (color_mode == "leiden") && has_leiden()

      if (use_leiden) {
        # Map islet_name -> leiden cluster
        lmap <- islet_leiden_map()
        if (!is.null(lmap)) {
          plot_df$cluster <- lmap[plot_df$islet_name]
          # Tissue cells with no islet get NA cluster
          plot_df$cluster[is.na(plot_df$cluster)] <- "tissue"
        } else {
          plot_df$cluster <- "tissue"
        }
      }

      # Phenotype filter (only in phenotype mode)
      selected_phenos <- input$pheno_filter
      if (!use_leiden && !is.null(selected_phenos) && length(selected_phenos) > 0) {
        plot_df <- plot_df[plot_df$phenotype %in% selected_phenos, , drop = FALSE]
      }

      # Split into foreground (highlighted) and background (dimmed) layers
      if (region_mode == "all") {
        # All cells colored; tissue background slightly dimmed
        fg <- plot_df[plot_df$cell_region %in% c("core", "peri"), , drop = FALSE]
        bg <- plot_df[plot_df$cell_region == "tissue", , drop = FALSE]
      } else if (region_mode == "core_peri") {
        fg <- plot_df[plot_df$cell_region %in% c("core", "peri"), , drop = FALSE]
        bg <- plot_df[plot_df$cell_region == "tissue", , drop = FALSE]
      } else {
        # Core only
        fg <- plot_df[plot_df$cell_region == "core", , drop = FALSE]
        bg <- plot_df[plot_df$cell_region != "core", , drop = FALSE]
      }

      p <- ggplot2::ggplot()

      color_bg <- isTRUE(input$color_background)

      # Background layer
      if (nrow(bg) > 0) {
        if (color_bg && !use_leiden) {
          # Color background cells by phenotype (dimmed)
          bg_phenos <- sort(unique(bg$phenotype))
          bg_pal <- palette()[bg_phenos]
          bg_pal[is.na(bg_pal)] <- "#CCCCCC"
          p <- p + ggplot2::geom_point(
            data = bg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = phenotype),
            size = 0.15, alpha = 0.25,
            inherit.aes = FALSE
          )
        } else {
          p <- p + ggplot2::geom_point(
            data = bg,
            ggplot2::aes(x = X_centroid, y = Y_centroid),
            color = "#d9d9d9", size = 0.15, alpha = 0.3,
            inherit.aes = FALSE
          )
        }
      }

      # Foreground layer: colored by phenotype or leiden
      if (nrow(fg) > 0) {
        if (use_leiden && "cluster" %in% colnames(fg)) {
          cluster_levels <- sort(unique(fg$cluster[fg$cluster != "tissue"]))
          fg$cluster <- factor(fg$cluster, levels = c(cluster_levels, "tissue"))
          n_clusters <- length(cluster_levels)
          pal_colors <- rep_len(leiden_palette, n_clusters)
          pal <- setNames(pal_colors, cluster_levels)
          pal["tissue"] <- "#d9d9d9"

          p <- p + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = cluster),
            size = 0.4, alpha = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(values = pal, name = "Cluster", na.value = "#d9d9d9",
                                        guide = ggplot2::guide_legend(override.aes = list(size = 4)))
        } else {
          # Phenotype coloring — legend in canonical palette order
          canonical_order <- names(palette())
          if (!is.null(selected_phenos) && length(selected_phenos) > 0) {
            present <- intersect(selected_phenos,
                                 unique(c(fg$phenotype, if (color_bg) bg$phenotype)))
          } else {
            present <- unique(c(fg$phenotype, if (color_bg) bg$phenotype))
          }
          legend_phenos <- intersect(canonical_order, present)
          legend_phenos <- c(legend_phenos, setdiff(present, legend_phenos))
          pal <- palette()[legend_phenos]
          pal[is.na(pal)] <- "#CCCCCC"

          p <- p + ggplot2::geom_point(
            data = fg,
            ggplot2::aes(x = X_centroid, y = Y_centroid, color = phenotype),
            size = 0.4, alpha = 0.6,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_color_manual(values = pal, breaks = legend_phenos,
                                        name = "Phenotype", na.value = "#CCCCCC",
                                        guide = ggplot2::guide_legend(override.aes = list(size = 4)))
        }
      }

      ds <- donor_status()

      # Zoom via coord_cartesian (clips display without dropping data)
      zoomed <- !is.null(scatter_zoom$xmin)
      if (zoomed) {
        # Brush coords are in display space (y already reversed by scale_y_reverse).
        # coord_cartesian ylim is in *data* space (before reversal), so swap y.
        p <- p + ggplot2::scale_y_reverse() +
          ggplot2::coord_cartesian(
            xlim = c(scatter_zoom$xmin, scatter_zoom$xmax),
            ylim = sort(c(scatter_zoom$ymin, scatter_zoom$ymax))
          )
      } else {
        p <- p + ggplot2::coord_fixed() +
          ggplot2::scale_y_reverse()
      }

      donor_id <- selected_donor()
      req(donor_id)
      p + ggplot2::labs(
          title = paste0("Donor ", donor_id, " (", ds, ")"),
          subtitle = paste0(nrow(fg), " foreground cells | ",
                            nrow(bg), " background cells"),
          x = expression(paste("X centroid (", mu, "m)")),
          y = expression(paste("Y centroid (", mu, "m)"))
        ) +
        theme_islet() +
        ggplot2::theme(
          legend.position = "right",
          legend.key.size = ggplot2::unit(0.7, "cm")
        )
    })

    output$tissue_scatter <- shiny::renderPlot({ spatial_tissue_scatter_plot() }, height = 800)
    output$dl_tissue_scatter <- shiny::downloadHandler(
      filename = function() {
        donor_id <- selected_donor()
        islet_png_filename(paste0("spatial_tissue_donor_", donor_id %||% "all"))
      },
      content = function(file) {
        p <- spatial_tissue_scatter_plot(); if (is.null(p)) return()
        ggplot2::ggsave(file, plot = p,
                        width = 12, height = 10,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units,
                        bg = "white")
      }
    )

    # ==== Card 3: Leiden UMAP (plotly, 1015 islets) ====
    output$leiden_not_available <- renderUI({
      if (!has_leiden()) {
        return(div(style = "font-style: italic; margin-bottom: 10px",
          "Islet-level Leiden clustering is not present in this dataset: the islet table carries ",
          "no ", tags$code("leiden_<resolution>"), " or ", tags$code("leiden_umap_1/2"),
          " columns. Every other card on this tab is unaffected."))
      }
      NULL
    })

    output$leiden_umap <- renderPlotly({
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp
      req("leiden_umap_1" %in% colnames(comp), "leiden_umap_2" %in% colnames(comp))

      leiden_col <- input$leiden_res %||% {
        nbr <- get_nbr_columns(colnames(comp))
        nbr$leiden[2] %||% nbr$leiden[1]
      }
      req(leiden_col %in% colnames(comp))

      # Filter by donor status
      plot_comp <- comp
      if (!is.null(active_groups()) && "Donor Status" %in% colnames(plot_comp)) {
        plot_comp <- plot_comp[plot_comp$`Donor Status` %in% active_groups(), , drop = FALSE]
      }

      umap1 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_1))
      umap2 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_2))
      cluster <- as.character(plot_comp[[leiden_col]])

      plot_df <- data.frame(
        umap1 = umap1, umap2 = umap2, cluster = cluster,
        donor_status = if ("Donor Status" %in% colnames(plot_comp)) as.character(plot_comp$`Donor Status`) else "",
        islet_key = if ("islet_key" %in% colnames(plot_comp)) as.character(plot_comp$islet_key) else "",
        stringsAsFactors = FALSE
      )
      plot_df <- plot_df[is.finite(plot_df$umap1) & is.finite(plot_df$umap2), , drop = FALSE]
      if (nrow(plot_df) == 0) return(plotly_empty() %>% layout(title = "No UMAP data"))

      # Sort cluster levels numerically
      cluster_levels <- sort(unique(plot_df$cluster))
      plot_df$cluster <- factor(plot_df$cluster, levels = cluster_levels)
      n_clusters <- length(cluster_levels)
      pal_colors <- rep_len(leiden_palette, n_clusters)
      pal <- setNames(pal_colors, cluster_levels)

      res_label <- gsub("^leiden_", "", leiden_col)

      fonts <- plotly_islet_fonts()
      plot_ly(plot_df, x = ~umap1, y = ~umap2, color = ~cluster,
              colors = pal,
              text = ~paste0("Cluster: ", cluster, "<br>",
                             "Status: ", donor_status, "<br>",
                             islet_key),
              hoverinfo = "text",
              type = "scatter", mode = "markers",
              marker = list(size = 3, opacity = 0.7)) %>%
        layout(
          font = fonts$global,
          title = list(text = paste0("Leiden ", res_label, " (", nrow(plot_df), " islets)"),
                       font = fonts$title),
          xaxis = list(title = list(text = "UMAP 1", font = fonts$axis_title),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE,
                       showline = FALSE, constrain = "domain"),
          yaxis = list(title = list(text = "UMAP 2", font = fonts$axis_title),
                       zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE,
                       showline = FALSE, scaleanchor = "x", scaleratio = 1,
                       constrain = "domain"),
          legend = list(title = list(text = "Cluster", font = fonts$legend_title),
                        font = fonts$legend_text,
                        itemsizing = "trace",
                        tracegroupgap = 2)
        ) %>%
        plotly_islet_config(paste0("leiden_umap_", res_label))
    })

    # ==== Donor Status UMAP (static ggplot, mirrors trajectory tab) ====
    spatial_umap_donor_plot <- reactive({
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp
      req("leiden_umap_1" %in% colnames(comp), "leiden_umap_2" %in% colnames(comp))

      plot_comp <- comp
      if (!is.null(active_groups()) && "Donor Status" %in% colnames(plot_comp)) {
        plot_comp <- plot_comp[plot_comp$`Donor Status` %in% active_groups(), , drop = FALSE]
      }

      umap1 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_1))
      umap2 <- suppressWarnings(as.numeric(plot_comp$leiden_umap_2))
      ds <- if ("Donor Status" %in% colnames(plot_comp)) as.character(plot_comp$`Donor Status`) else rep("?", nrow(plot_comp))

      plot_df <- data.frame(umap1 = umap1, umap2 = umap2, donor_status = ds, stringsAsFactors = FALSE)
      plot_df <- plot_df[is.finite(plot_df$umap1) & is.finite(plot_df$umap2), , drop = FALSE]
      if (nrow(plot_df) == 0) return(NULL)

      plot_df$donor_status <- factor(plot_df$donor_status, levels = c("ND", "Aab+", "T1D"))

      ggplot(plot_df, aes(x = umap1, y = umap2, color = donor_status)) +
        geom_point(alpha = 0.6, size = 1.0) +
        scale_color_manual(values = donor_colors_reactive()) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2", color = "Status") +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_islet() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank())
    })
    output$spatial_umap_donor <- renderPlot({ spatial_umap_donor_plot() })
    output$dl_spatial_umap_donor <- downloadHandler(
      filename = function() islet_png_filename("spatial_umap_donor"),
      content = function(file) {
        p <- spatial_umap_donor_plot(); if (is.null(p)) return()
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units,
                        bg = "white")
      }
    )

    # ==== Card 3 (bottom): Cluster Composition ====
    output$cluster_composition <- renderPlotly({
      req(has_leiden())
      pd <- prepared()
      comp <- pd$comp

      leiden_col <- input$leiden_res %||% {
        nbr <- get_nbr_columns(colnames(comp))
        nbr$leiden[2] %||% nbr$leiden[1]
      }
      req(leiden_col %in% colnames(comp))

      # Filter by donor status
      plot_comp <- comp
      if (!is.null(active_groups()) && "Donor Status" %in% colnames(plot_comp)) {
        plot_comp <- plot_comp[plot_comp$`Donor Status` %in% active_groups(), , drop = FALSE]
      }

      # Get phenotype proportion columns (prop_*), filtered by selected phenotypes
      prop_cols <- grep("^prop_", colnames(plot_comp), value = TRUE)
      if (length(prop_cols) == 0) {
        return(plotly_empty() %>% layout(title = "No phenotype proportion data"))
      }
      if (!is.null(input$pheno_filter) && length(input$pheno_filter) > 0) {
        keep_cols <- paste0("prop_", input$pheno_filter)
        prop_cols <- intersect(prop_cols, keep_cols)
        if (length(prop_cols) == 0) {
          return(plotly_empty() %>% layout(title = "No matching phenotypes selected"))
        }
      }

      cluster <- as.character(plot_comp[[leiden_col]])

      # Compute mean phenotype proportions per cluster
      cluster_levels <- sort(unique(cluster))
      mat <- matrix(NA_real_, nrow = length(prop_cols), ncol = length(cluster_levels),
                    dimnames = list(prop_cols, cluster_levels))
      for (cl in cluster_levels) {
        sub <- plot_comp[cluster == cl, prop_cols, drop = FALSE]
        if (nrow(sub) > 0) {
          mat[, cl] <- colMeans(sub, na.rm = TRUE)
        }
      }

      # Build stacked bar data
      bar_rows <- list()
      for (i in seq_along(prop_cols)) {
        for (j in seq_along(cluster_levels)) {
          bar_rows[[length(bar_rows) + 1]] <- data.frame(
            phenotype = gsub("^prop_", "", prop_cols[i]),
            cluster = cluster_levels[j],
            proportion = mat[i, j],
            stringsAsFactors = FALSE
          )
        }
      }
      bar_df <- do.call(rbind, bar_rows)

      # Remove phenotypes with zero proportion across all clusters
      max_per_pheno <- tapply(bar_df$proportion, bar_df$phenotype, max, na.rm = TRUE)
      keep_phenos <- names(max_per_pheno[max_per_pheno > 0])
      bar_df <- bar_df[bar_df$phenotype %in% keep_phenos, , drop = FALSE]

      # Use active phenotype palette where available
      pheno_present <- unique(bar_df$phenotype)
      pal <- palette()[pheno_present]
      pal[is.na(pal)] <- "#CCCCCC"

      canonical_order <- names(palette())
      ordered_phenos <- intersect(canonical_order, keep_phenos)
      ordered_phenos <- c(ordered_phenos, setdiff(keep_phenos, canonical_order))
      bar_df$phenotype <- factor(bar_df$phenotype, levels = rev(ordered_phenos))

      fonts <- plotly_islet_fonts()
      plot_ly(bar_df, x = ~cluster, y = ~proportion, color = ~phenotype,
              colors = pal,
              type = "bar") %>%
        layout(
          font = fonts$global,
          barmode = "stack",
          title = list(text = "Mean Phenotype Composition by Cluster", font = fonts$title),
          xaxis = list(title = list(text = "Cluster", font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(title = list(text = "Mean Proportion", font = fonts$axis_title),
                       tickfont = fonts$axis_tick, range = c(0, 1)),
          legend = list(font = list(size = max(ISLET_FONT_SIZES$dense_min,
                                                ISLET_FONT_SIZES$legend_text - 2)),
                        title = list(text = "", font = fonts$legend_title))
        ) %>%
        plotly_islet_config("cluster_composition")
    })

    # ---- Download handler ----
    output$download_spatial <- downloadHandler(
      filename = function() {
        paste0("spatial_neighborhood_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        pd <- prepared()
        df <- pd$comp
        if (is.null(df) || nrow(df) == 0) { df <- data.frame() }
        else {
          if (!is.null(active_groups()) && "Donor Status" %in% colnames(df))
            df <- df[df$`Donor Status` %in% active_groups(), , drop = FALSE]
          if ("total_cells_peri" %in% colnames(df))
            df <- df[!is.na(df$total_cells_peri) & df$total_cells_peri > 0, , drop = FALSE]
        }
        nbr <- get_nbr_columns(colnames(df))
        keep_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                       unlist(nbr, use.names = FALSE))
        keep_cols <- intersect(keep_cols, colnames(df))
        write.csv(df[, keep_cols, drop = FALSE], file, row.names = FALSE)
      }
    )

    # ======================================================================
    # NEIGHBORHOOD ANALYSIS CARDS (A: Infiltration, B: Enrichment, C: Proximity)
    # ======================================================================

    # ---- Shared reactive: filtered comp data for neighborhood cards ----
    nbr_comp <- reactive({
      pd <- prepared()
      req(pd$comp)
      comp <- pd$comp
      # The panel applies donor status, age, sex, AAb, min-cells and diameter in
      # one fixed order, and records which of them bound to a column here. The
      # hand-rolled block this replaces applied groups and min-cells only, so
      # the neighborhood cards were scoped differently from every other tab on
      # the same screen -- and its `nbr_min_cells` was a fourth unsynced copy of
      # a slider that also lived on Plot, Trajectory and Statistics.
      if (!is.null(cohort)) return(cohort$apply(comp, cells_expr = "core_peri"))

      if (!is.null(active_groups()) && "Donor Status" %in% colnames(comp))
        comp <- comp[comp$`Donor Status` %in% active_groups(), , drop = FALSE]
      min_cells <- input$nbr_min_cells %||% 1
      if (min_cells > 1) {
        tc_core <- if ("total_cells_core" %in% colnames(comp)) comp$total_cells_core else NA
        tc_peri <- if ("total_cells_peri" %in% colnames(comp)) comp$total_cells_peri else NA
        tc <- ifelse(is.finite(tc_core), tc_core, 0) + ifelse(is.finite(tc_peri), tc_peri, 0)
        comp <- comp[is.finite(tc) & tc >= min_cells, , drop = FALSE]
      }
      comp
    })

    # Phenotype options derived dynamically from prepared()$comp columns.
    # Each entry maps a phenotype label to its associated metric columns.
    # Replaces the previous hand-curated 7-immune-type list so users can
    # interrogate any phenotype the pipeline computed metrics for.
    phenotype_options <- reactive({
      comp <- nbr_comp()
      if (is.null(comp)) return(list())
      cn <- colnames(comp)
      ez_cols <- grep("^enrich_z_", cn, value = TRUE)
      if (length(ez_cols) == 0) return(list())
      build <- function(ez) {
        safe <- sub("^enrich_z_", "", ez)
        # Map safe-name back to the original phenotype label by checking
        # possible un-sanitised forms of the name against existing prop_* cols.
        candidates <- unique(c(
          gsub("plus", "+", gsub("_", " ", safe), fixed = TRUE),
          gsub("_", " ", safe),
          gsub("plus", "+", safe, fixed = TRUE),
          safe
        ))
        prop_core_col <- NA_character_
        label <- safe
        for (cand in candidates) {
          pc <- paste0("prop_", cand)
          if (pc %in% cn) { prop_core_col <- pc; label <- cand; break }
        }
        list(
          label      = label,
          safe       = safe,
          enrich     = ez,
          prop_core  = prop_core_col,
          peri_prop  = {p <- paste0("peri_prop_", safe);  if (p %in% cn) p else NA_character_},
          peri_count = {p <- paste0("peri_count_", safe); if (p %in% cn) p else NA_character_},
          min_dist   = {p <- paste0("min_dist_", safe);   if (p %in% cn) p else NA_character_}
        )
      }
      recs <- lapply(ez_cols, build)
      recs <- recs[order(sapply(recs, function(r) r$label))]
      names(recs) <- sapply(recs, function(r) r$label)
      recs
    })

    # Choice vector for selectInputs — names are user-facing labels, values are safe-names.
    phenotype_choices <- reactive({
      opts <- phenotype_options()
      if (length(opts) == 0) return(character(0))
      vec <- vapply(opts, function(o) o$safe, character(1))
      stats::setNames(vec, vapply(opts, function(o) o$label, character(1)))
    })

    # ---- Intermediate reactive for Card B: enrichment/proportion summary ----
    enrich_summary <- reactive({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp), "Case ID" %in% colnames(comp))

      clip <- isTRUE(input$enrich_clip)
      stat_mode <- input$enrich_stat %||% "Median"
      region <- input$enrich_region %||% "peri"

      statuses <- unique(comp$`Donor Status`)
      comp_cols <- colnames(comp)
      rows <- list()

      for (itm in phenotype_options()) {
        # Select column based on region
        if (region == "peri") {
          col <- itm$enrich      # enrichment z-score (peri)
        } else if (region == "core") {
          col <- itm$prop_core   # core proportion
        } else {
          col <- itm$peri_prop   # peri proportion
        }
        if (is.null(col) || is.na(col) || !(col %in% comp_cols)) next

        for (ds in statuses) {
          mask <- comp$`Donor Status` == ds
          vals <- comp[[col]][mask]
          donors <- comp$`Case ID`[mask]
          finite_mask <- is.finite(vals)
          vals <- vals[finite_mask]
          donors <- donors[finite_mask]
          if (region == "peri" && clip) {
            clip_mask <- vals >= -5 & vals <= 5
            vals <- vals[clip_mask]
            donors <- donors[clip_mask]
          }
          if (length(vals) == 0) next

          # Aggregate to donor-level means first (avoids pseudoreplication;
          # islet-level median is 0 for rare immune types because >50% of islets
          # have no cells of that type)
          donor_means <- tapply(vals, donors, mean, na.rm = TRUE)
          donor_means <- donor_means[is.finite(donor_means)]
          n <- length(donor_means)
          if (n == 0) next

          if (stat_mode == "Median") {
            z_summary <- median(donor_means)
            q <- quantile(donor_means, c(0.25, 0.75))
            z_lo <- q[1]
            z_hi <- q[2]
          } else {
            z_summary <- mean(donor_means)
            se <- if (n > 1) sd(donor_means) / sqrt(n) else 0
            z_lo <- z_summary - se
            z_hi <- z_summary + se
          }
          rows[[length(rows) + 1]] <- data.frame(
            col = col, cell_type = itm$label, donor_status = ds,
            z_summary = z_summary, z_lo = z_lo, z_hi = z_hi, n = n,
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(rows) == 0) return(NULL)
      do.call(rbind, rows)
    })

    # ---- Neighborhood cards UI (conditional on data availability) ----
    output$neighborhood_cards <- renderUI({
      # Cards A-C need the per-region composition columns. Returning NULL made three whole
      # sections disappear with nothing said, so a dataset that simply lacks the columns looked
      # identical to a dataset where they failed to load. Say which is missing instead.
      if (!has_neighborhood()) {
        return(div(class = "card",
          style = "padding: 15px; margin: 15px 0; border-left: 4px solid #CC6633;",
          tags$strong("Neighborhood analysis unavailable for this dataset."),
          tags$p(style = "margin: 6px 0 0",
            "Cards A-C read per-region composition columns (", tags$code("peri_prop_*"),
            ", ", tags$code("enrich_z_*"), ", ", tags$code("min_dist_*"),
            ") from the islet table, and none is present. They are written by ",
            tags$code("neighborhood_metrics.csv"), " and merged into ", tags$code("islets.parquet"),
            "; every other tab is unaffected.")
        ))
      }

      # section_heading helper (inline — cannot share with UI function directly)
      # Was a fourth inline copy of section_heading, with the comment
      # "inline -- cannot share with UI function directly". It can: the
      # helper is sourced from R/ui_modernist_bits.R like everything else.
      sec_heading <- ie_section_heading

      tagList(
        # ==== Display controls for the neighborhood cards ====
        # THREE CONTROLS LEFT, and the two that went are the point: "Exclude
        # outliers" and its z-threshold were non-namespaced MIRRORS of the
        # Plot and Trajectory copies, kept in step by three observers in app.R
        # that each had to name the other two. The cohort panel owns the rule
        # now. Min cells/islet went the same way -- it was a fourth unsynced
        # copy of a filter the panel applies to every tab.
        #
        # What is left is display only: point size and transparency change how
        # this card draws, not which islets it draws. Cards A (z-score) and C
        # (z on log1p(distance)) honour the panel's outlier rule; Card B's
        # "Clip z > 5" is display saturation and stays separate.
        fluidRow(column(12,
          div(class = "card", style = "padding: 12px 15px; margin-bottom: 10px; margin-top: 10px; display: flex; align-items: center; gap: 20px; overflow: visible",
            span(style = "font-weight: 600; font-size: 15px; white-space: nowrap", "Neighborhood analysis"),
            sliderInput(ns("nbr_pt_size"), "Point size", min = 2, max = 12, value = 5, step = 1, width = "150px"),
            sliderInput(ns("nbr_pt_alpha"), "Point transparency", min = 0, max = 0.95, value = 0.6, step = 0.05, width = "150px"),
            uiOutput(ns("nbr_islet_count"))
          )
        )),

        # ==== Card A: Immune Infiltration Overview ====
        fluidRow(column(12, sec_heading(
          "A", "Peri-Islet Phenotype Enrichment",
          "Per-phenotype peri-zone enrichment across disease stages, with peri vs core proportion scatter."
        ))),
        fluidRow(style = "display: flex; flex-wrap: wrap",
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box",
              div(style = "display: flex; gap: 15px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 8px",
                div(style = "flex: 1; min-width: 180px",
                  # Populated server-side from phenotype_choices().
                  selectInput(ns("infiltration_phenotype"), "Phenotype (peri enrichment z)",
                              choices = NULL, width = "100%"),
                  uiOutput(ns("infiltration_pheno_hint"))
                ),
                checkboxInput(ns("infiltration_show_outlier_table"),
                              "Show excluded-outlier table", value = FALSE),
                tags$small(style = "flex-basis: 100%",
                  "Outlier exclusion follows the global toggle above (Tukey 1.5×IQR per donor status on non-zero values, to tolerate phenotypes absent from many islets).")
              ),
              plotlyOutput(ns("infiltration_bars"), height = "400px"),
              uiOutput(ns("infiltration_outlier_info"))
            )
          ),
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box",
              div(style = "display: flex; gap: 15px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 8px",
                div(style = "flex: 1; min-width: 180px",
                  # Populated server-side from phenotype_choices().
                  selectInput(ns("scatter_phenotype"), "Phenotype (peri vs core)",
                              choices = NULL, width = "100%"),
                  uiOutput(ns("scatter_pheno_hint"))
                ),
                checkboxInput(ns("scatter_sqrt"), "Sqrt scale", value = TRUE),
                checkboxInput(ns("scatter_trend"), "Trend lines", value = FALSE)
              ),
              plotlyOutput(ns("infiltration_scatter"), height = "400px")
            )
          )
        ),
        fluidRow(column(12,
          div(style = "padding: 10px 15px; margin-bottom: 20px; font-size: 13px",
            tags$em("Immune infiltration quantifies the proportion of immune cells among all cells in the peri-islet zone (20\u00b5m expansion) and islet core. Higher fractions indicate increased immune surveillance or active infiltration.")
          )
        )),

        # ==== Card B: Immune Cell Enrichment by Type ====
        fluidRow(column(12, sec_heading(
          "B", "Phenotype Composition & Enrichment",
          "Which phenotypes are present in islet core vs peri-islet zones, and which are enriched vs tissue-wide background?"
        ))),
        fluidRow(style = "display: flex; flex-wrap: wrap",
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box",
              div(style = "display: flex; gap: 15px; align-items: center; flex-wrap: wrap; margin-bottom: 8px",
                radioButtons(ns("enrich_region"), "Region",
                             c("Peri-islet (enrichment z)" = "peri",
                               "Core (proportion)" = "core",
                               "Peri-islet (proportion)" = "peri_prop"),
                             selected = "peri", inline = TRUE),
                radioButtons(ns("enrich_stat"), "Summary", c("Median", "Mean"),
                             selected = "Median", inline = TRUE),
                conditionalPanel(
                  condition = paste0("input['", ns("enrich_region"), "'] == 'peri'"),
                  checkboxInput(ns("enrich_clip"), "Clip extreme z > 5", value = TRUE)
                )
              ),
              plotlyOutput(ns("enrichment_bars"), height = "420px")
            )
          ),
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%; box-sizing: border-box",
              h5("Heatmap", style = "font-size: 15px; margin-top: 0;"),
              plotlyOutput(ns("enrichment_heatmap"), height = "420px")
            )
          )
        ),
        fluidRow(column(12,
          div(style = "padding: 10px 15px; margin-bottom: 20px; font-size: 13px",
            tags$em("Enrichment z-scores (peri-islet mode) compare observed immune cell counts in the peri-islet zone vs expected counts based on tissue-wide proportions (Poisson model). z > 0 = enriched; z < 0 = depleted. Core and peri-islet proportion modes show raw cell type fractions (0\u20131) within each compartment.")
          )
        )),

        # ==== Card C: Immune Proximity to Islets ====
        fluidRow(column(12, sec_heading(
          "C", "Phenotype Proximity to Islet",
          "Minimum distance from islet core centroid to nearest peri-islet cells of each phenotype, plus signed-distance KDE for selected types."
        ))),
        fluidRow(style = "display: flex; flex-wrap: wrap",
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box",
              # Populated server-side: "Immune (all)" + every phenotype with min_dist_*.
              selectInput(ns("distance_metric"), "Distance to nearest",
                          choices = NULL, width = "100%"),
              uiOutput(ns("distance_pheno_hint")),
              checkboxInput(ns("distance_show_outlier_table"),
                            "Show excluded-outlier table", value = FALSE),
              tags$small(
                "Outlier exclusion follows the global toggle above (robust |z| from median/MAD on log(1+distance), per donor status — two-sided; threshold from the slider). Excluded points still appear as faint × markers."),
              plotlyOutput(ns("distance_boxplot"), height = "420px"),
              uiOutput(ns("distance_outlier_info"))
            )
          ),
          column(6,
            div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%; box-sizing: border-box",
              # Per-cell distance-to-islet, two views:
              #   * Density   — one population (a phenotype, a RESTORE-marker+ subset,
              #                 or all cells) as a signed-distance density, by status.
              #   * Composition — several RESTORE markers as grouped-stacked bars, to
              #                 see which markers sit at which distance from the islet.
              # This card has its OWN disease-status toggle (input$dist_status),
              # independent of the sidebar Donor Status filter.
              radioButtons(ns("dist_ptype"), "Plot",
                           choices = if (cell_distance_has_markers())
                                       c("Density curves" = "density", "Composition bars" = "composition")
                                     else c("Density curves" = "density"),
                           selected = "density", inline = TRUE),
              conditionalPanel(
                condition = "input.dist_ptype == 'density'", ns = ns,
                fluidRow(
                  column(5, radioButtons(ns("dist_mode"), "Group cells by",
                                         choices = if (cell_distance_has_markers())
                                                     c("Phenotype" = "phenotype", "Marker" = "marker")
                                                   else c("Phenotype" = "phenotype"),
                                         selected = "phenotype", inline = TRUE)),
                  column(7, selectInput(ns("dist_group"), "Show", choices = NULL, width = "100%"))
                )
              ),
              conditionalPanel(
                condition = "input.dist_ptype == 'composition'", ns = ns,
                if (cell_distance_has_markers())
                  checkboxGroupInput(ns("dist_comp_markers"), "Markers (RESTORE-positive)",
                                     choices = cell_distance_markers(),
                                     selected = intersect(c("INS", "CD3e", "CD20", "CD163", "CD31"),
                                                          cell_distance_markers()),
                                     inline = TRUE)
                else
                  tags$small(
                    "This dataset carries no per-marker positivity call, so composition bars are unavailable. Use the Density view, grouped by phenotype.")
              ),
              checkboxGroupInput(ns("dist_status"), "Disease status",
                                 choices = c("ND", "Aab+", "T1D"),
                                 selected = c("ND", "Aab+", "T1D"), inline = TRUE),
              uiOutput(ns("kde_pheno_hint")),
              plotlyOutput(ns("immune_distance_kde"), height = "520px")
            )
          )
        ),
        fluidRow(column(12,
          div(style = "padding: 10px 15px; margin-bottom: 20px; font-size: 13px",
            tags$em("Left: Distance metrics measure minimum Euclidean distance (\u00b5m) from islet core centroid to nearest peri-islet immune cells. NAs indicate no cells of that type in the peri-islet zone. Right: KDE of signed distance from islet boundary for individual immune cells. Negative = inside islet (core), positive = outside (peri-islet zone). The dashed line at zero marks the islet boundary.")
          )
        ))
      )
    })

    # ---- Islet count display ----
    output$nbr_islet_count <- renderUI({
      comp <- nbr_comp()
      n <- if (!is.null(comp)) nrow(comp) else 0
      span(style = "font-size: 13px; white-space: nowrap",
           paste0(formatC(n, big.mark = ","), " islets"))
    })

    # ---- Card A: excluded-outlier table ----
    # Mirrors the Plot tab's plot_outlier_info pattern: yellow warning banner
    # with header counts, body explanation, and a scrollable table. Three
    # display states — table, no-outliers note, or filter-OFF flagged note.
    output$infiltration_outlier_info <- renderUI({
      if (!isTRUE(input$infiltration_show_outlier_table)) return(NULL)
      tbl <- infiltration_outliers()
      filter_on <- isTRUE(remove_outliers())

      if (is.null(tbl) || nrow(tbl) == 0) {
        return(tags$div(
          style = "margin-top: 10px; padding: 8px 10px; font-size: 12px",
          "No outliers detected for the current phenotype."
        ))
      }

      header <- if (filter_on) {
        sprintf("⚠️ %d Outlier%s Excluded (Tukey 1.5×IQR per donor status, zero-aware)",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""))
      } else {
        sprintf("ⓘ %d Outlier%s Flagged (global filter OFF — still in summary)",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""))
      }
      body_text <- if (filter_on) {
        "The following islets were excluded from the donor-level enrichment summary because their value fell beyond the Tukey fences (Q1 − 1.5·IQR, Q3 + 1.5·IQR) of the non-zero distribution for their donor status:"
      } else {
        "These points would be excluded if the global outlier filter were turned on; they currently contribute to the donor-level summary:"
      }

      tags$div(
        style = "margin-top: 10px; padding: 10px",
        tags$h6(style = "margin-top: 0", header),
        tags$p(style = "font-size: 12px; margin-bottom: 10px", body_text),
        tags$div(style = "max-height: 220px; overflow-y: auto",
          renderTable(tbl, striped = TRUE, hover = TRUE, bordered = TRUE,
                      spacing = "xs", width = "100%")
        )
      )
    })

    # ---- Card C: excluded-outlier table ----
    output$distance_outlier_info <- renderUI({
      if (!isTRUE(input$distance_show_outlier_table)) return(NULL)
      tbl <- distance_outliers()
      filter_on <- isTRUE(remove_outliers())

      if (is.null(tbl) || nrow(tbl) == 0) {
        return(tags$div(
          style = "margin-top: 10px; padding: 8px 10px; font-size: 12px",
          "No outliers detected for the current distance metric."
        ))
      }

      thr <- suppressWarnings(as.numeric(outlier_threshold()))
      if (!is.finite(thr)) thr <- 3
      header <- if (filter_on) {
        sprintf("⚠️ %d Outlier%s Excluded (|z| > %g per donor status on log(1+distance))",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""), thr)
      } else {
        sprintf("ⓘ %d Outlier%s Flagged (global filter OFF — still in box plot)",
                nrow(tbl), ifelse(nrow(tbl) > 1, "s", ""))
      }
      body_text <- if (filter_on) {
        sprintf("The following islets were excluded from the box plot because their robust |z| on log(1 + distance) exceeded %g within their donor status (both high- and low-side tails). They still render on the plot as faint grey × markers:", thr)
      } else {
        "These points would be excluded if the global outlier filter were turned on; they currently sit in the box plot's distribution:"
      }

      tags$div(
        style = "margin-top: 10px; padding: 10px",
        tags$h6(style = "margin-top: 0", header),
        tags$p(style = "font-size: 12px; margin-bottom: 10px", body_text),
        tags$div(style = "max-height: 220px; overflow-y: auto",
          renderTable(tbl, striped = TRUE, hover = TRUE, bordered = TRUE,
                      spacing = "xs", width = "100%")
        )
      )
    })

    # ---- Populate dynamic phenotype dropdowns once data is available ----
    # Replaces the previous hand-curated 7-immune-type lists. Choices are
    # built from prepared()$comp columns (enrich_z_*, prop_*, peri_prop_*,
    # min_dist_*) so every phenotype the pipeline computed metrics for is
    # surfaced uniformly. Defaults pick a sensible immune phenotype when
    # available (Macrophage), else first.
    pick_default <- function(choices, prefer = c("Macrophage", "Immune", "T_cell", "CD8a_Tcell")) {
      if (length(choices) == 0) return(character(0))
      for (p in prefer) if (p %in% choices) return(p)
      choices[[1]]
    }

    observe({
      ch <- phenotype_choices()
      if (length(ch) == 0) return()
      # Card A-Left: enrichment z-score per phenotype
      cur_a <- isolate(input$infiltration_phenotype)
      sel_a <- if (!is.null(cur_a) && cur_a %in% ch) cur_a else pick_default(ch)
      updateSelectInput(session, "infiltration_phenotype", choices = ch, selected = sel_a)
      # Card A-Right: peri vs core proportion (only phenotypes with both prop_ + peri_prop_)
      opts <- phenotype_options()
      paired <- vapply(opts, function(o) !is.na(o$prop_core) && !is.na(o$peri_prop), logical(1))
      ch_b <- ch[paired]
      cur_b <- isolate(input$scatter_phenotype)
      sel_b <- if (!is.null(cur_b) && cur_b %in% ch_b) cur_b else pick_default(ch_b)
      updateSelectInput(session, "scatter_phenotype", choices = ch_b, selected = sel_b)
      # Card C-Left: distance metric — "Immune (all)" plus min_dist_<phenotype>
      dist_opts <- Filter(function(o) !is.na(o$min_dist), opts)
      dist_ch <- c("Immune (all)" = "min_dist_immune_mean",
                   stats::setNames(
                     vapply(dist_opts, function(o) o$min_dist, character(1)),
                     vapply(dist_opts, function(o) o$label, character(1))
                   ))
      cur_c <- isolate(input$distance_metric)
      sel_c <- if (!is.null(cur_c) && cur_c %in% dist_ch) cur_c else "min_dist_immune_mean"
      updateSelectInput(session, "distance_metric", choices = dist_ch, selected = sel_c)

      # Card C-Right: the per-cell distance plot's `dist_group` dropdown is
      # populated by its own observeEvent(input$dist_mode) below — phenotype vs
      # marker lists come from the DuckDB views, not phenotype_choices().
    })

    # Card C-Right (Density view): populate the `dist_group` dropdown, switching
    # between phenotype labels and the 10 RESTORE marker channels (both from the
    # cell_distance view) with an "All (all cells)" option at the top. observeEvent
    # on input$dist_mode only fires once the Card C renderUI has mounted the input
    # on the client, so it cannot race the mount.
    observeEvent(input$dist_mode, {
      mode <- input$dist_mode %||% "phenotype"
      items <- if (identical(mode, "marker")) cell_distance_markers() else cell_distance_phenotypes()
      ch <- c("All (all cells)" = "__all__", stats::setNames(items, items))
      cur <- isolate(input$dist_group)
      sel <- if (!is.null(cur) && cur %in% ch) cur else "__all__"
      updateSelectInput(session, "dist_group", choices = ch, selected = sel)
    }, ignoreNULL = FALSE)

    # ---- Inline marker-rule hints under the four phenotype selectors ----
    # Each looks up the selected phenotype's label via phenotype_options() and
    # renders phenotype_hint_ui(). "Immune (all)" aggregates have no single
    # rule — we display a description instead.
    .lookup_pheno_label <- function(safe_or_metric) {
      if (is.null(safe_or_metric) || !nzchar(safe_or_metric)) return(NULL)
      if (safe_or_metric == "min_dist_immune_mean" || safe_or_metric == "all") {
        return("Immune (all)")
      }
      opts <- phenotype_options()
      # Direct match by safe-name
      for (o in opts) if (identical(o$safe, safe_or_metric)) return(o$label)
      # Match by min_dist column
      for (o in opts) if (identical(o$min_dist, safe_or_metric)) return(o$label)
      # Fallback: assume it's already a label
      safe_or_metric
    }
    .hint_for <- function(label, inline_id) {
      if (is.null(label)) return(NULL)
      if (identical(label, "Immune (all)")) {
        return(shiny::tags$div(
          style = "margin-top: -6px; margin-bottom: 8px; padding: 6px 10px; font-size: 13px; line-height: 1.4",
          shiny::tags$strong("Markers: "),
          shiny::tags$em("Aggregate of all immune phenotypes (any of CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR positive)."),
          shiny::tags$span(" · "),
          shiny::actionLink(inline_id, "see all rules →",
                            style = "font-size: 12px")
        ))
      }
      phenotype_hint_ui(label, show_modal_link = inline_id)
    }
    output$infiltration_pheno_hint <- renderUI({
      .hint_for(.lookup_pheno_label(input$infiltration_phenotype),
                ns("show_phenotype_rules_inline_a"))
    })
    output$scatter_pheno_hint <- renderUI({
      .hint_for(.lookup_pheno_label(input$scatter_phenotype),
                ns("show_phenotype_rules_inline_scatter"))
    })
    output$distance_pheno_hint <- renderUI({
      .hint_for(.lookup_pheno_label(input$distance_metric),
                ns("show_phenotype_rules_inline_dist"))
    })
    output$kde_pheno_hint <- renderUI({
      box_style <- "margin-top: -6px; margin-bottom: 8px; padding: 6px 10px; font-size: 13px; line-height: 1.4"
      if (identical(input$dist_ptype %||% "density", "composition")) {
        return(shiny::tags$div(style = box_style,
          shiny::tags$strong("Composition: "),
          shiny::tags$em("stacked count of RESTORE-positive cells per marker at each distance bin (a cell may be positive for more than one marker).")))
      }
      mode <- input$dist_mode %||% "phenotype"
      grp  <- input$dist_group %||% "__all__"
      if (is.null(grp) || !nzchar(grp)) grp <- "__all__"
      if (identical(grp, "__all__")) {
        return(shiny::tags$div(style = box_style,
          shiny::tags$em("All cells (every phenotype) — signed distance to the nearest islet, split by disease status.")))
      }
      if (identical(mode, "marker")) {
        return(shiny::tags$div(style = box_style,
          shiny::tags$strong("Marker: "),
          shiny::tags$em(paste0(grp, "-positive cells by the RESTORE per-image threshold (calibrated on REDSEA-corrected MFI)."))))
      }
      # phenotype mode: the dropdown value IS the phenotype label from cell_distance
      tryCatch(phenotype_hint_ui(grp, show_modal_link = ns("show_phenotype_rules_inline_kde")),
               error = function(e) NULL)
    })

    # ==== Card A-Left: Per-phenotype peri-zone enrichment z-score (donor-level) ====
    output$infiltration_bars <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp), "Case ID" %in% colnames(comp))
      opts <- phenotype_options()
      sel_safe <- input$infiltration_phenotype
      req(sel_safe)
      itm <- Filter(function(o) o$safe == sel_safe, opts)
      req(length(itm) > 0)
      itm <- itm[[1]]
      metric <- itm$enrich
      req(metric %in% colnames(comp))
      metric_label <- paste0("Peri Enrichment z-score — ", itm$label)

      df <- data.frame(
        value = comp[[metric]],
        status = comp$`Donor Status`,
        case_id = comp$`Case ID`,
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$value), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = "No data"))

      # Outlier handling via shared helper. Phenotype enrichment z-scores are
      # zero-inflated for phenotypes absent from many islets, so we keep the
      # Tukey IQR rule with zero-aware fences (zeros never excluded).
      df <- flag_outliers(df, value_col = "value", group_col = "status",
                          method = "iqr", threshold = 1.5,
                          exclude_zeros = TRUE,
                          enabled = isTRUE(remove_outliers()))
      filter_summary <- summarize_outlier_filter(df)
      filter_summary$label <- if (filter_summary$enabled) "Tukey 1.5×IQR per donor status (zero-aware)" else filter_summary$label
      df_full <- df

      # Populate the per-card outlier table reactiveVal (consumed by the
      # `infiltration_outlier_info` renderUI below).
      flagged <- df_full[isTRUE(attr(df_full, "outlier_enabled")) & df_full$is_outlier, , drop = FALSE]
      if (!isTRUE(attr(df_full, "outlier_enabled"))) {
        flagged <- df_full[df_full$is_outlier, , drop = FALSE]
      }
      if (nrow(flagged) > 0) {
        infiltration_outliers(data.frame(
          Phenotype          = itm$label,
          Case_ID            = flagged$case_id,
          Donor_Status       = flagged$status,
          Value              = round(flagged$value, 3),
          Fence_Distance_IQR = round(flagged$outlier_score, 2),
          stringsAsFactors   = FALSE
        ))
      } else {
        infiltration_outliers(NULL)
      }

      if (isTRUE(remove_outliers())) {
        df <- df[!df$is_outlier | is.na(df$is_outlier), , drop = FALSE]
      }
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = "No data after outlier removal"))

      # Aggregate to donor-level means
      donor_df <- df %>%
        dplyr::group_by(case_id, status) %>%
        dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

      status_order <- c("ND", "Aab+", "T1D")
      donor_df$status <- factor(donor_df$status, levels = intersect(status_order, unique(donor_df$status)))

      # Donor-level summary: mean +/- SEM per group
      summary_df <- donor_df %>%
        dplyr::group_by(status) %>%
        dplyr::summarise(
          mean_val = mean(value, na.rm = TRUE),
          sd_val = sd(value, na.rm = TRUE),
          n = dplyr::n(),
          sem = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
          .groups = "drop"
        )

      dcols <- donor_colors_reactive()

      # Kruskal-Wallis on donor-level means (N=15)
      kw_p <- tryCatch({
        kt <- kruskal.test(value ~ status, data = donor_df)
        kt$p.value
      }, error = function(e) NA_real_)
      p_label <- if (is.finite(kw_p)) paste0("KW p = ", formatC(kw_p, format = "g", digits = 3)) else ""

      p <- plot_ly()
      for (s in levels(summary_df$status)) {
        row <- summary_df[summary_df$status == s, , drop = FALSE]
        if (nrow(row) == 0) next
        # Individual donor points
        donor_pts <- donor_df[donor_df$status == s, , drop = FALSE]
        p <- p %>% add_trace(
          x = s, y = donor_pts$value,
          type = "scatter", mode = "markers",
          marker = list(size = 8, opacity = 0.6, color = dcols[s]),
          name = s, showlegend = FALSE,
          hoverinfo = "text",
          text = paste0("Donor ", donor_pts$case_id, ": ", round(donor_pts$value, 4))
        )
        # Bar with SEM error bars
        p <- p %>% add_trace(
          x = s, y = row$mean_val,
          type = "bar",
          marker = list(color = paste0(dcols[s], "66"),
                        line = list(color = dcols[s], width = 1.5)),
          error_y = list(type = "data", array = row$sem, visible = TRUE,
                         color = dcols[s], thickness = 1.5, width = 6),
          name = s, showlegend = FALSE,
          hoverinfo = "text",
          text = paste0(s, "<br>Mean: ", round(row$mean_val, 4),
                       "<br>SEM: ", round(row$sem, 4),
                       "<br>N donors: ", row$n)
        )
      }
      excl_text <- if (filter_summary$enabled && filter_summary$n_excluded > 0) {
        sprintf(" | %d excluded (%s)", filter_summary$n_excluded, filter_summary$label)
      } else if (!filter_summary$enabled) {
        " | outlier filter OFF"
      } else {
        ""
      }
      fonts <- plotly_islet_fonts()
      p %>% layout(
        font = fonts$global,
        title = list(text = paste0(metric_label,
                                   "<br><sup>Donor-level means \u00b1 SEM (N=donors) | ", p_label, excl_text, "</sup>"),
                     font = fonts$title),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick,
                     categoryorder = "array",
                     categoryarray = c("ND", "Aab+", "T1D")),
        yaxis = list(title = list(text = metric_label, font = fonts$axis_title),
                     tickfont = fonts$axis_tick),
        barmode = "overlay",
        showlegend = FALSE
      ) %>%
        plotly_islet_config(paste0("spatial_infiltration_", itm$safe))
    })

    # ==== Card A-Right: Peri vs Core Proportion Scatter ====
    output$infiltration_scatter <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))

      opts <- phenotype_options()
      sel_safe <- input$scatter_phenotype
      req(sel_safe)
      itm <- Filter(function(o) o$safe == sel_safe, opts)
      req(length(itm) > 0)
      itm <- itm[[1]]
      cols <- list(peri = itm$peri_prop, core = itm$prop_core,
                   label = paste0(itm$label, " Proportion"))
      req(!is.na(cols$peri), !is.na(cols$core))
      req(cols$peri %in% colnames(comp), cols$core %in% colnames(comp))

      use_sqrt <- isTRUE(input$scatter_sqrt)

      df <- data.frame(
        peri = comp[[cols$peri]],
        core = comp[[cols$core]],
        status = comp$`Donor Status`,
        islet_key = if ("islet_key" %in% colnames(comp)) comp$islet_key else "",
        stringsAsFactors = FALSE
      )
      df <- df[is.finite(df$peri) & is.finite(df$core), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = "No data"))

      # Apply sqrt transform if requested
      if (use_sqrt) {
        df$peri_plot <- sqrt(df$peri)
        df$core_plot <- sqrt(df$core)
      } else {
        df$peri_plot <- df$peri
        df$core_plot <- df$core
      }

      status_order <- c("ND", "Aab+", "T1D")
      df$status <- factor(df$status, levels = intersect(status_order, unique(df$status)))
      dcols <- donor_colors_reactive()

      max_val <- max(c(df$peri_plot, df$core_plot), na.rm = TRUE) * 1.05

      pt_sz <- input$nbr_pt_size %||% 5
      pt_al <- 1 - (input$nbr_pt_alpha %||% 0.6)

      p <- plot_ly()
      for (s in levels(df$status)) {
        sub <- df[df$status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          data = sub, x = ~peri_plot, y = ~core_plot,
          text = ~paste0(islet_key, "<br>Status: ", status,
                        "<br>Peri: ", round(peri, 4),
                        "<br>Core: ", round(core, 4)),
          hoverinfo = "text", name = s,
          type = "scatter", mode = "markers",
          marker = list(size = pt_sz, opacity = pt_al, color = dcols[s])
        )
      }

      # Optional linear trend lines per status group
      if (isTRUE(input$scatter_trend)) {
        for (s in levels(df$status)) {
          sub <- df[df$status == s, , drop = FALSE]
          if (nrow(sub) < 3) next
          fit <- lm(core_plot ~ peri_plot, data = sub)
          pred_x <- seq(min(sub$peri_plot, na.rm = TRUE),
                        max(sub$peri_plot, na.rm = TRUE), length.out = 50)
          pred_y <- predict(fit, newdata = data.frame(peri_plot = pred_x))
          p <- p %>% add_trace(
            x = pred_x, y = pred_y,
            type = "scatter", mode = "lines",
            line = list(color = dcols[s], width = 2),
            name = paste0(s, " trend"),
            showlegend = FALSE,
            hoverinfo = "none"
          )
        }
      }

      # Axis labels
      scale_note <- if (use_sqrt) "\u221a" else ""
      x_title <- paste0(scale_note, cols$label, " (Peri)")
      y_title <- paste0(scale_note, cols$label, " (Core)")

      # Diagonal y=x reference
      fonts <- plotly_islet_fonts()
      p %>% add_trace(
        x = c(0, max_val), y = c(0, max_val),
        type = "scatter", mode = "lines",
        line = list(dash = "dash", color = "#999", width = 1),
        showlegend = FALSE, hoverinfo = "none"
      ) %>% layout(
        font = fonts$global,
        title = list(text = paste0("Peri vs Core: ", cols$label), font = fonts$title),
        xaxis = list(title = list(text = x_title, font = fonts$axis_title),
                     tickfont = fonts$axis_tick, range = c(0, max_val)),
        yaxis = list(title = list(text = y_title, font = fonts$axis_title),
                     tickfont = fonts$axis_tick, range = c(0, max_val)),
        legend = list(title = list(text = "Status", font = fonts$legend_title),
                      font = fonts$legend_text)
      ) %>%
        plotly_islet_config(paste0("spatial_infiltration_scatter_", itm$safe))
    })

    # ==== Card B-Left: Enrichment/Proportion Grouped Bar Chart ====
    output$enrichment_bars <- renderPlotly({
      es <- enrich_summary()
      req(es)

      status_order <- c("ND", "Aab+", "T1D")
      es$donor_status <- factor(es$donor_status, levels = intersect(status_order, unique(es$donor_status)))

      dcols <- donor_colors_reactive()

      stat_mode <- input$enrich_stat %||% "Median"
      region <- input$enrich_region %||% "peri"
      error_label <- if (stat_mode == "Median") "IQR" else "SEM"

      # Axis label and title depend on region
      if (region == "peri") {
        y_label <- "Enrichment z-score"
        title_prefix <- paste0(stat_mode, " Enrichment z-score")
      } else if (region == "core") {
        y_label <- "Proportion (core)"
        title_prefix <- paste0(stat_mode, " Core Proportion")
      } else {
        y_label <- "Proportion (peri-islet)"
        title_prefix <- paste0(stat_mode, " Peri-islet Proportion")
      }

      p <- plot_ly()
      for (s in levels(es$donor_status)) {
        sub <- es[es$donor_status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          x = sub$cell_type, y = sub$z_summary, name = s,
          type = "bar",
          marker = list(color = dcols[s]),
          error_y = list(
            type = "data",
            symmetric = FALSE,
            array = sub$z_hi - sub$z_summary,
            arrayminus = sub$z_summary - sub$z_lo,
            color = "#666", thickness = 1
          ),
          hovertemplate = paste0("%{x}<br>", s, "<br>", y_label, " = %{y:.3f}<extra></extra>")
        )
      }

      # Reference line at 0 for enrichment z-scores only
      shapes <- if (region == "peri") {
        list(list(type = "line", x0 = -0.5, x1 = length(unique(es$cell_type)) - 0.5, y0 = 0, y1 = 0,
                  line = list(color = "#999", dash = "dash", width = 1)))
      } else list()

      fonts <- plotly_islet_fonts()
      p %>% layout(
        font = fonts$global,
        barmode = "group",
        title = list(text = paste0(title_prefix, "<br><sup>Donor-level means, ", error_label, " (N=donors)</sup>"),
                     font = fonts$title),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, tickangle = -30,
                     categoryorder = "array",
                     categoryarray = unique(es$cell_type)),
        yaxis = list(title = list(text = y_label, font = fonts$axis_title),
                     tickfont = fonts$axis_tick),
        legend = list(title = list(text = "Status", font = fonts$legend_title),
                      font = fonts$legend_text),
        shapes = shapes
      ) %>%
        plotly_islet_config("spatial_enrichment_bars")
    })

    # ==== Card B-Right: Enrichment/Proportion Heatmap ====
    output$enrichment_heatmap <- renderPlotly({
      es <- enrich_summary()
      req(es)

      region <- input$enrich_region %||% "peri"

      # Pivot to matrix: cell_type x donor_status
      status_order <- c("ND", "Aab+", "T1D")
      cell_types <- unique(es$cell_type)
      statuses <- intersect(status_order, unique(es$donor_status))

      mat <- matrix(NA_real_, nrow = length(statuses), ncol = length(cell_types),
                    dimnames = list(statuses, cell_types))
      for (i in seq_len(nrow(es))) {
        s <- es$donor_status[i]
        ct <- es$cell_type[i]
        if (as.character(s) %in% rownames(mat) && ct %in% colnames(mat))
          mat[as.character(s), ct] <- es$z_summary[i]
      }

      # Annotation text: more decimals for proportions (small numbers)
      fmt <- if (region == "peri") "%.2f" else "%.4f"
      text_mat <- matrix(sprintf(fmt, mat), nrow = nrow(mat), ncol = ncol(mat))
      text_mat[is.na(mat)] <- ""

      if (region == "peri") {
        # Diverging for z-scores (symmetric around 0)
        z_abs_max <- max(abs(mat), na.rm = TRUE)
        if (!is.finite(z_abs_max) || z_abs_max == 0) z_abs_max <- 1
        colorscale <- list(c(0, "#2166AC"), c(0.5, "#FFFFFF"), c(1, "#B2182B"))
        zmin <- -z_abs_max; zmax <- z_abs_max
        cb_title <- "z-score"
        title_text <- paste0(input$enrich_stat %||% "Median", " Enrichment z-score")
      } else {
        # Sequential for proportions (0 to max)
        zmax <- max(mat, na.rm = TRUE)
        if (!is.finite(zmax) || zmax == 0) zmax <- 0.1
        colorscale <- list(c(0, "#FFFFFF"), c(0.5, "#FDB863"), c(1, "#B2182B"))
        zmin <- 0
        cb_title <- "Proportion"
        region_label <- if (region == "core") "Core" else "Peri-islet"
        title_text <- paste0(input$enrich_stat %||% "Median", " ", region_label, " Proportion")
      }

      fonts <- plotly_islet_fonts()
      plot_ly(
        x = colnames(mat), y = rownames(mat), z = mat,
        type = "heatmap",
        colorscale = colorscale,
        zmin = zmin, zmax = zmax,
        text = text_mat, texttemplate = "%{text}",
        textfont = list(size = ISLET_FONT_SIZES$dense_min),
        hovertemplate = paste0("%{x}<br>%{y}<br>", cb_title, " = %{z:.4f}<extra></extra>"),
        showscale = TRUE,
        colorbar = list(title = list(text = cb_title, font = fonts$colorbar_title),
                        tickfont = fonts$axis_tick)
      ) %>% layout(
        font = fonts$global,
        title = list(text = title_text, font = fonts$title),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, tickangle = -30),
        yaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, autorange = "reversed")
      ) %>%
        plotly_islet_config("spatial_enrichment_heatmap")
    })

    # ==== Card C-Left: Distance Box Plots ====
    output$distance_boxplot <- renderPlotly({
      comp <- nbr_comp()
      req(comp, "Donor Status" %in% colnames(comp))
      metric <- input$distance_metric %||% "min_dist_immune_mean"
      req(metric %in% colnames(comp))

      # Build label dynamically: aggregate "Immune (all)" or per-phenotype label
      if (metric == "min_dist_immune_mean") {
        metric_label <- "Min Distance to Immune (all)"
      } else {
        opts <- phenotype_options()
        match <- Filter(function(o) !is.na(o$min_dist) && o$min_dist == metric, opts)
        metric_label <- if (length(match) > 0)
          paste0("Min Distance to ", match[[1]]$label) else metric
      }

      df <- data.frame(
        value = comp[[metric]],
        status = comp$`Donor Status`,
        case_id = if ("Case ID" %in% colnames(comp)) comp$`Case ID` else NA,
        islet_key = if ("islet_key" %in% colnames(comp)) comp$islet_key else NA_character_,
        stringsAsFactors = FALSE
      )

      # Count non-NA per group before filtering
      status_order <- c("ND", "Aab+", "T1D")
      n_total <- tapply(df$value, df$status, length)
      n_valid <- tapply(df$value, df$status, function(x) sum(is.finite(x)))
      na_rate <- 1 - sum(is.finite(df$value)) / nrow(df)

      df <- df[is.finite(df$value), , drop = FALSE]
      if (nrow(df) == 0) return(plotly_empty() %>% layout(title = paste0(metric_label, " \u2014 all NA")))

      # Outlier handling via shared helper, applied to log1p(distance) to
      # tame the long right tail; ghost \u00d7 markers preserve visibility.
      df <- flag_outliers(df, value_col = "value", group_col = "status",
                          method = "zscore", threshold = outlier_threshold(),
                          transform = "log1p",
                          enabled = isTRUE(remove_outliers()))
      filter_summary <- summarize_outlier_filter(df)
      filter_summary$label <- if (filter_summary$enabled) {
        sprintf("|z| > %g per donor status on log(1+distance)",
                suppressWarnings(as.numeric(outlier_threshold())))
      } else filter_summary$label

      # Populate the per-card outlier table reactiveVal (consumed by the
      # `distance_outlier_info` renderUI below).
      flagged <- df[df$is_outlier, , drop = FALSE]
      if (nrow(flagged) > 0) {
        distance_outliers(data.frame(
          Phenotype     = metric_label,
          Case_ID       = flagged$case_id,
          Islet         = flagged$islet_key,
          Donor_Status  = flagged$status,
          Distance_um   = round(flagged$value, 1),
          Z_Score_log1p = round(flagged$outlier_score, 2),
          stringsAsFactors = FALSE
        ))
      } else {
        distance_outliers(NULL)
      }

      df_ghost <- df[df$is_outlier, , drop = FALSE]
      df_box   <- if (isTRUE(remove_outliers())) df[!df$is_outlier, , drop = FALSE] else df

      df_box$status <- factor(df_box$status, levels = intersect(status_order, unique(df_box$status)))
      if (nrow(df_ghost) > 0) {
        df_ghost$status <- factor(df_ghost$status, levels = levels(df_box$status))
      }
      dcols <- donor_colors_reactive()
      pt_sz <- input$nbr_pt_size %||% 5
      pt_al <- 1 - (input$nbr_pt_alpha %||% 0.6)

      # Subtitle with N per group
      n_labels <- sapply(levels(df_box$status), function(s) {
        nv <- n_valid[s]
        nt <- n_total[s]
        if (is.na(nv)) nv <- 0
        if (is.na(nt)) nt <- 0
        paste0(s, ": ", nv, "/", nt)
      })
      subtitle <- paste0("Non-NA: ", paste(n_labels, collapse = ", "),
                          " (", round(na_rate * 100, 0), "% NA overall)")
      excl_suffix <- if (filter_summary$enabled && filter_summary$n_excluded > 0) {
        sprintf(" | %d excluded (%s)", filter_summary$n_excluded, filter_summary$label)
      } else if (!filter_summary$enabled) {
        " | outlier filter OFF"
      } else {
        ""
      }

      p <- plot_ly()
      for (s in levels(df_box$status)) {
        sub <- df_box[df_box$status == s, , drop = FALSE]
        if (nrow(sub) == 0) next
        p <- p %>% add_trace(
          y = sub$value, x = s, name = s,
          type = "box",
          boxpoints = "all", jitter = 0.3, pointpos = -1.5,
          marker = list(size = pt_sz, opacity = pt_al, color = dcols[s]),
          line = list(color = dcols[s]),
          fillcolor = paste0(dcols[s], "44"),
          hoverinfo = "y"
        )
      }
      # Ghost layer: excluded outliers as faint grey \u00d7 markers per group
      if (nrow(df_ghost) > 0) {
        for (s in levels(df_box$status)) {
          gsub_ <- df_ghost[df_ghost$status == s, , drop = FALSE]
          if (nrow(gsub_) == 0) next
          p <- p %>% add_trace(
            y = gsub_$value, x = rep(s, nrow(gsub_)),
            type = "scatter", mode = "markers",
            marker = list(symbol = "x", size = pt_sz * 1.2,
                          color = "rgba(80,80,80,0.55)",
                          line = list(color = "rgba(50,50,50,0.7)", width = 0.8)),
            showlegend = FALSE,
            hoverinfo = "text",
            text = paste0("Excluded outlier \u00b7 z = ",
                          sprintf("%.2f", gsub_$outlier_score),
                          " \u00b7 distance = ", round(gsub_$value, 1), " \u00b5m"),
            name = paste0(s, " (excluded)")
          )
        }
      }
      fonts <- plotly_islet_fonts()
      p %>% layout(
        font = fonts$global,
        title = list(text = paste0(metric_label, "<br><sup>", subtitle, excl_suffix, "</sup>"),
                     font = fonts$title),
        xaxis = list(title = list(text = "", font = fonts$axis_title),
                     tickfont = fonts$axis_tick, categoryorder = "array",
                     categoryarray = c("ND", "Aab+", "T1D")),
        yaxis = list(title = list(text = paste0("Distance (\u00b5m)"), font = fonts$axis_title),
                     tickfont = fonts$axis_tick),
        showlegend = FALSE
      ) %>%
        plotly_islet_config(paste0("spatial_distance_box_", gsub("[^A-Za-z0-9_]", "_", metric)))
    })

    # ==== Card C-Right: per-cell distance-to-islet (DuckDB) ====
    # Two views over the `cell_distance` view (query_cell_distance_hist /
    # query_cell_distance_composition in spatial_helpers.R), aggregated in SQL so
    # only a tiny result crosses into R. Signed distance: negative = inside islet,
    # positive = outside. Marker positivity = RESTORE `_pos` (per-image gated).
    # Uses this card's OWN input$dist_status, NOT the sidebar.
    cell_dist_data <- reactive({
      query_cell_distance_hist(
        mode     = input$dist_mode %||% "phenotype",
        group    = input$dist_group %||% "__all__",
        statuses = input$dist_status %||% c("ND", "Aab+", "T1D")
      )
    })
    cell_comp_data <- reactive({
      query_cell_distance_composition(
        markers  = input$dist_comp_markers %||% character(0),
        statuses = input$dist_status %||% c("ND", "Aab+", "T1D")
      )
    })

    output$immune_distance_kde <- renderPlotly({
      fonts <- plotly_islet_fonts()
      empty_msg <- function(msg) plotly::plotly_empty(type = "scatter", mode = "markers") %>%
        layout(title = list(text = msg, font = fonts$title))
      if (!cell_distance_view_ready())
        return(empty_msg("Per-cell distance data not registered — run scripts/senior/build_cell_distance_parquet.py"))
      x_title <- "Signed distance (µm) — inside islet ← 0 → outside"

      # ---------- Composition: small-multiple stacked bars by status ----------
      # "Grouped stacked" = one clean stacked-marker bar panel PER disease status
      # (facets): stacked = markers, x = distance bin. plotly `offsetgroup`
      # side-by-side grouping renders muddy at these bar densities, so small
      # multiples (subplot rows) read far better and still compare across status.
      if (identical(input$dist_ptype %||% "density", "composition")) {
        res <- cell_comp_data()
        if (is.null(res) || is.null(res$comp) || nrow(res$comp) == 0)
          return(empty_msg("Pick at least one marker and one disease status."))
        comp <- res$comp
        markers <- res$markers
        status_order <- c("ND", "Aab+", "T1D")
        present <- intersect(status_order, unique(comp$donor_status))
        mcols <- CELL_DIST_MARKER_COLORS
        panels <- lapply(present, function(s) {
          ss <- comp[comp$donor_status == s, , drop = FALSE]
          pp <- plot_ly()
          for (m in markers) {
            sub <- ss[ss$marker == m, , drop = FALSE]
            sub <- sub[order(sub$binx), , drop = FALSE]
            pp <- pp %>% add_trace(
              x = sub$x, y = sub$n, type = "bar", textposition = "none",
              name = m, legendgroup = m, showlegend = identical(s, present[1]),
              marker = list(color = unname(mcols[m])),
              hoverinfo = "text",
              hovertext = paste0(s, " · ", m, "+<br>dist ≈ ", round(sub$x), " µm<br>",
                                 formatC(sub$n, big.mark = ","), " cells")
            )
          }
          # Islet-border line at 0 (a bar boundary sits exactly here, always shown).
          tot <- stats::aggregate(n ~ binx, ss, sum)
          ymax <- if (nrow(tot)) max(tot$n) else 1
          pp <- pp %>% add_trace(
            x = c(0, 0), y = c(0, ymax * 1.02), type = "scatter", mode = "lines",
            line = list(color = "#333", width = 1.4, dash = "dash"),
            showlegend = FALSE, hoverinfo = "none"
          )
          pp %>% layout(
            barmode = "stack",
            yaxis = list(title = list(text = s, font = fonts$axis_title), tickfont = fonts$axis_tick),
            xaxis = list(title = list(text = x_title, font = fonts$axis_title),
                         tickfont = fonts$axis_tick, zeroline = FALSE)
          )
        })
        fig <- if (length(panels) == 1) panels[[1]]
               else plotly::subplot(panels, nrows = length(panels), shareX = TRUE,
                                    titleY = TRUE, margin = 0.05)
        return(
          fig %>% layout(
            font = fonts$global, barmode = "stack",
            title = list(text = "Marker composition vs distance<br><sup>rows = disease status · y = RESTORE-positive cells</sup>",
                         font = fonts$title),
            legend = list(title = list(text = "Marker", font = fonts$legend_title),
                          font = fonts$legend_text)
          ) %>%
            plotly_islet_config("spatial_distance_composition")
        )
      }

      # ---------- Density: one population, split by status ----------
      res <- cell_dist_data()
      if (is.null(res)) return(empty_msg("Select at least one disease status."))
      mode  <- input$dist_mode %||% "phenotype"
      group <- input$dist_group %||% "__all__"
      if (!nzchar(group)) group <- "__all__"   # first render (before dropdown populates) sends ""
      type_label <- if (identical(group, "__all__")) "All cells"
        else if (identical(mode, "marker")) paste0(group, "+ cells")
        else group

      hist <- res$hist
      if (is.null(hist) || nrow(hist) == 0)
        return(empty_msg(paste0(type_label, " — no cells in this selection")))
      status_order <- c("ND", "Aab+", "T1D")
      present <- intersect(status_order, unique(as.character(hist$donor_status)))
      dcols <- donor_colors_reactive()
      y_max <- max(hist$dens, na.rm = TRUE)
      p <- plot_ly()
      for (s in present) {
        hs <- hist[hist$donor_status == s, , drop = FALSE]
        hs <- hs[order(hs$x), , drop = FALSE]
        n_cells <- res$n_by_status[[s]]
        p <- p %>% add_trace(
          x = hs$x, y = hs$dens, type = "scatter", mode = "lines",
          name = paste0(s, " (n=", formatC(n_cells, big.mark = ","), ")"),
          line = list(color = dcols[s], width = 2.5),
          fill = "tozeroy", fillcolor = paste0(dcols[s], "22"),
          hoverinfo = "text",
          text = paste0(s, "<br>dist = ", round(hs$real_x, 1), " µm<br>cells/µm = ", round(hs$dens, 2))
        )
      }
      p <- p %>% add_trace(
        x = c(0, 0), y = c(0, y_max * 1.05), type = "scatter", mode = "lines",
        line = list(color = "#333", width = 1.5, dash = "dash"),
        showlegend = FALSE, hoverinfo = "none"
      )
      ticks <- cd_symlog_ticks(res$lo, res$hi)   # real-µm labels at symlog positions
      p %>% layout(
        font = fonts$global,
        title = list(text = paste0(type_label, " — Distance from Islet"), font = fonts$title),
        xaxis = list(title = list(text = "Signed distance (µm, symlog) — inside islet ← 0 → outside",
                                  font = fonts$axis_title),
                     tickfont = fonts$axis_tick, zeroline = FALSE,
                     tickmode = "array", tickvals = ticks$vals, ticktext = ticks$text),
        yaxis = list(title = list(text = "Cells per µm of distance", font = fonts$axis_title),
                     tickfont = fonts$axis_tick),
        legend = list(title = list(text = "Status", font = fonts$legend_title),
                      font = fonts$legend_text)
      ) %>%
        plotly_islet_config(paste0("spatial_distance_", gsub("[^A-Za-z0-9_]", "_", type_label)))
    })

    # ---- Inline "see all rules" links under Card A/C phenotype hints ----
    for (id in c("show_phenotype_rules_inline_a",
                 "show_phenotype_rules_inline_scatter",
                 "show_phenotype_rules_inline_dist",
                 "show_phenotype_rules_inline_kde")) {
      local({
        nm <- id
        observeEvent(input[[nm]],
                     { showModal(phenotype_rules_modal_ui()) },
                     ignoreInit = TRUE)
      })
    }

  })
}
