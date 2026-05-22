trajectory_server <- function(id, prepared, selected_islet, forced_image, active_tab = reactive("Trajectory"),
                              donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                              remove_outliers = reactive(FALSE),
                              outlier_threshold = reactive(3.0)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---------- H5AD Path Resolution ----------
    resolve_traj_path <- function() {
      cand <- c(
        file.path("..", "..", "data", "adata_ins_root.h5ad"),
        file.path("..", "..", "scripts", "adata_ins_root.h5ad"),
        if (!is.null(project_root)) file.path(project_root, "data", "adata_ins_root.h5ad") else NULL,
        if (!is.null(project_root)) file.path(project_root, "scripts", "adata_ins_root.h5ad") else NULL,
        Sys.getenv("ADATA_PATH", unset = "")
      )
      cand <- cand[nzchar(cand)]
      for (p in cand) if (file.exists(p)) return(p)
      NULL
    }

    traj_path <- resolve_traj_path()

    # ---------- Reactive state ----------
    traj <- reactiveVal(NULL)
    traj_load_attempted <- reactiveVal(FALSE)
    traj_coord_lookup <- reactiveVal(NULL)
    traj_outliers <- reactiveVal(NULL)
    traj_outlier_summary <- reactiveVal(NULL)

    # ---------- Lazy H5AD loader ----------
    observe({
      if (traj_load_attempted()) return()

      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) return()

      traj_load_attempted(TRUE)

      if (!requireNamespace("anndata", quietly = TRUE)) {
        traj(list(error = "R package 'anndata' not installed.", error_detail = NULL))
        return()
      }

      ad <- NULL
      load_err <- NULL
      norm_path <- tryCatch(normalizePath(cur, mustWork = TRUE), error = function(e) cur)

      ad <- tryCatch(
        anndata::read_h5ad(norm_path),
        error = function(e) {
          load_err <<- conditionMessage(e)
          NULL
        }
      )

      if (is.null(ad)) {
        traj(list(error = sprintf("Failed to read H5AD: %s", norm_path), error_detail = load_err))
        return()
      }

      # Extract obs, obsm, uns
      obs  <- tryCatch(ad$obs,  error = function(e) NULL)
      obsm <- tryCatch(ad$obsm, error = function(e) NULL)
      uns  <- tryCatch(ad$uns,  error = function(e) NULL)
      var_names <- tryCatch(rownames(ad$var), error = function(e) NULL)

      if (is.null(obs)) {
        traj(list(error = "Could not extract obs data from AnnData", error_detail = NULL))
        return()
      }

      # Ensure obs is a data.frame
      if (!inherits(obs, "data.frame")) {
        obs <- tryCatch(as.data.frame(obs, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(obs)) {
          traj(list(error = "Obs coercion failed", error_detail = NULL))
          return()
        }
      }

      # Extract required columns. Both dpt_pseudotime (core) and
      # dpt_pseudotime_combined (core + peri) may be present; the user
      # toggles between them via input$pt_mode.
      cols <- intersect(
        c("combined_islet_id", "base_islet_id", "donor_status", "Case ID",
          "case_id", "imageid", "dpt_pseudotime", "dpt_pseudotime_combined",
          "age", "gender", "total_cells"),
        colnames(obs)
      )
      obs2 <- as.data.frame(obs[, cols, drop = FALSE])

      # Parse combined_islet_id for Case ID and islet_key
      if (!"Case ID" %in% names(obs2)) obs2$`Case ID` <- NA_character_

      if ("imageid" %in% names(obs2)) {
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, as.character(obs2$imageid))
      }

      if ("combined_islet_id" %in% names(obs2)) {
        cid <- stringr::str_extract(obs2$combined_islet_id, "^[0-9]{3,4}")
        islet_num <- stringr::str_extract(obs2$combined_islet_id, "(?<=_Islet_)[0-9]{1,3}")
        obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, cid)
        obs2$islet_key <- ifelse(!is.na(islet_num), paste0("Islet_", islet_num), NA_character_)
      }

      # Normalize donor_status
      if ("donor_status" %in% names(obs2)) {
        obs2$donor_status <- stringr::str_trim(obs2$donor_status)
        obs2$donor_status <- dplyr::case_when(
          tolower(obs2$donor_status) %in% c("nd", "non-diabetic", "nondiabetic") ~ "ND",
          tolower(obs2$donor_status) %in% c("aab+", "aab", "autoantibody") ~ "Aab+",
          tolower(obs2$donor_status) %in% c("t1d", "type1", "diabetic") ~ "T1D",
          TRUE ~ obs2$donor_status
        )
      }

      # Create generic pseudotime column (default to PAGA dpt_pseudotime)
      if ("dpt_pseudotime" %in% names(obs2)) {
        obs2$pseudotime <- obs2$dpt_pseudotime
      } else {
        obs2$pseudotime <- NA_real_
      }

      traj(list(obs = obs2, obsm = obsm, uns = uns, var_names = var_names, adata = ad))
    })

    # ---------- Pseudotime reactive ----------
    traj_with_pseudotime <- reactive({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error)) return(tr)

      obs <- tr$obs
      ad  <- tr$adata

      # Switch pseudotime mode based on user selection.
      # "core" uses dpt_pseudotime (core scVI features).
      # "combined" uses dpt_pseudotime_combined (core + peri-islet scVI features).
      # Falls back to whichever column is available if the requested one isn't.
      mode <- if (!is.null(input$pt_mode)) input$pt_mode else "core"
      if (mode == "combined" && "dpt_pseudotime_combined" %in% names(obs)) {
        obs$pseudotime <- obs$dpt_pseudotime_combined
      } else if ("dpt_pseudotime" %in% names(obs)) {
        obs$pseudotime <- obs$dpt_pseudotime
      } else if ("dpt_pseudotime_combined" %in% names(obs)) {
        obs$pseudotime <- obs$dpt_pseudotime_combined
      }

      list(obs = obs, obsm = tr$obsm, uns = tr$uns, var_names = tr$var_names, adata = ad)
    })

    # ---------- Feature selector UI ----------
    output$traj_feature_selector <- renderUI({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error)) {
        return(selectInput(ns("traj_feature"), "Select Feature:",
                           choices = c("No features available" = ""), selected = ""))
      }

      # Get available features from AnnData
      available_features <- NULL
      if (!is.null(tr$var_names)) {
        available_features <- tr$var_names
      } else {
        # Fallback - try to get from the loaded trajectory
        tryCatch({
          adata_path <- resolve_traj_path()
          if (!is.null(adata_path) && file.exists(adata_path)) {
            adata <- anndata::read_h5ad(adata_path)
            available_features <- rownames(adata$var)
          }
        }, error = function(e) NULL)
      }

      if (is.null(available_features) || length(available_features) == 0) {
        return(selectInput(ns("traj_feature"), "Select Feature:",
                           choices = c("No features available" = ""), selected = ""))
      }

      # Organize features into categories
      hormone_markers  <- intersect(c("INS", "GCG", "SST"), available_features)
      immune_markers   <- intersect(c("CD3e", "CD4", "CD8a", "CD68", "CD163", "CD20", "CD45", "HLADR"), available_features)
      vascular_markers <- intersect(c("CD31", "CD34", "SMA", "ColIV"), available_features)
      neural_markers   <- intersect(c("B3TUBB", "GAP43", "PGP9.5"), available_features)
      other_markers    <- setdiff(available_features, c(hormone_markers, immune_markers, vascular_markers, neural_markers))

      # Build choices list with categories.
      # Note: peri-islet immune/density columns are NOT exposed here as
      # Y-axis features. Users select pseudotime mode (Core / Core+Peri)
      # to incorporate the peri-zone and structural-density signals into
      # the X-axis (pseudotime) instead, then visualise the associated
      # markers (CD3e, CD8a, CD45, CD31, CD34, PDPN, etc.) on Y.
      choices <- list()
      if (length(hormone_markers)  > 0) choices[["Hormone Markers"]]  <- setNames(hormone_markers, hormone_markers)
      if (length(immune_markers)   > 0) choices[["Immune Markers"]]   <- setNames(immune_markers, immune_markers)
      if (length(vascular_markers) > 0) choices[["Vascular Markers"]] <- setNames(vascular_markers, vascular_markers)
      if (length(neural_markers)   > 0) choices[["Neural Markers"]]   <- setNames(neural_markers, neural_markers)
      if (length(other_markers)    > 0) choices[["Other Features"]]   <- setNames(other_markers, other_markers)

      # Default selection
      default_feature <- if ("INS" %in% available_features) "INS" else available_features[1]

      selectInput(ns("traj_feature"), "Select Feature:",
                  choices = choices,
                  selected = default_feature)
    })

    # ---------- Validate selected feature ----------
    traj_validate <- reactive({
      selected_feature <- input$traj_feature
      if (is.null(selected_feature) || !nzchar(selected_feature)) return(NULL)

      adata_path <- resolve_traj_path()
      if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)

      adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
      if (is.null(adata)) return(NULL)

      var_names <- rownames(adata$var)
      if (!selected_feature %in% var_names) return(NULL)

      list(feature = selected_feature, status = "ready")
    })

    # ---------- Main data pipeline ----------
    traj_data_clean <- reactive({
      ms <- traj_validate()
      if (is.null(ms)) return(NULL)

      selected_feature <- ms$feature
      if (is.null(selected_feature)) return(NULL)

      adata_path <- resolve_traj_path()
      if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)

      tr <- traj_with_pseudotime()
      if (is.null(tr) || !is.null(tr$error)) return(NULL)

      # Use cached adata object instead of re-reading
      adata <- tr$adata
      if (is.null(adata)) {
        adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
        if (is.null(adata)) return(NULL)
      }

      var_names <- rownames(adata$var)
      if (!selected_feature %in% var_names) return(NULL)

      obs_df <- tr$obs
      if (is.null(obs_df)) return(NULL)

      # Extract expression values from .X (markers only).
      expression_vals <- tryCatch({
        idx <- which(var_names == selected_feature)
        as.numeric(adata$X[, idx])
      }, error = function(e) {
        message("[traj_data_clean] Error extracting feature: ", e$message)
        rep(NA_real_, nrow(obs_df))
      })

      # Get donor_status with proper normalization
      donor_clean <- obs_df$donor_status
      if (is.null(donor_clean)) donor_clean <- rep(NA_character_, nrow(obs_df))

      # Extract UMAP coordinates
      umap_coords <- NULL
      if (!is.null(tr$obsm) && "X_umap" %in% names(tr$obsm)) {
        umap_coords <- tr$obsm$X_umap
      }

      # Get donor ID for coloring
      donor_ids <- as.character(obs_df$imageid)

      # Build result data.frame
      islet_keys <- as.character(obs_df$base_islet_id)
      islet_keys <- gsub("^Islet_Islet_", "Islet_", islet_keys)

      result_df <- data.frame(
        case_id = as.character(obs_df$imageid),
        islet_key = islet_keys,
        combined_islet_id = as.character(obs_df$combined_islet_id),
        pt = as.numeric(obs_df$pseudotime),
        donor_status = donor_clean,
        donor_id = donor_ids,
        value = as.numeric(expression_vals),
        feature_name = selected_feature,
        total_cells = as.integer(obs_df$total_cells),
        stringsAsFactors = FALSE
      )

      # Guard: fill missing total_cells with 1 (safe default weight)
      if (any(is.na(result_df$total_cells))) {
        result_df$total_cells[is.na(result_df$total_cells)] <- 1L
      }

      # Add UMAP coordinates
      if (!is.null(umap_coords) && nrow(umap_coords) == nrow(result_df)) {
        result_df$umap_1 <- as.numeric(umap_coords[, 1])
        result_df$umap_2 <- as.numeric(umap_coords[, 2])
      } else {
        result_df$umap_1 <- NA_real_
        result_df$umap_2 <- NA_real_
      }

      # Merge diameter and Leiden cluster cols from prepared() comp table.
      prep_data <- prepared()
      if (!is.null(prep_data) && !is.null(prep_data$comp)) {
        leiden_cols <- intersect(c("leiden_0.3", "leiden_0.5", "leiden_0.8", "leiden_1.0"),
                                 colnames(prep_data$comp))
        lookup_cols <- c("Case ID", "Donor Status", "islet_key", "islet_diam_um",
                         leiden_cols)
        size_lookup <- prep_data$comp[intersect(lookup_cols, colnames(prep_data$comp))]
        # comp$`Case ID` is integer (data_loading.R:346); match the type so merge succeeds.
        result_df$`Case ID` <- suppressWarnings(as.integer(result_df$case_id))
        result_df$`Donor Status` <- result_df$donor_status

        # Base R merge (faster for this use case)
        merged <- merge(result_df, size_lookup,
                        by = c("Case ID", "Donor Status", "islet_key"),
                        all.x = TRUE, sort = FALSE)

        match_idx <- match(
          paste(result_df$`Case ID`, result_df$`Donor Status`, result_df$islet_key),
          paste(merged$`Case ID`, merged$`Donor Status`, merged$islet_key)
        )

        if ("islet_diam_um" %in% colnames(merged)) {
          result_df$islet_diam_um <- merged$islet_diam_um[match_idx]
        }
        for (lc in leiden_cols) {
          if (lc %in% colnames(merged)) {
            result_df[[lc]] <- as.character(merged[[lc]][match_idx])
          }
        }
      }

      if (!"islet_diam_um" %in% colnames(result_df)) {
        result_df$islet_diam_um <- NA_real_
      }

      # Add spatial coordinates from segmentation data
      if (!is.null(islet_spatial_lookup)) {
        result_df$centroid_x_um <- NA_real_
        result_df$centroid_y_um <- NA_real_
        result_df$has_spatial <- FALSE

        for (i in seq_len(nrow(result_df))) {
          annot <- get_islet_annotations(result_df$case_id[i], result_df$islet_key[i])
          if (!is.null(annot)) {
            result_df$centroid_x_um[i] <- annot$centroid_x
            result_df$centroid_y_um[i] <- annot$centroid_y
            result_df$has_spatial[i] <- TRUE
          }
        }

        cat("[traj_data_clean] Added spatial coordinates for",
            sum(result_df$has_spatial), "islets\n")
      }

      # Filter valid rows
      valid_idx <- which(
        is.finite(result_df$pt) &
        !is.na(result_df$value) &
        !is.na(result_df$case_id) &
        !is.na(result_df$islet_key)
      )

      if (length(valid_idx) == 0) {
        message("[traj_data_clean] No valid rows after filtering")
        return(NULL)
      }

      result_df <- result_df[valid_idx, , drop = FALSE]

      # Apply min cells/islet filter
      min_cells <- input$min_cells %||% 1
      if (min_cells > 1 && "total_cells" %in% colnames(result_df)) {
        tc <- suppressWarnings(as.integer(result_df$total_cells))
        result_df <- result_df[is.finite(tc) & tc >= min_cells, , drop = FALSE]
        if (nrow(result_df) == 0) return(NULL)
      }

      # Rank-transform pseudotime for visualization (eliminates DPT compression artifacts)
      # Raw DPT compresses 95% of data between 0.25-0.70 due to diffusion manifold geometry.
      # Rank percentile preserves ordering while spreading points uniformly.
      result_df$pt_raw <- result_df$pt
      result_df$pt <- rank(result_df$pt, ties.method = "average") / nrow(result_df)

      message("[traj_data_clean] Successfully processed ", nrow(result_df),
              " observations for ", selected_feature)
      return(result_df)
    })

    # ---------- Status bar ----------
    output$traj_status <- renderUI({
      cur <- traj_path
      if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
      if (is.null(cur) || !file.exists(cur)) {
        return(tags$div(style = "color:#b00;",
          "AnnData not found. Place 'adata_ins_root.h5ad' under data/ or scripts/, or set ADATA_PATH."))
      }
      if (!requireNamespace("anndata", quietly = TRUE)) {
        return(tags$div(style = "color:#b00;",
          "R package 'anndata' not installed. Install with BiocManager::install('anndata')."))
      }
      tr <- traj_with_pseudotime()
      if (is.null(tr)) {
        return(tags$div(style = "color:#b00;",
          "AnnData object not available or H5AD failed to load."))
      }
      if (!is.null(tr$error)) {
        det <- if (!is.null(tr$error_detail) && nzchar(tr$error_detail)) paste0(" Details: ", tr$error_detail) else ""
        install_hint <- if (grepl("anndata", tr$error, ignore.case = TRUE)) " Install with: BiocManager::install('anndata')" else ""
        return(tags$div(style = "color:#b00;",
          sprintf("Trajectory load error: %s%s%s", tr$error, det, install_hint)))
      }
      ms <- traj_validate()
      jn <- traj_data_clean()
      nm <- if (is.null(ms)) 0 else nrow(ms)
      nj <- if (is.null(jn)) 0 else nrow(jn)
      ndist <- if (is.null(jn)) 0 else nrow(dplyr::distinct(jn, case_id, islet_key))
      src <- if (!is.null(traj_path) && file.exists(traj_path)) traj_path else cur

      tags$div(style = "color:#555;",
        sprintf("Loaded AnnData (%d obs). Joined %d rows (%d unique islets) of %d islet metrics. Source: %s",
                nrow(tr$obs), nj, ndist, nm, basename(src))
      )
    })

    # ---------- Trajectory scatter plot ----------
    output$traj_scatter <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) {
        return(plotly_empty() %>% layout(title = "No trajectory data available"))
      }

      # Flag outliers via the shared helper. Excluded points are drawn as
      # ghost × markers (not deleted) so users can see WHERE the exclusions
      # land along pseudotime; trend lines are drawn only from kept rows.
      df <- flag_outliers(df, value_col = "value", group_col = "donor_status",
                          method = "zscore", threshold = outlier_threshold(),
                          enabled = isTRUE(remove_outliers()))
      outlier_summary <- summarize_outlier_filter(df)
      outliers <- df[df$is_outlier, , drop = FALSE]

      if (nrow(outliers) > 0) {
        traj_outliers(data.frame(
          Feature      = outliers$feature_name,
          Case_ID      = outliers$case_id,
          Islet        = outliers$islet_key,
          Donor_Status = outliers$donor_status,
          Pseudotime   = round(outliers$pt, 3),
          Value        = round(outliers$value, 3),
          Z_Score      = round(outliers$outlier_score, 2),
          stringsAsFactors = FALSE
        ))
      } else {
        traj_outliers(NULL)
      }

      # Split: kept rows feed the scatter, color scales, and trend lines.
      # Outlier rows are rendered as a final ghost layer.
      df_keep    <- df[!df$is_outlier, , drop = FALSE]
      df_ghost   <- df[df$is_outlier, , drop = FALSE]

      cat(sprintf("[TRAJECTORY] %d outliers flagged (%s); %d kept, %d drawn as ghost ×\n",
                  outlier_summary$n_excluded, outlier_summary$label,
                  nrow(df_keep), nrow(df_ghost)))
      traj_outlier_summary(outlier_summary)

      df <- df_keep

      # Determine color mapping and aesthetics
      color_by   <- input$traj_color_by %||% "donor_status"
      size_by    <- input$traj_point_size %||% "uniform"
      trend_type <- input$traj_show_trend %||% "by_donor"
      alpha_val  <- 1 - (input$traj_alpha %||% 0.4)
      point_size <- input$traj_point_size_slider %||% 3.0

      # Always apply jitter — use the kept-rows range as the jitter scale
      # so ghost points land in the same visual frame.
      set.seed(42)
      jitter_amount_x <- diff(range(df$pt, na.rm = TRUE)) * 0.012
      jitter_amount_y <- diff(range(df$value, na.rm = TRUE)) * 0.02
      df$pt    <- pmax(0, pmin(1, df$pt + runif(nrow(df), -jitter_amount_x, jitter_amount_x)))
      df$value <- df$value + runif(nrow(df), -jitter_amount_y, jitter_amount_y)
      if (nrow(df_ghost) > 0) {
        df_ghost$pt    <- pmax(0, pmin(1, df_ghost$pt + runif(nrow(df_ghost), -jitter_amount_x, jitter_amount_x)))
        df_ghost$value <- df_ghost$value + runif(nrow(df_ghost), -jitter_amount_y, jitter_amount_y)
      }

      # Log-transform cell counts for LOESS weights (raw counts too skewed: median=9, max=1902)
      df$loess_weight <- log1p(df$total_cells)

      # Create base aesthetic mapping with cell count for hover
      df$hover_cells <- paste0("Cells: ", df$total_cells)
      df$click_key <- paste(df$case_id, df$islet_key, sep = "|")
      aes_mapping <- aes(x = pt, y = value, text = hover_cells, key = click_key)

      # Add color aesthetic
      if (color_by == "donor_status") {
        df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
        aes_mapping$colour <- as.name("donor_status")
        color_title <- "Donor Status"
      } else if (color_by == "donor_id") {
        df$donor_id <- factor(df$donor_id)
        aes_mapping$colour <- as.name("donor_id")
        color_title <- "Donor ID"
      } else if (startsWith(color_by, "leiden_") && color_by %in% colnames(df)) {
        # Sort factor levels numerically so legend reads 0,1,2,... instead of 0,1,10,11
        lvl_vals <- unique(stats::na.omit(df[[color_by]]))
        num_vals <- suppressWarnings(as.numeric(lvl_vals))
        sorted_lvls <- if (!any(is.na(num_vals))) lvl_vals[order(num_vals)] else sort(lvl_vals)
        df[[color_by]] <- factor(df[[color_by]], levels = sorted_lvls)
        aes_mapping$colour <- as.name(color_by)
        color_title <- paste("Leiden", sub("^leiden_", "res ", color_by))
      } else if (startsWith(color_by, "leiden_")) {
        # Column not present in df (e.g., Excel fallback) — fall back to donor_status
        df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
        aes_mapping$colour <- as.name("donor_status")
        color_title <- "Donor Status"
        color_by <- "donor_status"
      }

      # Add size aesthetic if not uniform
      if (size_by != "uniform") {
        if (size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
          if (!all(is.na(df$islet_diam_um))) {
            aes_mapping$size <- as.name("islet_diam_um")
          }
        } else if (size_by == "total_cells" && "total_cells" %in% colnames(df)) {
          df$size_cells <- sqrt(df$total_cells)
          aes_mapping$size <- as.name("size_cells")
        }
      }

      # Create plot
      g <- ggplot(df, aes_mapping)
      if (size_by == "uniform") {
        g <- g + geom_point(alpha = alpha_val, size = point_size, stroke = 0)
      } else {
        g <- g + geom_point(alpha = alpha_val, stroke = 0)
      }
      # Ghost layer: excluded outliers as faint grey × markers. Drawn before
      # trend lines so the LOESS curves visibly track the kept distribution.
      if (nrow(df_ghost) > 0) {
        df_ghost$hover_cells <- paste0("Excluded outlier · Cells: ", df_ghost$total_cells,
                                       sprintf(" · z = %.2f", df_ghost$outlier_score))
        df_ghost$click_key <- paste(df_ghost$case_id, df_ghost$islet_key, sep = "|")
        g <- g + geom_point(
          data = df_ghost,
          aes(x = pt, y = value, text = hover_cells, key = click_key),
          inherit.aes = FALSE,
          shape = 4,
          size = point_size * 1.1,
          stroke = 0.9,
          colour = "grey45",
          alpha = 0.55
        )
      }
      g <- g +
        theme_islet() +
        theme(
          legend.position = "none",
          plot.margin = unit(c(5, 5, 5, 5), "pt"),
          plot.title = element_text(margin = margin(b = 15))
        )

      # Apply color scales
      if (color_by == "donor_status") {
        g <- g + scale_color_manual(
          values = donor_colors_reactive(),
          name = color_title, drop = FALSE)
      } else if (color_by == "donor_id") {
        df$donor_id <- factor(df$donor_id)
        n_donors <- length(levels(df$donor_id))
        if (n_donors <= 12) {
          colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
          names(colors) <- levels(df$donor_id)[1:length(colors)]
        } else {
          colors <- rainbow(n_donors, s = 0.8, v = 0.8)
          names(colors) <- levels(df$donor_id)
        }
        g <- g + scale_color_manual(values = colors, name = color_title, drop = FALSE)
      } else if (startsWith(color_by, "leiden_") && color_by %in% colnames(df)) {
        clust_levels <- levels(df[[color_by]])
        n_clust <- length(clust_levels)
        clust_colors <- if (n_clust <= 12) {
          RColorBrewer::brewer.pal(max(3, min(12, n_clust)), "Set3")[seq_len(n_clust)]
        } else {
          scales::hue_pal()(n_clust)
        }
        names(clust_colors) <- clust_levels
        g <- g + scale_color_manual(values = clust_colors, name = color_title, drop = FALSE)
      }

      # Size scale
      if (size_by != "uniform" && size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
        if (!all(is.na(df$islet_diam_um))) {
          size_range <- c(0.5 * point_size, 2.5 * point_size)
          g <- g + scale_size_continuous(name = "Islet Diameter (um)", range = size_range, guide = "none")
        }
      } else if (size_by == "total_cells" && "size_cells" %in% colnames(df)) {
        size_range <- c(0.3 * point_size, 3.0 * point_size)
        g <- g + scale_size_continuous(name = "Cell Count", range = size_range, guide = "none")
      }

      # Trend lines
      pt_range <- range(df$pt, na.rm = TRUE)

      if (trend_type == "overall") {
        g <- g + geom_smooth(method = "loess", se = FALSE, color = "black",
                             size = 1, fullrange = FALSE, span = 0.75,
                             inherit.aes = FALSE,
                             aes(x = pt, y = value, weight = loess_weight))
      } else if (trend_type == "by_donor") {
        if ("donor_status" %in% names(df) && !all(is.na(df$donor_status))) {
          df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
          trend_colors <- donor_colors_reactive()

          if (color_by == "donor_id" || startsWith(color_by, "leiden_")) {
            # Use explicit color per-status loop — adding a second scale_color_manual
            # after a Leiden/donor_id cluster scale would override it and blank out the points.
            for (status in c("ND", "Aab+", "T1D")) {
              df_subset <- df[df$donor_status == status & !is.na(df$donor_status), ]
              if (nrow(df_subset) > 0) {
                g <- g + geom_smooth(data = df_subset,
                                     aes(x = pt, y = value, weight = loess_weight),
                                     inherit.aes = FALSE,
                                     method = "loess", se = FALSE,
                                     size = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
                                     color = trend_colors[status])
              }
            }
          } else {
            g <- g + geom_smooth(
              method = "loess", se = FALSE,
              size = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
              inherit.aes = FALSE,
              mapping = aes(x = pt, y = value, group = donor_status, color = donor_status, weight = loess_weight)) +
              scale_color_manual(values = trend_colors,
                                name = if (color_by == "donor_status") color_title else "Trend Lines",
                                drop = FALSE,
                                guide = if (color_by == "donor_status") "legend" else "none")
          }
        }
      }

      # Labels
      selected_feature <- input$traj_feature %||% "Selected feature"
      metric_label <- "— normalized expression (mean per islet, 0–1)"
      x_label <- "Pseudotime (rank-normalized)"

      g <- g + labs(
        x = x_label,
        y = paste(selected_feature, metric_label),
        title = NULL
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = range(df$value, na.rm = TRUE) * c(1, 1) + c(-0.05, 0.05) * diff(range(df$value, na.rm = TRUE)))

      # Convert to plotly with custom data for reliable click handling
      p <- ggplotly(g, tooltip = c("x", "y", "colour", "text"), source = ns("traj_scatter"))
      # Strip "colour: black" leak (from overall LOESS line's color="black"
      # parameter) and collapse duplicate donor_status entries in hover text.
      p <- strip_plotly_color_leak(p)

      # Store coordinate lookup as fallback for click handling
      if ("combined_islet_id" %in% colnames(df) && "case_id" %in% colnames(df)) {
        traj_coord_lookup(data.frame(
          pt                = df$pt,
          value             = df$value,
          case_id           = df$case_id,
          islet_key         = df$islet_key,
          combined_islet_id = df$combined_islet_id,
          donor_status      = df$donor_status,
          stringsAsFactors  = FALSE
        ))
      }

      # Layout and register click events
      fonts <- plotly_islet_fonts()
      p %>%
        layout(
          showlegend = FALSE,
          title = NULL,
          font = fonts$global,
          margin = list(l = 80, r = 20, t = 30, b = 70, pad = 5),
          xaxis = list(fixedrange = FALSE, automargin = TRUE, range = c(0, 1),
                       title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(fixedrange = FALSE, automargin = TRUE,
                       title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(font = fonts$legend_text,
                        title = list(font = fonts$legend_title))
        ) %>%
        plotly_islet_config("trajectory_scatter") %>%
        event_register("plotly_click")
    })

    # ---------- Click handler for trajectory scatter ----------
    observeEvent(event_data("plotly_click", source = ns("traj_scatter")), {
      click_data <- event_data("plotly_click", source = ns("traj_scatter"))
      if (is.null(click_data)) return()

      # Check if segmentation viewer is available
      if (!SF_AVAILABLE) {
        showNotification("Segmentation viewer requires the 'sf' package. Install with: install.packages('sf')",
                         type = "warning", duration = 5)
        return()
      }

      if (is.null(islet_spatial_lookup)) {
        showNotification("Islet spatial lookup data not available", type = "warning", duration = 5)
        return()
      }

      cat("[TRAJ CLICK] Received click event\n")

      # --- Priority 1: Parse key (format "case_id|islet_key") ---
      case_id   <- NULL
      islet_key <- NULL
      click_key <- click_data$key
      if (!is.null(click_key) && length(click_key) > 0 && !is.na(click_key[1]) && nzchar(click_key[1])) {
        key_parts <- strsplit(as.character(click_key[1]), "\\|")[[1]]
        if (length(key_parts) >= 2) {
          case_id   <- key_parts[1]
          islet_key <- key_parts[2]
          cat("[TRAJ CLICK] From key:", case_id, islet_key, "\n")
        }
      }

      # --- Priority 2: Fallback to customdata ---
      if (is.null(case_id) || is.null(islet_key)) {
        custom <- click_data$customdata
        if (!is.null(custom) && length(custom) >= 2) {
          case_id   <- custom[[1]]
          islet_key <- custom[[2]]
          cat("[TRAJ CLICK] From customdata:", case_id, islet_key, "\n")
        }
      }

      # --- Priority 3: Fallback to coordinate lookup ---
      if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
        lookup <- traj_coord_lookup()
        if (!is.null(lookup) && nrow(lookup) > 0) {
          tolerance <- 0.01
          matches <- which(
            abs(lookup$pt - click_data$x) < tolerance &
            abs(lookup$value - click_data$y) < tolerance
          )
          if (length(matches) > 0) {
            idx <- matches[1]
            case_id   <- lookup$case_id[idx]
            islet_key <- lookup$islet_key[idx]
            cat("[TRAJ CLICK] From coord lookup:", case_id, islet_key, "\n")
          }
        }
      }

      if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
        showNotification("Could not identify clicked islet", type = "warning", duration = 3)
        return()
      }

      # Look up centroid coordinates from spatial lookup
      spatial_match <- islet_spatial_lookup[
        islet_spatial_lookup$case_id == case_id & islet_spatial_lookup$islet_key == islet_key,
      ]

      # Try zero-padded variant (e.g. "112" -> "0112")
      if (nrow(spatial_match) == 0) {
        case_id_num <- suppressWarnings(as.integer(case_id))
        if (!is.na(case_id_num)) {
          case_id_padded <- sprintf("%04d", case_id_num)
          spatial_match <- islet_spatial_lookup[
            islet_spatial_lookup$case_id == case_id_padded & islet_spatial_lookup$islet_key == islet_key,
          ]
          if (nrow(spatial_match) > 0) case_id <- case_id_padded
        }
      }

      if (nrow(spatial_match) == 0) {
        showNotification(paste("Spatial data not found for", islet_key, "in case", case_id),
                         type = "warning", duration = 3)
        return()
      }

      centroid_x <- spatial_match$centroid_x_um[1]
      centroid_y <- spatial_match$centroid_y_um[1]

      cat("[TRAJ CLICK] Centroid: x=", centroid_x, ", y=", centroid_y, "\n")

      # Look up donor info from prepared data
      pd <- prepared()
      donor_status_val <- NULL
      donor_age <- NULL
      donor_gender <- NULL
      if (!is.null(pd$comp) && "Case ID" %in% colnames(pd$comp)) {
        donor_row <- pd$comp[pd$comp$`Case ID` == as.integer(case_id), , drop = FALSE]
        if (nrow(donor_row) > 0) {
          if ("Donor Status" %in% colnames(donor_row)) donor_status_val <- as.character(donor_row$`Donor Status`[1])
          if ("age" %in% colnames(donor_row)) donor_age <- donor_row$age[1]
          if ("gender" %in% colnames(donor_row)) donor_gender <- as.character(donor_row$gender[1])
        }
      }

      # Store selected islet info - triggers embedded viewer to update
      selected_islet(list(
        case_id      = case_id,
        islet_key    = islet_key,
        centroid_x   = centroid_x,
        centroid_y   = centroid_y,
        donor_status = donor_status_val,
        donor_age    = donor_age,
        donor_gender = donor_gender
      ))
    })

    # ---------- Segmentation viewer panel ----------
    output$segmentation_viewer_panel <- renderUI({
      # Only render when Trajectory tab is active to avoid duplicate non-namespaced output IDs
      if (active_tab() != "Trajectory") return(NULL)
      info <- selected_islet()
      if (is.null(info)) return(NULL)

      # Check if single-cell drill-down data exists for this islet
      has_cells <- drilldown_available() &&
        !is.null(load_islet_cells(info$case_id, info$islet_key))

      # Marker choices for color-by dropdown
      marker_choices <- c("phenotype" = "phenotype")
      if (has_cells) {
        cells <- load_islet_cells(info$case_id, info$islet_key)
        marker_cols <- setdiff(colnames(cells),
                               c("X_centroid", "Y_centroid", "phenotype", "cell_region",
                                 "Cell.Area", "Cell Area", "Nucleus.Area", "Nucleus Area"))
        if (length(marker_cols) > 0) {
          marker_choices <- c(marker_choices, setNames(marker_cols, marker_cols))
        }
      }

      # Build title: "Islet_N (case_id, status, age, gender)"
      title_detail <- info$case_id
      if (!is.null(info$donor_status)) title_detail <- paste0(title_detail, ", ", info$donor_status)
      if (!is.null(info$donor_age) && !is.na(info$donor_age)) title_detail <- paste0(title_detail, ", age ", round(info$donor_age))
      if (!is.null(info$donor_gender) && !is.na(info$donor_gender) && nzchar(info$donor_gender)) title_detail <- paste0(title_detail, ", ", info$donor_gender)

      div(class = "card", style = "padding: 20px; margin-bottom: 20px; border: 2px solid #0066CC;",
        # Header: title + close button
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px;",
          h4(paste0(info$islet_key, " (", title_detail, ")"),
             style = "margin: 0; color: #0066CC; font-size: 20px;"),
          actionButton(ns("clear_segmentation"), "Close", class = "btn btn-sm btn-outline-secondary")
        ),
        # Main content: controls sidebar (left) + plot (right)
        fluidRow(
          # Left sidebar: all controls
          column(2,
            div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
              if (has_cells) {
                tagList(
                  h5("Layers", style = "margin-top: 0; font-size: 15px;"),
                  checkboxInput("drilldown_show_cells", "Single Cells", value = TRUE),
                  checkboxInput("drilldown_show_peri_boundary", "Peri Boundary", value = TRUE),
                  checkboxInput("drilldown_show_structures", "Structures", value = TRUE),
                  hr(style = "margin: 8px 0;"),
                  conditionalPanel(
                    condition = "input.drilldown_show_cells",
                    selectInput("drilldown_color_by", "Color by",
                                choices = marker_choices, selected = "phenotype"),
                    selectInput("drilldown_palette", "Phenotype Palette",
                                choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                                selected = "High Contrast"),
                    checkboxInput("drilldown_show_peri", "Show peri-islet cells", value = TRUE)
                  )
                )
              } else {
                tagList(
                  h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-size: 15px;"),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #0066CC; margin-right: 8px;"),
                    span("Islet", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #00CCCC; margin-right: 8px;"),
                    span("+20\u00b5m", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #CC00CC; margin-right: 8px;"),
                    span("Nerve", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #CC0000; margin-right: 8px;"),
                    span("Capillary", style = "font-size: 14px;")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 24px; height: 3px; background-color: #00AA00; margin-right: 8px;"),
                    span("Lymphatic", style = "font-size: 14px;")
                  )
                )
              },
              hr(style = "margin: 8px 0;"),
              span(paste("Centroid:", round(info$centroid_x, 0), "/", round(info$centroid_y, 0), "\u00b5m"),
                   style = "color: #666; font-size: 13px;")
            )
          ),
          # Center: plot
          column(6,
            plotOutput("islet_segmentation_view", height = "550px"),
            div(style = "text-align: right; margin-top: 6px;",
              downloadButton("dl_islet_segmentation_view", "Download PNG",
                             class = "btn-sm btn-outline-primary"))
          ),
          # Right: composition sidebar
          column(4,
            if (has_cells) {
              conditionalPanel(
                condition = "input.drilldown_show_cells",
                div(style = "padding: 12px; background-color: #f8f9fa; border-radius: 5px;",
                  h5("Cell Composition", style = "margin-top: 0; font-size: 16px;"),
                  plotOutput("islet_drilldown_summary", height = "300px"),
                  div(style = "text-align: right; margin-top: 6px;",
                    downloadButton("dl_islet_drilldown_summary", "Download PNG",
                                   class = "btn-sm btn-outline-primary")),
                  tableOutput("islet_drilldown_table"),
                  uiOutput("islet_drilldown_rules")
                )
              )
            }
          )
        )
      )
    })

    # ---------- Clear segmentation viewer ----------
    observeEvent(input$clear_segmentation, {
      selected_islet(NULL)
    })

    # Note: islet_segmentation_view renderPlot is at root level in app.R
    # (shared between Plot modal and Trajectory embedded panel)

    # ---------- Outlier badge (above scatter) ----------
    output$traj_outlier_badge <- renderUI({
      summary <- traj_outlier_summary()
      if (is.null(summary)) return(NULL)
      suffix <- outlier_title_suffix(summary)
      if (!nzchar(suffix)) {
        if (!isTRUE(summary$enabled)) {
          return(tags$div(
            style = "background: #fdfbe7; border: 1px solid #f0d97a; border-radius: 4px; padding: 4px 10px; margin-bottom: 6px; font-size: 12px; color: #6e6300;",
            tags$span(title = GLOBAL_OUTLIER_TOOLTIP,
              "\u24d8 Outlier filter OFF \u2014 all points contribute to trend lines.")
          ))
        }
        return(NULL)
      }
      tags$div(
        style = "background: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; padding: 4px 10px; margin-bottom: 6px; font-size: 12px; color: #5a4400;",
        tags$span(title = GLOBAL_OUTLIER_TOOLTIP,
          sprintf("\u26a0 %d outlier%s excluded from trend lines (%s) \u2014 drawn as grey \u00d7 markers.",
                  summary$n_excluded,
                  if (summary$n_excluded == 1) "" else "s",
                  summary$label))
      )
    })

    # ---------- Outlier table (below scatter) ----------
    output$traj_outlier_info <- renderUI({
      show_table <- isTRUE(input$show_outlier_table)
      if (!show_table) return(NULL)
      outlier_data <- traj_outliers()
      summary <- traj_outlier_summary()
      filter_on <- isTRUE(summary$enabled)

      if (is.null(outlier_data) || nrow(outlier_data) == 0) {
        return(tagList(
          div(style = "background-color: #e8f4fd; border: 1px solid #b8daff; border-radius: 5px; padding: 8px 12px; margin-top: 10px; color: #004085; font-size: 12px;",
            "No outliers detected for the current selection.")
        ))
      }

      thr <- suppressWarnings(as.numeric(outlier_threshold()))
      if (!is.finite(thr)) thr <- 3
      header <- if (filter_on) {
        sprintf("\u26a0\ufe0f %d Outlier%s Excluded (|z| > %g within donor status)",
                nrow(outlier_data), ifelse(nrow(outlier_data) > 1, "s", ""), thr)
      } else {
        sprintf("\u24d8 %d Outlier%s Flagged (filter OFF \u2014 still in trend)",
                nrow(outlier_data), ifelse(nrow(outlier_data) > 1, "s", ""))
      }
      body_text <- if (filter_on) {
        sprintf("The following data points were excluded from trend lines and rendered as faint \u00d7 markers because their robust |z| (median/MAD modified z-score) exceeds %g \u2014 flagged on both the high and low side of the donor-status distribution:", thr)
      } else {
        "These points would be excluded if the global outlier filter were turned on; currently they are kept and contribute to trend lines:"
      }

      tagList(
        div(style = "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; padding: 10px; margin-top: 10px; color: #000000;",
          h5(style = "color: #000000; margin-top: 0;", header),
          p(style = "color: #000000; font-size: 12px; margin-bottom: 10px;", body_text),
          div(style = "max-height: 200px; overflow-y: auto;",
            renderTable({
              outlier_data
            }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "xs", width = "100%")
          )
        )
      )
    })

    # ---------- Legend for trajectory scatter ----------
    output$traj_legend <- renderUI({
      color_by <- input$traj_color_by %||% "donor_status"

      if (color_by == "donor_status") {
        dc <- donor_colors_reactive()
        tagList(
          div(style = "display: flex; align-items: center; margin-bottom: 8px;",
            div(style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 3px; margin-right: 8px;", dc[["ND"]])),
            span("ND", style = "font-size: 15px;")
          ),
          div(style = "display: flex; align-items: center; margin-bottom: 8px;",
            div(style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 3px; margin-right: 8px;", dc[["Aab+"]])),
            span("Aab+", style = "font-size: 15px;")
          ),
          div(style = "display: flex; align-items: center;",
            div(style = sprintf("width: 20px; height: 20px; background-color: %s; border-radius: 3px; margin-right: 8px;", dc[["T1D"]])),
            span("T1D", style = "font-size: 15px;")
          )
        )
      } else if (color_by == "donor_id") {
        df <- traj_data_clean()
        if (!is.null(df) && nrow(df) > 0) {
          donor_ids <- sort(unique(df$donor_id))
          n_donors <- length(donor_ids)

          # Generate colors matching the plot
          if (n_donors <= 12) {
            colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
          } else {
            colors <- rainbow(n_donors, s = 0.8, v = 0.8)
          }

          # Split into two columns
          mid_point <- ceiling(n_donors / 2)
          col1_ids <- donor_ids[1:mid_point]
          col2_ids <- if (n_donors > mid_point) donor_ids[(mid_point + 1):n_donors] else NULL

          col1_items <- lapply(seq_along(col1_ids), function(i) {
            div(style = "display: flex; align-items: center; margin-bottom: 6px;",
              div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[i])),
              span(col1_ids[i], style = "font-size: 14px;")
            )
          })

          col2_items <- if (!is.null(col2_ids)) {
            lapply(seq_along(col2_ids), function(i) {
              idx <- mid_point + i
              div(style = "display: flex; align-items: center; margin-bottom: 6px;",
                div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[idx])),
                span(col2_ids[i], style = "font-size: 14px;")
              )
            })
          } else {
            NULL
          }

          tagList(
            div(style = "display: flex; gap: 10px;",
              div(style = "flex: 1;", col1_items),
              if (!is.null(col2_items)) div(style = "flex: 1;", col2_items)
            )
          )
        } else {
          p("No data available", style = "font-size: 12px; color: #999;")
        }
      } else if (startsWith(color_by, "leiden_")) {
        df <- traj_data_clean()
        if (is.null(df) || !(color_by %in% colnames(df))) {
          return(p("Leiden clusters unavailable (rebuild islet_explorer.h5ad).",
                   style = "font-size: 12px; color: #999;"))
        }
        lvl_vals <- unique(stats::na.omit(df[[color_by]]))
        num_vals <- suppressWarnings(as.numeric(lvl_vals))
        clust_levels <- if (!any(is.na(num_vals))) lvl_vals[order(num_vals)] else sort(lvl_vals)
        n_clust <- length(clust_levels)
        clust_colors <- if (n_clust <= 12) {
          RColorBrewer::brewer.pal(max(3, min(12, n_clust)), "Set3")[seq_len(n_clust)]
        } else {
          scales::hue_pal()(n_clust)
        }
        mid_point <- ceiling(n_clust / 2)
        col1_idx <- seq_len(mid_point)
        col2_idx <- if (n_clust > mid_point) (mid_point + 1):n_clust else integer(0)
        make_row <- function(i) {
          div(style = "display: flex; align-items: center; margin-bottom: 6px;",
            div(style = sprintf(
              "width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;",
              clust_colors[i])),
            span(paste("Cluster", clust_levels[i]), style = "font-size: 14px;")
          )
        }
        tagList(
          div(style = "display: flex; gap: 10px;",
            div(style = "flex: 1;", lapply(col1_idx, make_row)),
            if (length(col2_idx) > 0) div(style = "flex: 1;", lapply(col2_idx, make_row))
          )
        )
      }
    })

    # ---------- UMAP: Donor Status ----------
    traj_umap_donor_plot <- reactive({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
          theme_void())
      }

      df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))

      ggplot(df, aes(x = umap_1, y = umap_2, color = donor_status)) +
        geom_point(alpha = 0.6, size = 1.0) +
        scale_color_manual(values = donor_colors_reactive()) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP: Donor Status") +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_islet() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank())
    })
    output$traj_umap_donor <- renderPlot({ traj_umap_donor_plot() })
    output$dl_umap_donor <- downloadHandler(
      filename = function() islet_png_filename("trajectory_umap_donor"),
      content = function(file) {
        p <- traj_umap_donor_plot(); if (is.null(p)) return()
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ---------- UMAP: Selected Feature ----------
    traj_umap_feature_plot <- reactive({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
          theme_void())
      }

      selected_feature <- input$traj_feature %||% "Selected feature"

      g <- ggplot(df, aes(x = umap_1, y = umap_2, color = value)) +
        geom_point(alpha = 0.6, size = 1.0) +
        scale_x_continuous(expand = expansion(mult = 0.02)) +
        scale_y_continuous(expand = expansion(mult = 0.02)) +
        coord_fixed() +
        labs(x = "UMAP 1", y = "UMAP 2",
             title = paste("UMAP:", selected_feature),
             color = "Expression") +
        theme_islet() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank(), axis.line = element_blank())

      # Continuous colormap scaled to data min/max
      val_range <- range(df$value, na.rm = TRUE)
      g <- g + scale_color_viridis_c(option = "inferno", limits = val_range,
                                      na.value = "#bbbbbb")

      g
    })
    output$traj_umap_feature <- renderPlot({ traj_umap_feature_plot() })
    output$dl_umap_feature <- downloadHandler(
      filename = function() islet_png_filename("trajectory_umap_feature"),
      content = function(file) {
        p <- traj_umap_feature_plot(); if (is.null(p)) return()
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                        dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
      }
    )

    # ---------- Heatmap: Dominant Status × Feature Intensity ----------
    # Hue per pseudotime bin = dominant donor_status, computed as the
    # argmax of per-status NORMALIZED density (count_in_bin / total_for_status).
    # Normalization corrects for global cohort imbalance
    # (ND ~2300, Aab+ ~1668, T1D ~1246 islets) so dominance reflects
    # enrichment relative to each status's pool, not raw count.
    # Intensity per bin = bin-mean of selected feature, normalized via
    # 5–95th percentile clip, blended in [INTENSITY_FLOOR, 1] toward the
    # donor color so even the dimmest bin is a clear tint (never white).
    # Rendered with ggplot + geom_tile to guarantee bars span [0, 1].
    INTENSITY_FLOOR <- 0.50  # minimum blend factor — bars never go pale-white
    output$traj_heatmap <- renderPlotly({
      df <- traj_data_clean()
      if (is.null(df) || nrow(df) == 0) return(NULL)

      feature_name <- df$feature_name[1]
      donor_cols   <- donor_colors_reactive()
      n_bins       <- 25
      brks         <- seq(0, 1, length.out = n_bins + 1)
      bin_w        <- diff(brks)[1]
      mids         <- head(brks, -1) + bin_w / 2

      df$pt_bin <- cut(pmax(0, pmin(1, df$pt)), breaks = brks,
                       include.lowest = TRUE, right = FALSE)
      bin_levels <- levels(df$pt_bin)

      status_levels <- c("ND", "Aab+", "T1D")
      status_totals <- sapply(status_levels, function(s) sum(df$donor_status == s))

      bin_rows <- lapply(seq_along(bin_levels), function(i) {
        bd <- df[!is.na(df$pt_bin) & df$pt_bin == bin_levels[i], , drop = FALSE]
        if (nrow(bd) == 0) {
          return(data.frame(idx = i, x = mids[i], n = 0L,
                            n_ND = 0L, n_Aab = 0L, n_T1D = 0L,
                            dens_ND = 0, dens_Aab = 0, dens_T1D = 0,
                            dominant = NA_character_, feat_mean = NA_real_,
                            stringsAsFactors = FALSE))
        }
        counts <- c(ND = sum(bd$donor_status == "ND"),
                    `Aab+` = sum(bd$donor_status == "Aab+"),
                    T1D = sum(bd$donor_status == "T1D"))
        dens <- ifelse(status_totals > 0, counts / status_totals, 0)
        names(dens) <- status_levels
        dominant <- names(dens)[which.max(dens)]
        data.frame(idx = i, x = mids[i], n = nrow(bd),
                   n_ND = counts[["ND"]], n_Aab = counts[["Aab+"]], n_T1D = counts[["T1D"]],
                   dens_ND = dens[["ND"]], dens_Aab = dens[["Aab+"]], dens_T1D = dens[["T1D"]],
                   dominant = dominant,
                   feat_mean = mean(bd$value, na.rm = TRUE),
                   stringsAsFactors = FALSE)
      })
      hm <- do.call(rbind, bin_rows)
      if (is.null(hm) || nrow(hm) == 0) return(NULL)

      # Intensity from feature value, robust 5–95th percentile clip
      q   <- stats::quantile(hm$feat_mean, c(0.05, 0.95), na.rm = TRUE)
      rng <- diff(q)
      raw_intensity <- if (is.finite(rng) && rng > 0) {
        pmax(0, pmin(1, (hm$feat_mean - q[1]) / rng))
      } else {
        rep(1, nrow(hm))
      }
      raw_intensity[!is.finite(raw_intensity)] <- INTENSITY_FLOOR
      # Rescale so the lowest bin is INTENSITY_FLOOR-saturated (clearly tinted),
      # the highest bin is fully saturated.
      hm$intensity <- INTENSITY_FLOOR + (1 - INTENSITY_FLOOR) * raw_intensity

      blend_to_white <- function(hex, t) {
        rgb_full <- grDevices::col2rgb(hex) / 255
        grDevices::rgb(red   = t * rgb_full[1] + (1 - t) * 1,
                       green = t * rgb_full[2] + (1 - t) * 1,
                       blue  = t * rgb_full[3] + (1 - t) * 1)
      }
      hm$fill <- mapply(function(dom, ity) {
        if (is.na(dom) || is.null(donor_cols[[dom]])) return("#BBBBBB")
        blend_to_white(donor_cols[[dom]], ity)
      }, hm$dominant, hm$intensity, USE.NAMES = FALSE)

      hm$hover <- sprintf(
        paste0(
          "PT: %.2f<br>Dominant: %s<br>",
          "n: %d (ND %d, Aab+ %d, T1D %d)<br>",
          "Normalized density: ND %.1f%%, Aab+ %.1f%%, T1D %.1f%%<br>",
          "%s mean: %.3f"
        ),
        hm$x, ifelse(is.na(hm$dominant), "n/a", hm$dominant),
        hm$n, hm$n_ND, hm$n_Aab, hm$n_T1D,
        100 * hm$dens_ND, 100 * hm$dens_Aab, 100 * hm$dens_T1D,
        feature_name, hm$feat_mean
      )

      g <- ggplot(hm, aes(x = x, y = 0.5, fill = fill, text = hover)) +
        geom_tile(width = bin_w, height = 1) +
        scale_fill_identity(guide = "none") +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        guides(fill = "none") +
        theme_void() +
        theme(legend.position = "none",
              plot.margin = margin(5, 20, 10, 60, "pt"))

      ggplotly(g, tooltip = "text") %>%
        layout(
          showlegend = FALSE,
          xaxis = list(fixedrange = TRUE, showticklabels = FALSE,
                       zeroline = FALSE, range = c(0, 1)),
          yaxis = list(fixedrange = TRUE, showticklabels = FALSE,
                       zeroline = FALSE, range = c(0, 1))
        ) %>%
        config(displayModeBar = FALSE)
    })

    # ---------- Multi-Feature Heatmap ----------

    # Marker selector UI
    output$multi_heatmap_marker_selector <- renderUI({
      tr <- traj()
      if (is.null(tr) || !is.null(tr$error) || is.null(tr$var_names)) return(NULL)
      checkboxGroupInput(ns("multi_heatmap_markers"), "Select markers:",
                         choices = tr$var_names, selected = tr$var_names, inline = TRUE)
    })

    # Multi-feature data reactive
    traj_multi_feature_data <- reactive({
      markers <- input$multi_heatmap_markers
      if (is.null(markers) || length(markers) == 0) return(NULL)
      tr <- traj_with_pseudotime()
      if (is.null(tr) || !is.null(tr$error) || is.null(tr$adata)) return(NULL)
      valid <- intersect(markers, tr$var_names)
      if (length(valid) == 0) return(NULL)
      idx <- match(valid, tr$var_names)
      expr <- tryCatch(as.matrix(tr$adata$X[, idx, drop = FALSE]), error = function(e) NULL)
      if (is.null(expr)) return(NULL)
      colnames(expr) <- valid
      result <- data.frame(
        pt = as.numeric(tr$obs$pseudotime),
        total_cells = as.integer(tr$obs$total_cells),
        stringsAsFactors = FALSE)
      result <- cbind(result, as.data.frame(expr))
      result <- result[is.finite(result$pt), , drop = FALSE]

      # Apply min cells/islet filter (matches traj_data_clean)
      min_cells <- input$min_cells %||% 1
      if (min_cells > 1) {
        result <- result[is.finite(result$total_cells) & result$total_cells >= min_cells, , drop = FALSE]
        if (nrow(result) == 0) return(NULL)
      }

      # Rank-normalize pseudotime (matches scatter x-axis)
      result$pt <- rank(result$pt, ties.method = "average") / nrow(result)
      result$total_cells <- NULL

      list(data = result, markers = valid)
    })

    # Multi-row heatmap renderPlot
    traj_multi_heatmap_plot <- reactive({
      mfd <- traj_multi_feature_data()
      if (is.null(mfd)) return(NULL)

      df <- mfd$data
      markers <- mfd$markers
      nb <- input$multi_heatmap_nbins %||% 20

      # Bin pseudotime
      pt_range <- range(df$pt, na.rm = TRUE)
      df$pt_norm <- (df$pt - pt_range[1]) / (pt_range[2] - pt_range[1])
      brks <- seq(0, 1, length.out = nb + 1)
      df$pt_bin <- cut(pmax(0, pmin(1, df$pt_norm)), breaks = brks, include.lowest = TRUE, right = FALSE)
      bin_levels <- levels(df$pt_bin)

      # Compute per-bin mean for each marker, then z-score
      all_mids <- head(brks, -1) + diff(brks)[1] / 2
      hm_rows <- list()

      for (mk in markers) {
        bin_means <- numeric(length(bin_levels))
        bin_counts <- integer(length(bin_levels))
        for (i in seq_along(bin_levels)) {
          bin_data <- df[!is.na(df$pt_bin) & df$pt_bin == bin_levels[i], mk]
          bin_counts[i] <- length(bin_data)
          bin_means[i] <- if (length(bin_data) >= 1) mean(bin_data, na.rm = TRUE) else NA_real_
        }
        # Z-score across bins
        valid_means <- bin_means[is.finite(bin_means)]
        if (length(valid_means) < 2) next
        mu <- mean(valid_means)
        sd_val <- sd(valid_means)
        if (!is.finite(sd_val) || sd_val == 0) sd_val <- 1
        z <- (bin_means - mu) / sd_val
        z <- pmax(-2.5, pmin(2.5, z))  # clamp

        for (i in seq_along(bin_levels)) {
          hm_rows[[length(hm_rows) + 1]] <- data.frame(
            marker = mk,
            x = all_mids[i] * (pt_range[2] - pt_range[1]) + pt_range[1],
            z = z[i],
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(hm_rows) == 0) return(NULL)
      hm <- do.call(rbind, hm_rows)

      # Order markers: hormones first, then immune, then others
      hormone_order <- intersect(c("INS", "GCG", "SST"), markers)
      immune_order <- intersect(c("CD3e", "CD4", "CD8a", "CD20", "CD45", "CD68", "CD163", "HLADR"), markers)
      other_order <- setdiff(markers, c(hormone_order, immune_order))
      marker_order <- c(hormone_order, immune_order, other_order)
      hm$marker <- factor(hm$marker, levels = rev(marker_order))

      n_markers <- length(marker_order)

      ggplot(hm, aes(x = x, y = marker, fill = z)) +
        geom_tile() +
        scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                             midpoint = 0, limits = c(-2.5, 2.5),
                             name = "Z-score", na.value = "#f7f7f7") +
        scale_x_continuous(limits = c(pt_range[1], pt_range[2]), expand = c(0, 0)) +
        labs(x = "Pseudotime", y = "", title = "Multi-Feature Expression Along Pseudotime") +
        theme_islet() +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size = max(ISLET_FONT_SIZES$dense_min, ISLET_FONT_SIZES$axis_text - 2))
        )
    })
    output$traj_multi_heatmap <- renderPlot({ traj_multi_heatmap_plot() },
      height = function() {
        n <- length(input$multi_heatmap_markers %||% character(0))
        max(220, 60 + n * 36)
      })
    output$dl_multi_heatmap <- downloadHandler(
      filename = function() islet_png_filename("trajectory_multi_heatmap"),
      content = function(file) {
        p <- traj_multi_heatmap_plot(); if (is.null(p)) return()
        n <- length(input$multi_heatmap_markers %||% character(0))
        h <- max(6, 1 + n * 0.4)
        ggplot2::ggsave(file, plot = p,
                        width = ISLET_HEATMAP_PNG_DIMS$width, height = h,
                        dpi = ISLET_HEATMAP_PNG_DIMS$dpi, units = ISLET_HEATMAP_PNG_DIMS$units)
      }
    )

  })
}
