#' Build a per-cell cluster-assignment data.frame for export.
#'
#' Pulls every `seurat_clusters_res_*` column plus key per-cell metadata
#' (cell_id / barcode / x_centroid / y_centroid / orig.ident if present)
#' into a tidy CSV-shaped data.frame.
build_cluster_csv <- function(xen) {
  if (is.null(xen)) {
    stop("No clustered Seurat object available — run the Cluster pipeline first.")
  }
  res_cols <- cluster_resolution_columns(xen)
  if (!length(res_cols)) {
    stop("No seurat_clusters_res_* columns on this object.")
  }
  md <- xen@meta.data
  keep_cols <- intersect(c("cell_id", "x_centroid", "y_centroid",
                           "orig.ident",
                           paste0("nCount_",   SeuratObject::DefaultAssay(xen)),
                           paste0("nFeature_", SeuratObject::DefaultAssay(xen))),
                         names(md))
  out <- data.frame(cell = colnames(xen),
                    stringsAsFactors = FALSE)
  for (k in keep_cols) out[[k]] <- md[[k]]
  for (k in res_cols)  out[[k]] <- as.character(md[[k]])
  out
}

#' Combine a list of marker tables (one per (object,column) key) into a
#' single long data.frame. Keys are split into source / group_col columns.
build_markers_csv <- function(markers_cache) {
  if (!length(markers_cache)) {
    stop("Markers cache is empty — compute markers in the Markers tab first.")
  }
  rows <- lapply(names(markers_cache), function(k) {
    parts <- strsplit(k, "::", fixed = TRUE)[[1]]
    df <- markers_cache[[k]]
    df$source_run_id <- parts[1]
    df$group_col     <- parts[2]
    df
  })
  do.call(rbind, rows)
}

#' Render the session-summary report as a self-contained HTML file.
#'
#' Avoids pandoc / quarto entirely — uses htmltools to build the page.
#' Sections are added only for state that exists (loaded dataset,
#' clustered object, subcluster stack, marker tables).
render_session_report <- function(file, panels, app_state,
                                  title = "Xenium Panel Explorer — session report",
                                  custom_label = NULL) {
  if (is.null(custom_label)) {
    custom_label <- custom_panel_label(app_state$custom_panel_status, panels)
  }
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  card <- function(title, ...) {
    htmltools::div(
      style = "border:1px solid #dee2e6;border-radius:6px;padding:12px 16px;margin:12px 0;",
      htmltools::h3(title), ...)
  }
  kv <- function(...) {
    pairs <- list(...)
    htmltools::tags$dl(class = "row",
      lapply(seq(1, length(pairs), by = 2), function(i) {
        htmltools::tagList(
          htmltools::tags$dt(style = "display:inline-block;width:200px;font-weight:600;",
                              pairs[[i]]),
          htmltools::tags$dd(style = "display:inline-block;",
                              pairs[[i + 1L]]))
      }))
  }
  df_to_html <- function(df, max_rows = 50L) {
    if (is.null(df) || !nrow(df)) return(htmltools::p(htmltools::em("(empty)")))
    df <- utils::head(df, max_rows)
    rows <- htmltools::tagList(lapply(seq_len(nrow(df)), function(i) {
      htmltools::tags$tr(lapply(df[i, ], function(v)
        htmltools::tags$td(format(v, scientific = FALSE))))
    }))
    htmltools::tags$table(
      class = "table table-sm table-striped",
      style = "border-collapse:collapse;",
      htmltools::tags$thead(
        htmltools::tags$tr(lapply(names(df),
          function(n) htmltools::tags$th(n)))),
      htmltools::tags$tbody(rows))
  }

  body <- htmltools::tagList()

  # 1. Audit overview
  s <- panels$summary
  custom_n <- if (!is.null(panels$custom)) nrow(panels$custom) else 0L
  body <- htmltools::tagAppendChildren(body, card("Audit overview",
    kv("Tissue",          panels$tissue$display_name %||% "(unknown)",
       "Generated",       ts,
       "Subpanels",       nrow(s),
       "5K reference genes", nrow(panels$reference_5k),
       sprintf("%s genes", custom_label), custom_n,
       "Excluded genes",  nrow(panels$excluded)),
    htmltools::p(htmltools::strong("Top by gene count: ")),
    df_to_html(s[order(-s$n_genes), c("subpanel","n_genes","description")],
               max_rows = 10L)
  ))

  # 2. Loaded dataset + validation
  if (!is.null(app_state$xen)) {
    rep <- panel_validate(rownames(app_state$xen), panels)
    body <- htmltools::tagAppendChildren(body, card("Loaded dataset",
      kv("Path",     app_state$xen_path %||% "(unknown)",
         "Cells",    ncol(app_state$xen),
         "Genes",    nrow(app_state$xen),
         "Validation", panel_validate_summary(rep,
                                              custom_label = custom_label))))
  }

  # 3. Cluster pipeline run log
  if (!is.null(app_state$xen_clustered)) {
    xc <- app_state$xen_clustered
    h  <- xc@misc$pipeline_history[[xc@misc$last_run_id]]
    if (!is.null(h)) {
      body <- htmltools::tagAppendChildren(body, card("Cluster pipeline",
        kv("Run id",        xc@misc$last_run_id,
           "When",          format(h$when, "%Y-%m-%d %H:%M:%S %Z"),
           "Cells out",     h$n_cells_out,
           "Features",      h$n_features,
           "PCs",           h$npcs,
           "Norm method",   h$norm_method,
           "Batch",         h$batch,
           "Algorithm",     h$algorithm,
           "Resolutions",   paste(round(h$resolutions, 3), collapse = ", "))))
    }
  }

  # 4. Subcluster stack
  stk <- app_state$cluster_stack %||% list()
  if (length(stk) > 1L) {
    rows <- do.call(rbind, lapply(seq_along(stk), function(i) {
      e <- stk[[i]]
      data.frame(level = i,
                 label = e$label %||% "(unnamed)",
                 cells = ncol(e$obj),
                 res_cols = length(cluster_resolution_columns(e$obj)),
                 stringsAsFactors = FALSE)
    }))
    body <- htmltools::tagAppendChildren(body, card("Subcluster stack",
      df_to_html(rows, max_rows = nrow(rows))))
  }

  # 5. Markers
  cache <- app_state$markers_cache %||% list()
  if (length(cache)) {
    last_key <- app_state$markers_last_key %||% names(cache)[1]
    m        <- cache[[last_key]]
    if (!is.null(m) && nrow(m)) {
      top <- top_markers(m, n = 10L, min_pct_in = 0.05, max_padj = 0.05)
      body <- htmltools::tagAppendChildren(body, card(
        sprintf("Top-10 markers per cluster — %s", last_key),
        df_to_html(top[, intersect(c("group", "feature", "auc", "logFC",
                                     "padj", "pct_in", "pct_out"),
                                   names(top))],
                   max_rows = 200L)))
    }
  }

  # 6. Session info
  body <- htmltools::tagAppendChildren(body, card("Session",
    htmltools::pre(paste(utils::capture.output(utils::sessionInfo()),
                          collapse = "\n"))))

  page <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title(title),
      htmltools::tags$style(htmltools::HTML("
        body { font-family: -apple-system, system-ui, sans-serif;
               max-width: 980px; margin: 1rem auto; padding: 0 1rem;
               color: #212529; }
        h1 { border-bottom: 1px solid #dee2e6; padding-bottom: 0.5rem; }
        h3 { margin-top: 0; color: #495057; }
        table { width: 100%; }
        td, th { padding: 4px 8px; border-bottom: 1px solid #f1f3f5;
                 font-size: 0.9rem; }
        pre { background: #f8f9fa; padding: 8px; border-radius: 4px;
              overflow-x: auto; font-size: 0.8rem; }
        dl { margin-bottom: 0; }
        dl dt, dl dd { padding: 2px 0; }
      "))),
    htmltools::tags$body(
      htmltools::h1(title),
      htmltools::p(htmltools::em("Generated ", ts)),
      body))

  htmltools::save_html(page, file, libdir = NULL)
  invisible(file)
}
