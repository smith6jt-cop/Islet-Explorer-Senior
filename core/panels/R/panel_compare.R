#' Per-gene detection percentage from a Seurat object.
#'
#' Returns a named numeric vector: percent of cells where the gene's raw
#' count is > 0. Computed once on `Matrix::rowMeans` of the boolean
#' counts matrix; cheap on Xenium-scale data (5K genes × 100k cells).
per_gene_detection_pct <- function(xen) {
  m <- SeuratObject::GetAssayData(xen, layer = "counts")
  100 * Matrix::rowMeans(m > 0)
}

#' Per-gene mean expression from raw counts.
per_gene_mean_expr <- function(xen) {
  m <- SeuratObject::GetAssayData(xen, layer = "counts")
  Matrix::rowMeans(m)
}

#' Cross-tabulate every subpanel against a loaded Xenium dataset.
#'
#' For every subpanel (and the residuals + custom + 5K reference) compute:
#'   * n_genes — declared in the subpanel
#'   * n_present — present in the dataset
#'   * n_passing — present AND detection ≥ `min_detection_pct`
#'   * frac_present — n_present / n_genes
#'   * frac_passing — n_passing / n_genes
#'   * mean_detection_pct — across present genes
#'   * median_detection_pct — across present genes
#'   * mean_expr — mean per-cell raw count, across present genes
#'   * median_expr — median per-cell raw count, across present genes
#'
#' @param panels output of [load_panels()]
#' @param xen Seurat object loaded via [load_xenium()]
#' @param min_detection_pct numeric ≥ 0, the per-gene detection cutoff
#'        applied to the loaded data
compute_subpanel_coverage <- function(panels, xen, min_detection_pct = 0,
                                      custom_label = NULL) {
  det <- per_gene_detection_pct(xen)
  mu  <- per_gene_mean_expr(xen)
  data_genes <- names(det)

  groups <- panels$subpanels
  if (!is.null(panels$custom) && nrow(panels$custom) > 0L) {
    lbl <- custom_label %||% custom_panel_label(NULL, panels)
    groups <- c(groups, stats::setNames(
      list(data.frame(gene = panels$custom$gene)), lbl))
  }
  # xenium5k_in_audit row reflects the constant 5K reference (not the
  # audit table, which may carry extra non-5K custom-panel genes).
  groups <- c(groups,
              list(xenium5k_in_audit = data.frame(gene = panels$reference_5k$gene)))

  rows <- lapply(names(groups), function(nm) {
    gset <- unique(as.character(groups[[nm]]$gene))
    present  <- intersect(gset, data_genes)
    pass_set <- present[det[present] >= min_detection_pct]

    data.frame(
      subpanel             = nm,
      n_genes              = length(gset),
      n_present            = length(present),
      n_passing            = length(pass_set),
      frac_present         = if (length(gset)) length(present) / length(gset)
                              else NA_real_,
      frac_passing         = if (length(gset)) length(pass_set) / length(gset)
                              else NA_real_,
      mean_detection_pct   = if (length(present)) mean(det[present])    else NA_real_,
      median_detection_pct = if (length(present)) stats::median(det[present]) else NA_real_,
      mean_expr            = if (length(present)) mean(mu[present])     else NA_real_,
      median_expr          = if (length(present)) stats::median(mu[present]) else NA_real_,
      stringsAsFactors     = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Top-N genes per subpanel ranked by detection percentage.
#'
#' Returns a long data.frame (subpanel, gene, detection_pct, mean_expr)
#' for use as a `geom_tile` heatmap source.
top_n_by_subpanel <- function(panels, xen, n = 10L,
                              min_detection_pct = 0) {
  det <- per_gene_detection_pct(xen)
  mu  <- per_gene_mean_expr(xen)
  data_genes <- names(det)

  rows <- lapply(names(panels$subpanels), function(nm) {
    gset    <- unique(as.character(panels$subpanels[[nm]]$gene))
    present <- intersect(gset, data_genes)
    if (!length(present)) return(NULL)
    pass    <- present[det[present] >= min_detection_pct]
    if (!length(pass)) return(NULL)
    o <- order(det[pass], decreasing = TRUE)
    pick <- pass[head(o, n)]
    data.frame(
      subpanel       = nm,
      gene           = pick,
      detection_pct  = unname(det[pick]),
      mean_expr      = unname(mu[pick]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}
