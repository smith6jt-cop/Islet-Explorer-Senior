#' Wilcoxon-AUC markers via presto.
#'
#' Calls `presto::wilcoxauc` directly on the data layer (genes × cells
#' matrix) and the group vector to avoid presto's Seurat dispatch, which
#' still passes the defunct `slot=` argument to `GetAssayData()` on
#' SeuratObject ≥ 5.0.
#'
#' @param xen   Seurat object (typically `app_state$xen_clustered`).
#' @param group character; either a column name in `xen@meta.data` or a
#'              vector of length `ncol(xen)` of group labels.
#' @param assay assay to pull data from (default = `DefaultAssay(xen)`).
#' @param layer slot/layer to use (default `"data"` after normalization;
#'              fall through to `"counts"` if `"data"` is empty).
#' @return a data.frame of one row per (group, gene) pair with the
#'         standard presto columns: feature, group, avgExpr, logFC,
#'         statistic, auc, pval, padj, pct_in, pct_out.
compute_markers <- function(xen, group,
                            assay = SeuratObject::DefaultAssay(xen),
                            layer = "data") {
  if (length(group) == 1L && group %in% names(xen@meta.data)) {
    y <- xen@meta.data[[group]]
  } else {
    if (length(group) != ncol(xen)) {
      stop("`group` must be a meta.data column name or a vector of length ncol(xen).")
    }
    y <- group
  }
  if (length(unique(y)) < 2L) {
    stop("Marker computation needs at least 2 groups; got: ",
         paste(unique(y), collapse = ", "))
  }

  X <- SeuratObject::GetAssayData(xen, assay = assay, layer = layer)
  if (!nrow(X) || sum(X@x %||% as.numeric(X)) == 0) {
    X <- SeuratObject::GetAssayData(xen, assay = assay, layer = "counts")
  }

  m <- presto::wilcoxauc(X, y)
  m[order(m$group, -m$auc), , drop = FALSE]
}

#' Top-N markers per group, ranked by AUC.
top_markers <- function(markers, n = 10L,
                        min_pct_in = 0.1,
                        max_padj = 0.05) {
  m <- markers
  if (!is.null(m$padj))   m <- m[!is.na(m$padj)   & m$padj   <= max_padj, ]
  if (!is.null(m$pct_in)) m <- m[!is.na(m$pct_in) & m$pct_in >= min_pct_in * 100, ]
  if (!nrow(m)) return(m)
  groups <- split(m, m$group)
  out <- lapply(groups, function(g) {
    g <- g[order(-g$auc), , drop = FALSE]
    head(g, n)
  })
  do.call(rbind, out)
}
