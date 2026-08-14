#' Stability summary for the resolution-sweep columns produced by
#' [run_cluster_pipeline()].
#'
#' For each resolution column (in ascending numeric order), reports
#' n_clusters and the agreement of cluster assignments against the
#' next-lower resolution. Agreement is reported with adjusted Rand index
#' (always available via the `mclust::adjustedRandIndex` reimplementation
#' below, with no extra dep) and — when `aricode` is installed — also
#' with normalized mutual information (NMI). When neither metric is
#' computable, the column is filled with NA but n_clusters is still
#' reported.
#'
#' @param xen a Seurat object with `seurat_clusters_res_*` columns.
#' @return data.frame with columns: resolution, n_clusters, ARI_vs_prev,
#'         NMI_vs_prev, frac_changed_vs_prev.
compute_stability_summary <- function(xen) {
  cols <- cluster_resolution_columns(xen)
  if (!length(cols)) {
    return(data.frame(resolution = numeric(),
                      n_clusters = integer(),
                      ARI_vs_prev = numeric(),
                      NMI_vs_prev = numeric(),
                      frac_changed_vs_prev = numeric()))
  }
  vals <- cluster_resolution_values(xen)
  ord  <- order(vals)
  cols <- cols[ord]; vals <- vals[ord]

  use_aricode <- requireNamespace("aricode", quietly = TRUE)

  n_clusters <- vapply(cols, function(c) {
    length(unique(xen@meta.data[[c]]))
  }, integer(1))

  ari <- rep(NA_real_, length(cols))
  nmi <- rep(NA_real_, length(cols))
  fch <- rep(NA_real_, length(cols))

  for (i in seq_along(cols)[-1L]) {
    a <- as.integer(factor(xen@meta.data[[cols[i - 1L]]]))
    b <- as.integer(factor(xen@meta.data[[cols[i]]]))
    ari[i] <- .adjusted_rand(a, b)
    if (use_aricode) nmi[i] <- aricode::NMI(a, b)
    # Best-match relabel of `b` against `a`, then count disagreements.
    fch[i] <- .frac_relabel_changed(a, b)
  }

  data.frame(resolution = vals,
             n_clusters = unname(n_clusters),
             ARI_vs_prev = ari,
             NMI_vs_prev = nmi,
             frac_changed_vs_prev = fch,
             stringsAsFactors = FALSE)
}

#' Adjusted Rand index — vectorised, no external deps.
.adjusted_rand <- function(a, b) {
  if (!length(a) || !length(b) || length(a) != length(b)) return(NA_real_)
  tab <- table(a, b)
  n   <- sum(tab)
  ai  <- rowSums(tab); bj <- colSums(tab)
  choose2 <- function(x) x * (x - 1) / 2
  index    <- sum(choose2(tab))
  expected <- sum(choose2(ai)) * sum(choose2(bj)) / choose2(n)
  max_idx  <- (sum(choose2(ai)) + sum(choose2(bj))) / 2
  if (max_idx == expected) return(NA_real_)
  (index - expected) / (max_idx - expected)
}

#' Best-match relabel of `b` to `a`, then return fraction of cells whose
#' assignment disagrees. Uses a greedy 1:1 cluster matching by cell count.
.frac_relabel_changed <- function(a, b) {
  tab <- table(a, b)
  used_a <- character()
  used_b <- character()
  map <- list()
  pairs <- arrayInd(order(-tab), dim(tab))
  for (k in seq_len(nrow(pairs))) {
    ra <- rownames(tab)[pairs[k, 1]]; cb <- colnames(tab)[pairs[k, 2]]
    if (ra %in% used_a || cb %in% used_b) next
    map[[cb]] <- ra
    used_a <- c(used_a, ra); used_b <- c(used_b, cb)
  }
  b_relabeled <- ifelse(as.character(b) %in% names(map),
                        unlist(map[as.character(b)]),
                        as.character(b))
  mean(b_relabeled != as.character(a))
}

#' Build the data.frame that backs the clustree edge-table download.
clustree_edge_table <- function(xen, prefix = "seurat_clusters_res_") {
  cols <- cluster_resolution_columns(xen)
  if (length(cols) < 2L) return(data.frame())
  ord <- order(cluster_resolution_values(xen))
  cols <- cols[ord]
  rows <- list()
  for (i in seq_along(cols)[-1L]) {
    a <- as.character(xen@meta.data[[cols[i - 1L]]])
    b <- as.character(xen@meta.data[[cols[i]]])
    tab <- as.data.frame(table(from = a, to = b),
                         stringsAsFactors = FALSE)
    tab <- tab[tab$Freq > 0, , drop = FALSE]
    if (!nrow(tab)) next
    src <- sub(prefix, "", cols[i - 1L])
    dst <- sub(prefix, "", cols[i])
    tab$from_resolution <- src
    tab$to_resolution   <- dst
    tab$in_prop <- tab$Freq / ave(tab$Freq, tab$to,   FUN = sum)
    rows[[i]] <- tab[, c("from_resolution","from","to_resolution","to",
                         "Freq","in_prop")]
  }
  do.call(rbind, rows)
}
