#' Demo-mode synthetic Xenium dataset.
#'
#' Mirrors the test fixture in `tests/testthat/helper-test-data.R` but
#' lives in `R/` so the running app can call it without test helpers.
#' Sampled gene names come from the loaded panels, so the panel
#' validation report on the Load Xenium tab shows realistic coverage.

#' Build a synthetic Seurat object for demo mode.
#'
#' @param panels output of [load_panels()]; supplies the gene pool.
#' @param n_cells number of synthetic cells (default 1500).
#' @param n_genes number of genes to sample (default 800).
#' @param seed RNG seed (default 1L) for reproducibility.
make_demo_seurat <- function(panels, n_cells = 1500L, n_genes = 800L,
                              seed = 1L) {
  pool <- unique(unlist(lapply(panels$subpanels, `[[`, "gene"),
                        use.names = FALSE))
  pool <- pool[!is.na(pool) & nzchar(pool)]
  if (!length(pool)) {
    stop("Demo data builder: no genes found in panels$subpanels.")
  }
  set.seed(seed)
  genes <- sample(pool, min(n_genes, length(pool)))

  # ~5% nonzero entries, Poisson(3)+1 counts.
  nnz <- as.integer(0.05 * length(genes) * n_cells)
  i <- sample.int(length(genes), nnz, replace = TRUE)
  j <- sample.int(n_cells, nnz, replace = TRUE)
  x <- stats::rpois(nnz, lambda = 3) + 1L
  M <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(genes), n_cells),
    dimnames = list(genes, sprintf("cell_%05d", seq_len(n_cells)))
  )

  # Two synthetic samples so cluster-by-sample plots have something to show.
  n_a <- floor(n_cells / 2)
  meta <- data.frame(
    cell_id           = colnames(M),
    x_centroid        = stats::runif(n_cells, 0, 1500),
    y_centroid        = stats::runif(n_cells, 0, 1500),
    transcript_counts = Matrix::colSums(M),
    orig.ident        = c(rep("demo_A", n_a),
                          rep("demo_B", n_cells - n_a)),
    stringsAsFactors  = FALSE,
    row.names         = colnames(M)
  )
  obj <- SeuratObject::CreateSeuratObject(
    counts    = M,
    assay     = "Xenium",
    meta.data = meta,
    project   = "demo"
  )
  obj@misc$is_demo <- TRUE
  obj
}

#' Sentinel value `chosen_path()` uses to mean "demo dataset requested".
demo_path_sentinel <- function() "<demo>"

#' Load (or build) the demo Seurat object, with a `qs2` cache.
#'
#' Cache filename is keyed by tissue so swapping tissues doesn't serve a
#' stale fixture sampled from the previous tissue's gene pool.
#'
#' @param panels output of [load_panels()]; passed to [make_demo_seurat()].
#' @param cache_dir cache directory (default `app_paths$cache`).
load_demo_xenium <- function(panels, cache_dir = app_paths$cache) {
  tid <- panels$tissue$id %||% "default"
  cache_file <- file.path(cache_dir, sprintf("demo_xenium_%s.qs2", tid))
  if (file.exists(cache_file)) {
    obj <- qs2::qs_read(cache_file)
    attr(obj, "load_xenium_cache") <- list(hit = TRUE, file = cache_file)
    return(obj)
  }
  obj <- make_demo_seurat(panels)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  qs2::qs_save(obj, cache_file)
  attr(obj, "load_xenium_cache") <- list(hit = FALSE, file = cache_file)
  obj
}
