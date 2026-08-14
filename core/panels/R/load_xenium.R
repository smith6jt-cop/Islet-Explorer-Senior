#' Ingest a Xenium dataset into a Seurat object.
#'
#' Three entry points (see CLAUDE.md §2):
#' \enumerate{
#'   \item Path to a Xenium output bundle directory containing
#'         `cell_feature_matrix.h5`, `cells.parquet` (or `cells.csv.gz`),
#'         and optionally `transcripts.parquet`.
#'   \item Path to a saved Seurat object (`.rds` or `.qs2`).
#'   \item Path to an AnnData export (`.h5ad`). Read with `zellkonverter`
#'         and converted to a Seurat object via `Seurat::as.Seurat`.
#'         Spatial coordinates are lifted from `reducedDim(sce, "spatial")`
#'         (the SpatialData / Xenium convention used by these h5ads) into
#'         `meta.data$x_centroid` / `meta.data$y_centroid`, which the rest
#'         of the app expects.
#' }
#'
#' Caches the resulting Seurat object under `cache/` keyed by the input
#' path's absolute path + mtime via `qs2::qs_save`. Reload of the same
#' input is a fast `qs_read`.
#'
#' @param path filesystem path (file or directory)
#' @param cache_dir cache directory (default `app_paths$cache`)
#' @param refresh logical; bypass cache (default FALSE)
#' @return a Seurat object (`xen`).
load_xenium <- function(path, cache_dir = app_paths$cache, refresh = FALSE) {
  if (!file.exists(path)) stop("input path does not exist: ", path)
  apath <- normalizePath(path, mustWork = TRUE)
  mtime <- as.numeric(file.mtime(apath))
  key   <- substr(rlang::hash(c(apath, mtime)), 1, 16)
  cache_file <- file.path(cache_dir, paste0("xen_", key, ".qs2"))

  if (!refresh && file.exists(cache_file)) {
    obj <- qs2::qs_read(cache_file)
    attr(obj, "load_xenium_cache") <- list(hit = TRUE, file = cache_file)
    return(obj)
  }

  obj <- if (dir.exists(apath)) {
    .load_xenium_bundle(apath)
  } else {
    ext <- tolower(tools::file_ext(apath))
    switch(ext,
           rds  = readRDS(apath),
           qs   = qs2::qs_read(apath),
           qs2  = qs2::qs_read(apath),
           h5ad = .load_xenium_h5ad(apath),
           stop("unsupported file extension: .", ext,
                " (expected .rds, .qs2, .h5ad, or a Xenium bundle directory)"))
  }

  if (!inherits(obj, "Seurat")) {
    stop("ingest did not return a Seurat object (got ", class(obj)[1], ")")
  }

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  qs2::qs_save(obj, cache_file)
  attr(obj, "load_xenium_cache") <- list(hit = FALSE, file = cache_file)
  obj
}

#' Read a 10x Xenium output-bundle directory into a Seurat object.
#'
#' Required files in the bundle:
#'   * `cell_feature_matrix.h5` — sparse counts (genes × cells).
#'   * `cells.parquet` (or `cells.csv.gz`) — per-cell metadata with
#'     centroid `x_centroid` / `y_centroid` columns and `cell_id`.
#' Optional:
#'   * `transcripts.parquet` — left on disk; not loaded into the object
#'     (we keep its path in `xen@misc$transcripts_parquet` for downstream
#'     modules that may want to stream it).
.load_xenium_bundle <- function(dir) {
  h5_path <- file.path(dir, "cell_feature_matrix.h5")
  if (!file.exists(h5_path)) {
    stop("Xenium bundle missing cell_feature_matrix.h5 in ", dir)
  }

  counts <- .read_cell_feature_matrix_h5(h5_path)

  cells <- .read_cells_table(dir)
  # Align cells to the count matrix's columns (barcodes).
  cell_ids <- as.character(cells$cell_id)
  if (!all(colnames(counts) %in% cell_ids)) {
    # Fall back to row order if cell_id values don't line up
    if (nrow(cells) == ncol(counts)) {
      cells$cell_id <- colnames(counts)
    } else {
      stop("cells.parquet rows (", nrow(cells), ") do not match h5 cells (",
           ncol(counts), ").")
    }
  } else {
    cells <- cells[match(colnames(counts), as.character(cells$cell_id)), ,
                   drop = FALSE]
  }

  obj <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay  = "Xenium",
    meta.data = as.data.frame(cells, stringsAsFactors = FALSE),
    project = basename(dir)
  )

  trans <- file.path(dir, "transcripts.parquet")
  obj@misc$xenium_bundle_dir <- dir
  if (file.exists(trans)) {
    obj@misc$transcripts_parquet <- trans
  }
  obj
}

#' Read a 10x cell_feature_matrix.h5 into a sparse dgCMatrix.
.read_cell_feature_matrix_h5 <- function(h5_path) {
  # Standard 10x layout under /matrix:
  #   data, indices, indptr (CSC), shape, barcodes, features/name (or id)
  on.exit(try(rhdf5::h5closeAll(), silent = TRUE), add = TRUE)
  rh <- function(name) rhdf5::h5read(h5_path, name)

  data    <- as.numeric(rh("matrix/data"))
  indices <- as.integer(rh("matrix/indices"))
  indptr  <- as.integer(rh("matrix/indptr"))
  shape   <- as.integer(rh("matrix/shape"))   # c(n_features, n_cells)
  barcodes <- as.character(rh("matrix/barcodes"))

  ls_root <- rhdf5::h5ls(h5_path, recursive = FALSE)
  feat_paths <- rhdf5::h5ls(h5_path)$name
  features <- if ("name" %in% rhdf5::h5ls(h5_path,
                                          datasetinfo = FALSE)$name) {
    as.character(rh("matrix/features/name"))
  } else if ("id" %in% rhdf5::h5ls(h5_path, datasetinfo = FALSE)$name) {
    as.character(rh("matrix/features/id"))
  } else {
    as.character(rh("matrix/features"))
  }
  if (length(features) != shape[1]) {
    stop("feature count (", length(features),
         ") doesn't match shape[1] (", shape[1], ") in ", h5_path)
  }

  Matrix::sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = data,
    dims = shape,
    dimnames = list(features, barcodes),
    repr = "C"
  )
}

#' Read an AnnData `.h5ad` Xenium export into a Seurat object.
#'
#' Uses `zellkonverter::readH5AD` to produce a `SingleCellExperiment`,
#' picks the `counts` assay (preferred) or the first available assay,
#' converts to Seurat, renames the assay to `Xenium` so the rest of the
#' app's `nCount_Xenium` / `nFeature_Xenium` references line up, and
#' copies the SpatialData `spatial` reducedDim into `meta.data` as
#' `x_centroid` / `y_centroid`. Pre-existing `X_pca` / `X_umap`
#' reductions are preserved as `pca` / `umap` so the user can inspect
#' them without rerunning the pipeline.
.load_xenium_h5ad <- function(h5ad_path) {
  for (pkg in c("zellkonverter", "SingleCellExperiment", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("h5ad ingest requires the `", pkg, "` package; install it first.")
    }
  }

  sce <- zellkonverter::readH5AD(h5ad_path, reader = "R", use_hdf5 = FALSE)
  assay_names <- SummarizedExperiment::assayNames(sce)
  counts_name <- if ("counts" %in% assay_names) "counts" else assay_names[1]

  counts_mat <- SummarizedExperiment::assay(sce, counts_name)
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- methods::as(counts_mat, "CsparseMatrix")
  }
  md <- as.data.frame(SingleCellExperiment::colData(sce), stringsAsFactors = FALSE)

  obj <- SeuratObject::CreateSeuratObject(
    counts    = counts_mat,
    assay     = "Xenium",
    meta.data = md,
    project   = sub("\\.h5ad$", "", basename(h5ad_path))
  )

  rd_names <- SingleCellExperiment::reducedDimNames(sce)
  if ("spatial" %in% rd_names) {
    spat <- SingleCellExperiment::reducedDim(sce, "spatial")
    if (ncol(spat) >= 2 && nrow(spat) == ncol(obj)) {
      obj@meta.data$x_centroid <- as.numeric(spat[, 1])
      obj@meta.data$y_centroid <- as.numeric(spat[, 2])
    }
  }

  obj@misc$h5ad_source   <- h5ad_path
  obj@misc$h5ad_assays   <- assay_names
  obj@misc$h5ad_redDims  <- rd_names
  obj
}

#' Read cells.parquet (or cells.csv.gz) from a Xenium bundle.
.read_cells_table <- function(dir) {
  parq <- file.path(dir, "cells.parquet")
  csv  <- file.path(dir, "cells.csv.gz")
  if (file.exists(parq)) {
    df <- arrow::read_parquet(parq)
  } else if (file.exists(csv)) {
    df <- data.table::fread(csv, data.table = FALSE)
  } else {
    stop("Xenium bundle missing cells.parquet or cells.csv.gz in ", dir)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!"cell_id" %in% names(df)) {
    if ("barcode" %in% names(df)) df$cell_id <- df$barcode
    else stop("cells table missing `cell_id` column (have: ",
              paste(names(df), collapse = ", "), ")")
  }
  df
}
