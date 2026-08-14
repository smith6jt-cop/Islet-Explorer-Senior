#' Default options accepted by [run_cluster_pipeline()].
default_cluster_opts <- function() {
  list(
    # Gene set
    subpanels         = character(),    # subpanel keys whose union -> features
    extra_genes       = character(),    # individually added genes
    use_all_features  = FALSE,          # if TRUE, ignore subpanel/extra and use all rownames
    # Cell filters
    nCount_min        = 0,
    nCount_max        = Inf,
    nFeature_min      = 0,
    nFeature_max      = Inf,
    meta_filter_col   = NULL,
    meta_filter_vals  = NULL,
    # Normalization
    norm_method       = "LogNormalize",  # "LogNormalize" / "SCT" / "skip"
    do_scale          = TRUE,
    sct_warn_threshold = 200000L,
    # PCA
    npcs              = 30,
    # Batch
    batch             = "None",          # "None" / "Harmony"
    batch_var         = "orig.ident",
    harmony_theta     = 2,
    # UMAP / neighborhood
    k_param           = 30,
    min_dist          = 0.3,
    n_neighbors       = 30,
    metric            = "cosine",
    # Clustering
    cluster_algorithm = "Louvain",       # "Louvain" / "Leiden"
    resolutions       = seq(0.1, 1.0, by = 0.1),
    # Repro
    seed              = 1L
  )
}

.merge_opts <- function(opts) {
  modifyList(default_cluster_opts(), opts %||% list())
}

# Forward an increment to the user-supplied progress handler, or no-op.
.tick <- function(progress, amount, message, detail = NULL) {
  if (is.function(progress)) {
    progress(amount = amount, message = message, detail = detail)
  }
}

#' Build the variable-feature set used by RunPCA.
#'
#' Union of `subpanels` ∪ `extra_genes`, intersected with the loaded data;
#' if both are empty (or `use_all_features=TRUE`), returns all rownames.
.resolve_features <- function(xen, panels, opts) {
  data_genes <- rownames(xen)
  if (isTRUE(opts$use_all_features)) return(data_genes)

  set <- character()
  if (length(opts$subpanels)) {
    sp <- intersect(opts$subpanels, names(panels$subpanels))
    set <- unique(unlist(lapply(panels$subpanels[sp], `[[`, "gene"),
                          use.names = FALSE))
  }
  if (length(opts$extra_genes)) {
    set <- unique(c(set, opts$extra_genes))
  }
  if (!length(set)) return(data_genes)
  intersect(set, data_genes)
}

#' Apply user cell filters and return the column index to keep.
.resolve_cells <- function(xen, opts) {
  md     <- xen@meta.data
  assay  <- SeuratObject::DefaultAssay(xen)
  ncount_col   <- paste0("nCount_",   assay)
  nfeature_col <- paste0("nFeature_", assay)
  keep <- rep(TRUE, ncol(xen))
  if (ncount_col %in% names(md)) {
    keep <- keep & md[[ncount_col]]   >= opts$nCount_min
    keep <- keep & md[[ncount_col]]   <= opts$nCount_max
  }
  if (nfeature_col %in% names(md)) {
    keep <- keep & md[[nfeature_col]] >= opts$nFeature_min
    keep <- keep & md[[nfeature_col]] <= opts$nFeature_max
  }
  if (!is.null(opts$meta_filter_col) && nzchar(opts$meta_filter_col) &&
      opts$meta_filter_col %in% names(md) &&
      length(opts$meta_filter_vals)) {
    keep <- keep & (md[[opts$meta_filter_col]] %in% opts$meta_filter_vals)
  }
  which(keep)
}

#' Run the M5 cluster pipeline on a Seurat object.
#'
#' Pure function: takes a Seurat object + opts, returns a new Seurat
#' object with PCA + (optional Harmony) + UMAP reductions and
#' `seurat_clusters_res_<r>` metadata columns for every resolution
#' requested. Records the run parameters under
#' `xen@misc$pipeline_history[[run_id]]`.
#'
#' @param xen a Seurat object (typically `app_state$xen`).
#' @param panels output of [load_panels()]; supplies subpanel gene sets.
#' @param opts named list, see [default_cluster_opts()].
#' @param progress optional `function(amount, message, detail = NULL)` invoked
#'   at each major step. Drive this from `shiny::withProgress` /
#'   `shiny::incProgress` to stream step-level progress to the UI; pass `NULL`
#'   (the default) to run silently.
#' @return a new Seurat object.
run_cluster_pipeline <- function(xen, panels, opts = list(), progress = NULL) {
  if (!inherits(xen, "Seurat")) stop("xen must be a Seurat object")
  opts <- .merge_opts(opts)

  set.seed(opts$seed)

  # 1. Cell + feature subset
  .tick(progress, 0.02, "Subsetting cells and features")
  keep_cells <- .resolve_cells(xen, opts)
  if (!length(keep_cells)) stop("Cell filter removed all cells.")
  obj <- xen[, keep_cells]

  feats <- .resolve_features(obj, panels, opts)
  if (length(feats) < 5L) {
    stop("Feature set has <5 genes after filtering; cannot run pipeline.")
  }
  npcs_use <- min(opts$npcs, length(feats) - 1L, ncol(obj) - 1L)

  # 2. Normalization
  used_assay <- SeuratObject::DefaultAssay(obj)
  norm <- opts$norm_method
  if (identical(norm, "SCT")) {
    if (ncol(obj) > opts$sct_warn_threshold) {
      warning(sprintf(
        "SCTransform requested on %d cells (> %d); falling back to LogNormalize.",
        ncol(obj), opts$sct_warn_threshold))
      norm <- "LogNormalize"
    }
  }
  .tick(progress, 0.10, "Normalizing", detail = norm)
  obj <- switch(norm,
    LogNormalize = Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                                         verbose = FALSE),
    SCT          = Seurat::SCTransform(obj, assay = used_assay,
                                       verbose = TRUE,
                                       return.only.var.genes = FALSE),
    skip         = obj,
    stop("unknown norm_method: ", norm)
  )
  active_assay <- if (identical(norm, "SCT")) "SCT" else used_assay
  SeuratObject::DefaultAssay(obj) <- active_assay

  # 3. Variable features + scaling
  SeuratObject::VariableFeatures(obj) <- feats
  if (isTRUE(opts$do_scale) && !identical(norm, "SCT")) {
    .tick(progress, 0.05, "Scaling")
    obj <- Seurat::ScaleData(obj, features = feats, verbose = FALSE)
  }

  # 4. PCA on the variable feature set
  .tick(progress, 0.10, "Running PCA",
        detail = sprintf("%d PCs over %d features", npcs_use, length(feats)))
  obj <- Seurat::RunPCA(obj, features = feats, npcs = npcs_use,
                        verbose = TRUE, seed.use = opts$seed)

  # 5. Optional Harmony batch correction
  reduction_for_neighbors <- "pca"
  if (identical(opts$batch, "Harmony")) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      warning("harmony not installed; skipping batch correction.")
    } else if (!opts$batch_var %in% names(obj@meta.data)) {
      warning(sprintf("batch_var '%s' not in meta.data; skipping Harmony.",
                      opts$batch_var))
    } else if (length(unique(obj@meta.data[[opts$batch_var]])) < 2L) {
      warning(sprintf(
        "batch_var '%s' has only one level; skipping Harmony.",
        opts$batch_var))
    } else {
      .tick(progress, 0.08, "Harmony batch correction",
            detail = sprintf("group.by = %s, theta = %.1f",
                             opts$batch_var, opts$harmony_theta))
      obj <- harmony::RunHarmony(obj,
                                 group.by.vars = opts$batch_var,
                                 theta         = opts$harmony_theta,
                                 verbose       = TRUE)
      reduction_for_neighbors <- "harmony"
    }
  }

  # 6. Neighbors + UMAP
  .tick(progress, 0.05, "Finding neighbors",
        detail = sprintf("k.param = %d", opts$k_param))
  obj <- Seurat::FindNeighbors(obj,
                               reduction = reduction_for_neighbors,
                               dims      = seq_len(npcs_use),
                               k.param   = opts$k_param,
                               verbose   = TRUE)
  .tick(progress, 0.15, "Computing UMAP",
        detail = sprintf("n.neighbors = %d, min.dist = %.2f",
                         opts$n_neighbors, opts$min_dist))
  obj <- Seurat::RunUMAP(obj,
                         reduction   = reduction_for_neighbors,
                         dims        = seq_len(npcs_use),
                         min.dist    = opts$min_dist,
                         n.neighbors = opts$n_neighbors,
                         metric      = opts$metric,
                         seed.use    = opts$seed,
                         verbose     = TRUE)

  # 7. Resolution sweep
  algo_id <- 1L  # Louvain
  if (identical(opts$cluster_algorithm, "Leiden")) {
    if (!requireNamespace("leidenAlg", quietly = TRUE)) {
      warning("leidenAlg not installed; falling back to Louvain.")
    } else {
      algo_id <- 4L
    }
  }

  used_so_far <- 0.02 + 0.10 +
    (if (isTRUE(opts$do_scale) && !identical(norm, "SCT")) 0.05 else 0) +
    0.10 +
    (if (identical(reduction_for_neighbors, "harmony")) 0.08 else 0) +
    0.05 + 0.15
  per_res <- max(0, (1 - used_so_far)) / max(1L, length(opts$resolutions))
  algo_label <- if (algo_id == 4L) "Leiden" else "Louvain"

  for (i in seq_along(opts$resolutions)) {
    r <- opts$resolutions[[i]]
    .tick(progress, per_res, "Clustering",
          detail = sprintf("%s — resolution %.2f (%d/%d)",
                           algo_label, r, i, length(opts$resolutions)))
    obj <- Seurat::FindClusters(obj,
                                resolution = r,
                                algorithm  = algo_id,
                                verbose    = TRUE,
                                random.seed = opts$seed)
    col <- sprintf("seurat_clusters_res_%g", r)
    obj@meta.data[[col]] <- obj@meta.data[["seurat_clusters"]]
  }
  obj@meta.data[["seurat_clusters"]] <- NULL  # drop intermediate

  # 8. Run log
  run_id <- format(Sys.time(), "run_%Y%m%dT%H%M%OS3")
  history <- obj@misc$pipeline_history %||% list()
  history[[run_id]] <- list(
    when         = Sys.time(),
    n_cells_in   = ncol(xen),
    n_cells_out  = ncol(obj),
    n_features   = length(feats),
    npcs         = npcs_use,
    norm_method  = norm,
    batch        = if (identical(reduction_for_neighbors, "harmony")) "Harmony"
                   else "None",
    algorithm    = if (algo_id == 4L) "Leiden" else "Louvain",
    resolutions  = opts$resolutions,
    opts         = opts
  )
  obj@misc$pipeline_history <- history
  obj@misc$last_run_id <- run_id
  obj
}

#' Resolution-column names produced by [run_cluster_pipeline()] in this object.
cluster_resolution_columns <- function(xen) {
  grep("^seurat_clusters_res_", names(xen@meta.data), value = TRUE)
}

#' Numeric resolution value parsed out of `seurat_clusters_res_<r>` names.
cluster_resolution_values <- function(xen) {
  cols <- cluster_resolution_columns(xen)
  as.numeric(sub("^seurat_clusters_res_", "", cols))
}
