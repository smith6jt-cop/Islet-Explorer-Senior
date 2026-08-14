# Synthetic Xenium-like fixtures used by every test that needs a Seurat
# object (M3 onward). Built once per testthat session in tempdir().
#
# Per CLAUDE.md §10: "Generate the synthetic test object inside
# tests/testthat/helper-test-data.R using the actual gene names from
# subpanel_summary_v2.csv (sample ~500 genes across subpanels) so
# panel-validation logic is exercised."

#' Resolve the absolute path to the project root.
test_project_root <- function() {
  rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))
}

#' Resolve the absolute path to the on-disk reference 5K directory.
test_reference_5k_dir <- function() {
  file.path(test_project_root(), "data", "reference_5k")
}

#' Resolve the absolute path to the on-disk shared subpanels directory.
test_shared_subpanels_dir <- function() {
  file.path(test_project_root(), "data", "subpanels_shared")
}

#' Resolve the absolute path to the on-disk tissues root.
test_tissues_root <- function() {
  file.path(test_project_root(), "data", "tissues")
}

#' Resolve the absolute path to the on-disk tissue with the given
#' `tissue_id`. The manifest's `tissue_id` is authoritative (folder names
#' are arbitrary), so this consults the manifest index, not the directory
#' name. Returns NA when the tissue is not registered.
test_tissue_dir <- function(tissue_id = "pancreas") {
  idx <- tissues_index(test_tissues_root())
  if (!(tissue_id %in% names(idx))) return(NA_character_)
  unname(idx[tissue_id])
}

#' Test predicate: is a tissue with this `tissue_id` present on disk?
test_tissue_present <- function(tissue_id = "pancreas") {
  d <- test_tissue_dir(tissue_id)
  !is.na(d) && dir.exists(d)
}

#' Load the on-disk pancreas tissue's panels with the legacy composer
#' shape. Skips the test if the tissue directory is missing.
test_load_panels <- function(tissue_id = "pancreas") {
  load_panels(
    tissue_id              = tissue_id,
    reference_5k_path      = test_reference_5k_dir(),
    subpanels_shared_path  = test_shared_subpanels_dir(),
    tissues_root           = test_tissues_root()
  )
}

#' Build a minimal tissue tree on the fly under `root`.
#'
#' Useful for tests that need a tissue layout without touching the real
#' data/tissues/ tree. Produces:
#'   <root>/reference_5k/xenium5k_genes.csv  (3 genes, biology cols)
#'   <root>/subpanels_shared/01_shared_demo.csv  (gene + category)
#'   <root>/tissues/<tissue>/manifest.yml
#'   <root>/tissues/<tissue>/subpanel_summary.csv
#'   <root>/tissues/<tissue>/subpanels/02_tissue_demo.csv
#'   <root>/tissues/<tissue>/audit/xenium5k_audit.csv
#'   <root>/tissues/<tissue>/custom_panel.csv  (when has_custom = TRUE)
make_test_tissue <- function(root        = tempfile("xen_tissue_"),
                              tissue      = "demo",
                              has_custom  = TRUE,
                              ref_runs    = c("run_alpha", "run_beta")) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  ref_dir <- file.path(root, "reference_5k")
  shared_dir <- file.path(root, "subpanels_shared")
  tdir <- file.path(root, "tissues", tissue)
  for (d in c(ref_dir, shared_dir,
              file.path(tdir, "subpanels"),
              file.path(tdir, "audit"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # Reference 5K (constant) — 6 genes, biology only.
  ref5k <- data.frame(
    gene             = c("AAA", "BBB", "CCC", "DDD", "EEE", "FFF"),
    gene_id          = sprintf("ENS%06d", seq_len(6)),
    full_name        = sprintf("gene %s full", c("AAA","BBB","CCC","DDD","EEE","FFF")),
    location         = "Cytoplasm",
    cell_type        = "test cell",
    cellchat_pathway = NA_character_,
    stringsAsFactors = FALSE
  )
  data.table::fwrite(ref5k, file.path(ref_dir, "xenium5k_genes.csv"), na = "")

  # One shared subpanel, biology-only.
  shared_df <- data.frame(
    gene                   = c("AAA", "BBB", "CCC"),
    category               = "shared_demo",
    rationale              = "test",
    single_gene_panel_note = NA_character_,
    stringsAsFactors       = FALSE
  )
  data.table::fwrite(shared_df, file.path(shared_dir, "01_shared_demo.csv"), na = "")

  # One tissue-specific subpanel, biology-only.
  tissue_df <- data.frame(
    gene                   = c("DDD", "EEE"),
    category               = sprintf("%s_demo", tissue),
    rationale              = "test",
    single_gene_panel_note = NA_character_,
    stringsAsFactors       = FALSE
  )
  data.table::fwrite(tissue_df,
                    file.path(tdir, "subpanels", "02_tissue_demo.csv"),
                    na = "")

  # Audit table — detection_pct columns for both reference runs.
  audit_df <- data.frame(
    gene = ref5k$gene,
    exclude_recommended = "no",
    stringsAsFactors = FALSE
  )
  for (run in ref_runs) {
    audit_df[[paste0("detection_pct_", run)]] <- c(10, 20, 30, 40, 50, 60)
    audit_df[[paste0("n_cells_", run)]]       <- c(100, 200, 300, 400, 500, 600)
  }
  audit_df$log2_detection_ratio <- log2(0.9)
  data.table::fwrite(audit_df,
                    file.path(tdir, "audit", "xenium5k_audit.csv"),
                    na = "")

  # Per-tissue subpanel summary
  summary_df <- data.frame(
    subpanel = c("01_shared_demo", "02_tissue_demo"),
    n_genes  = c(3L, 2L),
    description = c("shared demo subpanel", sprintf("%s demo subpanel", tissue)),
    pct_of_5K_kept_pool = c(0.06, 0.04),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(summary_df,
                    file.path(tdir, "subpanel_summary.csv"), na = "")

  custom_block <- if (has_custom) {
    custom_df <- data.frame(
      gene                   = c("AAA", "FFF"),
      category               = "custom_demo",
      rationale              = "test",
      single_gene_panel_note = NA_character_,
      stringsAsFactors       = FALSE
    )
    data.table::fwrite(custom_df,
                      file.path(tdir, "custom_panel.csv"), na = "")
    paste0(
      'custom_panel:\n',
      '  file: custom_panel.csv\n',
      sprintf('  display_name: "%s custom"\n', tissue))
  } else ''

  manifest <- paste0(
    sprintf('tissue_id: %s\n', tissue),
    sprintf('display_name: "%s test"\n', tissue),
    'reference_runs:\n',
    paste(sprintf('  - "%s"', ref_runs), collapse = "\n"), "\n",
    custom_block
  )
  writeLines(manifest, con = file.path(tdir, "manifest.yml"))

  invisible(list(root = root, tissue = tissue,
                 ref5k_dir = ref_dir, shared_dir = shared_dir,
                 tissue_dir = tdir))
}

#' Sample ~500 unique gene names spanning multiple subpanels.
make_test_genes <- function(panels = NULL, n = 500L, seed = 1L) {
  if (is.null(panels)) panels <- test_load_panels()
  set.seed(seed)
  pool <- unique(unlist(lapply(panels$subpanels, `[[`, "gene"),
                        use.names = FALSE))
  pool <- pool[!is.na(pool) & nzchar(pool)]
  sample(pool, min(n, length(pool)))
}

#' Build a synthetic Seurat object for testing.
#'
#' n_cells × n_genes Poisson counts; centroid coords on a 1000×1000 grid
#' jittered slightly. Genes are real audit names so panel_validate sees
#' realistic input.
make_test_seurat <- function(n_cells = 800L, n_genes = 500L, seed = 1L,
                              panels = NULL) {
  if (is.null(panels)) panels <- test_load_panels()
  genes <- make_test_genes(panels, n = n_genes, seed = seed)
  set.seed(seed)
  # ~5% non-zero entries
  nnz <- as.integer(0.05 * length(genes) * n_cells)
  i <- sample.int(length(genes), nnz, replace = TRUE)
  j <- sample.int(n_cells, nnz, replace = TRUE)
  x <- stats::rpois(nnz, lambda = 3) + 1L
  M <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(genes), n_cells),
    dimnames = list(genes, sprintf("cell_%05d", seq_len(n_cells)))
  )
  meta <- data.frame(
    cell_id    = colnames(M),
    x_centroid = stats::runif(n_cells, 0, 1000),
    y_centroid = stats::runif(n_cells, 0, 1000),
    transcript_counts = Matrix::colSums(M),
    stringsAsFactors = FALSE,
    row.names = colnames(M)
  )
  SeuratObject::CreateSeuratObject(
    counts    = M,
    assay     = "Xenium",
    meta.data = meta,
    project   = "synthetic_test"
  )
}

#' Save the synthetic Seurat object as a `.rds` and return the path.
make_test_rds <- function(dir = tempfile("xen_rds_"), ...) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  obj <- make_test_seurat(...)
  out <- file.path(dir, "synthetic_xenium.rds")
  saveRDS(obj, out)
  out
}

#' Write a synthetic Xenium output-bundle directory.
#'
#' Produces `cell_feature_matrix.h5` (10x layout) + `cells.parquet`
#' inside `dir`. Returns the directory path.
make_test_bundle <- function(dir = tempfile("xen_bundle_"), ...) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  obj <- make_test_seurat(...)

  counts <- SeuratObject::GetAssayData(obj, layer = "counts")
  features <- rownames(counts)
  barcodes <- colnames(counts)

  # 10x cell_feature_matrix.h5 expects CSC: indices = row indices, indptr
  # over columns. Matrix::sparseMatrix(repr="C") → dgCMatrix with @i (row
  # 0-based) and @p (column ptrs).
  csc <- methods::as(counts, "CsparseMatrix")
  h5  <- file.path(dir, "cell_feature_matrix.h5")
  if (file.exists(h5)) file.remove(h5)
  rhdf5::h5createFile(h5)
  rhdf5::h5createGroup(h5, "matrix")
  rhdf5::h5write(as.numeric(csc@x),       h5, "matrix/data")
  rhdf5::h5write(as.integer(csc@i),       h5, "matrix/indices")
  rhdf5::h5write(as.integer(csc@p),       h5, "matrix/indptr")
  rhdf5::h5write(as.integer(dim(csc)),    h5, "matrix/shape")
  rhdf5::h5write(as.character(barcodes),  h5, "matrix/barcodes")
  rhdf5::h5createGroup(h5, "matrix/features")
  rhdf5::h5write(as.character(features),  h5, "matrix/features/name")
  rhdf5::h5write(as.character(features),  h5, "matrix/features/id")
  rhdf5::h5closeAll()

  cells <- obj@meta.data
  cells$cell_id <- rownames(cells)
  arrow::write_parquet(cells, file.path(dir, "cells.parquet"))
  dir
}
