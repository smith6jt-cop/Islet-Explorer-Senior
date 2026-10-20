#!/usr/bin/env Rscript
# Arm C `tacit` -- the R half.  Driven by arms/tacit.py; never run by hand on cohort data.
#
#   Rscript tacit_run.R <scratch_dir>
#
# Reads   <scratch_dir>/input.parquet    object_id, in_subsample, <identity markers> (complete rows only,
#                                        already left-censored at 0 by the driver)
#         <scratch_dir>/signature.csv    TYPExMARKER 0/1 with the first column literally `cell_type`
#         <scratch_dir>/params.json      seed, p, r (null = calibrate), r_ladder, band_lo, band_hi
# Writes  <scratch_dir>/scores.parquet   object_id, in_subsample, cluster, tacit_label, score__<type>, pos__<type>
#         <scratch_dir>/learned.json     thresholds, r used and how, ladder table, checks, versions, timings
#         <scratch_dir>/tacit_stdout.log everything TACIT prints
#         <scratch_dir>/failure.json     only when the package cannot complete (exit status 2)
#
# WHY the shape is what it is
#   * TACIT(data_expression, r, p, Signature) learns one cut per type by (i) CLR -> scaled -> UMAP(p) ->
#     Louvain(r) microclusters over the cells it is given, (ii) the per-type relevance score
#     `ct = data_anb %*% t(Signature)` on the RAW input (the package does NOT CLR-normalise `ct`; CLR
#     feeds only the clustering path), (iii) a segmented fit over the ordered microcluster medians of
#     `ct` that splits clusters into a high and a low group, and (iv) the cut in [q05, q99] of the
#     per-cell score that minimises the two misclassification rates between those groups.  The cut is
#     therefore a property of the score's own bimodality in microcluster space -- no expected frequency
#     enters.  We call the package on the subsample and then apply the SAME expression for `ct` and the
#     SAME strict `>` comparison (threshold_function's last line) to every cell, so the whole-donor call
#     is TACIT's rule extended, not a re-derivation.
#   * `r` is calibrated on the package's own pre-processing.  The README says only "microcluster sizes
#     averaging 0.1-0.5 % of cells"; the resolution that achieves that depends on n and on the graph, so
#     we replicate TACIT's pre-processing verbatim (same calls, same fixed seeds: RunUMAP seed.use = 42,
#     FindClusters random.seed = 1), run FindClusters along a fixed ladder and take the FIRST ladder
#     value whose mean cluster size is <= band_hi.  Because the pre-processing is deterministic the
#     clusters TACIT then builds internally are the same (checked: `clusters_identical_to_prepass`).
#     If the ladder jumps over the band the choice is recorded as out_of_band; nothing else is tuned.
#   * Consistency check: the binaries we recompute for the subsample must equal TACIT's own
#     Group_threshold_data cell for cell (`n_binary_mismatch_vs_tacit` must be 0).
#   * Known package edge case (v1.0.0): when NO microcluster is mixed, or no cell is mixed, TACIT's
#     `for (i in 1:nrow(x))` loops index row 1 of an empty frame and the call errors.  That cannot be
#     fixed here without editing the package; it is reported as tacit_failed.
suppressPackageStartupMessages({
  library(Seurat)      # CreateSeuratObject/NormalizeData/.../FindClusters are called unqualified inside TACIT
  library(segmented)   # threshold_function calls segmented() unqualified via importFrom; attached for safety
  library(dplyr)       # data_frame (deprecated) is imported from dplyr by TACIT's NAMESPACE
  library(arrow)
  library(jsonlite)
  library(TACIT)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: Rscript tacit_run.R <scratch_dir>")
scratch <- args[[1]]
fail <- function(msg, phase) {
  write_json(list(error = as.character(msg), phase = phase), file.path(scratch, "failure.json"), auto_unbox = TRUE)
  message("tacit_run.R FAILED in ", phase, ": ", msg)
  quit(status = 2L, save = "no")
}
secs <- function(t0) as.numeric(Sys.time() - t0, units = "secs")

params <- fromJSON(file.path(scratch, "params.json"))
tab <- as.data.frame(read_parquet(file.path(scratch, "input.parquet")))
sig <- read.csv(file.path(scratch, "signature.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (colnames(sig)[1] != "cell_type") fail("signature first column must be `cell_type`", "input")
types <- as.character(sig$cell_type)
markers <- colnames(sig)[-1]
missing_cols <- setdiff(markers, colnames(tab))
if (length(missing_cols)) fail(paste("signature markers absent from input:", paste(missing_cols, collapse = ",")), "input")
X_all <- as.data.frame(tab[, markers, drop = FALSE])
in_sub <- as.logical(tab$in_subsample)
X_sub <- X_all[in_sub, , drop = FALSE]
n_all <- nrow(X_all); n_sub <- nrow(X_sub)
if (n_sub < 50L) fail(sprintf("only %d complete cells in the subsample", n_sub), "input")

# p cannot exceed the feature count (uwot refuses); on the real panel p = 10 < 27 so this is a no-op.
p_req <- as.integer(params$p)
p <- min(p_req, length(markers))
seed <- as.integer(params$seed)
timings <- list()

# ------------------------------------------------------------------ TACIT's pre-processing, verbatim
prepass <- function(data_anb, p) {
  orig_values <- t(as.matrix(data_anb)); colnames(orig_values) <- 1:nrow(data_anb)
  md <- data.frame(CellID = 1:nrow(data_anb)); row.names(md) <- md$CellID
  scfp <- CreateSeuratObject(counts = orig_values, meta.data = md)
  scfp <- NormalizeData(scfp, normalization.method = "CLR", margin = 2, verbose = FALSE)
  scfp <- FindVariableFeatures(scfp, selection.method = "disp", nfeatures = 500, verbose = FALSE)
  scfp <- ScaleData(scfp, features = rownames(scfp), verbose = FALSE)
  scfp <- RunPCA(scfp, features = VariableFeatures(object = scfp), verbose = FALSE)
  scfp <- RunUMAP(scfp, features = VariableFeatures(object = scfp), n.components = p, metric = "correlation",
                  umap.method = "uwot", n.neighbors = 15, learning.rate = 1, spread = 1, min.dist = 0.01,
                  set.op.mix.ratio = 1, local.connectivity = 1, verbose = FALSE)
  FindNeighbors(scfp, dims = 1:p, reduction = "umap", verbose = FALSE)
}

ladder_table <- NULL
r_source <- "given"
r_band_status <- NA_character_
if (is.null(params$r) || length(params$r) == 0L || is.na(params$r)) {
  r_source <- "calibrated_on_subsample_prepass"
  t0 <- Sys.time()
  sc <- tryCatch(suppressWarnings(suppressMessages(prepass(X_sub, p))), error = function(e) e)
  if (inherits(sc, "error")) fail(conditionMessage(sc), "prepass")
  timings$prepass_s <- secs(t0)
  band_lo <- as.numeric(params$band_lo); band_hi <- as.numeric(params$band_hi)
  ladder <- as.numeric(params$r_ladder)
  rows <- list(); r_use <- NA_real_
  for (rr in ladder) {
    t1 <- Sys.time()
    s2 <- FindClusters(sc, resolution = rr, random.seed = 1, verbose = FALSE)
    k <- length(unique(s2$seurat_clusters))
    frac <- (n_sub / k) / n_sub                       # mean cluster size as a fraction of the subsample
    rows[[length(rows) + 1L]] <- data.frame(r = rr, n_clusters = k, mean_size = n_sub / k,
                                            mean_size_frac = frac, seconds = secs(t1))
    if (frac <= band_hi) {                            # first ladder value that enters (or passes) the band
      r_use <- rr
      r_band_status <- if (frac >= band_lo) "in_band" else "out_of_band_below_lo"
      break
    }
  }
  ladder_table <- do.call(rbind, rows)
  if (is.na(r_use)) {                                 # the ladder never got fine enough: take its last value
    r_use <- ladder[length(ladder)]
    r_band_status <- "out_of_band_above_hi_ladder_exhausted"
  }
  prepass_clusters <- as.numeric(FindClusters(sc, resolution = r_use, random.seed = 1, verbose = FALSE)$seurat_clusters)
  rm(sc)
} else {
  r_use <- as.numeric(params$r)
  prepass_clusters <- NULL
}

# ------------------------------------------------------------------ the package call, on the subsample
set.seed(seed)                                        # only TACIT's runif jitter in its kNN step is unseeded
log_con <- file(file.path(scratch, "tacit_stdout.log"), open = "wt")
sink(log_con); sink(log_con, type = "message")
t0 <- Sys.time()
res <- tryCatch(TACIT(data_expression = X_sub, r = r_use, p = p, Signature = sig), error = function(e) e)
sink(type = "message"); sink(); close(log_con)
timings$tacit_s <- secs(t0)
if (inherits(res, "error")) fail(conditionMessage(res), "TACIT")

final_threshold <- res[[1]]                           # data.frame(value, Name)
gtd <- res[[2]]                                       # ID, clusters, <type> binaries (subsample rows)
mem <- res[[3]]                                       # mem, ID
clusters <- res[[4]]
thr <- setNames(as.numeric(final_threshold$value), as.character(final_threshold$Name))
if (!all(types %in% names(thr))) fail("TACIT returned no threshold for some type", "postprocess")
if (any(!is.finite(thr[types]))) fail(paste("non-finite threshold for", paste(types[!is.finite(thr[types])], collapse = ",")), "postprocess")

# ------------------------------------------------------------------ TACIT's `ct`, for EVERY cell
Signature_matrix <- as.matrix(sig[, -1, drop = FALSE])
ct_all <- as.matrix(X_all) %*% t(Signature_matrix)   # the identical expression to TACIT's `ct`
colnames(ct_all) <- types
pos_all <- sapply(types, function(t) as.integer(ct_all[, t] > thr[[t]]))   # threshold_function's strict `>`
if (is.null(dim(pos_all))) pos_all <- matrix(pos_all, ncol = length(types), dimnames = list(NULL, types))

# consistency: our binaries on the subsample rows == TACIT's Group_threshold_data
gtd_bin <- as.matrix(gtd[order(gtd$ID), types, drop = FALSE])
n_mismatch <- sum(pos_all[in_sub, , drop = FALSE] != gtd_bin)

tacit_label <- rep(NA_character_, n_all)
sub_rows <- which(in_sub)
mem_ok <- !is.na(mem$ID) & mem$ID >= 1 & mem$ID <= n_sub
tacit_label[sub_rows[mem$ID[mem_ok]]] <- as.character(mem$mem[mem_ok])
cluster_all <- rep(NA_integer_, n_all)
cluster_all[sub_rows] <- as.integer(clusters)

out <- data.frame(object_id = as.character(tab$object_id), in_subsample = in_sub, cluster = cluster_all,
                  tacit_label = tacit_label, stringsAsFactors = FALSE)
for (t in types) { out[[paste0("score__", t)]] <- as.numeric(ct_all[, t]); out[[paste0("pos__", t)]] <- as.integer(pos_all[, t]) }
write_parquet(out, file.path(scratch, "scores.parquet"))

desc <- packageDescription("TACIT")
learned <- list(
  r_used = r_use, r_source = r_source, r_band_status = r_band_status,
  r_ladder = if (is.null(ladder_table)) NULL else ladder_table,
  p_requested = p_req, p_used = p, seed = seed,
  n_cells_scored = n_all, n_subsample = n_sub,
  n_clusters = length(unique(clusters)), mean_cluster_size = n_sub / length(unique(clusters)),
  mean_cluster_size_frac = 1 / length(unique(clusters)),
  clusters_identical_to_prepass = if (is.null(prepass_clusters)) NA else identical(as.numeric(clusters), prepass_clusters),
  thresholds = as.list(thr[types]),
  n_binary_mismatch_vs_tacit = n_mismatch,
  n_tacit_labelled = sum(!is.na(tacit_label)), n_tacit_others = sum(tacit_label == "Others", na.rm = TRUE),
  tacit_internal_seeds = list(RunUMAP_seed.use = 42, FindClusters_random.seed = 1, algorithm = "Louvain (Seurat FindClusters algorithm 1)"),
  timings_s = timings,
  versions = list(R = R.version.string, TACIT = desc$Version, TACIT_remote_ref = desc$RemoteRef, TACIT_sha = desc$RemoteSha,
                  Seurat = as.character(packageVersion("Seurat")), SeuratObject = as.character(packageVersion("SeuratObject")),
                  uwot = as.character(packageVersion("uwot")), segmented = as.character(packageVersion("segmented")),
                  caret = as.character(packageVersion("caret")), dplyr = as.character(packageVersion("dplyr")),
                  arrow = as.character(packageVersion("arrow"))),
  cores_visible = parallel::detectCores()
)
write_json(learned, file.path(scratch, "learned.json"), auto_unbox = TRUE, digits = NA, null = "null", na = "null", pretty = TRUE)
message(sprintf("tacit_run.R ok: n=%d subsample=%d r=%g (%s) p=%d clusters=%d mismatch=%d", n_all, n_sub, r_use, r_source, p, length(unique(clusters)), n_mismatch))
