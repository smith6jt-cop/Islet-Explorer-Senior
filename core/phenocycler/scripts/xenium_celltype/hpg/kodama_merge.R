#!/usr/bin/env Rscript
# Merge the M shards of one (donor, arm) into the embedding a monolithic KODAMA.matrix
# would have produced, then visualise it.
#
# WHAT IS REPRODUCED, AND WHY IT IS SAFE TO REPRODUCE IT. KODAMA.matrix ends with two steps
# the shards cannot share: the feature-space kNN (which uses only `data`, so it is identical
# for every shard AND for both arms) and the dissimilarity reweighting, which needs the FULL
# M x n label matrix. Those ~20 lines are re-implemented here against the pinned KODAMA
# commit -- the same arithmetic, vectorised over each row's 500 neighbours and blocked over
# rows so a 1.1M-cell donor neither exhausts memory nor loses everything to a wall clock.
#
# Re-implementation is a real risk and it is met with a test, not a comment:
# `check_shard_equivalence.R` runs a monolithic KODAMA.matrix and this shard+merge path on
# the same input and asserts `identical()` on `res`, on the reweighted neighbour indices and
# distances, and on the embedding. The submit script runs it and stops if it fails, so the
# production run cannot start against a KODAMA version this code has not been proven against.
#
# TWO FIXES TO THE UPSTREAM VISUALISATION CALL, both about reproducibility:
#   * `KODAMA.visualization` sets config$n_sgd_threads <- n.cores whenever
#     config$define.n.cores is not TRUE. uwot documents multi-threaded SGD as
#     non-reproducible. A run whose embedding changes between invocations would make the
#     KMeans groups cut from it change too, so the config is pinned here with
#     define.n.cores = TRUE, n_sgd_threads = 1 and an explicit random_state; n_threads keeps
#     the parallel parts parallel.
#   * only the first `numap = min(3 * n_neighbors, n - 1, ncol(indices))` = 45 neighbour
#     columns are read by the UMAP, so the checkpointed reweighted graph is trimmed to 64
#     columns. That is provably the same input, at an eighth of the memory.
#
# Writes <out>/embedding_<arm>.csv, <out>/kodama_<arm>.json, and
#        <out>/dissim_<arm>/block_<NNNN>.rds  (row-block checkpoints; safe to delete after)

suppressWarnings(suppressMessages({
  library(Matrix); library(KODAMA); library(jsonlite); library(Rnanoflann)
}))

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  o <- list(input = NA, output = NA, arm = "nonspatial", M = 96L, m = 32L,
            landmarks = 2000L, n_cores = 32L, seed = 0L, block = 50000L,
            keep_cols = 64L, keep_blocks = FALSE)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]; v <- if (i < length(a)) a[[i + 1]] else NA
    switch(k,
      "--input" = {o$input <- v; i <- i + 2}, "--output" = {o$output <- v; i <- i + 2},
      "--arm" = {o$arm <- v; i <- i + 2}, "--M" = {o$M <- as.integer(v); i <- i + 2},
      "--m" = {o$m <- as.integer(v); i <- i + 2},
      "--landmarks" = {o$landmarks <- as.integer(v); i <- i + 2},
      "--n-cores" = {o$n_cores <- as.integer(v); i <- i + 2},
      "--seed" = {o$seed <- as.integer(v); i <- i + 2},
      "--block" = {o$block <- as.integer(v); i <- i + 2},
      "--keep-cols" = {o$keep_cols <- as.integer(v); i <- i + 2},
      "--keep-blocks" = {o$keep_blocks <- TRUE; i <- i + 1},
      {i <- i + 1})
  }
  if (is.na(o$input) || is.na(o$output)) stop("--input and --output are required")
  o
}
args <- parse_args()

# =========================================================================================
# 1. Collect the shards and prove they tile 1..M exactly.
# =========================================================================================
collect_res <- function(dir, arm, M, m) {
  files <- sort(list.files(dir, pattern = sprintf("^res_%s_shard[0-9]+\\.rds$", arm),
                           full.names = TRUE))
  if (!length(files)) stop(sprintf("no shards for arm=%s under %s", arm, dir))
  parts <- lapply(files, readRDS)
  ks <- sort(unlist(lapply(parts, `[[`, "k_global")))
  if (!identical(ks, seq_len(M))) {
    have <- paste(range(ks), collapse = "..")
    stop(sprintf(paste0("shards for arm=%s cover %d iterations (%s), not 1..%d. A merge over ",
                        "a partial or overlapping iteration set is a DIFFERENT computation ",
                        "from the one recorded in the manifest and is refused."),
                 arm, length(ks), have, M))
  }
  if (length(unique(vapply(parts, `[[`, integer(1), "m"))) != 1L) stop("shards disagree on m")
  ord <- order(vapply(parts, function(p) p$k_global[1], integer(1)))
  parts <- parts[ord]
  res <- do.call(rbind, lapply(parts, `[[`, "res"))
  storage.mode(res) <- "double"                       # exactly the upstream storage mode
  list(res = res,
       acc = unlist(lapply(parts, `[[`, "acc")),
       v = do.call(rbind, lapply(parts, `[[`, "v")),
       res_constrain = { rc <- do.call(rbind, lapply(parts, `[[`, "res_constrain"))
                         storage.mode(rc) <- "double"; rc },
       n_shards = length(parts), files = basename(files))
}

X   <- as.matrix(readMM(file.path(args$input, "pcs.mtx"))); storage.mode(X) <- "double"
ids <- read.csv(file.path(args$input, "cell_id.csv"))[[1]]
nsample <- nrow(X)
stopifnot(nsample == length(ids))

S <- collect_res(args$output, args$arm, args$M, args$m)
res <- S$res
stopifnot(ncol(res) == nsample, nrow(res) == args$M)
cat(sprintf("[merge] %s: %d shards -> res %d x %d\n", args$arm, S$n_shards,
            nrow(res), ncol(res)))

# =========================================================================================
# 2. Feature-space kNN. Identical to KODAMA.matrix's, and identical for BOTH arms, because
#    it reads `data` only. Cached on disk so the second arm does not pay for it again.
# =========================================================================================
neighbors <- max(1L, floor(min(c(args$landmarks, nsample * 0.75 - 1), 500)))
nn_call <- function(x, y, k, method, cores) {
  if (cores > 1L) {
    out <- try(Rnanoflann::nn(x, y, k, method = method, parallel = TRUE, cores = cores),
               silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  Rnanoflann::nn(x, y, k, method = method)
}
knn_cache <- file.path(args$input, sprintf("knn_%d.rds", neighbors))
if (file.exists(knn_cache)) {
  cat("[merge] reusing cached feature kNN\n"); knn <- readRDS(knn_cache)
} else {
  t0 <- Sys.time()
  knn <- nn_call(X, X, neighbors + 1L, "euclidean", args$n_cores)
  dn <- if (!is.null(knn$distances)) "distances" else "distance"
  if (is.null(knn[[dn]])) stop("Rnanoflann::nn did not return distances")
  knn$distances <- knn[[dn]][, -1, drop = FALSE]
  knn$indices   <- knn$indices[, -1, drop = FALSE]
  knn <- list(indices = knn$indices, distances = knn$distances)
  saveRDS(knn, knn_cache, compress = FALSE)
  cat(sprintf("[merge] feature kNN (k=%d) in %.1f s\n", neighbors,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopifnot(nrow(knn$indices) == nsample, ncol(knn$indices) == neighbors)

# =========================================================================================
# 3. Dissimilarity reweighting, row-blocked and checkpointed.
#    Upstream, per row k and each of its `neighbors` neighbours j:
#        kod = mean(res[, k] == res[, idx[j]], na.rm = TRUE)
#        d[j] = (1 + d[j]) / kod^2
#    then reorder by d. colMeans over the neighbour block is the same arithmetic with one
#    R-level call instead of 500.
# =========================================================================================
keep <- min(args$keep_cols, neighbors)
numap_needed <- min(3L * 15L, nsample - 1L)          # config.umap.default$n_neighbors = 15
if (keep < numap_needed)
  stop(sprintf("--keep-cols %d is below the %d columns the UMAP reads", keep, numap_needed))

blkdir <- file.path(args$output, sprintf("dissim_%s", args$arm))
dir.create(blkdir, showWarnings = FALSE, recursive = TRUE)
starts <- seq.int(1L, nsample, by = args$block)

# `Rnanoflann::nn` returns DOUBLE indices and KODAMA keeps them double, so the checkpoint
# keeps them double too. Storing the same numbers as int32 would halve the file and fail
# `identical()` against upstream for a reason that has nothing to do with the arithmetic --
# and the equivalence gate is only worth having if it is comparing like with like.
reweight_block <- function(rows) {
  idx <- matrix(0, length(rows), keep); dst <- matrix(0, length(rows), keep)
  out <- parallel::mclapply(seq_along(rows), function(i) {
    k <- rows[i]
    ii <- knn$indices[k, ]; dd <- knn$distances[k, ]
    kod <- colMeans(res[, ii, drop = FALSE] == res[, k], na.rm = TRUE)
    dd <- (1 + dd) / (kod^2)
    oo <- order(dd)
    list(i = ii[oo][seq_len(keep)], d = dd[oo][seq_len(keep)])
  }, mc.cores = args$n_cores, mc.preschedule = TRUE)
  for (i in seq_along(rows)) { idx[i, ] <- out[[i]]$i; dst[i, ] <- out[[i]]$d }
  list(indices = idx, distances = dst)
}

t0 <- Sys.time()
for (b in seq_along(starts)) {
  f <- file.path(blkdir, sprintf("block_%04d.rds", b))
  if (file.exists(f)) next                            # restart: finished blocks are not redone
  rows <- starts[b]:min(starts[b] + args$block - 1L, nsample)
  saveRDS(reweight_block(rows), f, compress = FALSE)
  cat(sprintf("[merge] %s dissimilarity block %d/%d (%d rows) %.0f s elapsed\n",
              args$arm, b, length(starts), length(rows),
              as.numeric(difftime(Sys.time(), t0, units = "secs")))); flush(stdout())
}
IDX <- matrix(0, nsample, keep); DST <- matrix(0, nsample, keep)
for (b in seq_along(starts)) {
  rows <- starts[b]:min(starts[b] + args$block - 1L, nsample)
  z <- readRDS(file.path(blkdir, sprintf("block_%04d.rds", b)))
  IDX[rows, ] <- z$indices; DST[rows, ] <- z$distances
}
dissim_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# =========================================================================================
# 4. Visualisation, with the config pinned so the embedding is reproducible.
# =========================================================================================
kk <- list(acc = S$acc, v = S$v, res = res,
           knn_Rnanoflann = list(indices = IDX, distances = DST, neighbors = neighbors),
           data = X, res_constrain = S$res_constrain, n.cores = args$n_cores)
cfg <- KODAMA::config.umap.default
cfg$define.n.cores <- TRUE
cfg$n_threads      <- args$n_cores
cfg$n_sgd_threads  <- 1L
cfg$random_state   <- args$seed
t1 <- Sys.time()
emb <- as.matrix(KODAMA.visualization(kk, method = "UMAP", config = cfg))
vis_s <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

# TWO copies of the embedding, on purpose. The CSV is what downstream reads; write.csv
# round-trips through text and loses ~5e-14, which is irrelevant to a KMeans over it and
# fatal to an exactness proof. The RDS carries the doubles unchanged, so
# check_shard_equivalence.R can assert bit-equality against a monolithic run rather than
# assert it within a tolerance chosen to make the test pass.
write.csv(data.frame(cell_id = as.character(ids), dim1 = emb[, 1], dim2 = emb[, 2]),
          file.path(args$output, sprintf("embedding_%s.csv", args$arm)), row.names = FALSE)
saveRDS(list(cell_id = as.character(ids), embedding = emb),
        file.path(args$output, sprintf("embedding_%s.rds", args$arm)), compress = TRUE)
write_json(list(complete = TRUE, arm = args$arm, n_cells = nsample, n_comps = ncol(X),
                M = args$M, m_shard = args$m, n_shards = S$n_shards, shard_files = S$files,
                landmarks = args$landmarks, neighbors = neighbors, keep_cols = keep,
                base_seed = args$seed, n_cores = args$n_cores,
                umap_config = list(n_sgd_threads = 1L, random_state = args$seed,
                                   define.n.cores = TRUE),
                kodama_version = as.character(packageVersion("KODAMA")),
                r_version = R.version.string,
                dissimilarity_seconds = dissim_s, visualization_seconds = vis_s,
                mean_iteration_accuracy = mean(S$acc, na.rm = TRUE)),
           file.path(args$output, sprintf("kodama_%s.json", args$arm)),
           auto_unbox = TRUE, pretty = TRUE)

if (!args$keep_blocks) unlink(blkdir, recursive = TRUE)
cat(sprintf("[merge] %s done: dissimilarity %.1f s, UMAP %.1f s -> embedding_%s.csv\n",
            args$arm, dissim_s, vis_s, args$arm))
