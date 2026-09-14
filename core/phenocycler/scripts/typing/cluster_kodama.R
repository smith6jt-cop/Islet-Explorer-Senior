#!/usr/bin/env Rscript
# The R half of scripts/typing/cluster_kodama.py: KODAMA's learned dissimilarity graph.
#
# It returns a GRAPH, not an embedding.  `KODAMA.matrix` already ends by rewriting its own
# k-nearest-neighbour arrays in place -- `(1 + euclidean_ij) / p_ij^2`, where p_ij is the
# fraction of the M Monte-Carlo iterations that gave cells i and j the same label -- and then
# re-sorting each row by that distance.  `KODAMA.visualization` only lays those same arrays out
# with UMAP.  Handing the first `graph_k` columns straight to the caller therefore skips the
# projection without skipping any of KODAMA: it is the representation the projection is made of.
#
# Reads  <scratch>/Z.parquet  (cells x markers, already standardised per marker by the caller)
# Writes <scratch>/knn.parquet (i01..iK int32 1-based, d01..dK double) and <scratch>/kodama.json
#
#   Rscript cluster_kodama.R --scratch <dir> --M 100 --landmarks 500 --graph-k 15 \
#                            --n-cores 5 --seed 0

suppressWarnings(suppressMessages({library(KODAMA); library(arrow); library(jsonlite)}))

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  o <- list(scratch = NA, M = 100L, landmarks = 500L, ncomp = 50L, graph_k = 15L,
            n_cores = 5L, seed = 0L)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]; v <- if (i < length(a)) a[[i + 1]] else NA
    switch(k,
      "--scratch"   = {o$scratch   <- v;              i <- i + 2},
      "--M"         = {o$M         <- as.integer(v);  i <- i + 2},
      "--landmarks" = {o$landmarks <- as.integer(v);  i <- i + 2},
      "--ncomp"     = {o$ncomp     <- as.integer(v);  i <- i + 2},
      "--graph-k"   = {o$graph_k   <- as.integer(v);  i <- i + 2},
      "--n-cores"   = {o$n_cores   <- as.integer(v);  i <- i + 2},
      "--seed"      = {o$seed      <- as.integer(v);  i <- i + 2},
      stop(sprintf("unknown argument %s", k)))
  }
  if (is.na(o$scratch)) stop("--scratch is required")
  o
}
args <- parse_args()

# KODAMA's M loop dispatches `mclapply(chunk, mc.cores = n.cores)` with
# chunk_size = max(1, min(16, ceiling(M / (4 * n.cores)))), so it cannot occupy more than
# min(16, ceiling(M / (4 * cores))) workers -- it needs M >= 4 * cores^2 before `cores` cores are
# all busy.  Refuse rather than burn the difference; the same gate is in
# scripts/xenium_celltype/hpg/kodama_shard.R and the reason is measured there.
occupancy <- function(m, cores) {
  if (cores <= 1L) return(1L)
  min(max(1L, min(16L, ceiling(m / (cores * 4L)))), cores)
}
occ <- occupancy(args$M, args$n_cores)
if (occ < args$n_cores) {
  stop(sprintf(paste0("KODAMA's M loop would occupy %d of the %d cores asked for; the chunk rule ",
                      "needs M >= 4 * cores^2, i.e. M >= %d here. Raise M or drop n_cores to %d."),
               occ, args$n_cores, 4L * args$n_cores^2, max(1L, as.integer(floor(sqrt(args$M / 4))))))
}

# Pin, checked before the hours rather than after them.  An unmodified fork of
# tkcaccia/KODAMAlegacy; cluster_kodama.py records the same two strings.
KODAMA_VERSION <- "3.3"
KODAMA_SHA <- "d2dbe30ee66509b82c084616d6961c5b292cd059"
desc <- read.dcf(system.file("DESCRIPTION", package = "KODAMA"))
got_ver <- as.character(packageVersion("KODAMA"))
got_sha <- if ("RemoteSha" %in% colnames(desc)) unname(desc[1, "RemoteSha"]) else NA_character_
if (!identical(got_ver, KODAMA_VERSION) || !identical(got_sha, KODAMA_SHA)) {
  stop(sprintf("installed KODAMA is %s / %s; this script is pinned to %s / %s",
               got_ver, got_sha, KODAMA_VERSION, KODAMA_SHA))
}

Z <- as.matrix(arrow::read_parquet(file.path(args$scratch, "Z.parquet")))
storage.mode(Z) <- "double"
if (anyNA(Z)) stop("Z.parquet contains NA; a missing measurement is unavailable, not a value")
if (args$graph_k < 2L) stop("--graph-k must be at least 2")

# KODAMA's own neighbour count, reproduced here so a graph too narrow for --graph-k fails in
# seconds instead of after the fit.  It is derived from the ORIGINAL --landmarks, before the
# `if (nsample <= landmarks)` adjustment inside the package, and is capped at 500.
neighbors <- max(1L, floor(min(c(args$landmarks, nrow(Z) * 0.75 - 1), 500)))
if (neighbors < args$graph_k) {
  stop(sprintf(paste0("KODAMA would build a %d-neighbour graph but --graph-k is %d; `neighbors` ",
                      "is min(landmarks, 0.75n - 1, 500), so raise --landmarks (now %d) or lower ",
                      "--graph-k. n = %d."),
               neighbors, args$graph_k, args$landmarks, nrow(Z)))
}

cat(sprintf("[kodama] %d cells x %d markers | M=%d landmarks=%d ncomp=%d cores=%d (occupancy %d) seed=%d\n",
            nrow(Z), ncol(Z), args$M, args$landmarks, min(args$ncomp, ncol(Z)),
            args$n_cores, occ, args$seed))

t0 <- Sys.time()
kk <- KODAMA.matrix(Z, M = args$M, ncomp = min(args$ncomp, ncol(Z)),
                    landmarks = args$landmarks, n.cores = args$n_cores, seed = args$seed)
t1 <- Sys.time()

idx <- kk$knn_Rnanoflann$indices
dst <- kk$knn_Rnanoflann$distances
if (ncol(idx) != neighbors || ncol(idx) < args$graph_k) {
  stop(sprintf("KODAMA returned a %d-neighbour graph; %d was predicted and %d asked for",
               ncol(idx), neighbors, args$graph_k))
}
K <- args$graph_k
idx <- idx[, seq_len(K), drop = FALSE]
dst <- dst[, seq_len(K), drop = FALSE]

# A non-finite kept distance means p_ij = 0: that cell was never once co-assigned with the
# neighbour, so M is too small to place it.  Silently dropping the edge would change the graph
# for exactly the cells the graph is least sure about, so this stops instead.
bad <- !is.finite(dst)
if (any(bad)) {
  stop(sprintf(paste0("%d of %d kept edges have p_ij = 0 (distance Inf) across %d cells: M = %d ",
                      "is too small for this data. Raise --M."),
               sum(bad), length(dst), sum(rowSums(bad) > 0), args$M))
}

storage.mode(idx) <- "integer"
out <- c(setNames(as.data.frame(idx), sprintf("i%02d", seq_len(K))),
         setNames(as.data.frame(dst), sprintf("d%02d", seq_len(K))))
arrow::write_parquet(as.data.frame(out), file.path(args$scratch, "knn.parquet"),
                     compression = "zstd")

write_json(list(complete = TRUE, n_cells = nrow(Z), n_markers = ncol(Z),
                M = args$M, landmarks = args$landmarks, ncomp = min(args$ncomp, ncol(Z)),
                graph_k = K, neighbors = kk$knn_Rnanoflann$neighbors,
                n_cores = args$n_cores, m_loop_occupancy = occ, seed = args$seed,
                mean_cv_accuracy = mean(kk$acc, na.rm = TRUE),
                kodama_version = got_ver, kodama_sha = got_sha,
                r_version = R.version.string,
                blas_threads = Sys.getenv("OPENBLAS_NUM_THREADS", ""),
                matrix_seconds = as.numeric(difftime(t1, t0, units = "secs")),
                total_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
           file.path(args$scratch, "kodama.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("[kodama] done in %.1f s -> knn.parquet (%d x %d)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(idx), K))
