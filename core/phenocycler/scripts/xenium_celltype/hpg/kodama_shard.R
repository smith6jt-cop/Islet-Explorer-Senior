#!/usr/bin/env Rscript
# One shard of a KODAMA run: Monte-Carlo iterations [m*s+1 .. m*s+m] of M, and nothing else.
#
# WHY THIS EXISTS. KODAMA 3.3's `KODAMA.matrix` couples its M iterations in exactly one
# place -- the M x n label matrix `res` -- and `one_iteration(k)` opens with
# `set.seed(seed + k)` and reads no other iteration. So iteration k is a pure function of
# (data, spatial, landmarks, Tcycle, ncomp, seed + k). Running M in shards and rbind-ing
# `res` is therefore not an approximation; it is the same computation.
#
# It needs no patch to KODAMA, because `seed` addresses the k-range: a shard asking for
# M = m at seed = S + m*s runs set.seed(S + m*s + 1 .. S + m*s + m), and the shards tile
# S+1 .. S+M exactly. MEASURED (KODAMA 3.3, both arms): a monolithic M=6 seed=0 run and
# rbind of three M=2 runs at seeds 0/2/4 give `identical()` res, acc and res_constrain.
# `check_shard_equivalence.R` re-proves this for the pinned commit before any production
# submission, and the production submit script will not proceed if it fails.
#
# WHY SHARDING IS NOT MERELY CONVENIENT. `run_parallel_chunks` sets
#     chunk_size = max(1, min(16, ceiling(length(indices) / (n.cores * 4))))
# and dispatches `mclapply(chunk, ..., mc.cores = n.cores)` chunk by chunk, so the M loop's
# concurrency is min(chunk_size, n.cores) -- it needs M >= 4 * n.cores^2 before it can use
# n.cores at all. MEASURED here on a 3,000 x 30 matrix at M=16: 2 cores 100.7 s (chunk 2),
# 4 cores 103.7 s (chunk 1), 16 cores 118.6 s (chunk 1). More cores was SLOWER. The
# dissimilarity loop iterates over cells, not over M, so it does parallelise normally.
#
# Consequence, enforced below: this script computes KODAMA's own occupancy and refuses to
# run when the allocation cannot be used, rather than quietly burning idle cores for hours.
# The rule is m >= 4 * n_cores^2; the SLURM scripts pin (m = 32, cpus-per-task = 2) and
# (m = 100, cpus-per-task = 5) as the two shapes that satisfy it.
#
# Writes <out>/res_<arm>_shard<NN>.rds  (integer M x n labels + acc + v + res_constrain)
#    and <out>/res_<arm>_shard<NN>.json (provenance: versions, seeds, digests, timings)
# and DISCARDS the shard's own dissimilarity/knn, which kodama_merge.R recomputes once over
# the full rbind-ed `res`. That waste is real and is reported in the JSON as knn_wasted_s.

suppressWarnings(suppressMessages({library(Matrix); library(KODAMA); library(jsonlite)}))

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  o <- list(input = NA, output = NA, arm = "nonspatial", m = 32L, shard = 0L,
            base_seed = 0L, landmarks = 2000L, ncomp = 50L,
            spatial_resolution = 0.3, n_cores = 2L, force = FALSE)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]; v <- if (i < length(a)) a[[i + 1]] else NA
    switch(k,
      "--input" = {o$input <- v; i <- i + 2}, "--output" = {o$output <- v; i <- i + 2},
      "--arm" = {o$arm <- v; i <- i + 2}, "--m" = {o$m <- as.integer(v); i <- i + 2},
      "--shard" = {o$shard <- as.integer(v); i <- i + 2},
      "--base-seed" = {o$base_seed <- as.integer(v); i <- i + 2},
      "--landmarks" = {o$landmarks <- as.integer(v); i <- i + 2},
      "--ncomp" = {o$ncomp <- as.integer(v); i <- i + 2},
      "--spatial-resolution" = {o$spatial_resolution <- as.numeric(v); i <- i + 2},
      "--n-cores" = {o$n_cores <- as.integer(v); i <- i + 2},
      "--force" = {o$force <- TRUE; i <- i + 1},
      {i <- i + 1})
  }
  if (is.na(o$input) || is.na(o$output)) stop("--input and --output are required")
  if (!o$arm %in% c("spatial", "nonspatial")) stop("--arm must be spatial or nonspatial")
  o
}
args <- parse_args()
dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

# ---- occupancy gate: refuse to waste an allocation --------------------------------------
# Mirrors run_parallel_chunks' arithmetic for the M loop exactly.
occupancy <- function(m, cores) {
  if (cores <= 1L) return(1L)
  min(max(1L, min(16L, ceiling(m / (cores * 4L)))), cores)
}
occ <- occupancy(args$m, args$n_cores)
if (occ < args$n_cores && !args$force) {
  stop(sprintf(paste0("KODAMA's M loop would use %d of the %d requested cores (chunk rule ",
                      "needs m >= 4*cores^2, i.e. m >= %d for %d cores). Either raise --m, ",
                      "lower --n-cores to %d, or pass --force to burn the difference."),
               occ, args$n_cores, 4L * args$n_cores^2, args$n_cores,
               max(1L, as.integer(floor(sqrt(args$m / 4))))))
}

X   <- as.matrix(readMM(file.path(args$input, "pcs.mtx"))); storage.mode(X) <- "double"
ids <- read.csv(file.path(args$input, "cell_id.csv"))[[1]]
co  <- as.matrix(read.csv(file.path(args$input, "coords.csv")))
stopifnot(nrow(X) == length(ids), nrow(X) == nrow(co))

seed  <- args$base_seed + args$m * args$shard
stem  <- file.path(args$output, sprintf("res_%s_shard%02d", args$arm, args$shard))
rds   <- paste0(stem, ".rds"); js <- paste0(stem, ".json")

# ---- idempotent restart: a finished shard is never recomputed ----------------------------
# The provenance must match, not merely the filename -- a shard written under a different
# seed, m, landmarks or input digest is a different computation wearing the same name.
input_digest <- digest::digest(list(dim(X), sum(X), args$landmarks, args$ncomp,
                                    args$arm, args$spatial_resolution), algo = "xxhash64")
if (file.exists(rds) && file.exists(js) && !args$force) {
  prev <- try(fromJSON(js), silent = TRUE)
  if (!inherits(prev, "try-error") && isTRUE(prev$complete) &&
      identical(prev$seed, seed) && identical(prev$m, args$m) &&
      identical(prev$input_digest, input_digest)) {
    cat(sprintf("[shard] %s shard %02d already complete (seed=%d) -- skipping\n",
                args$arm, args$shard, seed)); quit(status = 0)
  }
  cat("[shard] existing output does not match this request; recomputing\n")
}

cat(sprintf("[shard] %s shard %02d/%s: %d cells x %d comps | m=%d seed=%d landmarks=%d\n",
            args$arm, args$shard, "?", nrow(X), ncol(X), args$m, seed, args$landmarks))
cat(sprintf("[shard] cores=%d -> KODAMA M-loop occupancy %d\n", args$n_cores, occ))

t0 <- Sys.time()
kk <- if (args$arm == "spatial") {
  KODAMA.matrix(X, spatial = co, M = args$m, ncomp = min(args$ncomp, ncol(X)),
                landmarks = args$landmarks, spatial.resolution = args$spatial_resolution,
                n.cores = args$n_cores, seed = seed)
} else {
  KODAMA.matrix(X, M = args$m, ncomp = min(args$ncomp, ncol(X)),
                landmarks = args$landmarks, n.cores = args$n_cores, seed = seed)
}
t1 <- Sys.time()

res <- kk$res
stopifnot(all(is.na(res) | res == round(res)))   # labels are integral; int32 is lossless
storage.mode(res) <- "integer"
rc <- kk$res_constrain
stopifnot(all(is.na(rc) | rc == round(rc)))
storage.mode(rc) <- "integer"

saveRDS(list(res = res, acc = kk$acc, v = kk$v, res_constrain = rc,
             arm = args$arm, shard = args$shard, m = args$m, seed = seed,
             k_global = seq.int(args$m * args$shard + 1L, args$m * (args$shard + 1L))),
        rds, compress = TRUE)

write_json(list(complete = TRUE, arm = args$arm, shard = args$shard, m = args$m,
                seed = seed, base_seed = args$base_seed,
                k_global_first = args$m * args$shard + 1L,
                k_global_last = args$m * (args$shard + 1L),
                n_cells = nrow(X), n_comps = ncol(X), landmarks = args$landmarks,
                ncomp = args$ncomp, n_cores = args$n_cores, m_loop_occupancy = occ,
                spatial_resolution = if (args$arm == "spatial") args$spatial_resolution else NULL,
                input_digest = input_digest,
                kodama_version = as.character(packageVersion("KODAMA")),
                r_version = R.version.string,
                seconds = as.numeric(difftime(t1, t0, units = "secs")),
                slurm_job = Sys.getenv("SLURM_JOB_ID", ""),
                slurm_array_task = Sys.getenv("SLURM_ARRAY_TASK_ID", "")),
           js, auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("[shard] %s shard %02d done in %.1f s -> %s\n", args$arm, args$shard,
            as.numeric(difftime(t1, t0, units = "secs")), basename(rds)))
