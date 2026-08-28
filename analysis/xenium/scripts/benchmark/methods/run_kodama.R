#!/usr/bin/env Rscript
# KODAMA runner (tkcaccia/KODAMA v3, spatially-informed) -- spatial-domain track.
# Runs in the xenium_analysis R 4.5.3 with packages from R_libs.
#
# KODAMA.matrix(data, spatial=coords, ...) builds a spatially-constrained,
# cross-validated dissimilarity; KODAMA.visualization(kk, "UMAP") embeds it;
# we KMeans the embedding to --target-domains (granularity-matched).
#
# Reads the shared bridge written by prep_benchmark_input.py:
#   {s}_bench_{v}_expr.mtx  (ROWS=cells x COLS=genes, lognorm)
#   {s}_bench_{v}_coords.csv, _barcodes.csv
# Writes results/kodama_{v}_labels.csv (cell_id, domain).
#
# Usage:
#   Rscript scripts/benchmark/methods/run_kodama.R --sample 0041435 --variant crop25k [--n-cores 8]

suppressWarnings(suppressMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER", "/blue/maigan/smith6jt/R_libs"), .libPaths()))
  library(Matrix); library(KODAMA)
}))

# Set by run_benchmark.sh, which derives it from its own location.  No default:
# a wrong ROOT here writes results into someone else's tree, silently.
ROOT <- Sys.getenv("XENIUM_ANALYSIS_ROOT")
if (!nzchar(ROOT)) stop("XENIUM_ANALYSIS_ROOT is unset; run via run_benchmark.sh")
SEED <- 0L

# ---- tiny arg parser --------------------------------------------------------
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  # Tractable, DOCUMENTED defaults: KODAMA's full M=100/landmarks=10000 is
  # ~hours on 25-50k cells (landmarks=10000 half-samples the data). landmarks
  # is KODAMA's intended scalability knob; M is the Monte-Carlo robustness knob.
  # The benchmark reports these settings + flags KODAMA as the slowest method.
  out <- list(sample = NA, variant = "crop50k", n_cores = 1L,
              target_domains = 12L, M = 50L, ncomp = 50L, landmarks = 2000L)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]; v <- if (i < length(a)) a[[i + 1]] else NA
    switch(k,
      "--sample" = {out$sample <- v; i <- i + 2},
      "--variant" = {out$variant <- v; i <- i + 2},
      "--n-cores" = {out$n_cores <- as.integer(v); i <- i + 2},
      "--target-domains" = {out$target_domains <- as.integer(v); i <- i + 2},
      "--M" = {out$M <- as.integer(v); i <- i + 2},
      "--ncomp" = {out$ncomp <- as.integer(v); i <- i + 2},
      "--landmarks" = {out$landmarks <- as.integer(v); i <- i + 2},
      {i <- i + 1})
  }
  if (is.na(out$sample)) stop("--sample required")
  out
}
args <- parse_args()

bench <- file.path(ROOT, "data/processed/benchmark", args$sample)
stem  <- file.path(bench, sprintf("%s_bench_%s", args$sample, args$variant))
resdir <- file.path(bench, "results"); dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

expr   <- readMM(paste0(stem, "_expr.mtx"))          # cells x genes (sparse)
coords <- as.matrix(read.csv(paste0(stem, "_coords.csv")))   # cols x,y
ids    <- read.csv(paste0(stem, "_barcodes.csv"))[[1]]
data   <- as.matrix(expr)
storage.mode(data) <- "double"
cat(sprintf("  [kodama] %d cells x %d genes; n.cores=%d, M=%d, target=%d\n",
            nrow(data), ncol(data), args$n_cores, args$M, args$target_domains))
stopifnot(nrow(data) == nrow(coords), nrow(data) == length(ids))

t0 <- Sys.time()
kk <- KODAMA.matrix(data, spatial = coords, M = args$M,
                    ncomp = min(args$ncomp, ncol(data)),
                    landmarks = args$landmarks, spatial.resolution = 0.3,
                    n.cores = args$n_cores, seed = SEED)
vis <- KODAMA.visualization(kk, method = "UMAP")
emb <- as.matrix(vis)
cat(sprintf("  [kodama] embedding %d x %d in %.1f s\n",
            nrow(emb), ncol(emb), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

set.seed(SEED)
km  <- kmeans(emb, centers = args$target_domains, nstart = 10, iter.max = 300)
lab <- km$cluster

out <- data.frame(cell_id = as.character(ids), domain = as.character(lab))
outfile <- file.path(resdir, sprintf("kodama_%s_labels.csv", args$variant))
write.csv(out, outfile, row.names = FALSE)
cat(sprintf("  [kodama] wrote %d domains x %d cells -> %s\n",
            length(unique(lab)), nrow(out), basename(outfile)))
