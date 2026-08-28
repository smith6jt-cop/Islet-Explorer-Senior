#!/usr/bin/env Rscript
# MISTy runner (saezlab/mistyR) -- INTERACTION track (NOT a domain method).
# Runs in the xenium_analysis R 4.5.3 with packages from R_libs.
#
# Models each target marker from three spatial views (intra / juxta / para);
# outputs per-target intra.R2/multi.R2/gain.R2, view contributions, and
# predictor->target importances. Run on the MISTy marker subset (the tissue
# MARKER_DICT union) so per-target ranger cost is tractable on the 5k panel.
#
# Reads the shared bridge (expr.mtx cells x genes, genes.csv, coords.csv,
# misty_markers.txt). Writes results/misty_{v}/ (run_misty) + tidy CSVs
# results/misty_{v}_{improvements,importances,contributions,...}.csv.
#
# Usage:
#   Rscript scripts/benchmark/methods/run_misty.R --sample 0041435 --variant crop25k [--n-cores 8]

suppressWarnings(suppressMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER", "/blue/maigan/smith6jt/R_libs"), .libPaths()))
  library(Matrix); library(mistyR); library(future); library(dplyr)
}))

# Set by run_benchmark.sh, which derives it from its own location.  No default:
# a wrong ROOT here writes results into someone else's tree, silently.
ROOT <- Sys.getenv("XENIUM_ANALYSIS_ROOT")
if (!nzchar(ROOT)) stop("XENIUM_ANALYSIS_ROOT is unset; run via run_benchmark.sh")
SEED <- 0L

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(sample = NA, variant = "crop50k", n_cores = 1L,
              neighbor_thr = 30, l = 100)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]; v <- if (i < length(a)) a[[i + 1]] else NA
    switch(k,
      "--sample" = {out$sample <- v; i <- i + 2},
      "--variant" = {out$variant <- v; i <- i + 2},
      "--n-cores" = {out$n_cores <- as.integer(v); i <- i + 2},
      "--neighbor-thr" = {out$neighbor_thr <- as.numeric(v); i <- i + 2},
      "--l" = {out$l <- as.numeric(v); i <- i + 2},
      {i <- i + 1})
  }
  if (is.na(out$sample)) stop("--sample required")
  out
}
args <- parse_args()

bench <- file.path(ROOT, "data/processed/benchmark", args$sample)
stem  <- file.path(bench, sprintf("%s_bench_%s", args$sample, args$variant))
resdir <- file.path(bench, "results"); dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
misty_folder <- file.path(resdir, sprintf("misty_%s", args$variant))
unlink(misty_folder, recursive = TRUE)   # run_misty appends; start clean

expr   <- readMM(paste0(stem, "_expr.mtx"))                 # cells x genes
genes  <- read.csv(paste0(stem, "_genes.csv"))[[1]]
markers <- readLines(paste0(stem, "_misty_markers.txt"))
markers <- markers[nzchar(markers)]
coords <- read.csv(paste0(stem, "_coords.csv"))             # x, y
colnames(expr) <- genes
keep <- intersect(markers, genes)
em <- as.data.frame(as.matrix(expr[, keep, drop = FALSE]))
colnames(em) <- keep
cat(sprintf("  [misty] %d cells, %d target markers; juxta.thr=%g, para.l=%g, workers=%d\n",
            nrow(em), length(keep), args$neighbor_thr, args$l, args$n_cores))
stopifnot(nrow(em) == nrow(coords))

if (args$n_cores > 1) plan(multisession, workers = args$n_cores) else plan(sequential)

t0 <- Sys.time()
views <- create_initial_view(em) |>
  add_juxtaview(positions = coords, neighbor.thr = args$neighbor_thr) |>
  add_paraview(positions = coords, l = args$l, family = "gaussian")
run_misty(views, results.folder = misty_folder, seed = SEED)
res <- collect_results(misty_folder)
cat(sprintf("  [misty] modelled %d targets in %.1f s\n",
            length(keep), as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# tidy CSVs for the Python evaluator (write every data-frame element)
for (nm in names(res)) {
  if (is.data.frame(res[[nm]])) {
    f <- file.path(resdir, sprintf("misty_%s_%s.csv", args$variant, nm))
    write.csv(res[[nm]], f, row.names = FALSE)
  }
}
saveRDS(res, file.path(resdir, sprintf("misty_%s_results.rds", args$variant)))
cat(sprintf("  [misty] wrote tidy CSVs -> %s/misty_%s_*.csv\n", basename(resdir), args$variant))
