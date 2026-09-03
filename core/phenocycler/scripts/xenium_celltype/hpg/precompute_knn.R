#!/usr/bin/env Rscript
# Build the feature-space kNN once per donor and cache it beside pcs.mtx.
#
# KODAMA.matrix computes `nn_call(data, data, neighbors + 1, metrics)` from `data` alone --
# not from `spatial`, and not from anything the Monte-Carlo loop produces. So it is the same
# array for the spatial arm, the nonspatial arm and every shard of both, and at the cohort's
# largest donor (1,109,281 cells x 501 neighbours, two double arrays) it is ~8.8 GB and the
# dominant fixed cost of every one of those calls. Computing it here means kodama_merge.R
# reads it instead of rebuilding it, twice per donor.
#
# `neighbors` is reproduced exactly, including that KODAMA derives it from the ORIGINAL
# --landmarks before the `if (nsample <= landmarks)` adjustment further down, and that it is
# capped at 500 -- which is why raising landmarks above 500 buys a better per-iteration
# partition but not a wider neighbour graph.
#
#   Rscript precompute_knn.R --input <prep dir> --landmarks 2000 --n-cores 16

suppressWarnings(suppressMessages({library(Matrix); library(Rnanoflann); library(jsonlite)}))
a <- commandArgs(trailingOnly = TRUE)
getopt <- function(k, d) { i <- match(k, a); if (is.na(i)) d else a[[i + 1]] }
INPUT <- getopt("--input", NA); LAND <- as.integer(getopt("--landmarks", "2000"))
CORES <- as.integer(getopt("--n-cores", "16")); FORCE <- "--force" %in% a
if (is.na(INPUT)) stop("--input is required")

X <- as.matrix(readMM(file.path(INPUT, "pcs.mtx"))); storage.mode(X) <- "double"
nsample <- nrow(X)
neighbors <- max(1L, floor(min(c(LAND, nsample * 0.75 - 1), 500)))
out <- file.path(INPUT, sprintf("knn_%d.rds", neighbors))
if (file.exists(out) && !FORCE) {
  cat(sprintf("[knn] %s already present -- skipping\n", basename(out))); quit(status = 0)
}

cat(sprintf("[knn] %d cells x %d comps, k=%d, cores=%d\n", nsample, ncol(X), neighbors, CORES))
t0 <- Sys.time()
knn <- if (CORES > 1L) {
  z <- try(Rnanoflann::nn(X, X, neighbors + 1L, method = "euclidean",
                          parallel = TRUE, cores = CORES), silent = TRUE)
  if (inherits(z, "try-error")) Rnanoflann::nn(X, X, neighbors + 1L, method = "euclidean") else z
} else Rnanoflann::nn(X, X, neighbors + 1L, method = "euclidean")
dn <- if (!is.null(knn$distances)) "distances" else "distance"
if (is.null(knn[[dn]])) stop("Rnanoflann::nn did not return distances")
knn <- list(indices = knn$indices[, -1, drop = FALSE], distances = knn[[dn]][, -1, drop = FALSE])
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

saveRDS(knn, out, compress = FALSE)
write_json(list(complete = TRUE, n_cells = nsample, n_comps = ncol(X), landmarks = LAND,
                neighbors = neighbors, n_cores = CORES, seconds = secs,
                indices_typeof = typeof(knn$indices),
                bytes = as.numeric(file.size(out)),
                Rnanoflann = as.character(packageVersion("Rnanoflann"))),
           file.path(INPUT, sprintf("knn_%d.json", neighbors)), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("[knn] %.1f s -> %s (%.1f GiB)\n", secs, basename(out), file.size(out) / 2^30))
