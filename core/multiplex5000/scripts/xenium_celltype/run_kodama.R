#!/usr/bin/env Rscript
# KODAMA embedding for Xenium cell typing -- one donor, one arm.
#
# WHY PCs AND NOT GENES. Measured on donor 6539, identical n/landmarks/M (5000/500/5):
# 5,101 genes = 1583.2 s; 50 principal components = 44.7 s. A 35x difference, and the cost
# is the input width, not the cell count. Feeding the gene matrix directly is what makes
# KODAMA "the slowest method" in the sibling spatial-domain benchmark; it is not inherent.
# PCA is unsupervised and identity-agnostic, so it cannot smuggle the answer in -- but the
# component count IS a choice, and the driver records it and ships a stability check.
#
# TWO ARMS, DELIBERATELY. `--arm spatial` passes `spatial=coords`, which lets a cell's
# neighbourhood inform its own representation. That is the right thing for spatial DOMAINS
# and the wrong thing for a cell TYPE call in this study: the question is whether immune
# cells sit near islets, so a type that already knows where the cell is makes the answer
# partly circular. Both arms are run and compared; the difference between them is the
# measurement of how much location leaks into identity.
#
# Reads:  <in>/pcs.mtx (cells x components), <in>/coords.csv (x,y microns), <in>/cell_id.csv
# Writes: <out>/embedding_<arm>.csv  (cell_id, dim1, dim2), <out>/kodama_<arm>.json (timings)

suppressWarnings(suppressMessages({library(Matrix); library(KODAMA); library(jsonlite)}))

p <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  o <- list(input = NA, output = NA, arm = "nonspatial", M = 20L,
            landmarks = 2000L, spatial_resolution = 0.3, n_cores = 16L, seed = 0L)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]; v <- if (i < length(a)) a[[i + 1]] else NA
    switch(k,
      "--input"={o$input<-v; i<-i+2}, "--output"={o$output<-v; i<-i+2},
      "--arm"={o$arm<-v; i<-i+2}, "--M"={o$M<-as.integer(v); i<-i+2},
      "--landmarks"={o$landmarks<-as.integer(v); i<-i+2},
      "--spatial-resolution"={o$spatial_resolution<-as.numeric(v); i<-i+2},
      "--n-cores"={o$n_cores<-as.integer(v); i<-i+2},
      "--seed"={o$seed<-as.integer(v); i<-i+2}, {i<-i+1})
  }
  if (is.na(o$input) || is.na(o$output)) stop("--input and --output are required")
  if (!o$arm %in% c("spatial", "nonspatial")) stop("--arm must be spatial or nonspatial")
  o
}
args <- p()
dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

X <- as.matrix(readMM(file.path(args$input, "pcs.mtx"))); storage.mode(X) <- "double"
ids <- read.csv(file.path(args$input, "cell_id.csv"))[[1]]
co  <- as.matrix(read.csv(file.path(args$input, "coords.csv")))
stopifnot(nrow(X) == length(ids), nrow(X) == nrow(co))

cat(sprintf("[kodama] %s: %d cells x %d comps, M=%d landmarks=%d cores=%d\n",
            args$arm, nrow(X), ncol(X), args$M, args$landmarks, args$n_cores))
t0 <- Sys.time()
# Braces are load-bearing: at top level R parses `if (...) expr` as complete, so an `else`
# beginning the next line is a syntax error.
kk <- if (args$arm == "spatial") {
  KODAMA.matrix(X, spatial = co, M = args$M, ncomp = min(50L, ncol(X)),
                landmarks = args$landmarks, spatial.resolution = args$spatial_resolution,
                n.cores = args$n_cores, seed = args$seed)
} else {
  KODAMA.matrix(X, M = args$M, ncomp = min(50L, ncol(X)),
                landmarks = args$landmarks, n.cores = args$n_cores, seed = args$seed)
}
t1 <- Sys.time()
emb <- as.matrix(KODAMA.visualization(kk, method = "UMAP"))
t2 <- Sys.time()

out <- data.frame(cell_id = as.character(ids), dim1 = emb[, 1], dim2 = emb[, 2])
write.csv(out, file.path(args$output, sprintf("embedding_%s.csv", args$arm)), row.names = FALSE)
write_json(list(arm = args$arm, n_cells = nrow(X), n_comps = ncol(X), M = args$M,
                landmarks = args$landmarks,
                spatial_resolution = if (args$arm == "spatial") args$spatial_resolution else NULL,
                seed = args$seed, n_cores = args$n_cores,
                kodama_version = as.character(packageVersion("KODAMA")),
                matrix_seconds = as.numeric(difftime(t1, t0, units = "secs")),
                visualization_seconds = as.numeric(difftime(t2, t1, units = "secs"))),
           file.path(args$output, sprintf("kodama_%s.json", args$arm)),
           auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("[kodama] %s done in %.1f s -> %s\n", args$arm,
            as.numeric(difftime(t2, t0, units = "secs")), args$output))
