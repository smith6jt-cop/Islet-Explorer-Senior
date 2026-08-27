#!/usr/bin/env Rscript
# RELEASE GATE. Proves, for the KODAMA build actually installed here, that
#     kodama_shard.R x S  ->  kodama_merge.R
# computes what a single KODAMA.matrix(M = m*S) + KODAMA.visualization would have computed.
#
# It runs the production scripts themselves, not a copy of their logic, so it fails if
# either script drifts. It asserts on three things, in increasing distance from the source:
#   res              -- the shard claim: seed = base + m*shard tiles the iteration set
#   indices/distances-- the merge claim: the re-implemented dissimilarity reweighting
#   embedding        -- the whole path, end to end
#
# It is CHEAP (seconds to a couple of minutes at the default size), so submit_all.sh runs it
# every time and refuses to submit if it fails. An upgraded KODAMA that changes
# `one_iteration`'s seeding or `run_parallel_chunks`' arithmetic will be caught here rather
# than by a cohort of embeddings that are quietly not what the manifest says they are.
#
#   Rscript check_shard_equivalence.R [--n 1200] [--m 3] [--shards 3] [--arm both] [--input DIR]
# --input takes a real prep directory (pcs.mtx/coords.csv/cell_id.csv) and crops it, which is
# the version worth running before a production submission; with no --input it uses a
# synthetic matrix, which still exercises every line.

suppressWarnings(suppressMessages({library(Matrix); library(KODAMA); library(jsonlite)}))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                       value = TRUE)[1])))

a <- commandArgs(trailingOnly = TRUE)
getopt <- function(k, d) { i <- match(k, a); if (is.na(i)) d else a[[i + 1]] }
N       <- as.integer(getopt("--n", "1200"))
M_SHARD <- as.integer(getopt("--m", "3"))
NSHARD  <- as.integer(getopt("--shards", "3"))
ARMS    <- if (getopt("--arm", "both") == "both") c("nonspatial", "spatial") else getopt("--arm", "both")
INPUT   <- getopt("--input", NA_character_)
LAND    <- as.integer(getopt("--landmarks", "150"))
CORES   <- as.integer(getopt("--n-cores", "1"))
SEED    <- as.integer(getopt("--seed", "0"))
M_TOTAL <- M_SHARD * NSHARD

tmp <- file.path(tempdir(), sprintf("kodama_equiv_%d", Sys.getpid()))
inp <- file.path(tmp, "input"); out <- file.path(tmp, "out")
dir.create(inp, recursive = TRUE); dir.create(out, recursive = TRUE)

if (!is.na(INPUT)) {
  X <- as.matrix(readMM(file.path(INPUT, "pcs.mtx"))); storage.mode(X) <- "double"
  co <- as.matrix(read.csv(file.path(INPUT, "coords.csv")))
  ids <- read.csv(file.path(INPUT, "cell_id.csv"))[[1]]
  keep <- seq_len(min(N, nrow(X)))
  X <- X[keep, , drop = FALSE]; co <- co[keep, , drop = FALSE]; ids <- ids[keep]
} else {
  set.seed(11); p <- 20
  X <- matrix(rnorm(N * p), N, p)
  X[seq_len(N %/% 3), 1:5] <- X[seq_len(N %/% 3), 1:5] + 3
  co <- cbind(runif(N, 0, 400), runif(N, 0, 400)); ids <- sprintf("cell_%05d", seq_len(N))
}
invisible(writeMM(as(X, "CsparseMatrix"), file.path(inp, "pcs.mtx")))
write.csv(co, file.path(inp, "coords.csv"), row.names = FALSE)
write.csv(data.frame(cell_id = ids), file.path(inp, "cell_id.csv"), row.names = FALSE)
# writeMM round-trips through text; compare against what the scripts will actually read.
X <- as.matrix(readMM(file.path(inp, "pcs.mtx"))); storage.mode(X) <- "double"

Rbin <- file.path(R.home("bin"), "Rscript")
fails <- character(0)
for (arm in ARMS) {
  cat(sprintf("\n=== arm = %s ===\n", arm))
  for (s in seq_len(NSHARD) - 1L) {
    st <- system2(Rbin, c(file.path(HERE, "kodama_shard.R"),
      "--input", inp, "--output", out, "--arm", arm, "--m", M_SHARD, "--shard", s,
      "--base-seed", SEED, "--landmarks", LAND, "--ncomp", min(50L, ncol(X)),
      "--n-cores", CORES, "--force"), stdout = NULL, stderr = NULL)
    if (st != 0) stop(sprintf("kodama_shard.R exited %d (arm=%s shard=%d)", st, arm, s))
  }
  st <- system2(Rbin, c(file.path(HERE, "kodama_merge.R"),
    "--input", inp, "--output", out, "--arm", arm, "--M", M_TOTAL, "--m", M_SHARD,
    "--landmarks", LAND, "--n-cores", CORES, "--seed", SEED, "--block", 500,
    "--keep-cols", 64, "--keep-blocks"), stdout = NULL, stderr = NULL)
  if (st != 0) stop(sprintf("kodama_merge.R exited %d (arm=%s)", st, arm))

  # ---- the reference: one monolithic call, same settings ----
  invisible(capture.output(mono <- if (arm == "spatial")
      KODAMA.matrix(X, spatial = co, M = M_TOTAL, ncomp = min(50L, ncol(X)),
                    landmarks = LAND, spatial.resolution = 0.3, n.cores = CORES, seed = SEED)
    else
      KODAMA.matrix(X, M = M_TOTAL, ncomp = min(50L, ncol(X)), landmarks = LAND,
                    n.cores = CORES, seed = SEED)))
  cfg <- KODAMA::config.umap.default
  cfg$define.n.cores <- TRUE; cfg$n_threads <- CORES
  cfg$n_sgd_threads <- 1L; cfg$random_state <- SEED
  mono_emb <- as.matrix(KODAMA.visualization(mono, method = "UMAP", config = cfg))

  # ---- reassemble what the merge produced ----
  parts <- lapply(sort(list.files(out, sprintf("^res_%s_shard[0-9]+\\.rds$", arm),
                                  full.names = TRUE)), readRDS)
  ord <- order(vapply(parts, function(p) p$k_global[1], integer(1)))
  got_res <- do.call(rbind, lapply(parts[ord], `[[`, "res"))
  storage.mode(got_res) <- "double"

  blk <- sort(list.files(file.path(out, sprintf("dissim_%s", arm)), full.names = TRUE))
  zz <- lapply(blk, readRDS)
  got_idx <- do.call(rbind, lapply(zz, `[[`, "indices"))
  got_dst <- do.call(rbind, lapply(zz, `[[`, "distances"))
  got_emb <- readRDS(file.path(out, sprintf("embedding_%s.rds", arm)))$embedding
  csv_emb <- as.matrix(read.csv(file.path(out, sprintf("embedding_%s.csv", arm)))[, c("dim1", "dim2")])
  csv_err <- max(abs(unname(csv_emb) - unname(got_emb)))

  nu <- min(45L, nrow(X) - 1L, ncol(got_idx))     # the columns the UMAP actually reads
  chk <- list(
    res        = identical(mono$res, got_res),
    acc        = identical(mono$acc, unlist(lapply(parts[ord], `[[`, "acc"))),
    knn_idx    = identical(mono$knn_Rnanoflann$indices[, seq_len(nu), drop = FALSE],
                           got_idx[, seq_len(nu), drop = FALSE]),
    knn_dist   = identical(mono$knn_Rnanoflann$distances[, seq_len(nu), drop = FALSE],
                           got_dst[, seq_len(nu), drop = FALSE]),
    embedding  = isTRUE(all.equal(unname(mono_emb), unname(got_emb), tolerance = 0)))
  for (nm in names(chk))
    cat(sprintf("  %-10s %s\n", nm, if (chk[[nm]]) "IDENTICAL" else "*** DIFFERS ***"))
  # Reported, never gated on: the CSV is a text rendering of the RDS above, so its error is
  # a property of write.csv, not of the computation. Downstream reads the CSV, so the size
  # of that error belongs in the log where anyone can see it is ~1e-13 and not ~1e-3.
  cat(sprintf("  %-10s %.3g (write.csv text round trip; not gated)\n", "csv_err", csv_err))
  bad <- names(chk)[!unlist(chk)]
  if (length(bad)) fails <- c(fails, sprintf("%s: %s", arm, paste(bad, collapse = ", ")))
}

unlink(tmp, recursive = TRUE)
cat("\n")
if (length(fails)) {
  cat(sprintf("SHARD EQUIVALENCE FAILED (KODAMA %s, %s)\n",
              as.character(packageVersion("KODAMA")), R.version.string))
  for (f in fails) cat("  ", f, "\n")
  quit(status = 1)
}
cat(sprintf("SHARD EQUIVALENCE OK  n=%d M=%d (%d x %d) landmarks=%d arms=%s | KODAMA %s | %s\n",
            nrow(X), M_TOTAL, NSHARD, M_SHARD, LAND, paste(ARMS, collapse = "+"),
            as.character(packageVersion("KODAMA")), R.version.string))
