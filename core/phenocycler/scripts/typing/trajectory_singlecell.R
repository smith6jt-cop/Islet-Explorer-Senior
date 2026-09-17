#!/usr/bin/env Rscript
# Single-cell trajectory over islet and peri-islet cells, by the published algorithm.
#
# ALGORITHM, pinned from the paper's own Methods before this was written -- Steenbuck, Damond,
# Engler, Kusmartseva, Posgai, Drotar, Williams, de Souza, Brusko, Brusko, Wasserfall, Atkinson,
# Bodenmiller, bioRxiv 2025.03.05.641526 (v1, 11 March 2025), section "Pseudotime analysis":
#
#   "Pseudotime analysis of beta-cell expression profiles was performed using slingshot (v. 2.10.0).
#    We used default parameters to identify the global lineage structure, and fit smooth branching
#    curves to infer the pseudotime variables for each beta-cell. Notably, we did not provide start
#    and leaf nodes or cluster information and thus, slingshot inferred the global trajectory in an
#    unbiased manner. As input, we used diffusion map embedded single-cell beta-cell expression
#    profiles, which were computed using destiny (v. 3.16.0)."
#
# So: destiny diffusion map -> slingshot, DEFAULT PARAMETERS, no start/leaf nodes, no cluster
# labels. Nothing below chooses a start point or a lineage structure; both are inferred.
#
# Versions here are destiny 3.22.0 and slingshot 2.16.0 against the paper's 3.16.0 and 2.10.0 --
# recorded rather than pinned, because the older pair does not install on R 4.5.
#
# The cells and the markers are this study's, not the paper's; see the docstring of
# trajectory_singlecell_features.py for why, and what that changes.
#
# WHY THE ORDERING IS FIT ON A SUBSAMPLE, AND WHY NOTHING IS PROJECTED ONTO IT.
# destiny is roughly quadratic here -- 25,000 cells take 94 s, so the full 800,562 would be days.
# The first version of this script fit a subsample and extended the diffusion map to the rest with
# destiny's Nystrom predictor, but that forms an (out-of-sample x fit) kernel block: 37 GB resident
# for a 5,000-cell fit, and ~200 chunked passes at a realistic fit size. More to the point, the
# paper has no projection step -- it orders the cells it has. So this orders a proportional
# per-donor subsample and stops there. Cells outside the sample carry NO pseudotime, and the file
# says so; an imputed ordering would look identical to a measured one and is not one.
# Proportional per-donor quotas are this lane's own established sampling, from the coarse embedding.
#
#   Rscript core/phenocycler/scripts/typing/trajectory_singlecell.R [n_fit]

suppressMessages({
  library(arrow); library(destiny); library(slingshot)
})

set.seed(20260914)          # the subsample is reproducible; nothing else is random
args <- commandArgs(trailingOnly = TRUE)
N_FIT <- if (length(args) >= 1) as.integer(args[[1]]) else 100000L
N_DCS <- 20L                # diffusion components carried into slingshot

root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", "..", ".."))
DIR <- file.path(root, "data/development/typing/trajectory")

cat("reading features\n"); flush.console()
cells <- as.data.frame(read_parquet(file.path(DIR, "singlecell_features.parquet")))
meta_cols <- c("donor_id", "object_id", "islet_key", "region", "phenotype", "donor_status",
               "stain_group", "X_centroid", "Y_centroid", "signed_dist_islet_um")
markers <- setdiff(names(cells), meta_cols)
X <- as.matrix(cells[, markers])
cat(sprintf("  %s cells x %d markers\n", format(nrow(X), big.mark = ","), ncol(X)))
flush.console()

# --- Proportional per-donor quotas ----------------------------------------------------------
if (N_FIT < nrow(X)) {
  # as.numeric: N_FIT and the table counts are both integer, and their product overflows 32-bit
  # for any realistic fit size (50,000 x 122,951). R warns and returns NA, and sample() then
  # fails on an NA size -- loudly here, but the same arithmetic silently truncates elsewhere.
  quota <- round(as.numeric(N_FIT) * table(cells$donor_id) / nrow(cells))
  fit_idx <- sort(unlist(lapply(names(quota), function(d) {
    rows <- which(cells$donor_id == d)
    sample(rows, min(quota[[d]], length(rows)))
  })))
} else {
  fit_idx <- seq_len(nrow(X))
}
cat(sprintf("  ordering %s cells (%.1f%% of the population), proportional per donor\n",
            format(length(fit_idx), big.mark = ","), 100 * length(fit_idx) / nrow(X)))
print(table(cells$donor_id[fit_idx]))
flush.console()

cat("diffusion map (destiny, default parameters)\n"); flush.console()
t0 <- Sys.time()
dm <- DiffusionMap(X[fit_idx, , drop = FALSE], n_eigs = N_DCS, verbose = FALSE)
cat(sprintf("  %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()

rd <- as.matrix(eigenvectors(dm))[, seq_len(N_DCS), drop = FALSE]
colnames(rd) <- paste0("DC", seq_len(N_DCS))

cat("slingshot (default parameters, no start/leaf nodes, no cluster labels)\n"); flush.console()
t0 <- Sys.time()
sl <- slingshot(rd)
curves <- slingCurves(sl)
cat(sprintf("  %.1f min; lineages inferred: %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(curves)))
flush.console()

pt <- slingPseudotime(sl)
colnames(pt) <- paste0("pseudotime_lineage_", seq_len(ncol(pt)))

out <- cells[fit_idx, meta_cols]
out <- cbind(out, as.data.frame(pt), as.data.frame(rd))
out$pseudotime <- pt[, 1]

# --- Reported, never gated. A result is not a verification target. ---------------------------
finite <- is.finite(out$pseudotime)
cat("\nreported\n")
cat(sprintf("  cells ordered              %s (%s without a finite pseudotime)\n",
            format(sum(finite), big.mark = ","), format(sum(!finite), big.mark = ",")))
for (col in c("region", "donor_status")) {
  med <- tapply(out$pseudotime[finite], out[[col]][finite], median)
  cat(sprintf("  median pseudotime by %-13s %s\n", col,
              paste(sprintf("%s=%.3g", names(med), med), collapse = "  ")))
}
med_ph <- sort(tapply(out$pseudotime[finite], out$phenotype[finite], median))
cat("  median pseudotime by phenotype\n")
for (n in names(med_ph)) cat(sprintf("      %-14s %.4g\n", n, med_ph[[n]]))

eta2 <- function(v, g) {
  v <- v[is.finite(v)]; g <- g[is.finite(out$pseudotime)][seq_along(v)]
  grand <- mean(v)
  between <- sum(vapply(unique(g), function(u)
    sum(g == u) * (mean(v[g == u]) - grand)^2, numeric(1)))
  total <- sum((v - grand)^2)
  if (total > 0) between / total else NA_real_
}
cat(sprintf("  eta2 donor %.3f  |  eta2 region %.3f  |  eta2 donor_status %.3f  |  eta2 stain_group %.3f\n",
            eta2(out$pseudotime, out$donor_id), eta2(out$pseudotime, out$region),
            eta2(out$pseudotime, out$donor_status), eta2(out$pseudotime, out$stain_group)))

write_parquet(out, file.path(DIR, "singlecell_pseudotime.parquet"), compression = "zstd")
jsonlite::write_json(list(
  algorithm = "destiny DiffusionMap -> slingshot, default parameters, no start/leaf nodes, no cluster labels",
  source = "Steenbuck, Damond, ... Bodenmiller, bioRxiv 2025.03.05.641526, 'Pseudotime analysis'",
  versions = list(destiny = as.character(packageVersion("destiny")),
                  slingshot = as.character(packageVersion("slingshot")),
                  R = R.version.string,
                  paper_versions = list(destiny = "3.16.0", slingshot = "2.10.0")),
  population = "islet core and peri-islet (20 um) cells; Excluded removed",
  n_population = nrow(X), n_ordered = length(fit_idx),
  sampling = "proportional per-donor quotas, seed 20260914",
  not_ordered = "cells outside the subsample carry NO pseudotime; none is imputed",
  n_markers = length(markers), markers = markers, n_dcs = N_DCS,
  n_lineages = length(curves),
  eigenvalues = unname(eigenvalues(dm)[seq_len(N_DCS)])
), file.path(DIR, "singlecell_trajectory.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8)

cat(sprintf("\n%s cells, %d lineage(s) -> %s\n", format(nrow(out), big.mark = ","),
            length(curves), file.path(DIR, "singlecell_pseudotime.parquet")))
