#!/usr/bin/env Rscript
# Single-cell trajectory over islet and peri-islet cells, FIT WITHIN DONOR.
#
# ALGORITHM, pinned from the paper's own Methods -- Steenbuck, Damond, Engler, Kusmartseva, Posgai,
# Drotar, Williams, de Souza, Brusko, Brusko, Wasserfall, Atkinson, Bodenmiller, bioRxiv
# 2025.03.05.641526 (v1, 11 March 2025), section "Pseudotime analysis":
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
# ============================================================================================
# WITHIN DONOR, 2026-09-17, by maintainer decision. The previous build fit ONE diffusion map over
# all donors' cells together. Acquisition batch is aliased with disease status at the coarse level
# in this cohort -- every ND donor is in batch 1, every batch-2 donor is Aab+ or T1D -- so a
# pooled embedding measures acquisition, and whatever measures acquisition is then read as
# progression. The same shape was hit by the pLN clustering lane, which handled it by clustering
# within donor and relating across donors afterwards with a permutation null on every cross-donor
# quantity. This follows that.
#
# So: ONE DIFFUSION MAP AND ONE SLINGSHOT CURVE PER DONOR. No donor's cells ever enter another
# donor's kernel. There is no pooled embedding anywhere in this file.
#
# THREE CONSEQUENCES, none of them hidden.
#
#  1. EACH DONOR'S AXIS HAS ITS OWN SCALE. slingshot's pseudotime is arc length along that donor's
#     own curve through that donor's own diffusion space; the units are not shared. So the column
#     every consumer reads, `pseudotime`, is the WITHIN-DONOR QUANTILE -- rank / (n + 1) over that
#     donor's ordered cells, in (0, 1). The raw arc length is kept beside it as `pseudotime_raw`
#     for anyone who wants it, and is not what anything downstream reads.
#
#  2. EACH DONOR'S AXIS HAS ITS OWN DIRECTION. slingshot given no start node picks an endpoint
#     arbitrarily, so donor A's "late" can be donor B's "early" and the two are not comparable
#     until a direction is fixed. The direction is fixed by a DECLARED CONVENTION, applied
#     identically to every donor: flip the axis if the donor's pooled core cells sit ABOVE its
#     pooled peri cells, so that low pseudotime is the core end. That is a coordinate convention on
#     a sign-arbitrary axis, not a measurement -- and the cost is stated rather than buried: it
#     makes the COHORT-LEVEL SIGN of peri-minus-core a convention. What is NOT forced by it, and is
#     measured, is how far apart the two regions sit, how much islets within a donor differ, and
#     which islets run the other way. `orientation_flipped` records, per donor, whether the raw
#     axis already ran that way.
#
#  3. THE DONOR IS NO LONGER A SOURCE OF VARIANCE INSIDE THE FIT, so the ordering cannot be a donor
#     or batch axis by construction. Relating the donors to each other afterwards is a SEPARATE
#     step with its own null: `trajectory_crossdonor_relation.py`.
#
# WHY THE ORDERING IS STILL FIT ON A SUBSAMPLE, AND WHY NOTHING IS PROJECTED ONTO IT.
# destiny is roughly quadratic, so the population cannot be ordered whole. The first version of
# this script extended a subsample to the rest with destiny's Nystrom predictor, which forms an
# (out-of-sample x fit) kernel block -- 37 GB resident for a 5,000-cell fit. More to the point, the
# paper has no projection step -- it orders the cells it has. So this orders a per-donor subsample
# and stops there. Cells outside the sample carry NO pseudotime, and the file says so; an imputed
# ordering would look identical to a measured one and is not one.
#
# Fitting within donor is what makes the sample larger, not smaller: the cost is quadratic in the
# cells of ONE donor, so twenty separate fits of N cells cost 20*N^2 where one pooled fit of 20*N
# costs 400*N^2. The per-donor budget below therefore orders several times as many cells as the
# pooled build did, in less wall time.
#
#   Rscript core/multiplex5000/scripts/typing/trajectory_singlecell.R [n_fit_per_donor]

suppressMessages({
  library(arrow); library(destiny); library(slingshot)
})

# A REPLICATE PASS, and why one exists. A second argument is a seed offset: the same fit over an
# independently drawn subsample of the same donor's cells, written beside the first rather than
# over it. `trajectory_reproducibility.py` then measures, per donor, whether the two passes order
# the cells the same way. That is internal coherence -- does this build reproduce itself -- and not
# a biological check; it gates nothing, and a donor whose two passes disagree keeps its ordering
# and carries the number. It was added because on the first twenty-donor build two of three donors
# probed disagreed with their own replicate, which nothing in the artifact would otherwise say.
#
#   Rscript core/multiplex5000/scripts/typing/trajectory_singlecell.R 25000 99001

args <- commandArgs(trailingOnly = TRUE)
#: Cells ordered per donor. A donor with fewer eligible cells contributes all of them.
N_FIT_PER_DONOR <- if (length(args) >= 1) as.integer(args[[1]]) else 15000L
REPLICATE <- if (length(args) >= 2) as.integer(args[[2]]) else NA_integer_
SUFFIX <- if (is.na(REPLICATE)) "" else "_replicate"
N_DCS <- 20L                # diffusion components carried into slingshot
#: Below this a donor gets no ordering at all and carries the reason. A diffusion map over a few
#: hundred cells is not an ordering of anything; the value is declared, not tuned to an outcome.
MIN_CELLS_FIT <- 500L
BASE_SEED <- 20260917L

root <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", "..", ".."))
DIR <- file.path(root, "data/development/typing/trajectory")

cat("reading features\n"); flush.console()
cells <- as.data.frame(read_parquet(file.path(DIR, "singlecell_features.parquet")))
meta_cols <- c("donor_id", "object_id", "islet_key", "region", "phenotype", "donor_status",
               "stain_group", "X_centroid", "Y_centroid", "signed_dist_islet_um")
markers <- setdiff(names(cells), meta_cols)
cat(sprintf("  %s cells x %d markers over %d donors\n", format(nrow(cells), big.mark = ","),
            length(markers), length(unique(cells$donor_id))))
flush.console()

donors <- sort(unique(as.character(cells$donor_id)))
pieces <- list()
report <- list()

for (d in donors) {
  rows <- which(as.character(cells$donor_id) == d)
  Xd <- as.matrix(cells[rows, markers])
  # A missing measurement is unavailable, never zero: a cell without an aligned marker row cannot
  # enter a kernel that is a distance over those markers. It is dropped WITH ITS REASON, not
  # filled. These are the re-exported batch-2 cells whose object_id has no row in the aligned
  # values -- the join their rebuild reuses.
  complete <- which(stats::complete.cases(Xd))
  n_pop <- length(rows)
  n_drop <- n_pop - length(complete)

  entry <- list(donor_id = d,
                donor_status = as.character(cells$donor_status[rows[1]]),
                stain_group = as.character(cells$stain_group[rows[1]]),
                n_population = n_pop,
                n_eligible = length(complete),
                n_dropped_no_aligned_marker_row = n_drop,
                n_islets_population = length(unique(cells$islet_key[rows])))

  if (length(complete) < MIN_CELLS_FIT) {
    entry$n_ordered <- 0L
    entry$state <- "unavailable"
    entry$reason <- sprintf("below_min_cells_for_a_diffusion_map_%d(eligible=%d)",
                            MIN_CELLS_FIT, length(complete))
    cat(sprintf("  %-6s SKIPPED -- %s\n", d, entry$reason)); flush.console()
    report[[d]] <- entry
    next
  }

  seed <- BASE_SEED + as.integer(sub("^0+", "", gsub("[^0-9]", "", d)))
  if (!is.na(REPLICATE)) seed <- seed + REPLICATE
  set.seed(seed)
  take <- if (length(complete) > N_FIT_PER_DONOR) sort(sample(complete, N_FIT_PER_DONOR)) else complete
  idx <- rows[take]
  Xfit <- Xd[take, , drop = FALSE]
  # A marker with no variation inside this donor contributes nothing to a distance and can make the
  # kernel degenerate. Dropped per donor, counted, never silently.
  keep_m <- which(apply(Xfit, 2, function(v) stats::sd(v) > 0))
  entry$markers_used <- length(keep_m)
  entry$markers_constant_within_donor <- setdiff(markers, markers[keep_m])
  Xfit <- Xfit[, keep_m, drop = FALSE]

  t0 <- Sys.time()
  dm <- DiffusionMap(Xfit, n_eigs = N_DCS, verbose = FALSE)
  rd <- as.matrix(eigenvectors(dm))[, seq_len(N_DCS), drop = FALSE]
  colnames(rd) <- paste0("DC", seq_len(N_DCS))
  sl <- slingshot(rd)
  curves <- slingCurves(sl)
  pt <- slingPseudotime(sl)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  raw <- pt[, 1]
  region <- as.character(cells$region[idx])
  # --- the declared orientation convention; see consequence 2 in the header --------------------
  med_core <- stats::median(raw[region == "core"], na.rm = TRUE)
  med_peri <- stats::median(raw[region == "peri"], na.rm = TRUE)
  flipped <- is.finite(med_core) && is.finite(med_peri) && med_core > med_peri
  if (flipped) raw <- -raw
  finite <- is.finite(raw)
  # within-donor quantile in (0, 1); the only scale that is shared between donors, and shared only
  # because it is a rank, not because the axes are the same axis.
  q <- rep(NA_real_, length(raw))
  q[finite] <- rank(raw[finite], ties.method = "average") / (sum(finite) + 1)

  out <- cells[idx, meta_cols]
  out <- cbind(out, as.data.frame(rd))
  out$pseudotime_raw <- raw
  out$pseudotime <- q
  out$fit_scope <- "within_donor"
  pieces[[d]] <- out

  entry$n_ordered <- length(idx)
  entry$sampled_fraction <- length(idx) / length(complete)
  entry$n_islets_ordered <- length(unique(cells$islet_key[idx]))
  entry$n_lineages <- length(curves)
  entry$lineage_used <- 1L
  entry$orientation_flipped <- flipped
  entry$n_nonfinite_pseudotime <- sum(!finite)
  entry$eigenvalues <- unname(eigenvalues(dm)[seq_len(N_DCS)])
  entry$fit_minutes <- round(mins, 2)
  entry$seed <- seed
  entry$state <- "measured"
  entry$reason <- ""
  report[[d]] <- entry

  cat(sprintf("  %-6s %-4s  ordered %6s of %7s eligible (%5.1f%%)  islets %5d  lineages %d  %s  %.1f min\n",
              d, entry$donor_status, format(length(idx), big.mark = ","),
              format(length(complete), big.mark = ","), 100 * entry$sampled_fraction,
              entry$n_islets_ordered, entry$n_lineages,
              if (flipped) "flipped" else "as-fit", mins))
  flush.console()
}

stopifnot(length(pieces) > 0)
all_out <- do.call(rbind, pieces)
rownames(all_out) <- NULL

# --- Reported, never gated. A result is not a verification target. ---------------------------
cat("\nreported -- WITHIN DONOR, so none of this is a cross-donor claim\n")
n_ordered <- nrow(all_out)
n_population <- nrow(cells)
cat(sprintf("  cells ordered              %s of %s in the population (%.1f%%)\n",
            format(n_ordered, big.mark = ","), format(n_population, big.mark = ","),
            100 * n_ordered / n_population))
cat(sprintf("  donors ordered             %d of %d\n", length(pieces), length(donors)))

med_reg <- tapply(all_out$pseudotime, all_out$region, median)
cat(sprintf("  median within-donor quantile by region   %s\n",
            paste(sprintf("%s=%.3f", names(med_reg), med_reg), collapse = "  ")))
med_ph <- sort(tapply(all_out$pseudotime, all_out$phenotype, median))
cat("  median within-donor quantile by phenotype\n")
for (n in names(med_ph)) cat(sprintf("      %-16s %.3f\n", n, med_ph[[n]]))

eta2 <- function(v, g) {
  ok <- is.finite(v); v <- v[ok]; g <- as.character(g)[ok]
  grand <- mean(v)
  total <- sum((v - grand)^2)
  if (total <= 0) return(NA_real_)
  between <- sum(vapply(unique(g), function(u) sum(g == u) * (mean(v[g == u]) - grand)^2,
                        numeric(1)))
  between / total
}
cat(sprintf("  eta2 donor %.4f  |  eta2 region %.4f  |  eta2 donor_status %.4f  |  eta2 stain_group %.4f\n",
            eta2(all_out$pseudotime, all_out$donor_id), eta2(all_out$pseudotime, all_out$region),
            eta2(all_out$pseudotime, all_out$donor_status),
            eta2(all_out$pseudotime, all_out$stain_group)))
cat("  (donor eta2 on a within-donor quantile is near zero BY CONSTRUCTION -- every donor is\n")
cat("   rank-uniform on its own axis. It is not evidence that batch was removed. The cross-donor\n")
cat("   relation is measured separately, with a permutation null, by\n")
cat("   trajectory_crossdonor_relation.py.)\n")

write_parquet(all_out, file.path(DIR, sprintf("singlecell_pseudotime%s.parquet", SUFFIX)),
              compression = "zstd")
jsonlite::write_json(list(
  fit_scope = "within_donor",
  algorithm = paste("destiny DiffusionMap -> slingshot, default parameters, no start/leaf nodes,",
                    "no cluster labels; ONE FIT PER DONOR, no pooled embedding"),
  source = "Steenbuck, Damond, ... Bodenmiller, bioRxiv 2025.03.05.641526, 'Pseudotime analysis'",
  versions = list(destiny = as.character(packageVersion("destiny")),
                  slingshot = as.character(packageVersion("slingshot")),
                  R = R.version.string,
                  paper_versions = list(destiny = "3.16.0", slingshot = "2.10.0")),
  population = "islet core and peri-islet (20 um) cells; Excluded removed",
  n_population = n_population, n_ordered = n_ordered,
  n_donors_population = length(donors), n_donors_ordered = length(pieces),
  sampling = sprintf(paste("per-donor subsample, up to %d cells per donor, seed %d + donor id;",
                           "no cross-donor quota and no pooled embedding"),
                     N_FIT_PER_DONOR, BASE_SEED),
  not_ordered = "cells outside the subsample carry NO pseudotime; none is imputed",
  pseudotime_column = paste("`pseudotime` is the WITHIN-DONOR QUANTILE, rank/(n+1) in (0,1).",
                            "`pseudotime_raw` is the oriented slingshot arc length on that",
                            "donor's own axis and is NOT comparable between donors."),
  orientation = paste("declared convention: each donor's axis is flipped if its pooled core cells",
                      "sit above its pooled peri cells, so low pseudotime is the core end. The",
                      "cohort-level SIGN of peri-minus-core is therefore a convention, not a",
                      "measurement; its spread across islets is not."),
  min_cells_fit = MIN_CELLS_FIT,
  n_markers = length(markers), markers = markers, n_dcs = N_DCS,
  per_donor = unname(report),
  replicate_seed_offset = if (is.na(REPLICATE)) NULL else REPLICATE,
  is_replicate = !is.na(REPLICATE)
), file.path(DIR, sprintf("singlecell_trajectory%s.json", SUFFIX)),
   auto_unbox = TRUE, pretty = TRUE, digits = 8)

cat(sprintf("\n%s cells over %d donors -> %s\n", format(n_ordered, big.mark = ","),
            length(pieces), file.path(DIR, sprintf("singlecell_pseudotime%s.parquet", SUFFIX))))
