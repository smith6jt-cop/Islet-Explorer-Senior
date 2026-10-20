#!/usr/bin/env Rscript
# The R half of broadtype arm D `celesta`.  Driven by arms/celesta.py as a subprocess; the
# Python side owns the record, the reasons, and every join.  This script runs CELESTA with its
# DEFAULT parameters on one region at a time and writes what the package exposes.
#
#   Rscript celesta_run.R --scratch <dir> --seed <n>
#
#   <dir>/prior.csv             the 1/0/NA signature: col 1 cell type, col 2 lineage code
#                               ('round_parent_number', a flat round-1 set is 1_0_1 .. 1_0_k),
#                               cols 3+ one column per marker (names must match the region files)
#   <dir>/in/region_<r>.csv     object_id, X, Y, <markers...>  (X, Y literal: CELESTA reads
#                               imaging_data_file$X / $Y); one file per region, > 10 finite cells
#   <dir>/env.json              written first: versions, the CELESTA commit, and the package
#                               DEFAULTS actually in force (formals()), so Python can refuse to
#                               record a parameter the package no longer has
#   <dir>/out/region_<r>.csv    object_id, celesta_label ('Unknown' = unassigned), celesta_prob
#                               (the assigned type's entry of @current_cell_prob; NA if Unknown),
#                               celesta_anchor (assigned in the anchor step, before neighbours)
#   <dir>/out/region_<r>.status.json   status ok | too_small | activation_nan | error, timing, message
#   <dir>/celesta_stdout.log    CELESTA's own print()/cat() chatter (sink'ed), for the record
#
# Why one region per call inside one process: the per-region status file is written as soon as
# the region finishes, so when the Python side kills this process at the wall-clock cap the
# finished regions are kept and only the unfinished ones become Unavailable.
#
# Method notes that live here because they are about the package, not the arm:
#   * CreateCelestaObject ignores its `cofactor` argument and calls GetMarkerExpMatrix with the
#     literal 10, so the transform is asinh(x / 10) regardless.  Recorded as cofactor_effective.
#   * FitGmmModel fits Rmixmod 'Gaussian_p_Lk_Ck' with equal.proportions = TRUE -- a 50/50
#     mixture per marker.  The protocol flags that in-method prior; nothing here tunes around it.
#   * FitGmmModel thins exact zeros by dropping the FIRST n zero indices, so the fit depends on
#     row order.  The Python side sorts each region by object_id before writing it.
#   * FitGmmModel's set.seed(1) does NOT reach the fit: Rmixmod's mixmodCluster(seed = NULL)
#     seeds its C++ EM from its own source, and two identical calls under set.seed(1) return
#     different fits (verified 2026-09-05: likelihoods -825.0553 vs -825.0555 on one vector).
#     The package's NAMESPACE imports nothing, so its unqualified `mixmodCluster` resolves via
#     the global environment before the attached Rmixmod; the wrapper below supplies
#     `seed = <--seed>` and changes nothing else about the call.  The arm's seed therefore IS the
#     EM initialisation seed, and the screen's seed replicates measure exactly that.
#   * The pinned commit's NAMESPACE does not import `%<-%` from zeallot although every
#     constructor line uses it, so zeallot is attached here explicitly.
#   * No FilterCells (the artefact filter is off by default and the protocol says default).

suppressWarnings(suppressMessages({
  library(zeallot)   # `%<-%` -- see the NAMESPACE note above
  library(Rmixmod)
  library(spdep)
  library(CELESTA)
  library(jsonlite)
}))

args <- commandArgs(trailingOnly = TRUE)
scratch <- NULL
seed <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--scratch" && i < length(args)) { scratch <- args[i + 1]; i <- i + 2 }
  else if (args[i] == "--seed" && i < length(args)) { seed <- as.integer(args[i + 1]); i <- i + 2 }
  else { i <- i + 1 }
}
if (is.null(scratch) || !dir.exists(scratch)) stop("usage: Rscript celesta_run.R --scratch <existing dir> --seed <n>")
if (is.null(seed) || is.na(seed)) stop("--seed <integer> is required: Rmixmod's EM is otherwise seeded from its own clock")

# The seed wrapper (see the method note above).  Same arguments, one addition: seed.
mixmodCluster <- function(..., seed = NULL) Rmixmod::mixmodCluster(..., seed = getOption("celesta.rmixmod.seed"))
options(celesta.rmixmod.seed = seed)
in_dir <- file.path(scratch, "in")
out_dir <- file.path(scratch, "out")
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------------- env + defaults
desc <- read.dcf(system.file("DESCRIPTION", package = "CELESTA"))
one_of <- function(x) { u <- unique(as.numeric(x)); if (length(u) != 1) stop("threshold vector is not constant"); u }
f_assign <- formals(AssignCells)
f_create <- formals(CreateCelestaObject)
defaults <- list(
  # AssignCells -- evaluated from the package's own formals, never typed here
  max_iteration = eval(f_assign$max_iteration),
  cell_change_threshold = eval(f_assign$cell_change_threshold),
  min_diff = eval(f_assign$min_diff),
  min_probability = eval(f_assign$min_probability),
  high_expression_threshold_anchor = one_of(eval(f_assign$high_expression_threshold_anchor)),
  low_expression_threshold_anchor = one_of(eval(f_assign$low_expression_threshold_anchor)),
  high_expression_threshold_index = one_of(eval(f_assign$high_expression_threshold_index)),
  low_expression_threshold_index = one_of(eval(f_assign$low_expression_threshold_index)),
  scale_factor = eval(f_assign$scale_factor),
  bandwidth = eval(f_assign$bandwidth),
  # CreateCelestaObject
  cofactor_argument = eval(f_create$cofactor),
  cofactor_effective = 10,            # GetMarkerExpMatrix is called with the literal 10
  transform_type = eval(f_create$transform_type),
  number_of_neighbors = eval(f_create$number_of_neighbors),
  neighbor_bandwidth = eval(f_create$bandwidth)
)
env <- list(
  r_version = R.version.string,
  rscript = file.path(R.home("bin"), "Rscript"),
  celesta_version = as.character(packageVersion("CELESTA")),
  celesta_remote = paste0(desc[1, "RemoteUsername"], "/", desc[1, "RemoteRepo"]),
  celesta_sha = unname(desc[1, "RemoteSha"]),
  spdep_version = as.character(packageVersion("spdep")),
  rmixmod_version = as.character(packageVersion("Rmixmod")),
  zeallot_version = as.character(packageVersion("zeallot")),
  jsonlite_version = as.character(packageVersion("jsonlite")),
  zeallot_attached_explicitly = TRUE,
  rmixmod_seed = seed,
  rmixmod_seed_note = "supplied to mixmodCluster(seed=) by a global wrapper; the package's set.seed(1) does not reach Rmixmod",
  filter_cells = FALSE,
  defaults = defaults
)
write_json(env, file.path(scratch, "env.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

# ---------------------------------------------------------------------------- inputs
prior <- read.csv(file.path(scratch, "prior.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(prior) < 3) stop("prior.csv needs cell type, lineage code, and at least one marker column")
markers <- colnames(prior)[3:ncol(prior)]
for (m in markers) prior[[m]] <- as.integer(prior[[m]])   # 1 / 0 / NA; an all-NA column reads as logical

region_files <- sort(list.files(in_dir, pattern = "^region_.*\\.csv$", full.names = TRUE))
if (length(region_files) == 0) stop("no region_<r>.csv under ", in_dir)

log_con <- file(file.path(scratch, "celesta_stdout.log"), open = "wt")
write_status <- function(stem, status, ...) {
  write_json(c(list(region = stem, status = status), list(...)), file.path(out_dir, paste0(stem, ".status.json")),
             auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
}

# ---------------------------------------------------------------------------- per region
for (f in region_files) {
  stem <- sub("\\.csv$", "", basename(f))
  t0 <- Sys.time()
  cells <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(object_id = "character"))
  missing <- setdiff(c("object_id", "X", "Y", markers), colnames(cells))
  if (length(missing) > 0) {
    write_status(stem, "error", n_cells = nrow(cells), seconds = 0,
                 message = paste("region file lacks columns:", paste(missing, collapse = ",")))
    next
  }
  res <- tryCatch({
    sink(log_con, type = "output"); sink(log_con, type = "message")
    cat("\n==== ", stem, " (", nrow(cells), " cells) ====\n", sep = "")
    obj <- CreateCelestaObject(stem, prior, cells)                    # defaults: asinh(x/10), k = 5, bandwidth 100
    if (is.null(obj)) {
      list(status = "too_small", n_cells = nrow(cells), message = "CreateCelestaObject returned NULL (<= 10 cells)")
    } else {
      nan_markers <- colnames(obj@marker_exp_prob)[colSums(is.na(obj@marker_exp_prob)) > 0]
      if (length(nan_markers) > 0) {
        # A NaN activation column poisons every type score that references the marker (as 1 or
        # 0), which with shared negatives is every type.  Fail the region, name the markers.
        list(status = "activation_nan", n_cells = nrow(cells), markers = nan_markers,
             message = "sigmoid/GMM produced NaN activation probabilities")
      } else {
        obj <- AssignCells(obj, save_result = FALSE)                  # every threshold at the package default
        fin <- obj@final_cell_type_assignment
        label <- as.character(fin[, "Final cell type"])
        idx <- match(label, colnames(obj@current_cell_prob))          # NA for 'Unknown'
        prob <- rep(NA_real_, length(label))
        ok <- !is.na(idx)
        prob[ok] <- obj@current_cell_prob[cbind(which(ok), idx[ok])]
        anchor <- obj@anchor_cell_type_assignment[, 1] != 0
        out <- data.frame(object_id = cells$object_id, celesta_label = label, celesta_prob = prob,
                          celesta_anchor = anchor, stringsAsFactors = FALSE)
        write.csv(out, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)
        list(status = "ok", n_cells = nrow(cells), n_unknown = sum(label == "Unknown"), n_anchor = sum(anchor),
             counts = as.list(table(label)))
      }
    }
  }, error = function(e) list(status = "error", n_cells = nrow(cells), message = conditionMessage(e)),
  finally = { sink(type = "message"); sink(type = "output") })
  res$seconds <- as.numeric(Sys.time() - t0, units = "secs")
  do.call(write_status, c(list(stem = stem), res))
}
close(log_con)
cat("celesta_run.R: done\n")
