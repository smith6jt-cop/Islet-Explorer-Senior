#!/usr/bin/env Rscript
# ============================================================================
# Statistics-tab audit harness  (regression net for the Jul 2026 correctness pass)
# ============================================================================
#
# WHY THIS EXISTS
# ---------------
# The six Statistics defects fixed in `37ee4b7` were all invisible to code
# review: silent no-ops, size-dependent method switches, and cross-module
# coupling. They were found by driving the REAL modules headlessly with
# `shiny::testServer` against the live H5AD and measuring. That harness was
# throwaway, so the audit could not be repeated without rewriting it -- which is
# what this file fixes. Each check below pins one fixed defect.
#
# It is NOT in `tests/testthat/` because it needs the real `data/` artifacts
# (gitignored) plus `sf` / `reticulate` / `anndata`. CI runs the parse check and
# the data-independent unit tests; this is the manual/local layer.
#
#   cd apps/islet && Rscript tests/manual/stats_audit.R
#
# Exits 0 if every check passes, 1 otherwise. Skips cleanly (exit 0) when the
# data or packages are absent, so it is safe to wire into a broader script.
#
# THE TWO-STAGE SHAPE (the part that is easy to get wrong)
# --------------------------------------------------------
# `statistics_server` consumes `plot_returns$raw_df_prezero`, so ONE testServer
# call cannot drive both modules. Stage 1 drives `plot_server` and grabs
# `session$returned$raw_df_prezero()` as a plain data frame; stage 2 feeds that
# frame to `statistics_server` as a constant reactive. `stats_run()` is an
# `eventReactive(ignoreInit = TRUE)`, so it must be fired with
# `session$setInputs(run_tests = <increasing integer>)`.
# ============================================================================

suppressWarnings(suppressMessages({
  ok_pkgs <- all(vapply(c("shiny", "dplyr", "lmerTest", "emmeans"),
                        requireNamespace, logical(1), quietly = TRUE))
}))
if (!ok_pkgs) {
  cat("SKIP: required packages missing (shiny/dplyr/lmerTest/emmeans)\n")
  quit(save = "no", status = 0)
}

# Run from apps/islet so the "../../data" paths in 00_globals.R resolve.
if (!file.exists("app.R") && file.exists("apps/islet/app.R")) setwd("apps/islet")
if (!file.exists("app.R")) {
  cat("SKIP: run me from apps/islet (app.R not found)\n")
  quit(save = "no", status = 0)
}

# reticulate needs RETICULATE_PYTHON from .Renviron, or load_master_auto()
# silently falls back to Excel and the H5AD-only columns vanish.
if (file.exists(".Renviron")) readRenviron(".Renviron")

suppressWarnings(suppressMessages(library(shiny)))
for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  suppressWarnings(suppressMessages(source(f)))
}

if (!file.exists(h5ad_path) && !file.exists(master_path)) {
  cat("SKIP: no data/ artifacts present (H5AD and Excel both missing)\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Tiny assertion helpers (no testthat: this runs against live data, and a
# readable transcript matters more than a test-runner UI).
# ---------------------------------------------------------------------------
.fails <- 0L
chk <- function(label, cond, detail = "") {
  passed <- isTRUE(cond)
  if (!passed) .fails <<- .fails + 1L
  cat(sprintf("[%s] %s%s\n", if (passed) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  invisible(passed)
}
section <- function(x) cat("\n== ", x, " ", strrep("=", max(0, 60 - nchar(x))), "\n", sep = "")

# ---------------------------------------------------------------------------
# Load data once (~40 s). Cached to an RDS so repeat runs are fast.
# ---------------------------------------------------------------------------
# NOT tempdir(): that is per-R-session, so the cache never survived to a second
# run and every invocation paid the full ~40 s load. Stable path, overridable.
cache <- Sys.getenv("ISLET_AUDIT_CACHE",
                    file.path(dirname(tempdir()), "islet_prepared_cache.rds"))
PD <- if (file.exists(cache)) readRDS(cache) else {
  cat("Loading data (first run is slow; cached afterwards)...\n")
  p <- prep_data(load_master_auto(h5ad_path, master_path))
  saveRDS(p, cache); p
}
cat(sprintf("Data: %d comp rows, %d targets rows\n\n",
            nrow(PD$comp %||% data.frame()), nrow(PD$targets_all %||% data.frame())))

# ---------------------------------------------------------------------------
# Stage 1 -- drive plot_server, return its pre-zero per-islet frame.
# ---------------------------------------------------------------------------
#
# TRAP -- `aab_flags` MUST be set explicitly.
# `raw_df_base()` drops EVERY Aab+ islet when `input$aab_flags` is NULL or empty
# (`aabp <- aabp[0, ]`). The real UI ships a checkboxGroupInput with all five
# selected, so this never bites a user -- but testServer starts every unset
# input at NULL, so a harness that omits it silently tests a two-group cohort
# (ND + T1D) while appearing to request three. This is the same class of trap
# the audit itself was about: a hidden control still filtering upstream data.
ALL_AAB <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")

get_plot_frame <- function(mode, which, min_cells = 1) {
  out <- NULL
  testServer(plot_server,
             args = list(id = "plot", prepared = reactive(PD),
                         selected_islet = reactiveVal(NULL),
                         active_tab = reactive("Plot")), {
    session$setInputs(mode = mode, which = which, metric = "Percentage",
                      groups = c("ND", "Aab+", "T1D"), min_cells = min_cells,
                      aab_flags = ALL_AAB,
                      show_points = TRUE, exclude_zero_top = FALSE)
    out <<- tryCatch(as.data.frame(session$returned$raw_df_prezero()),
                     error = function(e) NULL)
  })
  out
}

# ---------------------------------------------------------------------------
# Stage 2 -- feed that frame to statistics_server and read stats_run().
# `extra` sets any additional stats_* inputs.
# ---------------------------------------------------------------------------
run_stats <- function(frame, extra = list(), want_html = NULL) {
  res <- NULL
  testServer(statistics_server,
             args = list(id = "stats", raw_df = reactive(frame),
                         get_selection_description = function() character(0)), {
    inputs <- utils::modifyList(
      list(stats_min_cells = 1, stats_exclude_zero = FALSE,
           stats_remove_outliers = FALSE, stats_outlier_z = 3.0),
      extra)
    do.call(session$setInputs, inputs)
    session$setInputs(run_tests = 1)
    st <- tryCatch(stats_run(), error = function(e) NULL)
    # renderUI/renderText outputs come back as character vectors or lists
    # depending on the renderer; flatten to one string so callers can grep it.
    html <- if (!is.null(want_html)) {
      tryCatch(paste(unlist(output[[want_html]]), collapse = " "),
               error = function(e) NA_character_)
    } else NA_character_
    res <<- list(st = st, html = html)
  })
  res
}

# ===========================================================================
section("Defect 1 -- 'Min cells/islet' must not be a silent no-op")
# The filter is guarded on `total_cells_core`/`total_cells_peri`. When
# prep_data() did not carry them onto targets_all, it SKIPPED instead of
# erroring: Capillary held at 5,380 islets and p = 0.3476 from min-cells 1
# through 200 while the CSV still reported the value set. Sweeping the knob and
# asserting the row count MOVES is the only reliable detection.
# ===========================================================================
micro <- get_plot_frame("Microenvironment", "Capillary|islet_core")
# Sanity gate FIRST: if the cohort silently lost a disease group, every check
# below is testing something weaker than it claims (see the aab_flags trap).
chk("all three donor groups reached the Statistics frame",
    !is.null(micro) && length(intersect(c("ND", "Aab+", "T1D"),
                                        unique(as.character(micro$donor_status)))) == 3,
    paste("present:", paste(sort(unique(as.character(micro$donor_status))), collapse = ", ")))
if (is.null(micro) || !nrow(micro)) {
  chk("Microenvironment frame available", FALSE, "plot_server returned nothing")
} else {
  chk("targets_all carries per-islet cell counts",
      all(c("total_cells_core", "total_cells_peri") %in% names(micro)),
      paste("cols:", paste(intersect(c("total_cells_core", "total_cells_peri"),
                                     names(micro)), collapse = ", ")))
  n_by_thresh <- vapply(c(1, 25, 100), function(mc) {
    st <- run_stats(micro, list(stats_min_cells = mc))$st
    if (is.null(st) || is.null(st$rdf)) NA_integer_ else nrow(st$rdf)
  }, integer(1))
  cat(sprintf("   islets at min_cells 1/25/100: %s\n",
              paste(n_by_thresh, collapse = " / ")))
  chk("min-cells actually filters in Microenvironment mode",
      !anyNA(n_by_thresh) && n_by_thresh[1] > n_by_thresh[3],
      "row count must strictly decrease as the threshold rises")
}

# ===========================================================================
section("Defect 2 -- p-values must not always be scientific")
# formatC(format='e') rendered p = 0.3476 as '3.48e-01', one misread order of
# magnitude from 0.0348. format_p_value() is the single rule: fixed to 4 dp
# down to 1e-4, scientific below.
# ===========================================================================
chk("format_p_value: mid-range p is fixed-point",
    format_p_value(0.3476) == "0.3476", format_p_value(0.3476))
chk("format_p_value: tiny p goes scientific",
    grepl("e-", format_p_value(1.6e-15), fixed = TRUE), format_p_value(1.6e-15))
chk("format_p_value: NA is not printed as a number",
    format_p_value(NA_real_) == "N/A", format_p_value(NA_real_))
if (!is.null(micro) && nrow(micro)) {
  h <- run_stats(micro, want_html = "primary_headline")$html
  h1 <- if (length(h) == 1) h else paste(h, collapse = " ")
  chk("rendered headline carries no bare 'e-01' style p",
      length(h1) != 1 || is.na(h1) || !grepl("[0-9]e-0[12]([^0-9]|$)", h1),
      "guards the exact regression that produced '3.48e-01'")
}

# ===========================================================================
section("Defect 3 -- emmeans must not silently switch d.f. method at n>3000")
# Above pbkrtest.limit/lmerTest.limit (3000) emmeans drops to asymptotic z
# (df = Inf), which is badly anti-conservative for a 15-donor design:
# p(ND-T1D) = 0.00007 at 5,212 islets vs 0.0027 at 2,874. The fix requests
# Satterthwaite explicitly and RECORDS the method, so a regression is visible.
# ===========================================================================
cells <- get_plot_frame("Cell Populations", "prop_Beta cell")
if (is.null(cells) || !nrow(cells)) {
  chk("Cell Populations frame available", FALSE)
} else {
  st <- run_stats(cells)$st
  chk("run produced a result", !is.null(st))
  if (!is.null(st)) {
    cat(sprintf("   n_islets = %s, d.f. method = %s\n",
                nrow(st$rdf), st$lmer_pairwise_df_method %||% "<unset>"))
    chk("d.f. method is recorded at all",
        !is.null(st$lmer_pairwise_df_method))
    chk("d.f. method is not asymptotic despite n > 3000",
        !identical(tolower(st$lmer_pairwise_df_method %||% ""), "asymptotic"),
        "must stay Satterthwaite above pbkrtest.limit")
    if (!is.null(st$lmer_pairwise) && nrow(st$lmer_pairwise)) {
      dfs <- suppressWarnings(as.numeric(st$lmer_pairwise$df))
      chk("pairwise contrasts have finite d.f. (not Inf = z-test)",
          all(is.finite(dfs)), paste("df:", paste(round(dfs, 1), collapse = ", ")))
    }
  }
}

# ===========================================================================
section("Defect 4 -- the autoantibody test must be donor-level")
# AAb count is a DONOR attribute, but the test ran across islets: p = 1.6e-15
# from 1,666 islets belonging to just 5 donors split 1/3/1. Donor-level: 0.67.
# ===========================================================================
if (!is.null(cells) && nrow(cells)) {
  st <- run_stats(cells)$st
  aab <- st$aab_kw
  if (is.null(aab)) {
    cat("   (no AAb result on this selection -- skipped)\n")
  } else {
    chk("AAb test reports donors per group, not islets",
        !is.null(aab$donors_per_group))
    chk("AAb test carries an underpowered flag",
        "underpowered" %in% names(aab))
    chk("AAb p is not an islet-level pseudoreplicated value",
        is.na(aab$p_value) || aab$p_value > 1e-10,
        sprintf("p = %s (islet-level gave 1.6e-15)", format_p_value(aab$p_value)))
  }
}

# ===========================================================================
section("Defect 5 -- AUC must obey THIS tab's filters, not the Plot sidebar's")
# AUC integrated the Plot module's summary_df(), so it ignored all four Advanced
# filters (identical at 5,212 vs 131 islets) and obeyed the Plot bin-width and
# Mean/Median controls. `summary_df` is now absent from the module signature.
# ===========================================================================
chk("summary_df is not a statistics_server argument",
    !("summary_df" %in% names(formals(statistics_server))),
    paste("formals:", paste(names(formals(statistics_server)), collapse = ", ")))
if (!is.null(cells) && nrow(cells)) {
  a_all <- run_stats(cells)$st
  a_cut <- run_stats(cells, list(stats_min_cells = 100))$st
  if (!is.null(a_all$auc_by_group) && !is.null(a_cut$auc_by_group)) {
    same <- isTRUE(all.equal(a_all$auc_by_group, a_cut$auc_by_group))
    chk("AUC changes when a Statistics filter changes", !same,
        "byte-identical AUC across filters was the original defect")
  }
}

# ===========================================================================
section("Defect 6 -- a cleared numericInput (NA) must not abort the run")
# `%||%` catches NULL but not NA; `if (NA > 1 && TRUE)` aborted the whole
# eventReactive with 'missing value where TRUE/FALSE needed', blanking the tab.
# ===========================================================================
if (!is.null(cells) && nrow(cells)) {
  for (bad in list(NA, NA_integer_, NULL)) {
    lbl <- if (is.null(bad)) "NULL" else "NA"
    st <- tryCatch(run_stats(cells, list(stats_min_cells = bad))$st,
                   error = function(e) structure(list(), err = conditionMessage(e)))
    chk(sprintf("min-cells = %s does not abort the run", lbl),
        !is.null(st) && is.null(attr(st, "err")),
        attr(st, "err") %||% "")
  }
}

# ===========================================================================
section("Regression guards for the 2026-08-14 pass")
# ===========================================================================
chk("feature_display_label resolves the prop_ encoding",
    feature_display_label("prop_Beta cell") == "Beta cell (core, % of cells)",
    feature_display_label("prop_Beta cell"))
chk("feature_display_label resolves the class|region encoding",
    feature_display_label("Capillary|islet_core") == "Capillary (islet core)",
    feature_display_label("Capillary|islet_core"))
chk("feature_display_label passes unknown names through",
    feature_display_label("Ins_any") == "Ins_any")
if (!is.null(micro) && nrow(micro)) {
  st <- run_stats(micro)$st
  chk("feature name shown to the user is not a raw column",
      !is.null(st) && !grepl("\\|", st$feature_name),
      st$feature_name %||% "<none>")
  chk("raw column is still retained for provenance",
      !is.null(st$feature_col) && grepl("\\|", st$feature_col),
      st$feature_col %||% "<none>")
  chk("KW effect size is labelled epsilon_sq, not eta_sq",
      is.null(st$eta_sq_label) || st$eta_sq_label %in% c("eta_sq", "epsilon_sq"),
      st$eta_sq_label %||% "<unset>")
}

# ---------------------------------------------------------------------------
cat(sprintf("\n%s\n%d check(s) failed.\n", strrep("=", 66), .fails))
quit(save = "no", status = if (.fails == 0L) 0L else 1L)
