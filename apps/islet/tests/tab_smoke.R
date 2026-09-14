#!/usr/bin/env Rscript
# Load every tab in a real browser and report what each one actually rendered.
#
# The testthat suite checks contracts -- which paths resolve, which columns must exist. It cannot
# catch a tab that resolves everything and then renders an empty panel, which is the failure this
# dataset kept producing: Spatial's cards returned NULL when a column was missing, the Trajectory
# tab showed nothing when its H5AD was absent, and the drill-down fell back to boundaries in
# silence. Those are only visible in a browser.
#
# Reports per tab: the R errors Shiny surfaced, the JavaScript console errors, how many plot and
# table outputs are present, and how many of them are EMPTY -- an output that exists but drew
# nothing is the thing worth knowing.
#
#   NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml Rscript tests/tab_smoke.R [outdir]
#
# NOT_CRAN is required: AppDriver$new() calls skip_on_cran() internally and, without it,
# exits with "Reason: On CRAN" -- which reads like a failure and is really a refusal to run.

suppressMessages(library(shinytest2))

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) out_dir <- file.path(tempdir(), "tab_smoke")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TABS <- c("Plot", "Trajectory", "Spatial", "Statistics", "Viewer")

cat("starting the app (data load can take a minute)\n"); flush.console()
# A tall viewport rather than a full-page capture. `get_screenshot(selector = "html")` sounds
# like the fix for a sidebar that scrolls past the fold, but it captures the html element's own
# width and clips the right-hand cards instead -- trading a vertical crop for a horizontal one.
# Sizing the window is what actually gets more of the page into one image.
app <- AppDriver$new(".", name = "tab-smoke", seed = 20260914,
                     load_timeout = 300000, timeout = 120000,
                     height = 2000, width = 1600,
                     options = list(shiny.autoreload = FALSE))
on.exit(try(app$stop(), silent = TRUE), add = TRUE)

report <- list()
for (tab in TABS) {
  cat(sprintf("\n== %s\n", tab)); flush.console()
  ok <- TRUE
  tryCatch({
    app$set_inputs(tabs = tab, timeout_ = 120000)
    app$wait_for_idle(duration = 4000, timeout = 120000)
  }, error = function(e) {
    cat("  navigation failed:", conditionMessage(e), "\n"); ok <<- FALSE
  })

  html <- tryCatch(app$get_html("body"), error = function(e) "")
  # Shiny renders an R error into .shiny-output-error with the message as text.
  r_errors <- regmatches(html, gregexpr("shiny-output-error[^>]*>[^<]*", html))[[1]]
  js <- tryCatch(app$get_logs(), error = function(e) NULL)
  js_errors <- if (!is.null(js)) {
    lines <- as.data.frame(js)
    sum(grepl("error|SEVERE", paste(lines$level, lines$message), ignore.case = TRUE))
  } else NA_integer_

  n_plots <- lengths(regmatches(html, gregexpr("class=\"[^\"]*shiny-plot-output", html)))
  n_plotly <- lengths(regmatches(html, gregexpr("plotly html-widget", html)))
  n_empty <- lengths(regmatches(html, gregexpr("shiny-html-output[^>]*>\\s*</div>", html)))

  file <- file.path(out_dir, sprintf("%s.png", tolower(tab)))
  tryCatch(app$get_screenshot(file), error = function(e) cat("  screenshot failed\n"))

  cat(sprintf("  R output errors %d | JS errors %s | plot outputs %d | plotly %d | empty html outputs %d\n",
              length(r_errors), js_errors, n_plots, n_plotly, n_empty))
  if (length(r_errors)) for (e in head(r_errors, 4)) cat("    ", substr(e, 1, 160), "\n")
  report[[tab]] <- list(ok = ok, r_errors = length(r_errors), js_errors = js_errors,
                        plots = n_plots, plotly = n_plotly, empty = n_empty)
}

cat("\n=== summary ===\n")
failed <- names(Filter(function(r) !r$ok || r$r_errors > 0, report))
for (tab in names(report)) {
  r <- report[[tab]]
  cat(sprintf("  %-12s %s  (R errors %d, plots %d, plotly %d)\n", tab,
              if (!r$ok || r$r_errors > 0) "FAIL" else "ok", r$r_errors, r$plots, r$plotly))
}
cat(sprintf("\nscreenshots -> %s\n", out_dir))
if (length(failed)) quit(status = 1)
