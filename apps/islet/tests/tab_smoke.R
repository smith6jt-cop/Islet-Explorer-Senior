#!/usr/bin/env Rscript
# Drive every tab in Chrome and inspect its visible DOM, not hidden outputs from
# other tabs. For registry cards, require a loaded plot/table or a named guard
# note. A heading or empty output is not proof that the card's server ran.
#
# NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml \
#   Rscript tests/tab_smoke.R /tmp/tabshots
#
# Saves tab/card screenshots, tab HTML, browser/app logs and report.json.
# NOT_CRAN is required by AppDriver$new(); otherwise it skips without running.

suppressMessages(library(shinytest2))

TABS <- c("Plot", "Trajectory", "Spatial", "Neighborhoods", "Statistics",
          "CrossModal", "LymphNode", "Viewer")

# Consume the actual registries, so a newly registered card is checked too.
registries <- new.env(parent = asNamespace("shiny"))
sys.source("R/mod_crossmodal_ui.R", envir = registries)
sys.source("R/mod_lymphnode_ui.R", envir = registries)
CARD_IDS <- setNames(list(
  vapply(registries$IE_CROSSMODAL_CARDS, `[[`, "", "id"),
  vapply(registries$IE_LYMPHNODE_CARDS, `[[`, "", "id")
), c(registries$IE_CROSSMODAL_TAB, registries$IE_LYMPHNODE_TAB))
NOT_BUILT <- c("ledger", "n1_sections")
stopifnot(all(names(CARD_IDS) %in% TABS))

run_smoke <- function(out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  cat("starting the app (data load can take a minute)\n"); flush.console()
  log_cursor <- Sys.time()
  app <- AppDriver$new(".", name = "tab-smoke", seed = 20260914,
                       load_timeout = 300000, timeout = 120000,
                       height = 2000, width = 1600,
                       options = list(shiny.autoreload = FALSE))
  on.exit(try(app$stop(), silent = TRUE), add = TRUE)
  # Keep startup diagnostics even if navigation or DOM inspection throws later.
  on.exit(try(writeLines(capture.output(print(app$get_logs())),
                         file.path(out_dir, "app.log")), silent = TRUE), add = TRUE, after = FALSE)

  app$run_js("window.tabSmokeInspect = function(tab) {
    const visible = el => el.getClientRects().length > 0 &&
      getComputedStyle(el).visibility !== 'hidden';
    const query = (root, selector) => Array.from(root.querySelectorAll(selector)).filter(visible);
    const pane = Array.from(document.querySelectorAll('.tab-pane')).find(el =>
      el.dataset.value === tab && el.classList.contains('active') && visible(el));
    if (!pane) return {visible: false, cards: [], errors: ['Requested tab is not visible'],
      plots: 0, drawn_plots: 0, plotly: 0, drawn_plotly: 0, tables: 0, drawn_tables: 0, empty_html: 0};
    const drawing = root => {
      const plots = query(root, '.shiny-plot-output');
      const plotly = query(root, '.plotly.html-widget');
      const tables = query(root, '.shiny-table-output, .datatables.html-widget');
      const plotDrawn = el => Array.from(el.querySelectorAll('img')).some(img =>
        visible(img) && img.complete && img.naturalWidth > 0 && img.naturalHeight > 0);
      const plotlyDrawn = el => !!el.querySelector('.main-svg, canvas');
      const tableDrawn = el => query(el, 'tbody tr').length > 0;
      return {
        plots: plots.length, drawn_plots: plots.filter(plotDrawn).length,
        plotly: plotly.length, drawn_plotly: plotly.filter(plotlyDrawn).length,
        tables: tables.length, drawn_tables: tables.filter(tableDrawn).length,
        empty_html: query(root, '.shiny-html-output').filter(el => !el.innerHTML.trim()).length
      };
    };
    const cards = query(pane, '.ie-card-stack > .ie-card').map(el => {
      const id = el.dataset.cardId || '';
      const notes = query(el, '.ie-card-note').map(n => n.innerText.trim()).filter(text =>
        text.includes('[' + id + ']') || text.includes('[guard_' + id + ']') ||
        text.includes('[' + id + '_') ||
        text.includes('unavailable(not_built): no card_' + id + '_ui()'));
      return Object.assign({id: id, title: el.querySelector('.ie-card-title')?.innerText || '',
        placeholder: el.classList.contains('ie-card-placeholder'), notes: notes}, drawing(el));
    });
    return Object.assign({visible: true, cards: cards, html: pane.outerHTML,
      errors: query(pane, '.shiny-output-error').map(el =>
        ({id: el.id, message: el.innerText.trim(), validation: el.classList.contains('shiny-output-error-validation')}))
    }, drawing(pane));
  };")

  report <- list()
  for (tab in TABS) {
    cat(sprintf("\n== %s\n", tab)); flush.console()
    ok <- TRUE
    tryCatch({
      # Selecting the initial tab, or a static placeholder tab, may produce no
      # output-change event. Wait for Shiny/DOM rendering explicitly below.
      app$set_inputs(tabs = tab, wait_ = FALSE)
      app$wait_for_idle(duration = 4000, timeout = 120000)
      app$wait_for_js("Array.from(document.querySelectorAll('.tab-pane.active .shiny-plot-output img')).every(img => img.complete && img.naturalWidth > 0)")
    }, error = function(e) {
      cat("  navigation/render wait failed:", conditionMessage(e), "\n"); ok <<- FALSE
    })

    dom <- app$get_js(sprintf("window.tabSmokeInspect(%s)", jsonlite::toJSON(tab, auto_unbox = TRUE)))
    if (!isTRUE(dom$visible)) ok <- FALSE
    if (!is.null(dom$html)) writeLines(dom$html, file.path(out_dir, paste0(tolower(tab), ".html")))
    dom$html <- NULL
    logs <- as.data.frame(app$get_logs())
    # Browser errors since the last tab (including startup on the first tab).
    # Matching 'error' in arbitrary Shiny/debug messages produced false counts.
    browser_errors <- logs[!is.na(logs$timestamp) & logs$timestamp >= log_cursor &
                             logs$location == "chromote" &
                             grepl("error|throw|SEVERE", logs$level, ignore.case = TRUE), , drop = FALSE]
    log_cursor <- Sys.time()
    js_errors <- nrow(browser_errors)
    writeLines(capture.output(print(browser_errors)),
               file.path(out_dir, paste0(tolower(tab), "-browser.log")))

    screenshot <- function(name, selector = NULL) {
      tryCatch({
        file <- file.path(out_dir, paste0(name, ".png"))
        # AppDriver refuses to overwrite an existing file; a repeat run should work.
        if (file.exists(file)) unlink(file)
        if (is.null(selector)) app$get_screenshot(file) else app$get_screenshot(file, selector = selector)
      }, error = function(e) {
        cat("  screenshot failed:", name, conditionMessage(e), "\n"); ok <<- FALSE
      })
    }
    screenshot(tolower(tab))

    expected <- CARD_IDS[[tab]]
    if (length(expected)) {
      actual <- vapply(dom$cards, `[[`, "", "id")
      if (!identical(actual, expected)) {
        cat("  FAIL registry cards: expected", paste(expected, collapse = ", "),
            "but saw", paste(actual, collapse = ", "), "\n")
        ok <- FALSE
      }
      for (card in dom$cards) {
        drawn <- card$drawn_plots + card$drawn_plotly + card$drawn_tables
        empty_drawings <- card$plots + card$plotly + card$tables - drawn
        named_note <- length(card$notes) > 0L
        valid <- nzchar(card$id) && nzchar(card$title) && empty_drawings == 0L &&
          if (card$id %in% NOT_BUILT) {
            isTRUE(card$placeholder) && named_note &&
              any(grepl("unavailable(not_built)", unlist(card$notes), fixed = TRUE))
          } else {
            !isTRUE(card$placeholder) && (drawn > 0L || named_note)
          }
        if (!valid) ok <- FALSE
        cat(sprintf("  card %-22s %s | plots %d | plotly %d | tables %d | named notes %d\n",
                    card$id, if (valid) "ok" else "FAIL", card$drawn_plots,
                    card$drawn_plotly, card$drawn_tables, length(card$notes)))
        for (note in card$notes) cat("    ", note, "\n", sep = "")
        if (nzchar(card$id)) screenshot(paste(tolower(tab), card$id, sep = "-"),
          sprintf('.tab-pane.active [data-card-id="%s"]', card$id))
      }
    }
    r_errors <- length(dom$errors)
    if (r_errors > 0L || js_errors > 0L) ok <- FALSE
    cat(sprintf(paste0("  R output errors %d | JS errors %d | plot outputs %d (drawn %d) | ",
                       "plotly %d (drawn %d) | tables %d (drawn %d) | empty html outputs %d\n"),
                r_errors, js_errors, dom$plots, dom$drawn_plots, dom$plotly,
                dom$drawn_plotly, dom$tables, dom$drawn_tables, dom$empty_html))
    for (e in dom$errors) cat("    ", jsonlite::toJSON(e, auto_unbox = TRUE), "\n", sep = "")
    report[[tab]] <- c(list(ok = ok, r_errors = r_errors, js_errors = js_errors), dom)
    jsonlite::write_json(report, file.path(out_dir, "report.json"), pretty = TRUE,
                         auto_unbox = TRUE, null = "null")
    flush.console()
  }

  cat("\n=== summary ===\n")
  for (tab in names(report)) {
    r <- report[[tab]]
    cat(sprintf("  %-14s %s  (R errors %d, JS errors %d, drawn plots %d, plotly %d, tables %d)\n",
                tab, if (r$ok) "ok" else "FAIL", r$r_errors, r$js_errors,
                r$drawn_plots, r$drawn_plotly, r$drawn_tables))
  }
  cat(sprintf("\nscreenshots and report -> %s\n", out_dir))
  report
}

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) out_dir <- file.path(tempdir(), "tab_smoke")
report <- run_smoke(out_dir)
if (any(!vapply(report, `[[`, logical(1), "ok"))) quit(status = 1)
