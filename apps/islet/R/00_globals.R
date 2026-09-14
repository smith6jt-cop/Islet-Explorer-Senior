# Debug modes are a deployment error in production. Keep this check before
# package loading so a misconfigured worker fails closed.
DEBUG_CREDS <- islet_flag_enabled(Sys.getenv("DEBUG_CREDENTIALS", "0"))
VIEWER_DEBUG_ENABLED <- islet_flag_enabled(Sys.getenv("VIEWER_DEBUG", "0"))
if (islet_is_production() && (DEBUG_CREDS || VIEWER_DEBUG_ENABLED)) {
  stop("Production does not permit credential or viewer debug modes.", call. = FALSE)
}

library(shiny)
library(shinyjs)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(broom)
library(jsonlite)

# ---- Scaling feature flags ------------------------------------------------
# Toggle these to fall back to the legacy in-memory code paths if something
# regresses. They are read by data_loading.R, spatial_helpers.R, and the
# stats/spatial modules.
USE_DUCKDB <- islet_flag_enabled(Sys.getenv("ISLET_USE_DUCKDB", "TRUE"))
USE_MIRAI  <- islet_flag_enabled(Sys.getenv("ISLET_USE_MIRAI", "TRUE"))
# AI Assistant panel: OFF by default so public deployments don't expose the UF
# Navigator API key (KEY/BASE). Set ISLET_ENABLE_AI=TRUE (e.g. local dev) to show it.
AI_ENABLED <- islet_flag_enabled(Sys.getenv("ISLET_ENABLE_AI", "FALSE"))

# Display geometry and grouping palettes come from the declared DOMAIN file.
# They are UI configuration, never cell-identity or threshold authorities.
PIXEL_SIZE_UM <- domain_pixel_size()
DONOR_PALETTES <- domain_group_palettes()
DONOR_COLORS <- domain_group_colors("Paul Tol (default)")
DONOR_COLORS_BRIGHT <- domain_group_colors("Bright")
DONOR_COLORS_COLORBLIND <- domain_group_colors("Okabe-Ito")
DONOR_COLORS_CLASSIC <- domain_group_colors("Diverging")

get_donor_color_palette <- function(name = "Bright") {
  DONOR_PALETTES[[name]] %||% DONOR_COLORS_BRIGHT
}

# Load sf package for GeoJSON/spatial operations (islet segmentation viewer)
SF_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (requireNamespace("sf", quietly = TRUE)) {
    library(sf)
    SF_AVAILABLE <- TRUE
    message("[SF] Package sf loaded for segmentation viewer")
  } else {
    message("[SF] Package 'sf' not found. Islet segmentation viewer will be disabled.")
    message("[SF] Install with: install.packages('sf')")
  }
})

# Check for httr2 package (required for AI assistant)
httr2_available <- requireNamespace("httr2", quietly = TRUE)
if (!httr2_available) {
  message("Package 'httr2' not found. The AI assistant will be disabled until you install it (install.packages('httr2')).")
}

# Load anndata package for trajectory analysis (install with BiocManager::install('anndata'))
suppressPackageStartupMessages({
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("Package 'anndata' not found. Install with: BiocManager::install('anndata')")
  }
})

# Load vitessceR for embedded image viewer (GitHub v0.1.0 required)
VITESSCE_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (!requireNamespace("vitessceR", quietly = TRUE)) {
    message("Package 'vitessceR' not found. Install with: remotes::install_github('vitessce/vitessceR')")
    message("Version required: >= 0.1.0")
  } else {
    # Verify we have the GitHub version with required features
    tryCatch({
      pkg_version <- packageVersion("vitessceR")
      if (pkg_version < "0.1.0") {
        warning("vitessceR version ", pkg_version, " detected. Please upgrade to >= 0.1.0 from GitHub")
      }
      test_vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.16", name = "Test")
      message("[VITESSCE] Package verified working (version ", pkg_version, ")")
      VITESSCE_AVAILABLE <<- TRUE
    }, error = function(e) {
      warning("[VITESSCE] Package installed but failed its capability check")
    })

    # CRITICAL FIX: Ensure all child routes from MultiImageWrapper are registered
    # This fixes segmentation layers not loading properly and appearing at wrong sizes
    if (VITESSCE_AVAILABLE) {
      tryCatch({
        vitessceR::MultiImageWrapper$set("public", "get_routes", function() {
          routes <- list()
          if (length(self$image_wrappers) == 0L) return(routes)
          for (w in self$image_wrappers) {
            wr <- w$get_routes()
            if (wr > 0L) routes <- c(routes, wr)
          }
          routes
        }, overwrite = TRUE)
        message("[VITESSCE] Applied MultiImageWrapper get_routes fix for segmentation layer routing")
      }, error = function(e) {
        warning("[VITESSCE] Failed to apply the segmentation routing fix")
      })
    }
  }
})

## for base64 encoding channel_config payload to Avivator
suppressPackageStartupMessages({
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Run scripts/install_shiny_deps.R to install dependencies.")
  }
})

# ---- Runtime diagnostics ----
APP_VERSION <- "modular-v1"
try({
  message("[APP LOAD] Version=", APP_VERSION, " time=", Sys.time())
  flag_path <- Sys.getenv("ISLET_APP_LOAD_FLAG", unset = "")
  if (nzchar(flag_path)) {
    temp_root <- normalizePath(tempdir(), mustWork = TRUE)
    flag_path <- normalizePath(flag_path, mustWork = FALSE)
    is_below_temp <- identical(flag_path, temp_root) ||
      startsWith(flag_path, paste0(temp_root, .Platform$file.sep))
    if (!is_below_temp) {
      warning("ISLET_APP_LOAD_FLAG must resolve below tempdir(); flag not written.")
    } else {
      if (!dir.exists(dirname(flag_path))) {
        dir.create(dirname(flag_path), recursive = TRUE, showWarnings = FALSE)
      }
      writeLines(paste0("loaded ", Sys.time(), " version=", APP_VERSION), flag_path)
    }
  }
}, silent = TRUE)

# Path constants. A DOMAIN may name its own data (see domain_data_root / domain_data_file in
# ieCore); a config without a `data:` section falls back to data/app_data, so the CODEX build
# is unaffected. This is what lets the PhenoCycler run cohort load as a SEPARATE dataset rather
# than being merged into the CODEX one -- the two share five donors of twenty and do not share
# a phenotype vocabulary, so a merged view would compare methods rather than biology.
master_path <- domain_data_file("master", "master_results.xlsx")
h5ad_path   <- domain_data_file("h5ad", "islet_explorer.h5ad")
parquet_dir <- file.path(domain_data_root(), "parquet")
project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# ---- DuckDB-backed data backend (Phase 1) ---------------------------------
# A single read-only DuckDB connection is created once at app startup and
# shared across all sessions. Parquet files produced by
# `scripts/convert_h5ad_to_parquet.py` are exposed as virtual views so each
# module can write lazy dbplyr pipelines (collect() only when an R-only
# function needs materialised data). When Parquet files are missing or the
# `USE_DUCKDB` flag is FALSE the app falls back to the existing
# H5AD/Excel in-memory loader.
con <- NULL
if (isTRUE(USE_DUCKDB)) {
  have_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
                 requireNamespace("DBI", quietly = TRUE)
  if (!have_duckdb) {
    message("[DUCKDB] Package 'duckdb' not installed; falling back to H5AD loader.")
    USE_DUCKDB <- FALSE
  } else {
    islets_parquet <- file.path(parquet_dir, "islets.parquet")
    tissue_root    <- file.path(parquet_dir, "tissue")
    if (!file.exists(islets_parquet)) {
      message("[DUCKDB] ", islets_parquet, " not found; run scripts/convert_h5ad_to_parquet.py. ",
              "Falling back to H5AD loader for this session.")
      USE_DUCKDB <- FALSE
    } else {
      con <- tryCatch(
        DBI::dbConnect(duckdb::duckdb(), read_only = TRUE),
        error = function(e) { message("event=duckdb_connect_failed"); NULL }
      )
      if (is.null(con)) {
        USE_DUCKDB <- FALSE
      } else {
        # Optional spatial extension — used for point-in-bbox click queries on
        # the tissue scatter. Swallow install errors so the app still boots
        # offline (extension lives in DuckDB's CDN by default).
        try(DBI::dbExecute(con, "INSTALL spatial; LOAD spatial;"), silent = TRUE)

        try({
          DBI::dbExecute(con, sprintf(
            "CREATE OR REPLACE VIEW islets AS SELECT * FROM parquet_scan('%s')",
            normalizePath(islets_parquet, winslash = "/", mustWork = TRUE)
          ))
          message("[DUCKDB] Registered view: islets -> ", islets_parquet)
        }, silent = FALSE)

        # Groovy tables — optional; only register if the converter emitted them.
        for (sheet in c("targets", "markers", "composition")) {
          pq <- file.path(parquet_dir, sprintf("uns_%s.parquet", sheet))
          if (file.exists(pq)) {
            try({
              DBI::dbExecute(con, sprintf(
                "CREATE OR REPLACE VIEW uns_%s AS SELECT * FROM parquet_scan('%s')",
                sheet, normalizePath(pq, winslash = "/", mustWork = TRUE)
              ))
            }, silent = TRUE)
          }
        }

        # Partitioned tissue dataset — register only if it exists so the
        # Spatial tab can fall back to CSVs otherwise.
        if (dir.exists(tissue_root)) {
          tissue_glob <- normalizePath(file.path(tissue_root, "**", "*.parquet"),
                                       winslash = "/", mustWork = FALSE)
          try({
            DBI::dbExecute(con, sprintf(
              "CREATE OR REPLACE VIEW tissue AS SELECT * FROM parquet_scan('%s', hive_partitioning=1)",
              tissue_glob
            ))
            message("[DUCKDB] Registered tissue view")
          }, silent = TRUE)
        }

        # Per-cell signed-distance dataset for the Spatial Card C right plot:
        # donor_status, cell_region, phenotype, dist_islet + 10 RESTORE `{m}_pos`
        # booleans, built by scripts/senior/build_cell_distance_parquet.py.
        # Partition-pruned; registered only if present so the app degrades
        # gracefully (the plot shows a "data missing" message) when absent.
        cell_distance_root <- file.path(parquet_dir, "cell_distance")
        if (dir.exists(cell_distance_root)) {
          try({
            DBI::dbExecute(con, sprintf(
              "CREATE OR REPLACE VIEW cell_distance AS SELECT * FROM parquet_scan('%s', hive_partitioning=1)",
              normalizePath(file.path(cell_distance_root, "**", "*.parquet"), winslash = "/", mustWork = FALSE)
            ))
            message("[DUCKDB] Registered cell-distance view")
          }, silent = TRUE)
        }

        # Ensure the connection is closed when the R process exits.
        # (shiny::onStop is only valid inside a running app, so register a
        # reg.finalizer on a harmless environment instead.)
        reg.finalizer(globalenv(), function(e) {
          try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
        }, onexit = TRUE)
      }
    }
  }
}

#' Helper: is the DuckDB backend active for this session?
duckdb_active <- function() isTRUE(USE_DUCKDB) && !is.null(con)

# ---- mirai worker pool for async tasks (Phase 4) --------------------------
# Long-running work (statistical batteries, heatmap resampling) is dispatched
# to a persistent mirai daemon pool so the UI stays responsive. Daemons are
# started once at app load and torn down on exit.
if (isTRUE(USE_MIRAI)) {
  have_mirai <- requireNamespace("mirai", quietly = TRUE)
  if (!have_mirai) {
    message("[MIRAI] Package 'mirai' not installed; async stats disabled.")
    USE_MIRAI <- FALSE
  } else {
    n_workers <- suppressWarnings(as.integer(Sys.getenv("ISLET_MIRAI_WORKERS", "4")))
    if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 4L
    try({
      mirai::daemons(n = n_workers, dispatcher = TRUE)
      message("[MIRAI] Started ", n_workers, " daemon workers")
    }, silent = TRUE)
    reg.finalizer(globalenv(), function(e) {
      try(mirai::daemons(0), silent = TRUE)
    }, onexit = TRUE)
  }
}

# Markers hidden from every app selector/plot (redundant epithelial keratins +
# nuclear stain + failed marker; mirrors scripts/senior/marker_taxonomy.py EXCLUDED).
# Pan_Cytokeratin is deliberately RETAINED (epithelial identity marker).
EXCLUDED_MARKERS <- c("DAPI", "IAPP", "Ker8_18", "Keratin_5")

# --- Make scaling globals visible to helpers that look them up in globalenv() ---
# Shiny auto-sources the R/ directory into a *support environment*, NOT the true
# globalenv(). Several helpers (data_loading.R, spatial_helpers.R,
# mod_statistics_server.R, rdeck_helpers.R) resolve these via
# get0(envir = globalenv()); without this mirror the running app silently falls
# back to the H5AD loader (and rdeck/mirai/partition-pruning never engage).
# Export the finalised values into the real global environment.
local({
  ge <- globalenv()
  for (nm in c("USE_DUCKDB", "USE_MIRAI", "con", "DONOR_COLORS", "EXCLUDED_MARKERS")) {
    if (exists(nm, inherits = TRUE)) assign(nm, get(nm, inherits = TRUE), envir = ge)
  }
})

# Ensure reticulate/anndata can discover a Python binary when RETICULATE_PYTHON
# isn't set (shiny-server environments often lack PATH).
#
# The interpreter MUST have the python `anndata` module, otherwise read_h5ad()
# throws and load_master_auto() silently falls back to the Excel source — which
# lacks the prop_*/peri_prop_* columns, gutting the Plot tab's cell-type
# selector. We therefore can't just grab the first python on PATH: the deploying
# user's venv (which ships anndata) may live under a different home than the one
# running the app. shiny-server runs as `shiny` (home /home/shiny), so the old
# path.expand("~") lookup found nothing and fell through to a bare system
# python3 with no anndata. Search several candidate locations and pick the first
# interpreter that can actually `import anndata`.
if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  .venv_rel <- file.path(".local", "share", "islet-explorer-py", "bin", "python")
  .py_candidates <- unique(c(
    file.path(path.expand("~"), .venv_rel),               # current user's venv
    Sys.glob(file.path("/home", "*", .venv_rel)),         # any user's venv (deployment)
    Sys.glob(file.path("/srv", "*", .venv_rel)),
    Sys.getenv("PYTHON", ""),
    Sys.which("python3"),
    Sys.which("python")
  ))
  .py_candidates <- .py_candidates[nzchar(.py_candidates) & file.exists(.py_candidates)]

  .py_has_anndata <- function(py) {
    status <- suppressWarnings(tryCatch(
      system2(py, c("-c", shQuote("import anndata")),
              stdout = FALSE, stderr = FALSE),
      error = function(e) 1L))
    identical(as.integer(status), 0L)
  }

  candidate_python <- ""
  for (.py in .py_candidates) {
    if (.py_has_anndata(.py)) { candidate_python <- .py; break }
  }
  # If none advertised anndata (e.g. the capability probe couldn't run), fall
  # back to the first existing candidate so behaviour is no worse than before.
  if (!nzchar(candidate_python) && length(.py_candidates)) {
    candidate_python <- .py_candidates[[1]]
    message("[PYTHON] WARNING: no candidate interpreter could import anndata; ",
            "H5AD loading may fail and fall back to Excel. Set RETICULATE_PYTHON ",
            "to a python that has the anndata module.")
  }

  if (nzchar(candidate_python)) {
    Sys.setenv(RETICULATE_PYTHON = candidate_python)
    message("[PYTHON] Compatible interpreter configured")
  } else {
    message("[PYTHON] No python executable found. Set RETICULATE_PYTHON to a valid interpreter.")
  }
}

# Feature flags
#
# Inflammation tab: HIDDEN by default (Aug 2026, at the user's request). The
# module and its docs are untouched — this only controls whether the tab is
# mounted, so re-enabling is a one-line env change rather than a code revert.
# Both the tabPanel and the inflammation_server() call in app.R are gated on it;
# leaving the server running without its UI would be harmless (Shiny reactives
# are lazy) but would still be dead wiring, so both are gated together.
#   SHOW_INFLAMMATION_TAB=1 in .Renviron  ->  tab reappears
SHOW_INFLAMMATION_TAB <- islet_flag_enabled(
  Sys.getenv("SHOW_INFLAMMATION_TAB", "0")
)

# ============================================================================
# The rail's destinations (Modernist shell)
# ============================================================================
# Each `key` is the VALUE the hidden tabset receives, and Senior's keys are the
# EXISTING tab labels rather than new lowercase ones. Junior's port renamed them
# and then had to chase every `active_tab() != "Trajectory"`, every
# `condition = "input.tabs == 'Plot'"` and viewer_server's own comparison. Here
# the navset keeps the id `tabs` and the same values, so the rail is a new way
# to reach the same tabs and nothing downstream of `input$tabs` changes.
#
# ORDER IS SENIOR'S OWN. Junior's redesign also reordered the tabs into its
# reading of the workflow and split Spatial in two. That is information
# architecture, not the look, and it was not what this port was asked for.
#
# The FIRST entry is what opens, and it is also what `navset_hidden` selects by
# default -- app.R's `nav_panel_hidden` order is kept in step so the two agree
# rather than one silently correcting the other on the first rail click.
#
# No cohort facts in `sub`: these strings are built once at startup, before any
# dataset is loaded, so a channel count or a donor count written here would be a
# claim about whichever domain happened to be configured first.
IE_PIPELINE <- c(
  list(
    list(key = "Plot",       num = "01", label = "Plot",
         sub = "Islet size \u00b7 composition"),
    list(key = "Trajectory", num = "02", label = "Trajectory",
         sub = "Pseudotime \u00b7 UMAP"),
    list(key = "Statistics", num = "03", label = "Statistics",
         sub = "Mixed-effects tests"),
    list(key = "Spatial",    num = "04", label = "Spatial",
         sub = "Tissue map \u00b7 neighbourhoods")
  ),
  # Gated exactly like the tabPanel it addresses, so the rail can never offer a
  # destination the tabset does not mount.
  if (SHOW_INFLAMMATION_TAB) list(
    list(key = "Inflammation", num = "05", label = "Inflammation",
         sub = "Cell types \u00b7 immune localisation")
  ) else list()
)
IE_REFERENCE <- list(
  list(key = "Viewer", label = "Viewer",
       sub = "The images at full resolution")
)

# How to read each tab. These used to be blue strips across the top of Plot and
# Trajectory, costing the first ~120px of the panel on every visit for text a
# reader needs once. The rail's Help entry shows the ACTIVE tab's tips on
# demand, so the cards start at the top of the panel.
#
# A tab with no entry gets a plain "nothing specific" note rather than an empty
# modal, so the Help row never looks broken.
IE_TAB_TIPS <- list(
  Plot = list(
    "Enable \"Show individual points\", then click any point on the Islet Size Distribution scatter to view that islet's segmentation map below.",
    "Enable \"Show individual points\" and set \"Color points by\" to Donor ID, then click a donor in the legend to show or hide it.",
    "Use browser zoom (Ctrl +/\u2212) to resize the cards. Scroll within a card to reach controls below the fold."
  ),
  Trajectory = list(
    "Click any point on the pseudotime scatter to view the segmentation map for that islet."
  ),
  Spatial = list(
    "Drag on the tissue scatter to zoom; double-click to reset."
  )
)

#' The Help > "Tips for this tab" modal, for whichever tab is active.
ie_tab_tips_modal <- function(tab) {
  tips <- IE_TAB_TIPS[[tab %||% ""]]
  shiny::modalDialog(
    title = paste0(tab %||% "This tab", " \u2014 how to read it"),
    easyClose = TRUE, footer = shiny::modalButton("Close"),
    if (length(tips)) {
      shiny::tags$ul(lapply(tips, function(t) shiny::tags$li(t)))
    } else {
      shiny::tags$p("Nothing specific to this tab. The phenotype rules are in the entry above.")
    }
  )
}
IE_TAB_KEYS <- c(vapply(IE_PIPELINE,  `[[`, "", "key"),
                 vapply(IE_REFERENCE, `[[`, "", "key"))

#' Keep N selectInputs and one reactiveVal in step, whichever the user touches.
#'
#' Replaces two hand-written observer sets in app.R -- a pair for the phenotype
#' palette and a triad for the donor palette -- that each restated the whole
#' fan-out per writer. Adding the rail's copy of each selector would have made
#' them a triad and a quartet of near-identical blocks, four places to forget
#' one id. The guard against a feedback loop is the same as before: a writer
#' that already agrees with the reactiveVal returns without updating anyone.
#'
#' @param ids character vector of NON-namespaced input ids. Every one is both a
#'   writer and a follower.
ie_sync_selects <- function(input, session, rv, ids) {
  for (this_id in ids) {
    local({
      writer <- this_id
      shiny::observeEvent(input[[writer]], {
        value <- input[[writer]]
        if (is.null(value) || identical(value, rv())) return()
        rv(value)
        for (other in setdiff(ids, writer)) {
          shiny::updateSelectInput(session, other, selected = value)
        }
      }, ignoreInit = TRUE)
    })
  }
  invisible(NULL)
}

# Robust click bridge for plotly scatters used by the islet drill-down.
# On macOS, plotly's own click/drag discrimination can suppress `plotly_click`
# entirely, so the server-side observer never fires (the symptom: clicking a
# point does nothing, no notification). This attaches a native mousedown/mouseup
# listener — which always fires regardless of platform — and, on a click-like
# gesture (movement < 6px), reads the hovered point's `key` ("case_id|islet_key",
# from the ggplot key aes) and pushes it to `input_id` via Shiny.setInputValue.
# The legacy plotly_click handler is kept too, so working platforms are
# unaffected. Pipe after the plotly build:
#   p %>% islet_click_bridge(session$ns("plot_native_click"))
islet_click_bridge <- function(p, input_id) {
  js <- sprintf("
function(el, x){
  // el is the plotly graph div; fall back to a child .js-plotly-plot if not.
  var gd = (el && typeof el.on === 'function') ? el
           : (el && el.querySelector ? el.querySelector('.js-plotly-plot') : el);
  if (!gd) gd = el;
  // onRender re-runs on every plot re-render; attach listeners only once per
  // graph div (plotly reuses the same node) so a click doesn't fire N times.
  if (gd.__isletBridge) return;
  gd.__isletBridge = true;
  var lastKey = null, downKey = null, downX = 0, downY = 0;
  function keyOf(pt){
    try {
      if (pt && pt.data && pt.data.key && pt.data.key[pt.pointNumber] != null)
        return pt.data.key[pt.pointNumber];
    } catch(e){}
    return null;
  }
  // plotly's live hover state — robust even if the plotly_hover event timing
  // is off, since at mousedown the cursor is over the target point.
  function hoverKey(){
    try { if (gd && gd._hoverdata && gd._hoverdata.length) return keyOf(gd._hoverdata[0]); } catch(e){}
    return null;
  }
  if (gd && typeof gd.on === 'function') {
    gd.on('plotly_hover', function(d){
      if (d && d.points && d.points.length){ var k = keyOf(d.points[0]); if (k != null) lastKey = k; }
    });
    gd.on('plotly_unhover', function(){ lastKey = null; });
  }
  (gd || el).addEventListener('mousedown', function(e){ downX = e.clientX; downY = e.clientY; downKey = lastKey || hoverKey(); }, true);
  (gd || el).addEventListener('mouseup', function(e){
    if (Math.abs(e.clientX - downX) > 5 || Math.abs(e.clientY - downY) > 5) return;
    var k = downKey || lastKey || hoverKey();
    if (k == null) return;
    Shiny.setInputValue('%s', {key: k, nonce: (Date.now() + Math.random())}, {priority: 'event'});
  }, true);
}", input_id)
  htmlwidgets::onRender(p, js)
}

# ============================================================================
# CSV export helper
# ============================================================================
# Every CSV the app emits should carry a `#`-commented provenance block naming
# the filters that produced it — otherwise two files with different row sets are
# indistinguishable after download. Writing to an open connection (rather than
# writeLines() + write.table(append = TRUE)) also avoids R's "appending column
# names to file" warning.
#
# @param file    Path or connection from a downloadHandler `content` function.
# @param df      Data frame to write. NULL/empty writes the header alone.
# @param header  Character vector of provenance lines (the "# " is added here).
write_islet_csv <- function(file, df, header = character(0)) {
  con <- base::file(file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  if (length(header) > 0) {
    writeLines(paste("#", header), con)
    writeLines("", con)
  }
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    writeLines("# No data for the current selection", con)
  } else {
    utils::write.csv(df, con, row.names = FALSE)
  }
  invisible(TRUE)
}

# Standard first lines of a provenance block. Callers append their own.
islet_csv_header <- function(title, ...) {
  c(title, paste("Generated:", Sys.time()), ...)
}

# Empty-state placeholder with the standard title font. `plotly_empty()` on its
# own renders its title at plotly's default size, so the nine "No data" messages
# across the app read noticeably smaller than everything around them.
plotly_islet_empty <- function(msg) {
  plotly::layout(plotly::plotly_empty(),
                 title = list(text = msg,
                              font = list(size = ISLET_FONT_SIZES$title)))
}

# ggplotly() silently DROPS `labs(subtitle = ...)`. Several plots put their
# p-value / correction annotations there, so those lines simply never rendered.
# Build a plotly title string that carries the subtitle as a `<sup>` second line
# and pass it to `layout(title = list(text = ...))`.
plotly_title_with_sub <- function(title, subtitle = NULL) {
  if (is.null(subtitle) || !nzchar(trimws(subtitle))) return(title)
  paste0(title, "<br><sup>", subtitle, "</sup>")
}

# Restyle the layout ANNOTATIONS that ggplotly() creates and that no `layout()`
# overlay can reach.
#
# For an unfaceted plot ggplotly keeps axis titles in `layout.xaxis.title`, so
# the usual `layout(xaxis = list(title = list(font = ...)))` works. Add a facet
# and it relocates BOTH the strip labels and the axis titles into
# `layout.annotations` — where that overlay is silently ignored and they render
# at ggplot's converted default (~14.6px). Verified with plotly_build().
#
# Strip labels have no `annotationType`; axis titles carry `"axis"`.
plotly_islet_strip_fonts <- function(gg,
                                     strip_size = ISLET_FONT_SIZES$strip,
                                     axis_title_size = ISLET_FONT_SIZES$axis_title) {
  ann <- gg$x$layout$annotations
  if (is.null(ann) || length(ann) == 0) return(gg)
  for (i in seq_along(ann)) {
    gg$x$layout$annotations[[i]]$font$size <-
      if (identical(ann[[i]]$annotationType, "axis")) axis_title_size else strip_size
  }
  gg
}

# Standard PNG download dimensions for ggplot ggsave() in downloadHandlers.
# Every such call also passes bg = "white" ON PURPOSE. theme_modernist() draws on a
# transparent background so figures sit directly on the Modernist card; a transparent
# PNG then inverts on a dark slide, so the saved file deliberately differs from the screen.
# Default: 10x7 in at 300 dpi -> 3000x2100 px. Heatmaps use wider 12x8.
ISLET_PNG_DIMS  <- list(width = 10, height = 7, dpi = 300, units = "in")
ISLET_HEATMAP_PNG_DIMS <- list(width = 12, height = 8, dpi = 300, units = "in")

islet_png_filename <- function(stem) {
  paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
}

# ============================================================================
# Phenotype rules (May 2026) — marker -> celltype definitions surfaced in UI
# ============================================================================
# Source of truth: data/phenotype_rules.csv (parsed once at startup).
# Column-1 is the parent rule, column-2 the phenotype, columns 3+ are markers
# with cell values:
#   "pos"    -> required positive
#   "anypos" -> at least one of the marker(s) in the same `anypos` group
#               (the parser groups all `anypos` cells per phenotype together)
#   "allneg" -> required negative
#
# PHENOTYPE_COLORS (drilldown_helpers.R) uses 21 names; a few don't have rows
# in the CSV (residuals from the upstream single-cell pipeline). Those get
# hand-written descriptions in PHENOTYPE_RESIDUALS below. The alias map
# normalises PHENOTYPE_COLORS keys -> CSV keys for the four known mismatches.

PHENOTYPE_RULES_PATH <- file.path(domain_data_root(), "phenotype_rules.csv")

# PHENOTYPE_COLORS-form -> CSV-form (and back via reverse lookup). Used to
# resolve a UI label to its rule row.
PHENOTYPE_ALIASES <- domain_phenotype_aliases()

# Phenotypes with no CSV row — residuals from the upstream pipeline.
# Each entry: list(parent, pos, neg, anypos_groups, description).
# This is a CATALOGUE, not the active set: .PHENOTYPE_RESIDUAL_CATALOGUE is filtered below to
# the phenotypes the active domain actually has. Unfiltered, the run cohort's rules modal listed
# four CODEX phenotypes (ECAD+, APCs, Structural, Unknown) that its typing never produces.
.PHENOTYPE_RESIDUAL_CATALOGUE <- list(
  `ECAD+` = list(
    parent = "Epithelial", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Epithelial cells (E-cadherin+); residual ductal/exocrine bucket from the upstream phenotyping pipeline. Not gated by the CSV rules."
  ),
  APCs = list(
    parent = "Immune", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Antigen-presenting cells aggregate (HLADR+ myeloid lineage). Used as an aggregate immune category in neighbourhood metrics; not a CSV-rule phenotype."
  ),
  Structural = list(
    parent = "all", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Residual structural cells negative for all endocrine/immune markers; bucket for fibroblast/stroma-like cells."
  ),
  Unknown = list(
    parent = "all", pos = character(), neg = character(),
    anypos_groups = list(),
    description = "Default fallback for cells that match no rule above."
  )
)

# The active set: the catalogue entries this domain's vocabulary contains, plus anything the
# domain config declares itself. A residual the domain never displays is not a residual here.
PHENOTYPE_RESIDUALS <- local({
  vocabulary <- tryCatch(domain_phenotype_order(), error = function(e) character())
  kept <- .PHENOTYPE_RESIDUAL_CATALOGUE[names(.PHENOTYPE_RESIDUAL_CATALOGUE) %in% vocabulary]
  declared <- tryCatch(domain_phenotype_residuals(), error = function(e) list())
  for (entry in declared) {
    name <- entry[["name"]]
    if (is.null(name) || name %in% names(kept)) next
    kept[[name]] <- list(parent = entry[["parent"]] %||% "all", pos = character(),
                         neg = character(), anypos_groups = list(),
                         description = entry[["description"]] %||% NA_character_)
  }
  kept
})

# Parse the CSV into a named list of rule records keyed by phenotype name.
# Returns NULL if the file is missing — callers should guard.
.parse_phenotype_rules_csv <- function(path) {
  if (!file.exists(path)) {
    message("[PHENOTYPE] Rules unavailable; modal and hints disabled")
    return(NULL)
  }
  raw <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"), check.names = FALSE)
  # Row 1 is the header; columns 3..N are marker names. An optional trailing column
  # headed `description` carries the rule in prose — needed where the call is not a flat
  # marker gate (the run cohort's classification asks a different question per context and
  # compares against a reference level, which "CD3e+, CD68+" would misstate as a plain OR).
  header <- as.character(raw[1, ])
  desc_col <- which(tolower(trimws(header)) %in% c("description", "notes"))
  desc_col <- if (length(desc_col)) desc_col[[1]] else NA_integer_
  marker_cols <- setdiff(seq.int(3, ncol(raw)), desc_col)
  marker_cols <- marker_cols[!is.na(header[marker_cols]) & nzchar(header[marker_cols])]
  markers <- header[marker_cols]
  data_rows <- raw[-1, , drop = FALSE]
  rules <- list()
  for (i in seq_len(nrow(data_rows))) {
    pheno <- as.character(data_rows[i, 2])
    if (is.na(pheno) || !nzchar(pheno)) next
    parent <- as.character(data_rows[i, 1])
    parent <- if (is.na(parent) || !nzchar(parent)) "all" else parent
    description <- if (is.na(desc_col)) NA_character_
                   else as.character(data_rows[i, desc_col])
    if (!is.na(description) && !nzchar(description)) description <- NA_character_
    vals <- as.character(data_rows[i, marker_cols])
    names(vals) <- markers
    pos <- names(vals)[!is.na(vals) & vals == "pos"]
    neg <- names(vals)[!is.na(vals) & vals == "allneg"]
    anypos_markers <- names(vals)[!is.na(vals) & vals == "anypos"]
    # CSV puts all `anypos` markers on one row as a single group.
    anypos_groups <- if (length(anypos_markers) > 0) list(anypos_markers) else list()
    rules[[pheno]] <- list(
      parent        = parent,
      pos           = pos,
      neg           = neg,
      anypos_groups = anypos_groups,
      description   = description
    )
  }
  rules
}

PHENOTYPE_RULES_CSV_RAW <- .parse_phenotype_rules_csv(PHENOTYPE_RULES_PATH)

# The markers that file actually declares, for the modal footer. Empty when no file.
PHENOTYPE_RULES_MARKERS <- if (is.null(PHENOTYPE_RULES_CSV_RAW)) {
  character()
} else {
  unique(unlist(lapply(PHENOTYPE_RULES_CSV_RAW,
                       function(r) c(r$pos, r$neg, unlist(r$anypos_groups)))))
}

# Compose hierarchical rules: a sub-rule (e.g. "CD8 T cell") inherits the
# positive/negative requirements of its parent chain ("T cell" -> "Immune").
# We walk the parent chain until we hit "all" (root) or a missing parent.
.compose_phenotype_rules <- function(raw) {
  if (is.null(raw)) return(NULL)
  composed <- list()
  compose <- function(name, seen = character()) {
    if (!is.null(composed[[name]])) return(composed[[name]])
    if (name %in% seen) return(raw[[name]])  # cycle guard
    if (!name %in% names(raw)) return(NULL)
    r <- raw[[name]]
    if (!is.null(r$parent) && nzchar(r$parent) && r$parent != "all"
        && r$parent %in% names(raw)) {
      p <- compose(r$parent, c(seen, name))
      if (!is.null(p)) {
        r$pos           <- unique(c(p$pos, r$pos))
        r$neg           <- unique(c(p$neg, r$neg))
        r$anypos_groups <- c(p$anypos_groups, r$anypos_groups)
      }
    }
    # Cleanup: drop anypos groups already satisfied by a required positive.
    r$anypos_groups <- Filter(function(grp) length(intersect(grp, r$pos)) == 0,
                              r$anypos_groups)
    # Deduplicate identical groups.
    if (length(r$anypos_groups) > 1) {
      keys <- vapply(r$anypos_groups, function(g) paste(sort(g), collapse = "|"), character(1))
      r$anypos_groups <- r$anypos_groups[!duplicated(keys)]
    }
    # Drop any group that is a strict SUPERSET of another group in the same rule
    # (the smaller group is the stricter, sufficient constraint).
    if (length(r$anypos_groups) > 1) {
      keep <- rep(TRUE, length(r$anypos_groups))
      for (i in seq_along(r$anypos_groups)) {
        for (j in seq_along(r$anypos_groups)) {
          if (i == j) next
          a <- r$anypos_groups[[i]]; b <- r$anypos_groups[[j]]
          if (length(b) < length(a) && all(b %in% a)) { keep[i] <- FALSE; break }
        }
      }
      r$anypos_groups <- r$anypos_groups[keep]
    }
    r$neg <- setdiff(r$neg, r$pos)
    composed[[name]] <<- r
    r
  }
  for (nm in names(raw)) compose(nm)
  composed
}

PHENOTYPE_RULES_CSV <- .compose_phenotype_rules(PHENOTYPE_RULES_CSV_RAW)

# Resolve a phenotype name (in either CSV or PHENOTYPE_COLORS form) to its rule.
# Returns NULL if the name has no rule and no residual entry.
get_phenotype_rule <- function(name) {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(NULL)
  # Try residuals first (small set, includes "Unknown" fallback).
  if (name %in% names(PHENOTYPE_RESIDUALS)) return(PHENOTYPE_RESIDUALS[[name]])
  # Apply alias (PHENOTYPE_COLORS -> CSV).
  csv_name <- if (name %in% names(PHENOTYPE_ALIASES)) unname(PHENOTYPE_ALIASES[[name]]) else name
  if (!is.null(PHENOTYPE_RULES_CSV) && csv_name %in% names(PHENOTYPE_RULES_CSV)) {
    return(PHENOTYPE_RULES_CSV[[csv_name]])
  }
  NULL
}

# Format the rule for display. Styles:
#   "inline" -> plain text, e.g. "CD3e+, CD8a+  —  CD4−, CD20−"
#   "html"   -> HTML with strong positives and muted negatives
#   "plain"  -> minimal text, safe for plot subtitles
format_phenotype_rule <- function(name, style = c("inline", "html", "plain")) {
  style <- match.arg(style)
  rule <- get_phenotype_rule(name)
  if (is.null(rule)) {
    msg <- paste0("No marker rule recorded for “", name, "”.")
    if (style == "html") return(htmltools::HTML(paste0("<em>", msg, "</em>")))
    return(msg)
  }
  if (!is.na(rule$description) && nzchar(rule$description)) {
    # Residual phenotype — just return the description.
    if (style == "html") return(htmltools::HTML(paste0("<em>", htmltools::htmlEscape(rule$description), "</em>")))
    return(rule$description)
  }
  # Build positive-side string
  pos_parts <- character()
  for (m in rule$pos) pos_parts <- c(pos_parts, paste0(m, "+"))
  for (grp in rule$anypos_groups) {
    if (length(grp) == 1) {
      pos_parts <- c(pos_parts, paste0(grp, "+"))
    } else {
      pos_parts <- c(pos_parts, paste0("any of ", paste(paste0(grp, "+"), collapse = " / ")))
    }
  }
  neg_parts <- if (length(rule$neg) > 0) paste0(rule$neg, "−") else character()
  if (style == "html") {
    pos_html <- if (length(pos_parts) > 0)
      paste0("<strong>", paste(htmltools::htmlEscape(pos_parts), collapse = ", "), "</strong>")
    else "<span style='color:#888;'>(no positive requirement)</span>"
    neg_html <- if (length(neg_parts) > 0)
      paste0("<span style='color:#666;'>", paste(htmltools::htmlEscape(neg_parts), collapse = ", "), "</span>")
    else ""
    sep <- if (nzchar(neg_html)) "  —  " else ""
    return(htmltools::HTML(paste0(pos_html, sep, neg_html)))
  }
  pos_txt <- if (length(pos_parts) > 0) paste(pos_parts, collapse = ", ") else "(no positive requirement)"
  neg_txt <- if (length(neg_parts) > 0) paste(neg_parts, collapse = ", ") else ""
  if (nzchar(neg_txt)) paste0(pos_txt, "  —  ", neg_txt) else pos_txt
}

# Fallback order for the reference modal — endocrine first, then immune (T cell subtypes
# grouped), then structural/endothelial, residuals last. The ACTIVE domain's declared
# `phenotypes.display_order` takes precedence; this list only covers names it omits, and is
# the whole order for a domain that declares none. Unqualified, this CODEX-shaped literal
# shares four spellings with the run cohort and sorted its other four alphabetically, so the
# modal opened on Exocrine with the endocrine classes scattered below.
PHENOTYPE_DISPLAY_ORDER <- c(
  "Beta cell", "Alpha cell", "Delta cell", "Endocrine",
  "Acinar", "Ductal", "Exocrine",
  "T cell", "CD8 T cell", "CD4 T cell", "B cell",
  "Macrophage", "Dendritic cell", "Immune", "APCs",
  "Lymphatic", "Blood vessel", "Endothelial",
  "Smooth Muscle", "Neural", "ECAD+", "Structural",
  "Unknown"
)

# Build a tidy data.frame for the modal table.
# Returns columns: Phenotype, Category, Positive markers, Negative markers, Notes.
phenotype_rules_table_df <- function() {
  csv_names <- if (is.null(PHENOTYPE_RULES_CSV)) character() else names(PHENOTYPE_RULES_CSV)
  all_names <- c(csv_names, names(PHENOTYPE_RESIDUALS))
  all_names <- unique(all_names)
  # Order: known order first, then anything else alphabetical.
  preferred <- tryCatch(domain_phenotype_order(), error = function(e) character())
  preferred <- c(preferred, setdiff(PHENOTYPE_DISPLAY_ORDER, preferred))
  ordered <- c(intersect(preferred, all_names), sort(setdiff(all_names, preferred)))
  rows <- lapply(ordered, function(nm) {
    rule <- get_phenotype_rule(nm)
    if (is.null(rule)) {
      return(data.frame(Phenotype = nm, Category = "?", Positive = "",
                        Negative = "", Notes = "",
                        stringsAsFactors = FALSE))
    }
    pos_parts <- character()
    for (m in rule$pos) pos_parts <- c(pos_parts, paste0(m, "+"))
    for (grp in rule$anypos_groups) {
      pos_parts <- c(pos_parts,
                     if (length(grp) == 1) paste0(grp, "+")
                     else paste0("any of ", paste(paste0(grp, "+"), collapse = " / ")))
    }
    neg_parts <- if (length(rule$neg) > 0) paste0(rule$neg, "−") else character()
    data.frame(
      Phenotype = nm,
      Category  = rule$parent,
      Positive  = if (length(pos_parts) > 0) paste(pos_parts, collapse = ", ") else "—",
      Negative  = if (length(neg_parts) > 0) paste(neg_parts, collapse = ", ") else "—",
      Notes     = if (!is.na(rule$description)) rule$description else "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Modal builder. Pure shiny — no DT dependency. Pass to showModal() from any tab's
# Phenotype-Rules button handler.
phenotype_rules_modal_ui <- function() {
  df <- phenotype_rules_table_df()
  header_row <- shiny::tags$tr(
    shiny::tags$th("Phenotype",          style = "text-align:left; padding:6px 10px"),
    shiny::tags$th("Category",           style = "text-align:left; padding:6px 10px"),
    shiny::tags$th("Positive markers",   style = "text-align:left; padding:6px 10px"),
    shiny::tags$th("Negative markers",   style = "text-align:left; padding:6px 10px"),
    shiny::tags$th("Notes",              style = "text-align:left; padding:6px 10px")
  )
  body_rows <- lapply(seq_len(nrow(df)), function(i) {
    shiny::tags$tr(
      shiny::tags$td(shiny::tags$strong(df$Phenotype[i]), style = "padding:6px 10px"),
      shiny::tags$td(df$Category[i],                       style = "padding:6px 10px"),
      shiny::tags$td(df$Positive[i],                       style = "padding:6px 10px"),
      shiny::tags$td(df$Negative[i],                       style = "padding:6px 10px"),
      shiny::tags$td(df$Notes[i],                          style = "padding:6px 10px; font-size: 13px")
    )
  })
  shiny::modalDialog(
    title = shiny::tagList(shiny::icon("dna"), "Cell Phenotype Rules — Marker Definitions"),
    size = "l", easyClose = TRUE, fade = TRUE,
    footer = shiny::modalButton("Close"),
    shiny::tags$p(
      "Every cell in the app is classified into one of these phenotypes by the upstream ",
      shiny::tags$code(PHENOTYPE_RULES_PATH), " rules. A ",
      shiny::tags$code("+"),
      " means the marker is required to be positive; a ",
      shiny::tags$code("−"),
      " means it must be negative. “any of X+ / Y+” means at least one of the listed markers must be positive.",
      " Where a phenotype's Notes column is filled, the call is not a flat marker gate and the note",
      " is the rule — read it rather than the marker columns."
    ),
    shiny::tags$div(style = "max-height: 60vh; overflow-y: auto",
      shiny::tags$table(
        class = "table table-striped",
        style = "width: 100%; font-size: 14px; margin: 0",
        shiny::tags$thead(header_row),
        shiny::tags$tbody(body_rows)
      )
    ),
    shiny::tags$p(style = "margin-top: 10px; font-size: 12px",
      # The marker list is read from the file that was actually parsed: the two domains carry
      # different panels, and a hard-coded list names the wrong markers for one of them.
      "Rule source: ", shiny::tags$code(PHENOTYPE_RULES_PATH),
      if (length(PHENOTYPE_RULES_MARKERS))
        paste0(" (", length(PHENOTYPE_RULES_MARKERS), " markers: ",
               paste(PHENOTYPE_RULES_MARKERS, collapse = ", "), ").")
      else "."
    )
  )
}

# Resolve a sanitised data-column name (e.g. "prop_CD8a_Tcell", "peri_prop_Beta_cell",
# "min_dist_Macrophage", "enrich_z_B_cell") back to a phenotype label. Returns
# NULL/the input unchanged if no known prefix is recognised.
extract_phenotype_from_column <- function(col) {
  if (is.null(col) || is.na(col) || !nzchar(col)) return(NULL)
  for (pre in c("combined_prop_", "peri_prop_", "peri_count_",
                "enrich_z_", "min_dist_", "prop_")) {
    if (startsWith(col, pre)) {
      stem <- sub(paste0("^", pre), "", col)
      return(gsub("_", " ", gsub("plus", "+", stem, fixed = TRUE), fixed = TRUE))
    }
  }
  col
}

# Human-readable label for a raw data-column / selection encoding, for use in
# plot axes, card titles and CSV headers.
#
# The Statistics tab took its feature name straight from the Plot module's
# `selection_used` attribute, so raw encodings leaked onto the page as if they
# were the feature's name: "prop_Beta cell", "peri_prop_CD8a_Tcell",
# "Capillary|islet_core". This resolves the two encodings the app actually uses:
#   - `<prefix>_<phenotype>`  (prop_ / peri_prop_ / combined_prop_ / ...)
#   - `<class>|<region>`      (Microenvironment targets)
# Anything unrecognised is returned unchanged, so this is always safe to apply.
feature_display_label <- function(col) {
  if (is.null(col) || length(col) != 1 || is.na(col) || !nzchar(col)) return(col)

  # Microenvironment encoding: "Capillary|islet_core"
  if (grepl("|", col, fixed = TRUE)) {
    parts  <- strsplit(col, "|", fixed = TRUE)[[1]]
    region <- switch(parts[2],
                     islet_core  = "islet core",
                     islet_band  = "peri-islet",
                     islet_union = "core+peri",
                     gsub("_", " ", parts[2], fixed = TRUE))
    return(paste0(parts[1], " (", region, ")"))
  }

  region <- if (startsWith(col, "combined_prop_")) "core+peri"
            else if (startsWith(col, "peri_prop_") || startsWith(col, "peri_count_")) "peri-islet"
            else if (startsWith(col, "prop_") || startsWith(col, "count_")) "core"
            else NULL

  measure <- if (grepl("^(combined_prop_|peri_prop_|prop_)", col)) "% of cells"
             else if (grepl("^(peri_count_|count_)", col)) "cells"
             else if (startsWith(col, "enrich_z_")) "enrichment z"
             else if (startsWith(col, "min_dist_")) "min distance, um"
             else NULL

  pheno <- extract_phenotype_from_column(col)
  if (is.null(pheno) || identical(pheno, col)) return(col)   # no prefix matched

  qual <- paste(c(region, measure), collapse = ", ")
  if (nzchar(qual)) paste0(pheno, " (", qual, ")") else pheno
}

# UI button that opens the Phenotype Rules modal. The caller chooses the input
# id (one per tab to avoid duplicate-DOM-id collisions across simultaneously
# rendered tabPanels) and registers a matching `observeEvent` that calls
# `showModal(phenotype_rules_modal_ui())`. Style mirrors the existing "About
# this tab" info banner.
phenotype_rules_button <- function(input_id,
                                   label = "Phenotype rules — marker definitions") {
  shiny::tags$div(
    style = "margin-bottom: 12px; padding: 8px 14px; display: inline-block",
    shiny::actionLink(input_id,
      shiny::tagList(shiny::icon("dna"),
                     shiny::tags$span(style = "margin-left:6px; font-weight:600",
                                      label)),
      style = "text-decoration: none")
  )
}

# UI helper: small muted hint block describing the currently selected phenotype.
# Pass a literal phenotype name. The output is a tagList suitable for embedding
# directly under a selectInput. `show_modal_link` controls the trailing
# "see all rules" actionLink; set to NULL to omit.
phenotype_hint_ui <- function(name, prefix = "Markers:",
                              show_modal_link = "show_phenotype_rules_inline") {
  if (is.null(name) || is.na(name) || !nzchar(name)) return(NULL)
  rule_html <- format_phenotype_rule(name, "html")
  link <- if (!is.null(show_modal_link)) {
    shiny::tagList(
      shiny::tags$span(" · "),
      shiny::actionLink(show_modal_link, "see all rules →",
                        style = "font-size: 12px")
    )
  }
  shiny::tags$div(
    style = "margin-top: -6px; margin-bottom: 8px; padding: 6px 10px; font-size: 13px; line-height: 1.4",
    shiny::tags$strong(paste0(prefix, " ")),
    shiny::tags$span(rule_html),
    link
  )
}
