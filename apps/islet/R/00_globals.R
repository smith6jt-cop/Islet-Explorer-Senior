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
USE_DUCKDB <- isTRUE(as.logical(Sys.getenv("ISLET_USE_DUCKDB", "TRUE")))
USE_MIRAI  <- isTRUE(as.logical(Sys.getenv("ISLET_USE_MIRAI",  "TRUE")))
# AI Assistant panel: OFF by default so public deployments don't expose the UF
# Navigator API key (KEY/BASE). Set ISLET_ENABLE_AI=TRUE (e.g. local dev) to show it.
AI_ENABLED <- isTRUE(as.logical(Sys.getenv("ISLET_ENABLE_AI", "FALSE")))

# Pixel size constant for coordinate conversion (micrometers per pixel)
# GeoJSON polygons use pixel coordinates, while islet_spatial_lookup.csv uses micrometers
# The QuPath structure GeoJSONs (results/*__<donor>_..._resolution__1.geojson) are
# exported at a ~2x-downsampled pyramid level: empirically 0.5078 µm/px (median over
# 2,355 islets, 3 donors; full-res base is ~0.2539). islet_spatial_lookup.csv is in µm,
# so this constant maps lookup µm -> GeoJSON pixel space. (Junior CODEX was 0.3774.)
PIXEL_SIZE_UM <- 0.5078  # micrometers per pixel (IO60panc2nd GeoJSON resolution #1)

# Donor status color palettes
DONOR_COLORS <- c("ND" = "#4477AA", "Aab+" = "#CC6633", "T1D" = "#228833")

DONOR_COLORS_BRIGHT <- c("ND" = "#3366CC", "Aab+" = "#FF6600", "T1D" = "#109618")

DONOR_COLORS_COLORBLIND <- c("ND" = "#0072B2", "Aab+" = "#D55E00", "T1D" = "#009E73")

DONOR_COLORS_CLASSIC <- c("ND" = "#2166AC", "Aab+" = "#B2182B", "T1D" = "#1B7837")

DONOR_PALETTES <- list(
  "Paul Tol (default)" = DONOR_COLORS,
  "Bright"             = DONOR_COLORS_BRIGHT,
  "Okabe-Ito"          = DONOR_COLORS_COLORBLIND,
  "Diverging"          = DONOR_COLORS_CLASSIC
)

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
      warning("[VITESSCE] Package installed but failed basic test: ", conditionMessage(e))
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
        warning("[VITESSCE] Failed to apply get_routes fix: ", conditionMessage(e))
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
  message("[APP LOAD] getwd=", tryCatch(getwd(), error = function(e) NA))
  message("[APP LOAD] list.files(.) contains app.R? ", any(grepl('^app.R$', list.files('.'))))
  flag_path <- file.path(dirname(sys.frame(1)$ofile %||% 'app.R'), 'APP_LOADED_FLAG')
  writeLines(paste0('loaded ', Sys.time(), ' version=', APP_VERSION), flag_path)
}, silent = TRUE)

# Path constants
master_path <- file.path("..", "..", "data", "app_data", "master_results.xlsx")
h5ad_path   <- file.path("..", "..", "data", "app_data", "islet_explorer.h5ad")
parquet_dir <- file.path("..", "..", "data", "app_data", "parquet")
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
        error = function(e) { message("[DUCKDB] connect failed: ", conditionMessage(e)); NULL }
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
            message("[DUCKDB] Registered view: tissue -> ", tissue_root)
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
            message("[DUCKDB] Registered view: cell_distance -> ", cell_distance_root)
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
    message("[PYTHON] RETICULATE_PYTHON=", candidate_python)
  } else {
    message("[PYTHON] No python executable found. Set RETICULATE_PYTHON to a valid interpreter.")
  }
}

# Debug flags
DEBUG_CREDS <- isTRUE(as.logical(Sys.getenv("DEBUG_CREDENTIALS", "0")))
VIEWER_DEBUG_ENABLED <- identical(Sys.getenv("VIEWER_DEBUG", "0"), "1")

# ============================================================================
# Plot font + theme helpers (May 2026)
# ============================================================================
# Centralised font sizing so every plot in the app reads clearly on screen and
# downloads as a presentation/publication-ready PNG. "Presentation-first"
# baseline: base 16pt, 18pt titles, 14pt legends, 12pt minimum for dense rows.

ISLET_FONT_SIZES <- list(
  base         = 16,
  title        = 18,
  subtitle     = 14,
  axis_title   = 16,
  axis_text    = 14,
  legend_title = 16,
  legend_text  = 14,
  strip        = 14,
  dense_min    = 12   # heatmap row labels, dense tick labels
)

# Standard ggplot theme. Use directly: `theme_islet()` or extend with
# additional `theme(...)` calls for plot-specific overrides.
theme_islet <- function(base_size = ISLET_FONT_SIZES$base) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$title, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ISLET_FONT_SIZES$subtitle, color = "#555"),
      axis.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_title),
      axis.text     = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_text),
      legend.title  = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_title, face = "bold"),
      legend.text   = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_text),
      strip.text    = ggplot2::element_text(size = ISLET_FONT_SIZES$strip, face = "bold")
    )
}

# Named list of plotly `list(size = N)` font specs ready to splice into
# layout() calls. Keys: global, title, axis_title, axis_tick, legend_title,
# legend_text, colorbar_title.
plotly_islet_fonts <- function() {
  list(
    global         = list(size = ISLET_FONT_SIZES$base),
    title          = list(size = ISLET_FONT_SIZES$title),
    axis_title     = list(size = ISLET_FONT_SIZES$axis_title),
    axis_tick      = list(size = ISLET_FONT_SIZES$axis_text),
    legend_title   = list(size = ISLET_FONT_SIZES$legend_title),
    legend_text    = list(size = ISLET_FONT_SIZES$legend_text),
    colorbar_title = list(size = ISLET_FONT_SIZES$legend_title)
  )
}

# Apply publication-quality PNG download config to a plotly object.
# Resulting PNG is `width * scale` x `height * scale` px (default 2800x1800).
# Use as `p %>% plotly_islet_config("filename_root")`.
plotly_islet_config <- function(p, filename_root = "islet_explorer_plot",
                                width = 1400, height = 900, scale = 2) {
  plotly::config(p,
    toImageButtonOptions = list(
      format   = "png",
      filename = filename_root,
      width    = width,
      height   = height,
      scale    = scale
    ),
    displaylogo = FALSE
  )
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
    gd.on('plotly_click', function(){ console.log('[ISLET] plotly_click fired'); });
  }
  (gd || el).addEventListener('mousedown', function(e){ downX = e.clientX; downY = e.clientY; downKey = lastKey || hoverKey(); }, true);
  (gd || el).addEventListener('mouseup', function(e){
    if (Math.abs(e.clientX - downX) > 5 || Math.abs(e.clientY - downY) > 5) return;
    var k = downKey || lastKey || hoverKey();
    if (k == null) return;
    console.log('[ISLET] native click -> ' + k);
    Shiny.setInputValue('%s', {key: k, nonce: (Date.now() + Math.random())}, {priority: 'event'});
  }, true);
}", input_id)
  htmlwidgets::onRender(p, js)
}

# Standard PNG download dimensions for ggplot ggsave() in downloadHandlers.
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

PHENOTYPE_RULES_PATH <- file.path("..", "..", "data", "app_data", "phenotype_rules.csv")

# PHENOTYPE_COLORS-form -> CSV-form (and back via reverse lookup). Used to
# resolve a UI label to its rule row.
PHENOTYPE_ALIASES <- c(
  "CD8a Tcell"   = "CD8 T cell",
  "CD4 Tcell"    = "CD4 T cell",
  "Blood Vessel" = "Blood vessel",
  "SMA+"         = "Smooth Muscle"
)

# Phenotypes with no CSV row — residuals from the upstream pipeline.
# Each entry: list(parent, pos, neg, anypos_groups, description).
PHENOTYPE_RESIDUALS <- list(
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

# Parse the CSV into a named list of rule records keyed by phenotype name.
# Returns NULL if the file is missing — callers should guard.
.parse_phenotype_rules_csv <- function(path) {
  if (!file.exists(path)) {
    message("[PHENOTYPE] Rules CSV not found at ", path, " — modal/hints disabled")
    return(NULL)
  }
  raw <- utils::read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"), check.names = FALSE)
  # Row 1 is the header; columns 3..N are marker names.
  markers <- as.character(raw[1, 3:ncol(raw)])
  markers <- markers[!is.na(markers) & nzchar(markers)]
  data_rows <- raw[-1, , drop = FALSE]
  rules <- list()
  for (i in seq_len(nrow(data_rows))) {
    pheno <- as.character(data_rows[i, 2])
    if (is.na(pheno) || !nzchar(pheno)) next
    parent <- as.character(data_rows[i, 1])
    parent <- if (is.na(parent) || !nzchar(parent)) "all" else parent
    vals <- as.character(data_rows[i, 3:(2 + length(markers))])
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
      description   = NA_character_
    )
  }
  rules
}

PHENOTYPE_RULES_CSV_RAW <- .parse_phenotype_rules_csv(PHENOTYPE_RULES_PATH)

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

# Order in which phenotypes are shown in the reference modal — endocrine first,
# then immune (T cell subtypes grouped), then structural/endothelial, residuals last.
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
  ordered <- c(intersect(PHENOTYPE_DISPLAY_ORDER, all_names),
               sort(setdiff(all_names, PHENOTYPE_DISPLAY_ORDER)))
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
    shiny::tags$th("Phenotype",          style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Category",           style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Positive markers",   style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Negative markers",   style = "text-align:left; padding:6px 10px;"),
    shiny::tags$th("Notes",              style = "text-align:left; padding:6px 10px;")
  )
  body_rows <- lapply(seq_len(nrow(df)), function(i) {
    shiny::tags$tr(
      shiny::tags$td(shiny::tags$strong(df$Phenotype[i]), style = "padding:6px 10px;"),
      shiny::tags$td(df$Category[i],                       style = "padding:6px 10px; color:#555;"),
      shiny::tags$td(df$Positive[i],                       style = "padding:6px 10px;"),
      shiny::tags$td(df$Negative[i],                       style = "padding:6px 10px; color:#666;"),
      shiny::tags$td(df$Notes[i],                          style = "padding:6px 10px; color:#444; font-size: 13px;")
    )
  })
  shiny::modalDialog(
    title = shiny::tagList(shiny::icon("dna"), "Cell Phenotype Rules — Marker Definitions"),
    size = "l", easyClose = TRUE, fade = TRUE,
    footer = shiny::modalButton("Close"),
    shiny::tags$p(
      "Every cell in the app is classified into one of these phenotypes by the upstream ",
      shiny::tags$code("data/phenotype_rules.csv"), " gating rules. A ",
      shiny::tags$code("+"),
      " means the marker is required to be positive; a ",
      shiny::tags$code("−"),
      " means it must be negative. “any of X+ / Y+” means at least one of the listed markers must be positive."
    ),
    shiny::tags$div(style = "max-height: 60vh; overflow-y: auto; border-top: 1px solid #ddd;",
      shiny::tags$table(
        class = "table table-striped",
        style = "width: 100%; font-size: 14px; margin: 0;",
        shiny::tags$thead(header_row),
        shiny::tags$tbody(body_rows)
      )
    ),
    shiny::tags$p(style = "margin-top: 10px; font-size: 12px; color: #777;",
      "Rule source: ",
      shiny::tags$code("data/phenotype_rules.csv"),
      " (18 markers: CD99, INS, GCG, SST, SMA, CD20, CD3e, CD68, CD163, CD8a, CD4, HLADR, PDPN, CD34, CD31, B3TUBB, CD56)."
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

# UI button that opens the Phenotype Rules modal. The caller chooses the input
# id (one per tab to avoid duplicate-DOM-id collisions across simultaneously
# rendered tabPanels) and registers a matching `observeEvent` that calls
# `showModal(phenotype_rules_modal_ui())`. Style mirrors the existing "About
# this tab" info banner.
phenotype_rules_button <- function(input_id,
                                   label = "Phenotype rules — marker definitions") {
  shiny::tags$div(
    style = "margin-bottom: 12px; padding: 8px 14px; background-color: #eef4fb; border: 1px solid #c8d8e9; border-radius: 6px; display: inline-block;",
    shiny::actionLink(input_id,
      shiny::tagList(shiny::icon("dna", style = "color:#4477AA;"),
                     shiny::tags$span(style = "margin-left:6px; color:#2a4d75; font-weight:600;",
                                      label)),
      style = "text-decoration: none;")
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
      shiny::tags$span(style = "color:#888;", " · "),
      shiny::actionLink(show_modal_link, "see all rules →",
                        style = "font-size: 12px; color: #4477AA;")
    )
  }
  shiny::tags$div(
    style = "margin-top: -6px; margin-bottom: 8px; padding: 6px 10px; background-color: #f7f9fc; border-left: 3px solid #4477AA; border-radius: 3px; font-size: 13px; color: #333; line-height: 1.4;",
    shiny::tags$strong(style = "color:#4477AA;", paste0(prefix, " ")),
    shiny::tags$span(rule_html),
    link
  )
}
