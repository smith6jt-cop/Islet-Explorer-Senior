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

# Pixel size constant for coordinate conversion (micrometers per pixel)
# GeoJSON polygons use pixel coordinates, while islet_spatial_lookup.csv uses micrometers
PIXEL_SIZE_UM <- 0.3774  # micrometers per pixel

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
master_path <- file.path("..", "..", "data", "master_results.xlsx")
h5ad_path   <- file.path("..", "..", "data", "islet_explorer.h5ad")
project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# Ensure reticulate/anndata can discover a Python binary when RETICULATE_PYTHON
# isn't set (shiny-server environments often lack PATH).
if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  candidate_python <- NULL
  venv_path <- file.path(path.expand("~"), ".local", "share", "islet-explorer-py", "bin", "python")
  if (file.exists(venv_path)) candidate_python <- venv_path
  if (is.null(candidate_python) || !nzchar(candidate_python)) {
    candidate_python <- Sys.getenv("PYTHON", Sys.which("python3"))
    if (!nzchar(candidate_python)) candidate_python <- Sys.which("python")
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

# Standard PNG download dimensions for ggplot ggsave() in downloadHandlers.
# Default: 10x7 in at 300 dpi -> 3000x2100 px. Heatmaps use wider 12x8.
ISLET_PNG_DIMS  <- list(width = 10, height = 7, dpi = 300, units = "in")
ISLET_HEATMAP_PNG_DIMS <- list(width = 12, height = 8, dpi = 300, units = "in")

islet_png_filename <- function(stem) {
  paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
}
