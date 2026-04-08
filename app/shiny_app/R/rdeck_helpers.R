# ---------- WebGL rendering helpers (Phase 3) ----------
#
# The legacy Spatial and Trajectory tabs render with ggplot2 (tissue scatter)
# and Plotly (UMAPs) respectively. Both perform fine at the current 5,214
# islet / ~177K cell-per-donor scale, but neither survives the jump to
# millions of cells. This module adds an opt-in WebGL backend powered by the
# `rdeck` package (deck.gl bindings for R), while leaving the existing code
# paths untouched for backward compatibility.
#
# Activation:
#   - Install: remotes::install_github("qfes/rdeck")
#   - Set the environment variable ISLET_USE_RDECK=TRUE before launching.
#   - When both are true and point counts exceed RDECK_POINT_THRESHOLD,
#     modules swap their scatter output for an rdeck widget.
#
# Boundaries:
#   - This file must not depend on any module's reactive values.
#   - Pure helpers only; modules own the `reactive` wiring so they can
#     collect() the right subset of data before calling us.
#
# Design notes:
#   - Tissue scatter expects a data.frame with X_centroid, Y_centroid,
#     phenotype, cell_region columns — identical to `load_donor_tissue()`.
#   - UMAP scatter expects a data.frame with umap_1, umap_2, and a
#     categorical color column.
#   - Color palettes are passed in so the module can reuse PHENOTYPE_COLORS
#     or DONOR_COLORS without this helper having to know about globals.

USE_RDECK <- isTRUE(as.logical(Sys.getenv("ISLET_USE_RDECK", "FALSE")))

# Only switch to rdeck for plots bigger than this. Below the threshold the
# ggplot path is usually crisper and animates faster for legends + brushes.
RDECK_POINT_THRESHOLD <- suppressWarnings(
  as.integer(Sys.getenv("ISLET_RDECK_THRESHOLD", "50000"))
)
if (!is.finite(RDECK_POINT_THRESHOLD) || RDECK_POINT_THRESHOLD < 1) {
  RDECK_POINT_THRESHOLD <- 50000L
}

#' Is rdeck installed and enabled for this session?
is_rdeck_available <- function() {
  isTRUE(USE_RDECK) && requireNamespace("rdeck", quietly = TRUE)
}

#' Decide whether a given scatter should use the rdeck WebGL backend.
#'
#' @param n_points Number of points to plot.
#' @return TRUE if rdeck should render this scatter.
should_use_rdeck <- function(n_points) {
  if (!is_rdeck_available()) return(FALSE)
  is.numeric(n_points) && length(n_points) == 1 && n_points >= RDECK_POINT_THRESHOLD
}

#' Hex palette lookup helper that tolerates unknown levels.
#'
#' rdeck's `scale_color_category` consumes a *named* vector mapping level to
#' hex. `fallback` is used for any level the palette doesn't know about so we
#' never emit NA fills.
resolve_color_vector <- function(levels, palette, fallback = "#999999") {
  levels <- unique(as.character(levels))
  out <- palette[levels]
  out[is.na(out)] <- fallback
  setNames(out, levels)
}

#' Build the rdeck widget for the Spatial-tab tissue scatter.
#'
#' @param cells data.frame with X_centroid, Y_centroid, phenotype, cell_region.
#' @param palette Named vector of phenotype hex colors (e.g. PHENOTYPE_COLORS).
#' @param title Optional title drawn above the canvas (plain text).
#' @return An rdeck htmlwidget (or NULL if rdeck is not installed).
render_tissue_deck <- function(cells, palette, title = NULL) {
  if (!is_rdeck_available()) return(NULL)
  if (is.null(cells) || nrow(cells) == 0) return(NULL)

  # Bounding box is just a sanity clamp; rdeck auto-fits if we omit it.
  xr <- range(cells$X_centroid, na.rm = TRUE)
  yr <- range(cells$Y_centroid, na.rm = TRUE)
  if (!all(is.finite(c(xr, yr)))) return(NULL)

  pal <- resolve_color_vector(cells$phenotype, palette)

  rdeck::rdeck(
    map_style = NULL,                 # No basemap — tissue um coords, not geographic
    initial_view_state = rdeck::view_state(
      center = c((xr[1] + xr[2]) / 2, (yr[1] + yr[2]) / 2),
      zoom = 0, pitch = 0, bearing = 0
    )
  ) |>
    rdeck::add_scatterplot_layer(
      data = cells,
      get_position = ~ cbind(X_centroid, Y_centroid),
      get_fill_color = rdeck::scale_color_category(
        col = phenotype,
        palette = pal
      ),
      get_radius = 3,
      radius_min_pixels = 1,
      radius_max_pixels = 5,
      pickable = TRUE,
      auto_highlight = TRUE
    )
}

#' Build the rdeck widget for islet-level UMAPs (trajectory + spatial tabs).
#'
#' @param umap data.frame with umap_1, umap_2 and the color column.
#' @param color_col Name of the categorical color column (default "donor_status").
#' @param palette Named hex vector for the category levels.
#' @return An rdeck htmlwidget (or NULL if rdeck is not installed).
render_umap_deck <- function(umap, color_col = "donor_status", palette = NULL) {
  if (!is_rdeck_available()) return(NULL)
  if (is.null(umap) || nrow(umap) == 0) return(NULL)
  if (!all(c("umap_1", "umap_2", color_col) %in% colnames(umap))) return(NULL)

  if (is.null(palette)) {
    # Default to the donor status palette from 00_globals.R.
    palette <- get0("DONOR_COLORS", envir = globalenv(),
                    ifnotfound = c("ND" = "#4477AA", "Aab+" = "#CC6633", "T1D" = "#228833"))
  }
  pal <- resolve_color_vector(umap[[color_col]], palette)

  xr <- range(umap$umap_1, na.rm = TRUE)
  yr <- range(umap$umap_2, na.rm = TRUE)
  if (!all(is.finite(c(xr, yr)))) return(NULL)

  rdeck::rdeck(
    map_style = NULL,
    initial_view_state = rdeck::view_state(
      center = c((xr[1] + xr[2]) / 2, (yr[1] + yr[2]) / 2),
      zoom = 0, pitch = 0, bearing = 0
    )
  ) |>
    rdeck::add_scatterplot_layer(
      data = umap,
      get_position = ~ cbind(umap_1, umap_2),
      get_fill_color = rdeck::scale_color_category(
        col = !!rlang::sym(color_col),
        palette = pal
      ),
      get_radius = 5,
      radius_min_pixels = 2,
      radius_max_pixels = 8,
      pickable = TRUE,
      auto_highlight = TRUE
    )
}

#' Output wrapper that returns the correct Shiny output slot.
#'
#' Modules use this in their UI (or in a `renderUI`) so the same code path
#' picks between rdeck and the legacy plot/plotly output. When rdeck is
#' disabled the `fallback_output` function is called with `id` and returned
#' as-is.
#'
#' @param id Namespaced output id (e.g. `ns("tissue_scatter")`).
#' @param fallback_output A function like `plotOutput` or `plotlyOutput` that
#'   will be called with `id` and any additional `...` args when rdeck is not
#'   available.
#' @param ... Passed through to `fallback_output`.
rdeck_or_output <- function(id, fallback_output, ...) {
  if (is_rdeck_available()) {
    return(rdeck::rdeckOutput(id, ...))
  }
  fallback_output(id, ...)
}
