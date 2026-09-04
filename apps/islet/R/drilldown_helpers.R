# ===== Single-Cell Drill-Down Helper Functions =====
# Provides cell-level visualization for the segmentation viewer.
# Loads per-islet CSV files from data/cells/ and renders cell scatter plots
# over the GeoJSON segmentation base.

# ── Phenotype display vocabulary ───────────────────────────────────────
# DOMAIN owns these UI labels/colors. It does not define upstream identities.
PHENOTYPE_COLORS <- domain_phenotype_colors()

# ── High Contrast palette: immune spread across purple/cyan/magenta/indigo ──
PHENOTYPE_COLORS_HIGH_CONTRAST <- c(
  "Beta cell"     = "#D32F2F",  # deep red
  "Alpha cell"    = "#FF8F00",  # amber
  "Delta cell"    = "#FDD835",  # bright yellow
  "Endocrine"     = "#8D6E63",  # brown
  "Acinar"        = "#00ACC1",  # cyan
  "Ductal"        = "#43A047",  # green
  "ECAD+"         = "#7CB342",  # lime green
  "Structural"    = "#9E9E9E",  # gray
  "SMA+"          = "#F06292",  # pink
  "Blood Vessel"  = "#E91E63",  # magenta-red
  "Endothelial"   = "#FF7043",  # deep orange
  "Lymphatic"     = "#26A69A",  # teal
  "Neural"        = "#AB47BC",  # purple
  "Macrophage"    = "#5C6BC0",  # indigo
  "APCs"          = "#26C6DA",  # light cyan
  "CD8a Tcell"    = "#7B1FA2",  # deep purple
  "CD4 Tcell"     = "#1565C0",  # strong blue
  "T cell"        = "#0097A7",  # dark cyan
  "B cell"        = "#EC407A",  # hot pink
  "Immune"        = "#283593",  # dark indigo
  "Unknown"       = "#BDBDBD"   # light gray
)

# ── Colorblind Safe palette: adapted from Paul Tol qualitative schemes ──
PHENOTYPE_COLORS_COLORBLIND <- c(
  "Beta cell"     = "#CC6677",  # rose
  "Alpha cell"    = "#DDCC77",  # sand
  "Delta cell"    = "#999933",  # olive
  "Endocrine"     = "#AA4499",  # purple
  "Acinar"        = "#44AA99",  # teal
  "Ductal"        = "#332288",  # indigo
  "ECAD+"         = "#117733",  # green
  "Structural"    = "#888888",  # gray
  "SMA+"          = "#882255",  # wine
  "Blood Vessel"  = "#EE8866",  # peach
  "Endothelial"   = "#EEDD88",  # light yellow
  "Lymphatic"     = "#77AADD",  # light blue
  "Neural"        = "#FFAABB",  # pink
  "Macrophage"    = "#44BB99",  # mint
  "APCs"          = "#BBCC33",  # pear
  "CD8a Tcell"    = "#AA4455",  # dark rose
  "CD4 Tcell"     = "#4477AA",  # steel blue
  "T cell"        = "#66CCEE",  # cyan
  "B cell"        = "#EE6677",  # coral
  "Immune"        = "#225566",  # dark teal
  "Unknown"       = "#BBBBBB"   # light gray
)

# ── Maximum Distinction palette: adapted from Kelly's 22 maximally distinct ──
PHENOTYPE_COLORS_MAX_DISTINCT <- c(
  "Beta cell"     = "#BE0032",  # vivid red
  "Alpha cell"    = "#F3C300",  # vivid yellow
  "Delta cell"    = "#875692",  # strong purple
  "Endocrine"     = "#F38400",  # vivid orange
  "Acinar"        = "#A1CAF1",  # very light blue
  "Ductal"        = "#008856",  # vivid green
  "ECAD+"         = "#C2B280",  # dark yellowish brown
  "Structural"    = "#848482",  # medium gray
  "SMA+"          = "#E68FAC",  # purplish pink
  "Blood Vessel"  = "#F99379",  # vivid orange-yellow
  "Endothelial"   = "#604E97",  # strong violet
  "Lymphatic"     = "#0067A5",  # strong blue
  "Neural"        = "#DCD300",  # vivid greenish yellow
  "Macrophage"    = "#7F180D",  # vivid reddish brown
  "APCs"          = "#B3446C",  # deep purplish pink
  "CD8a Tcell"    = "#F6A600",  # vivid yellowish brown
  "CD4 Tcell"     = "#2B3D26",  # dark olive green
  "T cell"        = "#8DB600",  # vivid yellowish green
  "B cell"        = "#654522",  # deep yellowish brown
  "Immune"        = "#E25822",  # vivid reddish orange
  "Unknown"       = "#C0C0C0"   # silver gray
)

# ── Palette registry and accessor ────────────────────────────────────────
complete_phenotype_palette <- function(palette, fallback = PHENOTYPE_COLORS) {
  missing <- setdiff(names(fallback), names(palette))
  c(palette, fallback[missing])
}

PHENOTYPE_PALETTES <- list(
  "Original"             = PHENOTYPE_COLORS,
  "High Contrast"        = complete_phenotype_palette(PHENOTYPE_COLORS_HIGH_CONTRAST),
  "Colorblind Safe"      = complete_phenotype_palette(PHENOTYPE_COLORS_COLORBLIND),
  "Maximum Distinction"  = complete_phenotype_palette(PHENOTYPE_COLORS_MAX_DISTINCT)
)

get_phenotype_palette <- function(name = "Original") {
  PHENOTYPE_PALETTES[[name]] %||% PHENOTYPE_COLORS
}

# ── Cell CSV cache (env-based, same pattern as geojson_cache) ──────────
drilldown_cache <- new.env()

# Path to per-islet cell CSVs
CELLS_DIR <- file.path("..", "..", "data", "app_data", "cells")

#' Check if single-cell drill-down data is available
drilldown_available <- function() {
  dir.exists(CELLS_DIR) && length(list.files(CELLS_DIR, pattern = "\\.csv$")) > 0
}

#' Load single-cell data for a specific islet
#' @param imageid Case/image ID (e.g., "6505")
#' @param islet_key Islet key (e.g., "Islet_284")
#' @return data.frame with columns: X_centroid, Y_centroid, phenotype, cell_region, + markers
load_islet_cells <- function(imageid, islet_key) {
  combined_id <- paste0(imageid, "_", islet_key)
  cache_key <- combined_id

  if (exists(cache_key, envir = drilldown_cache)) {
    return(get(cache_key, envir = drilldown_cache))
  }

  csv_path <- file.path(CELLS_DIR, paste0(combined_id, ".csv"))
  if (!file.exists(csv_path)) {
    return(NULL)
  }

  cells <- tryCatch({
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(NULL)
    df
  }, error = function(e) {
    message("[DRILLDOWN] Error loading cells for ", combined_id, ": ", e$message)
    NULL
  })

  if (!is.null(cells)) {
    assign(cache_key, cells, envir = drilldown_cache)
    message("[DRILLDOWN] Cached cells for ", combined_id, " (", nrow(cells), " cells)")
  }

  cells
}

#' Render single-cell drill-down plot over GeoJSON base
#' @param info List with case_id, islet_key, centroid_x, centroid_y
#' @param cells data.frame from load_islet_cells()
#' @param color_by "phenotype" or a marker column name
#' @param show_peri Logical; if FALSE, only show core cells
#' @return ggplot object
render_islet_drilldown_plot <- function(info, cells, color_by = "phenotype", show_peri = TRUE,
                                        show_peri_boundary = TRUE, show_structures = TRUE,
                                        palette = PHENOTYPE_COLORS) {
  if (is.null(info) || is.null(cells) || nrow(cells) == 0) return(NULL)

  # Build base plot: GeoJSON boundaries (islet core always drawn)
  base_plot <- build_segmentation_base_plot(info, show_peri_boundary = show_peri_boundary,
                                             show_structures = show_structures)
  if (is.null(base_plot)) base_plot <- ggplot2::ggplot() + ggplot2::theme_void()

  # Filter cells by region
  if (!show_peri && "cell_region" %in% colnames(cells)) {
    cells <- cells[cells$cell_region == "core", , drop = FALSE]
  }
  if (nrow(cells) == 0) return(base_plot)

  # Convert cell centroids from µm to pixels (matching GeoJSON coordinate space)
  cells$x_px <- cells$X_centroid / PIXEL_SIZE_UM
  cells$y_px <- cells$Y_centroid / PIXEL_SIZE_UM

  if (color_by == "phenotype" && "phenotype" %in% colnames(cells)) {
    # Categorical coloring by phenotype (using fill to avoid colour conflict with structures)
    pheno_present <- sort(unique(cells$phenotype))
    pal <- palette[pheno_present]
    # Fill missing phenotypes with gray
    pal[is.na(pal)] <- "#CCCCCC"

    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px, fill = phenotype),
        shape = 21, colour = "grey30", stroke = 0.3, size = 3.0, alpha = 0.9,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(values = pal, name = "Phenotype") +
      ggplot2::labs(title = paste(info$islet_key, "- Single Cells"))
  } else if (color_by %in% colnames(cells)) {
    # Continuous coloring by marker expression (using fill)
    cells$marker_val <- suppressWarnings(as.numeric(cells[[color_by]]))

    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px, fill = marker_val),
        shape = 21, colour = "grey30", stroke = 0.3, size = 3.0, alpha = 0.9,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_viridis_c(option = "inferno", name = color_by, na.value = "gray80") +
      ggplot2::labs(title = paste(info$islet_key, "-", color_by))
  } else {
    # Fallback: just plot cells in gray
    p <- base_plot +
      ggplot2::geom_point(
        data = cells,
        ggplot2::aes(x = x_px, y = y_px),
        color = "gray40", size = 2.5, alpha = 0.7,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(title = paste(info$islet_key, "- Single Cells"))
  }

  # theme_islet() FIRST, then the local overrides. Without it this plot inherited
  # whatever base_plot carried — and on the `theme_void()` fallback path (no
  # GeoJSON) that meant ggplot's 11pt default for the entire phenotype legend,
  # the densest text in the drill-down.
  p +
    theme_islet() +
    ggplot2::theme(
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "cm"),
      legend.text     = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_text),
      legend.title    = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_title,
                                              face = "bold"),
      # Spatial map: axes are pixel coordinates and carry no meaning.
      axis.title      = ggplot2::element_blank(),
      axis.text       = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank()
    )
}

#' Render phenotype composition summary bar chart
#' @param cells data.frame from load_islet_cells()
#' @return ggplot object
render_drilldown_summary <- function(cells, palette = PHENOTYPE_COLORS) {
  if (is.null(cells) || nrow(cells) == 0 || !"phenotype" %in% colnames(cells)) {
    return(NULL)
  }

  # Count phenotypes, split by region
  if ("cell_region" %in% colnames(cells)) {
    counts <- as.data.frame(table(cells$phenotype, cells$cell_region),
                            stringsAsFactors = FALSE)
    colnames(counts) <- c("phenotype", "region", "count")
    counts$region <- factor(counts$region, levels = c("core", "peri"))
  } else {
    counts <- as.data.frame(table(cells$phenotype), stringsAsFactors = FALSE)
    colnames(counts) <- c("phenotype", "count")
    counts$region <- "all"
  }

  # Order by total count descending
  total_by_pheno <- tapply(counts$count, counts$phenotype, sum)
  counts$phenotype <- factor(counts$phenotype,
                              levels = names(sort(total_by_pheno, decreasing = TRUE)))

  # Colors
  pheno_present <- levels(counts$phenotype)
  pal <- palette[pheno_present]
  pal[is.na(pal)] <- "#CCCCCC"

  if ("cell_region" %in% colnames(cells) && length(unique(counts$region)) > 1) {
    ggplot2::ggplot(counts, ggplot2::aes(x = phenotype, y = count, fill = phenotype, colour = region)) +
      ggplot2::geom_col(position = "stack", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::scale_colour_manual(values = c("core" = "grey20", "peri" = "white"), name = "Region") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Cell count") +
      theme_islet() +
      ggplot2::theme(legend.position = "bottom")
  } else {
    ggplot2::ggplot(counts, ggplot2::aes(x = phenotype, y = count, fill = phenotype)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Cell count") +
      theme_islet()
  }
}
