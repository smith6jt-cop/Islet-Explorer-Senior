# Shared presentation primitives. Data palettes come from DOMAIN; these helpers
# define typography only and do not encode scientific groups or calls.

ISLET_FONT_SIZES <- list(
  base = 18,
  title = 18,
  subtitle = 15,
  axis_title = 16,
  axis_text = 16,
  legend_title = 16,
  legend_text = 16,
  strip = 16,
  dense_min = 14,
  hover = 15,
  annotation = 15
)

ISLET_ANNOTATE_SIZE <- ISLET_FONT_SIZES$base / ggplot2::.pt

theme_islet <- function(base_size = ISLET_FONT_SIZES$base) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ISLET_FONT_SIZES$title, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ISLET_FONT_SIZES$subtitle, color = "#555"),
      axis.title = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_title),
      axis.text = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_text),
      legend.title = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_title, face = "bold"),
      legend.text = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_text),
      strip.text = ggplot2::element_text(size = ISLET_FONT_SIZES$strip, face = "bold")
    )
}

plotly_islet_fonts <- function() {
  list(
    global = list(size = ISLET_FONT_SIZES$base),
    title = list(size = ISLET_FONT_SIZES$title),
    axis_title = list(size = ISLET_FONT_SIZES$axis_title),
    axis_tick = list(size = ISLET_FONT_SIZES$axis_text),
    legend_title = list(size = ISLET_FONT_SIZES$legend_title),
    legend_text = list(size = ISLET_FONT_SIZES$legend_text),
    colorbar_title = list(size = ISLET_FONT_SIZES$legend_title)
  )
}

plotly_islet_config <- function(p, filename_root = "islet_explorer_plot",
                                width = 1400, height = 900, scale = 2,
                                uniform_text_mode = "hide") {
  p <- plotly::layout(
    p,
    hoverlabel = list(font = list(size = ISLET_FONT_SIZES$hover)),
    uniformtext = list(minsize = ISLET_FONT_SIZES$dense_min, mode = uniform_text_mode)
  )
  plotly::config(
    p,
    toImageButtonOptions = list(
      format = "png", filename = filename_root, width = width,
      height = height, scale = scale
    ),
    displaylogo = FALSE
  )
}
