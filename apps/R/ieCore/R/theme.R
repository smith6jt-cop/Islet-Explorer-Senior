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

#' @param family Plot font family. Defaults to `getOption("islet.plot.family")`, which is unset
#'   -- so ggplot's own default -- for every caller that does not opt in. The islet app sets it
#'   to the Archivo it vendors, in R/theme_modernist.R; the lymph fixture and any other consumer
#'   are unchanged.
#' @param transparent Draw on no background at all, from `getOption("islet.plot.transparent")`.
#'   The Modernist design puts figures directly on a card, so an opaque white plot background --
#'   which `theme_minimal()` sets, and which therefore survives any device or CSS setting --
#'   would paint a white rectangle inside a grey one. This is a SCREEN choice: a transparent PNG
#'   inverts on a dark slide, so the download handlers pass `bg = "white"` to ggsave and the
#'   saved file deliberately differs from the screen.
#'
#'   Both are OPTIONS rather than changed defaults so there is ONE definition of this theme. The
#'   app cannot simply redefine `theme_islet`: ieCore is `sys.source`d into the global
#'   environment rather than attached as a package, so an override would have to capture the
#'   original under another name, and a second source pass would then capture its own override
#'   and recurse.
theme_islet <- function(base_size = ISLET_FONT_SIZES$base,
                        family = getOption("islet.plot.family", ""),
                        transparent = getOption("islet.plot.transparent", FALSE)) {
  base <- ggplot2::theme_minimal(base_size = base_size, base_family = family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ISLET_FONT_SIZES$title, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = ISLET_FONT_SIZES$subtitle, color = "#555"),
      axis.title = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_title),
      axis.text = ggplot2::element_text(size = ISLET_FONT_SIZES$axis_text),
      legend.title = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_title, face = "bold"),
      legend.text = ggplot2::element_text(size = ISLET_FONT_SIZES$legend_text),
      strip.text = ggplot2::element_text(size = ISLET_FONT_SIZES$strip, face = "bold")
    )
  if (!isTRUE(transparent)) return(base)
  base + ggplot2::theme(
    plot.background   = ggplot2::element_rect(fill = "transparent", colour = NA),
    panel.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.key        = ggplot2::element_rect(fill = "transparent", colour = NA)
  )
}

#' Plotly font specs, ready to splice into `layout()`.
#'
#' @param family A CSS font stack, from `getOption("islet.plotly.family")`. Plotly renders in
#'   the browser, so this is a CSS stack and not the graphics-device family `theme_islet()`
#'   takes -- the ggplot figures and the plotly figures reach the same typeface by two entirely
#'   different routes, and a fallback list is mandatory here because a face the browser has not
#'   loaded is silently substituted. Unset by default, which is plotly's own default font.
plotly_islet_fonts <- function(family = getOption("islet.plotly.family", "")) {
  size_and_family <- function(size) {
    if (!nzchar(family)) list(size = size) else list(size = size, family = family)
  }
  list(
    global = size_and_family(ISLET_FONT_SIZES$base),
    title = size_and_family(ISLET_FONT_SIZES$title),
    axis_title = size_and_family(ISLET_FONT_SIZES$axis_title),
    axis_tick = size_and_family(ISLET_FONT_SIZES$axis_text),
    legend_title = size_and_family(ISLET_FONT_SIZES$legend_title),
    legend_text = size_and_family(ISLET_FONT_SIZES$legend_text),
    colorbar_title = size_and_family(ISLET_FONT_SIZES$legend_title)
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
