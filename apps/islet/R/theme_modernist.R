# ---------------------------------------------------------------------------
# Modernist plot constants + scale helpers
#
# Ported from design/delivered/shiny/R/theme_modernist.R. Read design/README.md
# for what was deliberately not carried across.
#
# SOURCE ORDER: Shiny auto-sources R/*.R alphabetically, so this file loads
# after 00_globals.R (ISLET_FONT_SIZES, ISLET_PLOT_FAMILY, DONOR_PALETTES) and
# after drilldown_helpers.R (PHENOTYPE_COLORS, get_phenotype_palette). The
# aliases below depend on that. Do NOT rename this file to a `0*_` prefix.
#
# THE RULE THIS FILE EXISTS TO ENFORCE: Modernist governs CHROME. It does not
# govern data. The app's donor and phenotype palettes are colourblind-safe and
# scientifically load-bearing -- four selectable variants each, synced through
# active_donor_colors() / active_palette() -- and the redesign does not touch
# them. Everything here is either a chrome constant or a thin wrapper over the
# existing palettes.
# ---------------------------------------------------------------------------

# ---- Chrome constants ------------------------------------------------------
# Used for plot ink, rules and muted text. NOT for data encodings; the fence in
# scripts/check_ui_style_fences.sh keeps IE_CRIMSON / IE_GOLD out of the module
# servers and the helpers. Measured for THIS app (CIEDE2000, sRGB/D65) over the 120
# distinct colours its data palettes declare, crimson has five collisions inside
# deltaE 10 -- 1.2 from #BE0032 (Beta cell, Max Distinction), 2.4 from #B2182B (Aab+,
# "Diverging"), 7.8 from #BC4749 (Blood vessel), 8.9 from #AA4455, 9.2 from #D32F2F --
# and gold is 9.3 from #CCBB44, which is PAL_LEIDEN slot 4, twelve lines below. Under
# deltaE 10 two colours are not reliably told apart, so a crimson rule beside a
# Beta-cell point reads as a Beta-cell key.
IE_INK     <- "#201e1d"
IE_GROUND  <- "#f3f2f2"
IE_RULE    <- "#c9c5c2"
IE_MUTED   <- "#6f6a67"
IE_CRIMSON <- "#b71533"
IE_GOLD    <- "#c8a02e"
# The 700 ramp step. Modernist's own rule is that accent text at body size uses
# a deep step rather than the accent itself (the accent-to-ground pair is tuned
# to 3:1, enough for icons and chrome, not for copy). deltaE 10.4 from #B2182B
# rather than crimson's 2.4 -- four times further, which is why it is the text step.
# It is NOT clean, though, and is not fenced: deltaE 8.0 from #7F180D (Macrophage,
# Max Distinction). Ink and copy only; never a fill.
IE_ACCENT_700 <- "#800c21"

# Disease-stage levels, in progression order. Replaces a dozen repeats of the
# literal c("ND", "Aab+", "T1D") across the module servers.
# NOTE: Junior's copy defines domain_group_levels() <- c("ND", "Aab+", "T1D") here and uses it as the
# `limits` and `breaks` of both donor scales. That literal is NOT ported. The two datasets declare
# their own levels, and the run cohort has no Aab+ donor at all -- a fixed `limits` with
# drop = FALSE would paint an empty Aab+ entry into every donor legend. The scales resolve
# domain_group_levels() at CALL time, not load time, so one running app can serve either domain.

# ---- Palette aliases -------------------------------------------------------
# The delivered file hardcoded its own palettes. These point at the canonical
# ones instead, so a retune in 00_globals.R / drilldown_helpers.R propagates.
#
# PAL_DONOR is byte-identical to DONOR_COLORS_BRIGHT today (#3366CC / #FF6600 /
# #109618) -- verified -- so the alias is a no-op now and stays correct later.
PAL_DONOR <- DONOR_COLORS_BRIGHT

# The delivered PAL_PHENO was a 12-entry DEMO set (`Beta`, `CD8 T`, `Vessel`, ...).
# The app has 21 real phenotypes (`Beta cell`, `CD8a Tcell`, `Blood Vessel`,
# `APCs`, `ECAD+`, `Unknown`, ...); only 5 of the 12 keys matched anything, and
# with scale_*_manual(drop = FALSE) the other 16 levels would render NA.
PAL_PHENO <- PHENOTYPE_COLORS

# Leiden runs at four resolutions and exceeds 6 clusters at the higher ones, so
# this is a recycling source rather than a fixed vector. Paul Tol bright, then
# muted, for colourblind safety at higher cluster counts.
PAL_LEIDEN <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377",
                "#BBBBBB", "#332288", "#88CCEE", "#44AA99", "#999933", "#882255")

# ---- Scale helpers ---------------------------------------------------------

#' Donor-status colour scale.
#'
#' @param values NAMED palette. **Always pass the reactive from the call site**
#'   -- `scale_colour_donor(donor_colors_reactive())`. The default exists for
#'   headless and test contexts only and matches get_donor_color_palette()'s own
#'   default. Omitting it at a live call site silently freezes the user's
#'   four-way palette selector, which is why this note is shouty.
scale_colour_donor <- function(values = get_donor_color_palette(),
                               name = "Donor status", drop = FALSE, ...) {
  ggplot2::scale_colour_manual(values = values, limits = domain_group_levels(),
                               breaks = domain_group_levels(), name = name, drop = drop, ...)
}
scale_fill_donor <- function(values = get_donor_color_palette(),
                             name = "Donor status", drop = FALSE, ...) {
  ggplot2::scale_fill_manual(values = values, limits = domain_group_levels(),
                             breaks = domain_group_levels(), name = name, drop = drop, ...)
}

#' Phenotype colour / fill scales.
#'
#' Both exist because the drill-down deliberately encodes phenotype with `fill`
#' (shape 21) -- `colour` is reserved for structure boundaries so the
#' "Structures" legend can coexist with the phenotype legend. The delivered file
#' shipped only the colour variant.
scale_colour_pheno <- function(values = get_phenotype_palette(),
                               name = "Phenotype", na.value = "#d9d9d9",
                               drop = FALSE, ...) {
  ggplot2::scale_colour_manual(values = values, name = name,
                               na.value = na.value, drop = drop, ...)
}
scale_fill_pheno <- function(values = get_phenotype_palette(),
                             name = "Phenotype", na.value = "#d9d9d9",
                             drop = FALSE, ...) {
  ggplot2::scale_fill_manual(values = values, name = name,
                             na.value = na.value, drop = drop, ...)
}

#' Leiden cluster scale. `n` recycles the palette to the cluster count.
scale_colour_leiden <- function(n = NULL, name = "Cluster",
                                na.value = "#d9d9d9", ...) {
  v <- if (is.null(n)) PAL_LEIDEN else rep_len(PAL_LEIDEN, n)
  ggplot2::scale_colour_manual(values = v, name = name, na.value = na.value, ...)
}

#' DEPRECATED. Significance is marked with a STAR, never with colour.
#'
#' Red/green reads as stop/go and bad/good; a p-value crossing a threshold is
#' neither -- it is a statement about evidence, not about whether the result is
#' welcome. Every site that used this now draws in one ink colour and appends a
#' star. Kept as an identity scale so any straggler call site degrades to "no
#' colour encoding" rather than erroring.
scale_colour_significance <- function(...) {
  warning("scale_colour_significance() is deprecated: mark significance with a star, not colour.",
          call. = FALSE)
  ggplot2::scale_colour_manual(values = c(`TRUE` = IE_INK, `FALSE` = IE_INK),
                               guide = "none", drop = FALSE)
}

# ---- theme_modernist -------------------------------------------------------

#' Modernist ggplot theme.
#'
#' A thin alias over `theme_islet()`, which stays the single source of truth for
#' plot styling (app/shiny_app/CLAUDE.md, "Plot Fonts & Publication-Quality
#' Downloads"). Exists so delivered/demo code calling `theme_modernist()` works
#' unchanged.
#'
#' Four things in the delivered theme are deliberately NOT adopted:
#'   * `base_size = 12`, which put axis text and titles at 10pt. This app's
#'     scale is base 18 / axis 16, raised on purpose in Jul 2026 for legibility,
#'     across 23 call sites.
#'   * `legend.title = element_blank()`. Four legends are named -- "Structures",
#'     "Phenotype", "Region", "Cluster" -- and the segmentation plot reserves the
#'     `colour` aesthetic for structures precisely so that legend can be titled.
#'   * `legend.position = "bottom"`. A 21-phenotype legend consumes four or five
#'     rows under a coord_fixed() plot, and 16 sites already set position
#'     per-plot.
#'   * `plot.background = IE_GROUND`. That is a SCREEN ground; baking it into
#'     the theme would put grey behind all 8 ggsave() publication PNGs. It
#'     belongs in CSS, on the plot container.
theme_modernist <- function(base_size = ISLET_FONT_SIZES$base,
                            grid = c("both", "x", "y", "none"),
                            family = ISLET_PLOT_FAMILY) {
  theme_islet(base_size = base_size, family = family, transparent = TRUE)
}

#' Deprecated. Font registration moved to `ie_register_plot_font()` in
#' 00_globals.R, which runs at source time and uses systemfonts rather than
#' showtext. Kept so the delivered app.R runs unchanged; returns the family that
#' was actually registered.
#' Register Archivo for ggplot from the copy vendored in www/fonts.
#'
#' Junior resolved the family by NAME, which only works where Archivo is installed system-wide;
#' `fc-list` finds none on this host, so a name lookup would silently fall back to the default
#' sans and every figure would disagree with the page around it. These faces are vendored in the
#' repo, so they are registered by PATH and the result is CHECKED rather than assumed. On failure
#' the family is "" -- ggplot's default, a plain figure, never an error.
ie_register_plot_font <- function() {
  faces <- file.path("www", "fonts", c("Archivo-Regular.ttf", "Archivo-ExtraBold.ttf"))
  if (!all(file.exists(faces)) || !requireNamespace("systemfonts", quietly = TRUE)) return("")
  ok <- tryCatch({
    systemfonts::register_font(name = "Archivo", plain = faces[[1]], bold = faces[[2]])
    "Archivo" %in% systemfonts::registry_fonts()$family
  }, error = function(e) FALSE)
  if (isTRUE(ok)) "Archivo" else ""
}

ISLET_PLOT_FAMILY <- ie_register_plot_font()

# Make the registered family and the card ground the app-wide default for every
# figure. There are 22 `theme_islet()` call sites across the module servers and
# the helpers; each also adds its own `+ theme(...)` afterwards, which still
# wins, so this changes the BASE only. Options rather than an override -- see
# theme_islet()'s own note in apps/R/ieCore/R/theme.R for why an override would
# recurse. Both domains of this app get the design, which is right: it is the
# same app skin, not a statement about either cohort's data.
options(islet.plot.family = ISLET_PLOT_FAMILY,
        islet.plot.transparent = TRUE)

# The plotly half of the same job, by a different route. plotly renders in the
# browser, so it needs a CSS stack rather than a registered device family, and
# the fallbacks are mandatory: a browser that has not loaded Archivo substitutes
# silently, and the figure would then disagree with the page around it. Archivo
# is named first only when it actually registered -- ISLET_PLOT_FAMILY is "" if
# the vendored faces are missing.
options(islet.plotly.family = paste0(
  if (nzchar(ISLET_PLOT_FAMILY)) paste0("'", ISLET_PLOT_FAMILY, "', ") else "",
  "'Liberation Sans', 'Helvetica Neue', Helvetica, Arial, system-ui, sans-serif"
))

#' Kept so the delivered app.R runs unchanged; returns the family actually registered.
register_modernist_fonts <- function(...) invisible(ISLET_PLOT_FAMILY)
