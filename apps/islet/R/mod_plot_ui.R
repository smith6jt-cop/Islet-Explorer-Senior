# ---------- Plot module UI ----------
# Exports: plot_sidebar_ui(id), plot_main_ui(id)
# All Plot-related inputs are namespaced by this module.

#' Sidebar UI for the Plot tab
#' Contains data filters, autoantibody filters, statistic selector, and export controls.
plot_sidebar_ui <- function(id) {
  ns <- NS(id)
  # No `.sidebar` class any more: this block now renders inside the shell's
  # frozen control panel (`.ie-aside` in app.R), which owns the surface. The
  # old class carried a dark navy ground of its own and would fight it.
  div(
    div(class = "ie-aside-head", "Data & Filters"),
    selectInput(ns("mode"), "Data Layer",
                choices = c("Microenvironment", "Cell Populations"),
                selected = "Cell Populations"),
    uiOutput(ns("dynamic_selector")),
    uiOutput(ns("metric_selector")),
    hr(),
    h5("Autoantibody filter (Aab+ donors only)"),
    uiOutput(ns("aab_warning")),
    checkboxGroupInput(ns("aab_flags"), NULL,
                       choices = c("GADA" = "AAb_GADA",
                                   "IA2A" = "AAb_IA2A",
                                   "ZnT8A" = "AAb_ZnT8A",
                                   "IAA" = "AAb_IAA",
                                   "mIAA" = "AAb_mIAA"),
                       selected = c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")),
    checkboxGroupInput(ns("groups"), "Donor Status",
                       choices = c("ND", "Aab+", "T1D"),
                       selected = c("ND", "Aab+", "T1D"),
                       inline = TRUE),
    conditionalPanel(
      condition = "input.tabs != 'Statistics'", ns = NS(NULL),
      numericInput(ns("min_cells"), "Min cells/islet",
                   value = 1, min = 1, max = 200, step = 1)
    ),
    uiOutput(ns("age_filter_ui")),
    uiOutput(ns("gender_filter_ui")),
    conditionalPanel(
      condition = "input.tabs != 'Statistics'", ns = NS(NULL),
      radioButtons(ns("stat"), "Statistic",
                   choices = c("Mean\u00b1SE" = "mean_se",
                               "Mean\u00b1SD" = "mean_sd",
                               "Median + IQR" = "median_iqr"),
                   selected = "mean_se"),
      # Global outlier toggle + z-threshold slider (non-namespaced so they
      # can be synced with the mirrored controls on the Trajectory and Spatial
      # tabs). Defaults to OFF; the threshold slider value (default 3, range
      # 0.5–10, step 0.5) feeds the z-score `flag_outliers()` calls across all
      # three tabs. Title attribute carries the global tooltip.
      tags$div(
        title = GLOBAL_OUTLIER_TOOLTIP,
        checkboxInput("sidebar_remove_outliers",
                      "Exclude outliers",
                      value = FALSE),
        sliderInput("sidebar_outlier_z", "Outlier z-threshold",
                    min = 0.5, max = 10, value = 3, step = 0.5)
      ),
      checkboxInput(ns("show_outlier_table"),
                    "Show excluded-outlier table",
                    value = FALSE),
      uiOutput(ns("plot_outlier_info"))
    ),
    conditionalPanel(
      condition = "input.tabs != 'Statistics'", ns = NS(NULL),
      hr(),
      h5("Color Palettes"),
      # Non-namespaced: global donor palette choice synced via app.R observers
      selectInput("sidebar_donor_palette", "Donor Status Colors",
                  choices = names(DONOR_PALETTES),
                  selected = "Bright")
    ),
    conditionalPanel(
      condition = "input.tabs != 'Statistics'", ns = NS(NULL),
      hr(),
      h5("Export"),
      downloadButton(ns("dl_summary"), "Download summary CSV")
    )
  )
}

#' Main content UI for the Plot tab
#' Left card: scatter plot with controls. Right card: distribution comparison.
plot_main_ui <- function(id, extra_panel = NULL) {
  ns <- NS(id)
  tagList(
    # The phenotype-rules pill and the three blue hint strips that sat here are
    # in the rail now -- Help > "Phenotype rules" and Help > "Tips for this tab".
    # They took a full row above the cards on every visit; the cards start at the
    # top of the panel instead.
    # Left card: Main scatter plot. 7/5 rather than 6/6: this one carries two
    # slider columns and a checkbox column beneath its figure, and at 6/6 the
    # "Color points by" select wrapped onto three lines.
    column(7,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px",
        h5("Islet Size Distribution", style = "margin-top: 0"),
        div(style = "flex: 1; min-height: 0; display: flex",
          plotlyOutput(ns("plt"), height = "100%")
        ),
        div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px; overflow-y: auto; max-height: 35%; padding-top: 10px",
          div(style = "display: flex; gap: 18px; align-items: stretch",
            # Slider column 1
            div(style = "flex: 1; min-width: 0; display: flex; flex-direction: column",
              sliderInput(ns("binwidth"), "Diameter bin width (\u00b5m)",
                          min = 1, max = 75, value = 50, step = 1),
              sliderInput(ns("diam_range"), "Islet diameter range (\u00b5m)",
                          min = 0, max = 500, value = c(0, 350), step = 10),
              selectInput(ns("plot_color_by"), "Color points by:",
                          choices = c("Donor Status" = "donor_status",
                                      "Donor ID" = "donor_id"),
                          selected = "donor_status")
            ),
            # Slider column 2
            div(style = "flex: 1; min-width: 0; display: flex; flex-direction: column",
              sliderInput(ns("pt_size"), "Point size",
                          min = 0.3, max = 4.0, value = 3.0, step = 0.1),
              sliderInput(ns("pt_alpha"), "Point transparency",
                          min = 0, max = 0.95, value = 0.4, step = 0.05),
              radioButtons(ns("add_smooth"), "Trend line",
                           choices = c("None", "LOESS"),
                           selected = "None", inline = TRUE)
            ),
            # Checkbox column 3 \u2014 vertical stack spanning the full height of the sliders
            div(style = "flex: 0 0 auto; min-width: 165px; display: flex; flex-direction: column; justify-content: space-between; padding-left: 16px",
              checkboxInput(ns("exclude_zero_top"), "Exclude zero values", value = FALSE),
              checkboxInput(ns("show_points"), "Show individual points", value = FALSE),
              checkboxInput(ns("log_scale"), "Log scale y-axis", value = FALSE),
              checkboxInput(ns("log_scale_x"), "Log2 scale x-axis", value = FALSE)
            )
          )
        )
      )
    ),
    # Right card: Distribution comparison
    # Widen to fill the row when no extra panel (e.g. AI chat) is present.
    column(if (is.null(extra_panel)) 5 else 3,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px",
        h5("Distribution Comparison", style = "margin-top: 0"),
        div(style = "flex: 1; min-height: 0; display: flex",
          plotlyOutput(ns("dist"), height = "100%")
        ),
        div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px; overflow-y: auto; max-height: 35%; padding-top: 10px",
          uiOutput(ns("dist_ui"))
        )
      )
    ),
    # Extra panel (e.g., AI chat) injected between cards and full-width panels
    extra_panel,
    # Full-width card below both plot cards: segmentation/drill-down viewer
    column(12,
      uiOutput(ns("segmentation_viewer_panel"))
    )
  )
}
