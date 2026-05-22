# ---------- Plot module UI ----------
# Exports: plot_sidebar_ui(id), plot_main_ui(id)
# All Plot-related inputs are namespaced by this module.

#' Sidebar UI for the Plot tab
#' Contains data filters, autoantibody filters, statistic selector, and export controls.
plot_sidebar_ui <- function(id) {
  ns <- NS(id)
  div(class = "sidebar", style = "height: 100%;",
    h4("Data & Filters"),
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
    column(10,
      phenotype_rules_button("show_phenotype_rules_plot"),
      div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
        div(style = "flex: 1; background-color: #e8f4fd; border: 1px solid #b8daff; border-radius: 5px; padding: 10px 12px; color: #004085; font-size: 13px;",
          "Enable \"Show individual points\" then click any point on the ",
          tags$strong("Islet Size Distribution"), " scatter to view its segmentation map below."
        ),
        div(style = "flex: 1; background-color: #e8f4fd; border: 1px solid #b8daff; border-radius: 5px; padding: 10px 12px; color: #004085; font-size: 13px;",
          "Enable \"Show individual points\" and set \"Color points by\" to Donor ID, then click a donor in the legend to show/hide it."
        ),
        div(style = "flex: 1; background-color: #e8f4fd; border: 1px solid #b8daff; border-radius: 5px; padding: 10px 12px; color: #004085; font-size: 13px;",
          "Use browser zoom (Ctrl +/\u2212) to resize cards. Scroll within cards to reveal hidden controls."
        )
      )
    ),
    # Left card: Main scatter plot
    column(6,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px;",
        h5("Islet Size Distribution", style = "margin-top: 0; color: #333;"),
        div(style = "flex: 1; min-height: 0; display: flex;",
          plotlyOutput(ns("plt"), height = "100%")
        ),
        div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px; overflow-y: auto; max-height: 35%; border-top: 1px solid #eee; padding-top: 10px;",
          fluidRow(
            column(6,
              sliderInput(ns("binwidth"), "Diameter bin width (\u00b5m)",
                          min = 1, max = 75, value = 50, step = 1),
              sliderInput(ns("diam_range"), "Islet diameter range (\u00b5m)",
                          min = 0, max = 500, value = c(0, 350), step = 10),
              selectInput(ns("plot_color_by"), "Color points by:",
                          choices = c("Donor Status" = "donor_status",
                                      "Donor ID" = "donor_id"),
                          selected = "donor_status")
            ),
            column(6,
              sliderInput(ns("pt_size"), "Point size",
                          min = 0.3, max = 4.0, value = 3.0, step = 0.1),
              sliderInput(ns("pt_alpha"), "Point transparency",
                          min = 0, max = 0.95, value = 0.4, step = 0.05),
              radioButtons(ns("add_smooth"), "Trend line",
                           choices = c("None", "LOESS"),
                           selected = "None", inline = TRUE)
            )
          ),
          div(style = "display: flex; flex-wrap: wrap; gap: 15px; justify-content: flex-end;",
            checkboxInput(ns("exclude_zero_top"), "Exclude zero values", value = FALSE),
            checkboxInput(ns("show_points"), "Show individual points", value = FALSE),
            checkboxInput(ns("log_scale"), "Log scale y-axis", value = FALSE),
            checkboxInput(ns("log_scale_x"), "Log2 scale x-axis", value = FALSE)
          )
        )
      )
    ),
    # Right card: Distribution comparison
    column(4,
      div(class = "card",
          style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px;",
        h5("Distribution Comparison", style = "margin-top: 0; color: #333;"),
        div(style = "flex: 1; min-height: 0; display: flex;",
          plotlyOutput(ns("dist"), height = "100%")
        ),
        div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px; overflow-y: auto; max-height: 35%; border-top: 1px solid #eee; padding-top: 10px;",
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
