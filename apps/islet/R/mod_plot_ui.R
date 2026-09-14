# ---------- Plot module UI ----------
# Exports: plot_sidebar_ui(id), plot_main_ui(id)
# All Plot-related inputs are namespaced by this module.

#' Sidebar UI for the Plot tab
#' Contains data filters, autoantibody filters, statistic selector, and export controls.
#' Selection block for the shell's frozen control panel.
#'
#' WHAT IT KEEPS is what SELECTS the data this tab plots: the data layer, the
#' feature within it, the metric, the summary statistic, and the outlier table.
#' WHAT IT LOST went to the shared cohort filter in the top bar, which now owns
#' the same controls for every tab:
#'
#'   groups / age_filter_ui / gender_filter_ui / aab_flags / min_cells
#'   sidebar_remove_outliers / sidebar_outlier_z
#'
#' `sidebar_donor_palette` moved to the rail footer with the phenotype palette.
#' Three of the four conditionalPanels are gone with the controls they hid on
#' the Statistics tab -- Statistics shares this block, because it has always
#' consumed the Plot module's frame.
plot_sidebar_ui <- function(id) {
  ns <- NS(id)
  # No `.sidebar` class any more: this block renders inside the shell's frozen
  # control panel (`.ie-aside` in app.R), which owns the surface. The old class
  # carried a dark navy ground of its own and would fight it.
  div(
    div(class = "ie-aside-head", "This tab"),
    selectInput(ns("mode"), "Data Layer",
                choices = c("Microenvironment", "Cell Populations"),
                selected = "Cell Populations", width = "100%"),
    uiOutput(ns("dynamic_selector")),
    uiOutput(ns("metric_selector")),
    # The AAb warning stays: it is about THIS tab's selection being unavailable
    # for Aab+ donors, not about the filter, which is in the top bar.
    uiOutput(ns("aab_warning")),
    conditionalPanel(
      condition = "input.tabs != 'Statistics'", ns = NS(NULL),
      radioButtons(ns("stat"), "Statistic",
                   choices = c("Mean\u00b1SE" = "mean_se",
                               "Mean\u00b1SD" = "mean_sd",
                               "Median + IQR" = "median_iqr"),
                   selected = "mean_se"),
      checkboxInput(ns("show_outlier_table"),
                    "Show excluded-outlier table",
                    value = FALSE),
      uiOutput(ns("plot_outlier_info")),
      downloadButton(ns("dl_summary"), "Download summary CSV", class = "btn-block")
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
              # `diam_range` moved to the cohort panel. Leaving it here would
              # have been worse than a duplicate: with the panel mounted, this
              # module routes its row filtering through cohort$apply(), so this
              # slider would have moved nothing at all while still looking like
              # a filter -- and its fixed 0-500 default clipped at 350, hiding
              # every islet above that in a cohort whose largest is 457 um.
              # The panel derives its bounds from the data instead.
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
