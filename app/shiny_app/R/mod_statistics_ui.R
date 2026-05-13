# ---------- Statistics module UI ----------
# Exports: statistics_ui(id)
# 5-section narrative layout: Configure, Primary Results, Size-Dependent,
#                             Confounders & Deeper, Methods Reference

statistics_ui <- function(id) {

  ns <- NS(id)

  # Helper: numbered section heading with subtitle

  section_heading <- function(step, title, subtitle) {
    div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px; border-bottom: 2px solid #d0e0f0;",
      div(style = "display: flex; align-items: baseline; gap: 10px;",
        span(step,
             style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                            " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                            " border-radius: 12px; min-width: 24px; text-align: center;")),
        span(title, style = "font-weight: 700; font-size: 17px; color: #333;")
      ),
      p(subtitle, style = "margin: 4px 0 0 0; color: #777; font-size: 13px;")
    )
  }

  tagList(

    # ===== SECTION 1: Configure Analysis =====
    section_heading("1", "Configure Analysis",
                    "Select your test parameters, then click Run Statistics."),

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
          fluidRow(
            column(8,
              uiOutput(ns("overview_banner"))
            ),
            column(4, style = "text-align: right; padding-top: 5px;",
              actionButton(ns("run_tests"), "Run Statistics",
                           class = "btn-primary", style = "margin-right: 8px;"),
              downloadButton(ns("download_stats"), "Download CSV",
                             style = "margin-top: 0;")
            )
          ),
          hr(style = "margin: 10px 0;"),
          # Row 1: Analysis settings
          fluidRow(
            column(2,
              radioButtons(ns("test_type"), "Test Type",
                           choices = c("Parametric", "Non-parametric"),
                           selected = "Parametric", inline = FALSE),
              tags$small(style = "color: #888; display: block; margin-top: -8px; line-height: 1.3;",
                tags$strong("Parametric"), " (ANOVA / t-test): assumes roughly normal data within groups.",
                tags$br(),
                tags$strong("Non-parametric"), " (Kruskal-Wallis / Wilcoxon): no normality assumption; based on ranks."
              )
            ),
            column(1,
              selectInput(ns("alpha"), "\u03b1",
                          choices = c("0.05", "0.01", "0.001"), selected = "0.05")
            ),
            column(2,
              # Outlier toggle uses the synced global control in the Plot
              # sidebar (`sidebar_remove_outliers`). The text below mirrors
              # the rule applied so users don't need to switch tabs to confirm.
              uiOutput(ns("stats_outlier_status")),
              numericInput(ns("stats_min_cells"), "Min cells/islet",
                           value = 1, min = 1, max = 200, step = 1),
              checkboxInput(ns("stats_no_binning"), "Skip size stratification", value = FALSE)
            ),
            column(3,
              sliderInput(ns("bin_width"), "Per-bin width (\u00b5m)",
                          min = 5, max = 150, value = 50, step = 5)
            ),
            column(3,
              sliderInput(ns("stats_diam_range"), "Diameter range (\u00b5m)",
                          min = 0, max = 500, value = c(0, 350), step = 10)
            )
          ),
          fluidRow(
            column(12,
              uiOutput(ns("normality_result"))
            )
          )
        )
      )
    ),

    # ===== SECTION 2: Primary Results =====
    section_heading("2", "Primary Results",
                    "Donor-level tests (N=15) to avoid pseudoreplication, with mixed-effects sensitivity analysis."),

    uiOutput(ns("pseudorep_banner")),

    fluidRow(style = "display: flex; flex-wrap: wrap;",
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Hypothesis Testing"),
          tableOutput(ns("test_results_table"))
        )
      ),
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Effect Size Forest Plot"),
          plotlyOutput(ns("forest_plot"), height = "350px")
        )
      ),
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Area Under Curve by Donor Group"),
          plotlyOutput(ns("auc_plot"), height = "280px"),
          hr(style = "margin: 10px 0;"),
          tableOutput(ns("auc_table")),
          uiOutput(ns("auc_interpretation"))
        )
      )
    ),

    # ===== SECTION 3: Size-Dependent Patterns =====
    conditionalPanel(
      condition = sprintf("!input['%s']", ns("stats_no_binning")),

    section_heading("3", "Size-Dependent Patterns",
                    "Does the group effect depend on islet size? Each diameter bin is tested independently (donor-level means per bin)."),

    fluidRow(style = "display: flex; flex-wrap: wrap;",
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Stratified Tests by Islet Diameter"),
          tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
            "Tests use donor-level means within each size bin. ",
            "P-values are BH-corrected across bins. ",
            "Bins with <2 donors per group are greyed out (untestable)."),
          plotlyOutput(ns("bin_heatmap"), height = "400px")
        )
      ),
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Trend Analysis (Kendall \u03c4)"),
          tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
            "Each point is one diameter bin. Kendall \u03c4 measures the ",
            "correlation between disease stage (ND=0, Aab+=1, T1D=2) and the feature value. ",
            tags$strong("\u03c4 > 0:"), " feature increases ND \u2192 Aab+ \u2192 T1D. ",
            tags$strong("\u03c4 < 0:"), " feature decreases. ",
            tags$strong("\u03c4 \u2248 0:"), " no monotonic trend. ",
            "Bins with <3 donors are skipped."),
          plotlyOutput(ns("trend_plot"), height = "400px")
        )
      )
    )
    ), # end conditionalPanel for Section 3

    # ===== SECTION 4 (or 3 if binning skipped): Confounders & Deeper Analysis =====
    uiOutput(ns("section_confounders_heading")),

    fluidRow(
      column(12,
        uiOutput(ns("demographics_card"))
      )
    ),

    # ===== SECTION 5 (or 4 if binning skipped): Methods Reference =====
    uiOutput(ns("section_methods_heading")),

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; background: #f8f9fa;",
          uiOutput(ns("stats_explanation"))
        )
      )
    )
  )
}
