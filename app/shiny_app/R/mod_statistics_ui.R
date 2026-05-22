# ---------- Statistics module UI ----------
# Exports: statistics_ui(id)
# 5-section narrative layout: Configure, Primary Results, Size-Dependent,
#                             Confounders & Deeper, Methods Reference

statistics_ui <- function(id) {

  ns <- NS(id)

  # Helper: numbered section heading with subtitle and optional tooltip
  # `tip` (NULL or character): if supplied, renders a small "(?)" icon next to
  # the title with a native browser tooltip describing the section's contents.

  section_heading <- function(step, title, subtitle, tip = NULL) {
    tip_icon <- if (!is.null(tip) && nzchar(tip)) {
      span("?", title = tip,
           style = paste0("display: inline-block; margin-left: 6px; width: 18px; height: 18px;",
                          " line-height: 18px; text-align: center; font-size: 12px; font-weight: 700;",
                          " color: #4477AA; background: #eaf2fa; border: 1px solid #c0d4e8;",
                          " border-radius: 50%; cursor: help;"))
    } else NULL
    div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px; border-bottom: 2px solid #d0e0f0;",
      div(style = "display: flex; align-items: baseline; gap: 10px;",
        span(step,
             style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                            " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                            " border-radius: 12px; min-width: 24px; text-align: center;")),
        span(title, style = "font-weight: 700; font-size: 17px; color: #333;"),
        tip_icon
      ),
      p(subtitle, style = "margin: 4px 0 0 0; color: #777; font-size: 13px;")
    )
  }

  tagList(

    phenotype_rules_button("show_phenotype_rules_stats"),

    # ===== Top: General-use description banner =====
    div(class = "card",
        style = paste0("padding: 14px 18px; margin-bottom: 18px; ",
                       "background: linear-gradient(135deg, #f4f8fc 0%, #eaf2fa 100%); ",
                       "border-left: 4px solid #4477AA;"),
      h4("About this tab",
         style = "margin: 0 0 6px 0; font-size: 16px; color: #2c4760;"),
      p(style = "margin: 0 0 6px 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("What it does."), " Compares the selected feature across ",
        tags$strong("ND, Aab+, and T1D"), " donors (N=15 total) using a donor-aware workflow ",
        "that avoids pseudoreplication. The primary test is a ",
        tags$strong("linear mixed-effects model"), " with donor as a random intercept; ",
        "a donor-level ANOVA / Kruskal-Wallis on per-donor means is reported as a conservative cross-check."),
      p(style = "margin: 0 0 6px 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("What you get."), " Hypothesis testing table, effect-size forest plot, ",
        "trapezoidal AUC, size-stratified per-bin tests, Kendall τ trend across disease stage, ",
        "and confounder checks (age, sex, AAb count, covariate-adjusted model)."),
      p(style = "margin: 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("How to use it."), " Pick a feature in the sidebar (shared with the Plot tab). ",
        "Choose Parametric or Non-parametric and an α, optionally toggle ",
        tags$em("Exclude zero-valued islets"), " for zero-inflated features, then click ",
        tags$strong("Run Statistics"), ". Hover the ", tags$em("(?)"), " icon on each section for a one-line description.")
    ),

    # ===== SECTION 1: Configure Analysis =====
    section_heading("1", "Configure Analysis",
                    "Select your test parameters, then click Run Statistics.",
                    tip = paste0("Pick test type, α, min-cells, diameter range, and optional ",
                                 "zero exclusion. Click Run Statistics to rebuild every output below.")),

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
          # Row 1: Controls, evenly spaced across the card.
          # Each control gets its own column; long captions sit under each input.
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
            column(2,
              selectInput(ns("alpha"), "\u03b1",
                          choices = c("0.05", "0.01", "0.001"), selected = "0.05"),
              numericInput(ns("stats_min_cells"), "Min cells/islet",
                           value = 1, min = 1, max = 200, step = 1)
            ),
            column(2,
              # Local outlier toggle + threshold slider (independent of the
              # Plot / Trajectory / Spatial global controls). Lets users
              # compare with and without outlier removal without switching
              # tabs.
              checkboxInput(ns("stats_remove_outliers"),
                            "Exclude outliers",
                            value = FALSE),
              sliderInput(ns("stats_outlier_z"), "Outlier z-threshold",
                          min = 0.5, max = 10, value = 3, step = 0.5),
              tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                "Local to Statistics: ", tags$strong("z above the slider value"),
                " per donor status. Excluded points are dropped from every test below."),
              uiOutput(ns("stats_outlier_count"))
            ),
            column(2,
              checkboxInput(ns("stats_no_binning"),
                            "Skip size stratification",
                            value = FALSE),
              tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                "Hide the per-bin (size-dependent) section.")
            ),
            column(2,
              checkboxInput(ns("stats_exclude_zero"),
                            "Exclude zero-valued islets",
                            value = FALSE),
              tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                "Filters value == 0 before all tests. ",
                "Useful for zero-inflated features. ",
                "AUC reads the Plot summary and is not affected.")
            ),
            column(2,
              # Reserved for future analytics summary / placeholder; keeps the
              # row at 12 columns and prevents the last control from feeling
              # orphaned next to a wide blank.
              tags$div()
            )
          ),
          # Row 2: Sliders (need width)
          fluidRow(
            column(6,
              sliderInput(ns("bin_width"), "Per-bin width (\u00b5m)",
                          min = 5, max = 150, value = 50, step = 5,
                          width = "100%")
            ),
            column(6,
              sliderInput(ns("stats_diam_range"), "Diameter range (\u00b5m)",
                          min = 0, max = 500, value = c(0, 350), step = 10,
                          width = "100%")
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
                    "Donor-level tests (N=15) to avoid pseudoreplication, with mixed-effects sensitivity analysis.",
                    tip = paste0("Mixed-effects model (donor random intercept) is the primary test. ",
                                 "Donor-level ANOVA / Kruskal on per-donor means is the conservative ",
                                 "cross-check. Forest plot shows pairwise Cohen's d with 95% CI. ",
                                 "AUC integrates the binned per-status means.")),

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
                    "Does the group effect depend on islet size? Each diameter bin is tested independently (donor-level means per bin).",
                    tip = paste0("Per-bin ANOVA and Kendall τ on donor means within each diameter bin. ",
                                 "BH-corrected across bins. A bin is tested only if at least 2 donor-status ",
                                 "groups each have >= 2 donors.")),

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
            "A bin is tested only if at least 2 donor-status groups each have >= 2 donors ",
            "(same rule as the heatmap)."),
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
