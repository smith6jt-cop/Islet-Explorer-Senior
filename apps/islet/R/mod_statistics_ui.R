# ---------- Statistics module UI ----------
# Exports: statistics_ui(id)
# 5-section narrative layout: Configure, Primary Result, Exploratory & Size-Dependent,
#                             Confounders & Demographics, Methods Reference
#
# Redesign (Jul 2026) — streamlined to a single principled path to resist p-hacking:
#   - Section 1 exposes only Run / Download on the primary path; the outlier, zero-value,
#     min-cells, and diameter-range knobs live inside a collapsed "Advanced / sensitivity
#     settings" <details> block so tuning them is a deliberate, off-path act.
#   - Test type is AUTO-selected from the Shapiro-Wilk normality result (no manual toggle);
#     alpha is FIXED at 0.05 (no selector). Both were classic "move until significant" levers.
#   - Section 2 shows one mixed-effects omnibus p + the effect-size forest plot; the full
#     hypothesis table and donor-level sensitivity test are collapsed under "Sensitivity check".
#   - AUC + size-stratified per-bin analysis are consolidated into Section 3 (fixed 50 µm bin
#     width), framed as exploratory / hypothesis-generating, not confirmatory.

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
        tags$strong("ND, Aab+, and T1D"), " donors (N=15 total) with a single, pre-specified ",
        "workflow that avoids pseudoreplication. The primary test is a ",
        tags$strong("linear mixed-effects model"), " with donor as a random intercept. ",
        "A conservative donor-level cross-check is reported as a sensitivity analysis, and its ",
        "parametric vs. non-parametric form is ", tags$strong("chosen automatically"),
        " from the normality of the donor means — not left as a knob to flip."),
      p(style = "margin: 0 0 6px 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("What you get."), " One headline result (mixed-effects p + ICC) and an ",
        "effect-size forest plot up front. Secondary tests, size-stratified per-bin patterns, ",
        "trapezoidal AUC, and confounder checks (age, sex, AAb count, covariate-adjusted model) ",
        "follow as clearly-labelled ", tags$em("exploratory"), " material."),
      p(style = "margin: 0; font-size: 13px; line-height: 1.5; color: #333;",
        tags$strong("How to use it."), " Pick a feature in the sidebar (shared with the Plot tab), ",
        "then click ", tags$strong("Run Statistics"), ". Significance is judged at a fixed ",
        tags$strong("α = 0.05"), ". Data-subset knobs (outliers, zero-values, min-cells, diameter ",
        "range) live under ", tags$em("Advanced / sensitivity settings"),
        " — use them to probe robustness, not to hunt for a p-value. Hover the ",
        tags$em("(?)"), " icon on each section for a one-line description.")
    ),

    # ===== SECTION 1: Configure Analysis =====
    section_heading("1", "Configure Analysis",
                    "Pick a feature in the sidebar, then click Run Statistics.",
                    tip = paste0("The primary test and α = 0.05 are fixed; the sensitivity ",
                                 "test is auto-selected from normality. Optional data-subset ",
                                 "knobs live under Advanced / sensitivity settings.")),

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
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

          # Normality readout — drives the AUTO test-type selection (read-only).
          uiOutput(ns("normality_result")),

          # ---- Advanced / sensitivity settings (collapsed by default) ----
          # These are data-subset knobs that each change which islets are tested.
          # Deliberately off the primary path: opening this block is an explicit
          # act of probing robustness, not part of the default analysis. Same
          # input IDs as before, so the server reads them unchanged.
          tags$details(
            style = paste0("margin-top: 12px; border: 1px solid #e0e0e0; border-radius: 6px; ",
                           "padding: 8px 12px; background: #fbfcfd;"),
            tags$summary(
              style = "cursor: pointer; font-weight: 600; font-size: 13px; color: #4477AA;",
              "Advanced / sensitivity settings ",
              tags$span(style = "font-weight: 400; color: #999; font-size: 12px;",
                        "— optional data-subset filters for robustness checks")),
            tags$p(style = "font-size: 12px; color: #888; margin: 8px 0 12px 0; line-height: 1.4;",
              "Each filter changes which islets enter every test below. Use them to check ",
              "whether a result holds on a cleaned or restricted subset — not to search for ",
              "significance. Whatever you set here is recorded in the downloaded CSV header."),
            fluidRow(
              column(3,
                checkboxInput(ns("stats_remove_outliers"),
                              "Exclude outliers",
                              value = FALSE),
                sliderInput(ns("stats_outlier_z"), "Outlier z-threshold",
                            min = 0.5, max = 10, value = 3, step = 0.5),
                tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                  "Robust median/MAD |z| per donor status. Excluded points drop from every test."),
                uiOutput(ns("stats_outlier_count"))
              ),
              column(3,
                checkboxInput(ns("stats_exclude_zero"),
                              "Exclude zero-valued islets",
                              value = FALSE),
                tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                  "Filters value == 0 before all tests. ",
                  "Useful for zero-inflated features. ",
                  "AUC is recomputed on the same subset, so it moves too.")
              ),
              column(3,
                numericInput(ns("stats_min_cells"), "Min cells/islet",
                             value = 1, min = 1, max = 200, step = 1),
                tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                  "Drops islets with fewer than this many core+peri cells. ",
                  "Applies in every mode, including Microenvironment. If the cell ",
                  "counts are unavailable the overview banner says so rather than ",
                  "ignoring the value silently.")
              ),
              column(3,
                sliderInput(ns("stats_diam_range"), "Diameter range (µm)",
                            min = 0, max = 500, value = c(0, 350), step = 10,
                            width = "100%"),
                tags$small(style = "color: #888; display: block; margin-top: -6px; line-height: 1.3;",
                  "Restrict islets to a diameter window before testing. Applied on ",
                  "top of the Plot sidebar's own diameter range, so it can only ",
                  "narrow further \u2014 widening past that range has no effect. The ",
                  "inherited range is listed in the overview banner above.")
              )
            )
          )
        )
      )
    ),

    # ===== SECTION 2: Primary Result =====
    section_heading("2", "Primary Result",
                    "One pre-specified test: a mixed-effects model with donor random intercept. Effect sizes are foregrounded over p-values.",
                    tip = paste0("The headline is the mixed-effects omnibus p (fixed α = 0.05) with ICC. ",
                                 "The forest plot shows pairwise Cohen's d with 95% CI. Secondary and ",
                                 "donor-level tests are collapsed under the Sensitivity check.")),

    uiOutput(ns("pseudorep_banner")),

    fluidRow(style = "display: flex; flex-wrap: wrap;",
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Primary Test — Mixed-Effects Model"),
          uiOutput(ns("primary_headline"))
        )
      ),
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          h5("Effect Size Forest Plot"),
          tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
            "Pairwise difference between donor means, in the feature's own units, with a ",
            "Welch 95% CI. Magnitude and its uncertainty carry more information than a ",
            "p-value alone, especially with N=5 per group. The unstandardised difference ",
            "leads because a pooled-SD effect size (Cohen's d) needs the groups to have ",
            "comparable spread, and they often do not here \u2014 d and the Glass's \u0394 range it ",
            "sits inside are in the Sensitivity check table."),
          plotlyOutput(ns("forest_plot"), height = "430px")
        )
      )
    ),

    # ---- Sensitivity check (collapsed): full hypothesis table + donor-level tests ----
    fluidRow(
      column(12,
        tags$details(
          style = paste0("margin-bottom: 15px; border: 1px solid #e0e0e0; border-radius: 6px; ",
                         "padding: 8px 14px; background: #fbfcfd;"),
          tags$summary(
            style = "cursor: pointer; font-weight: 600; font-size: 13px; color: #4477AA;",
            "Sensitivity check — secondary & donor-level tests ",
            tags$span(style = "font-weight: 400; color: #999; font-size: 12px;",
                      "— cross-checks for the primary result; interpret alongside, not instead of, it")),
          tags$p(style = "font-size: 12px; color: #888; margin: 8px 0 10px 0; line-height: 1.4;",
            "These tests re-examine the same data from more conservative angles (donor means, ",
            "mixed-effects pairwise contrasts). They exist to confirm the primary result is ",
            "robust — not as an alternative menu to pick the smallest p-value from."),
          h5("Hypothesis Testing"),
          tableOutput(ns("test_results_table"))
        )
      )
    ),

    # ===== SECTION 3: Exploratory & Size-Dependent Patterns =====
    section_heading("3", "Exploratory & Size-Dependent Patterns",
                    "Descriptive, hypothesis-generating views: integrated AUC and per-size-bin tests. Not confirmatory.",
                    tip = paste0("AUC integrates the binned per-status means. Per-bin ANOVA and ",
                                 "Kendall τ use donor means within fixed 50 µm diameter bins, ",
                                 "BH-corrected across bins. These are exploratory and are NOT corrected ",
                                 "across features — treat as hypothesis-generating.")),

    div(style = paste0("background: #fff8e1; border-left: 4px solid #f0ad4e; border-radius: 4px; ",
                       "padding: 10px 14px; margin-bottom: 12px; font-size: 13px; line-height: 1.5; color: #6a5200;"),
      tags$strong("Exploratory — not confirmatory. "),
      "The views below scan many diameter bins and integrate a summary curve. P-values are ",
      "BH-corrected within each per-bin family but ", tags$strong("not"),
      " across the many features you might try. Use them to generate hypotheses, then confirm ",
      "with the Section 2 primary test on a pre-specified feature."),

    checkboxInput(ns("stats_show_exploratory"),
                  "Show exploratory analyses", value = TRUE),

    conditionalPanel(
      condition = sprintf("input['%s']", ns("stats_show_exploratory")),

      fluidRow(style = "display: flex; flex-wrap: wrap;",
        column(4,
          div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
            h5("Area Under Curve by Donor Group"),
            tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
              "Trapezoidal integral of the binned mean feature value across islet diameters, ",
              "per group. Descriptive summary — no inferential test."),
            plotlyOutput(ns("auc_plot"), height = "380px"),
            hr(style = "margin: 10px 0;"),
            tableOutput(ns("auc_table")),
            uiOutput(ns("auc_interpretation"))
          )
        ),
        column(4,
          div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
            h5("Stratified Tests by Islet Diameter"),
            tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
              "Donor-level means within each fixed 50 µm size bin, BH-corrected across bins. ",
              "Hover a cell for the per-group coverage behind that specific test. ",
              "A bin needs \u2265 2 groups with \u2265 2 donors each; the same rule gates the ANOVA, ",
              "Kendall \u03c4 and the pairwise rows, so no pairwise contrast appears without an ",
              "omnibus behind it. Counts per bin, and any bins that failed the rule, are ",
              "listed below the plot."),
            plotlyOutput(ns("bin_heatmap"), height = "520px"),
            uiOutput(ns("bin_coverage_note"))
          )
        ),
        column(4,
          div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
            h5("Trend Analysis (Kendall τ)"),
            tags$small(style = "color: #888; display: block; margin-bottom: 8px;",
              "Each point is one 50 µm diameter bin. Kendall τ measures the ",
              "correlation between disease stage (ND=0, Aab+=1, T1D=2) and the feature value. ",
              tags$strong("τ > 0:"), " increases ND → Aab+ → T1D. ",
              tags$strong("τ < 0:"), " decreases. ",
              tags$strong("τ ≈ 0:"), " no monotonic trend. ",
              "Same coverage rule as the heatmap."),
            plotlyOutput(ns("trend_plot"), height = "400px")
          )
        )
      )
    ),

    # ===== SECTION 4: Confounders & Demographics =====
    section_heading("4", "Confounders & Demographics",
                    "Do demographics explain the effect? Age, sex, autoantibody count, and a covariate-adjusted donor-level model.",
                    tip = paste0("Demographic confounders (age, sex, AAb count) and a covariate-adjusted ",
                                 "donor-level linear model. Sex-stratified pairwise tests are BH-corrected ",
                                 "across the 6 contrasts (3 pairs x 2 sexes). Exploratory.")),

    fluidRow(
      column(12,
        uiOutput(ns("demographics_card"))
      )
    ),

    # ===== SECTION 5: Methods Reference =====
    section_heading("5", "Methods Reference",
                    "Technical details about the primary test, auto-selected sensitivity test, corrections, and guidelines.",
                    tip = paste0("Full methods text. Updates dynamically with the auto-selected sensitivity ",
                                 "test, outlier rule, zero-exclusion state, and the fixed α / bin width.")),

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; background: #f8f9fa;",
          uiOutput(ns("stats_explanation"))
        )
      )
    )
  )
}
