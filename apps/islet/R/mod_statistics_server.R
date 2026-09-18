# ---------- Statistics module server ----------
# Exports: statistics_server(id, raw_df, get_selection_description)
#
# NOTE: the module no longer takes the Plot module's `summary_df`. AUC used to be
# integrated from it, which made that one card obey the Plot sidebar's bin-width
# and Mean/Median controls (hidden on this tab) while ignoring every filter set
# here. AUC is now computed from this tab's own filtered frame at the fixed
# STATS_BIN_WIDTH_UM. Do not reintroduce a cross-tab data dependency.
#
# Phase 16: Donor-level aggregation to fix pseudoreplication.
# Tests use donor-level means (N~15) instead of islet-level (N~5,214).
# Mixed-effects model (lmerTest) as sensitivity analysis.
#
# Dependencies (auto-sourced):
#   summary_stats, per_bin_anova, per_bin_kendall   -- utils_stats.R
#   aggregate_to_donor, lmer_test_donor              -- utils_stats.R
#   per_bin_donor_anova, per_bin_donor_kendall       -- utils_stats.R
#   normality_tests                                  -- utils_stats.R
#   cohens_d, eta_squared                            -- utils_stats.R
#   bin_islet_sizes                                  -- data_loading.R

# Fixed per-bin width for the Section 3 size-stratified analysis. Intentionally a
# constant, not a slider: a tunable bin width is a multiple-comparison p-hacking
# surface (retune the bins until one is significant). Referenced by stats_run() and
# the methods text / CSV header.
STATS_BIN_WIDTH_UM <- 50L

#' @param cohort the shared cohort filter from cohortFilterServer(), or NULL.
#'   NULL keeps this module's own inputs, which is what the testthat suite and
#'   any standalone use get; app.R always passes one.
statistics_server <- function(id, raw_df, get_selection_description,
                              donor_colors = reactive(DONOR_COLORS),
                              remove_outliers = reactive(TRUE),
                              cohort = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # (Section headings are now static in mod_statistics_ui.R after the Jul 2026
    #  redesign — no server-side heading helper or dynamic section offset needed.
    #  Section numbering is a fixed 1–5.)

    # Compact caption beneath the local outlier toggle: shows how many islets
    # the most recent Run excluded, or an "OFF" notice when the toggle is off.
    # The toggle itself (`stats_remove_outliers`) lives in the UI and is local
    # to the Statistics tab — independent of the Plot / Trajectory / Spatial
    # synced control.
    output$stats_outlier_count <- renderUI({
      filter_on <- isTRUE(input$stats_remove_outliers)
      if (!filter_on) {
        return(tags$div(
          style = paste0("background: #fdfbe7; border: 1px solid #f0d97a;",
                         " border-radius: 4px; padding: 4px 6px; margin-top: 6px;",
                         " font-size: 11px; color: #6e6300;"),
          "Filter OFF — extreme islets contribute to every test."))
      }
      st <- tryCatch(stats_run(), error = function(e) NULL)
      n_excl <- if (!is.null(st) && !is.null(st$outlier_summary)) st$outlier_summary$n_excluded else NA_integer_
      label <- if (!is.na(n_excl)) sprintf("%d islet%s excluded last Run.",
                                            n_excl, if (n_excl == 1) "" else "s")
               else "Run Statistics to see how many islets are excluded."
      tags$div(
        style = paste0("margin-top: 6px; font-size: 11px; color: #2e7d32;",
                       " font-style: italic;"),
        label)
    })

    # ===========================================================================
    # Core computation: eventReactive triggered by Run button
    # ===========================================================================

    # Why the last Run produced no result. `stats_run()` returns NULL when the
    # filters leave nothing to test, and every card then falls back to its
    # pre-run placeholder ("click Run Statistics") -- indistinguishable from
    # never having clicked Run. This carries the reason so the banner can say
    # which filter emptied the frame. Cleared on every successful run.
    stats_empty_reason <- reactiveVal(NULL)

    stats_run <- eventReactive(input$run_tests, {
      empty <- function(reason) {
        stats_empty_reason(reason)
        showNotification(paste("No results:", reason), type = "warning", duration = 8)
        NULL
      }
      stats_empty_reason(NULL)

      rdf <- raw_df()
      if (is.null(rdf) || !nrow(rdf)) {
        return(empty(paste("the Plot sidebar's filters returned no islets, so this tab",
                           "received an empty data set. Widen them and Run again.")))
      }

      # Alpha is FIXED at 0.05 — a user-movable threshold is a p-hacking lever
      # (lower it post hoc to reframe a null, or to move a goalpost). Documented in
      # the methods text and CSV header.
      alpha_num <- 0.05
      # Test type (parametric vs non-parametric) is AUTO-selected from the donor-level
      # normality result — see `is_parametric` defined right after normality_tests()
      # below. The manual Parametric/Non-parametric toggle was removed because flipping
      # it to chase significance is a classic garden-of-forking-paths move.

      # ---- Zero-value exclusion (sensitivity for zero-inflated features) ----
      # Removes rows whose feature value == 0 before outlier flagging, donor-level
      # aggregation, mixed-effects modelling, demographics, and per-bin tests.
      # AUC is recomputed from the same filtered frame, so it IS affected.
      #
      # ORDER MATTERS: this runs BEFORE flag_outliers() so the robust centre and
      # scale (median / MAD) are computed on the values that actually enter the
      # tests. Flagging first would derive the MAD from a zero-inflated vector
      # and then discard the zeros — a different, and wrong, outlier set. This
      # also matches the Plot tab, which filters zeros in raw_df() before
      # plot_df() calls flag_outliers().
      #
      # `raw_df` here is the Plot module's `raw_df_prezero` (see app.R), so this
      # toggle is the ONLY zero filter applied to the statistics.
      exclude_zero <- isTRUE(input$stats_exclude_zero)
      zero_dropped <- 0L
      if (exclude_zero) {
        n_before <- nrow(rdf)
        rdf <- rdf[is.finite(rdf$value) & rdf$value != 0, , drop = FALSE]
        zero_dropped <- as.integer(n_before - nrow(rdf))
      }
      if (nrow(rdf) == 0) {
        return(empty(sprintf(paste("excluding zero-valued islets removed all %s of them.",
                                   "Untick 'Exclude zero values' in Advanced / sensitivity",
                                   "settings."), format(zero_dropped, big.mark = ","))))
      }

      # ---- Outlier removal via shared helper (controlled by LOCAL toggle) ----
      # `stats_remove_outliers` lives in the Statistics tab Section 1 controls
      # and is independent of the Plot / Trajectory / Spatial synced toggle.
      # The helper appends is_outlier + outlier_score columns; flagged rows
      # are dropped when the local filter is ON. The count is saved on the
      # result list so the methods text and download CSV can cite it.
      stats_filter_on <- isTRUE(input$stats_remove_outliers)
      stats_z_threshold <- {
        v <- suppressWarnings(as.numeric(input$stats_outlier_z))
        if (is.finite(v) && v >= 0.5) v else 3.0
      }
      rdf <- flag_outliers(rdf, value_col = "value", group_col = "donor_status",
                           method = "zscore", threshold = stats_z_threshold,
                           enabled = stats_filter_on)
      # NOTE the ordering. Flagging stays HERE, on the pre-min-cells/pre-diameter
      # population, so the z-scores (median/MAD per donor_status) are unchanged
      # by those filters -- moving it would silently change which rows are
      # flagged. But the *drop* and the *count* are deferred until after them:
      # `outlier_summary` used to be measured at this point and then reported as
      # "N islets excluded last Run", counting exclusions against a superset of
      # what was actually tested. Summarising post-filter makes the reported
      # count refer to the analysed set. See `outlier_summary` below.

      # ---- Min-cells filter (like spatial tab) ----
      # `%||%` only catches NULL; clearing the numericInput sends NA, and
      # `if (NA > 1 && TRUE)` aborts the whole eventReactive with "missing value
      # where TRUE/FALSE needed", blanking every card on the tab. Coerce to a
      # usable number first.
      min_cells <- {
        v <- suppressWarnings(as.numeric(input$stats_min_cells %||% 1))
        if (length(v) == 1 && is.finite(v) && v >= 1) v else 1
      }
      # `total_cells_*` are joined onto BOTH comp and targets_all by prep_data(),
      # so this filter now applies in Microenvironment mode too. The guard stays
      # for the Excel fallback path, where the neighborhood sheet is absent —
      # but it is no longer silent: see `min_cells_applied` below.
      min_cells_applied <- !(min_cells > 1) ||
        all(c("total_cells_core", "total_cells_peri") %in% colnames(rdf))
      if (min_cells > 1 &&
          "total_cells_core" %in% colnames(rdf) &&
          "total_cells_peri" %in% colnames(rdf)) {
        tc_core <- suppressWarnings(as.numeric(rdf$total_cells_core))
        tc_peri <- suppressWarnings(as.numeric(rdf$total_cells_peri))
        tc_core <- ifelse(is.finite(tc_core), tc_core, 0)
        tc_peri <- ifelse(is.finite(tc_peri), tc_peri, 0)
        rdf <- rdf[is.finite(tc_core + tc_peri) & (tc_core + tc_peri) >= min_cells, , drop = FALSE]
      }

      # ---- Diameter range filter (from Stats inline controls) ----
      # NOTE: this is `stats_diam_range`, NOT the Plot sidebar's `diam_range`.
      # Recorded on the result list so the CSV header can cite the range that
      # was actually applied (see the download handler).
      diam_range_used <- NULL
      if (!is.null(input$stats_diam_range) && length(input$stats_diam_range) == 2) {
        dmin <- as.numeric(input$stats_diam_range[1])
        dmax <- as.numeric(input$stats_diam_range[2])
        if (is.finite(dmin) && is.finite(dmax)) {
          diam_range_used <- c(dmin, dmax)
          rdf <- rdf[is.finite(rdf$islet_diam_um) &
                      rdf$islet_diam_um >= dmin &
                      rdf$islet_diam_um <= dmax, , drop = FALSE]
        }
      }

      # Check emptiness BEFORE the outlier drop so the message can name the
      # filter that actually emptied the frame, and so summarize_outlier_filter()
      # is never handed a zero-row data frame.
      if (nrow(rdf) == 0) {
        bits <- c(if (min_cells > 1) sprintf("min cells/islet = %g", min_cells),
                  if (!is.null(diam_range_used))
                    sprintf("diameter %g-%g um", diam_range_used[1], diam_range_used[2]))
        return(empty(paste0("no islets passed the Advanced / sensitivity filters",
                            if (length(bits)) paste0(" (", paste(bits, collapse = "; "), ")"),
                            ". Relax them and Run again.")))
      }

      # ---- Apply the deferred outlier drop (see the flagging block above) ----
      # Summarise first, so the count describes the rows that survived every
      # other filter -- i.e. the population actually tested.
      outlier_summary_pre <- summarize_outlier_filter(rdf)
      if (stats_filter_on) {
        rdf <- rdf[!rdf$is_outlier | is.na(rdf$is_outlier), , drop = FALSE]
      }

      if (nrow(rdf) == 0) {
        return(empty(sprintf(paste("the outlier filter (z > %s) removed every remaining",
                                   "islet. Raise the threshold or untick 'Remove outliers'."),
                             format(stats_z_threshold))))
      }

      rdf$donor_status <- factor(rdf$donor_status, levels = c("ND", "Aab+", "T1D"))
      n_islets_filtered <- nrow(rdf)

      # Detect feature name from raw_df attributes.
      # `selection_used` is the RAW column encoding ("prop_Beta cell",
      # "Capillary|islet_core"), which used to reach plot axes and card titles
      # verbatim. `feature_display_label()` resolves it to a readable name and
      # returns anything it does not recognise unchanged. The raw form is kept
      # as `feature_col` for the CSV, where the exact column matters.
      feature_col  <- attr(raw_df(), "selection_used") %||% "Unknown"
      feature_name <- feature_display_label(feature_col)

      # ---- Donor-level aggregation (conservative sensitivity) ----
      ddf <- aggregate_to_donor(rdf, agg_fn = "mean")
      ddf$donor_status <- factor(ddf$donor_status, levels = c("ND", "Aab+", "T1D"))
      n_donors <- nrow(ddf)

      # ---- Normality test on donor-level data ----
      norm_result <- normality_tests(ddf)

      # ---- AUTO test-type selection (replaces the removed Parametric/Non-parametric
      #      toggle) ----
      # Parametric only when every group with a computable Shapiro-Wilk test passes;
      # otherwise fall back to non-parametric. NA groups (too few donors) are ignored,
      # and an all-NA result defaults to non-parametric — the conservative choice when
      # normality cannot be established. This is the single source of truth for
      # is_parametric; the Section 1 normality readout reports the same decision.
      is_parametric <- {
        norm_flags <- norm_result$is_normal
        if (all(is.na(norm_flags))) FALSE else isTRUE(all(norm_flags, na.rm = TRUE))
      }

      # ---- PRIMARY: Mixed-effects model (preserves islet variation, accounts for donor clustering) ----
      lmer_result <- tryCatch(lmer_test_donor(rdf), error = function(e) {
        list(p_value = NA_real_, icc = NA_real_, fit = NULL,
             error_msg = paste("Unexpected error:", e$message))
      })
      lmer_p <- if (!is.null(lmer_result)) lmer_result$p_value else NA_real_
      icc <- if (!is.null(lmer_result)) lmer_result$icc else NA_real_
      lmer_error <- lmer_result$error_msg  # NULL if no error

      # Pairwise from mixed model via emmeans
      #
      # CRITICAL: emmeans caps its degrees-of-freedom machinery at
      # `pbkrtest.limit` / `lmerTest.limit` = 3000 observations. Above that it
      # silently drops to ASYMPTOTIC inference (z-tests, df = Inf), which for a
      # 15-donor design is badly anti-conservative — the same feature gave
      # p(ND-T1D) = 0.00007 at 5,212 islets (z) versus 0.0027 at 2,874 islets
      # (t, df ~ 12). The method therefore flipped on dataset size alone, with
      # nothing on screen to say so.
      #
      # Fix: request Satterthwaite explicitly with the limit raised past the row
      # count. Satterthwaite is also what lmerTest::anova() uses for the omnibus
      # p in the headline, so the pairwise contrasts are now consistent with it
      # (Kenward-Roger, emmeans' default, is far more expensive at this N).
      # `df_method` is recorded so the methods text can state it and so a silent
      # regression to df = Inf is visible rather than invisible.
      lmer_pairwise <- NULL
      lmer_pairwise_df_method <- NA_character_
      if (!is.null(lmer_result) && !is.null(lmer_result$fit)) {
        n_obs <- nrow(rdf)
        lmer_pairwise <- tryCatch({
          emm <- emmeans::emmeans(lmer_result$fit, pairwise ~ donor_status,
                                  lmer.df = "satterthwaite",
                                  lmerTest.limit = n_obs + 1L,
                                  pbkrtest.limit = n_obs + 1L)
          as.data.frame(emm$contrasts)
        }, error = function(e) {
          # Fall back to emmeans' own default rather than losing the contrasts
          tryCatch(as.data.frame(
            emmeans::emmeans(lmer_result$fit, pairwise ~ donor_status)$contrasts),
            error = function(e2) NULL)
        })
        if (!is.null(lmer_pairwise) && "df" %in% colnames(lmer_pairwise)) {
          lmer_pairwise_df_method <- if (all(!is.finite(lmer_pairwise$df))) {
            "asymptotic (z-test, df = Inf)"
          } else "Satterthwaite"
        }
      }

      # ---- SENSITIVITY: Donor-level test (conservative — collapses size variation) ----
      donor_test <- if (is_parametric) "ANOVA" else "Kruskal-Wallis"
      donor_p <- NA_real_
      eta_sq <- NA_real_
      # The two branches do NOT compute the same statistic. ANOVA gives eta^2
      # (SS_between / SS_total); Kruskal-Wallis gives epsilon^2 (H / (n-1)),
      # a different quantity on a different scale. Both used to be reported
      # under the single name "eta_sq" on screen and in the CSV, which misnames
      # the non-parametric result. The value stays in `eta_sq` (so nothing
      # reading the result list positionally breaks), but every place that
      # DISPLAYS or EXPORTS it now labels it with `eta_sq_label`.
      eta_sq_label <- if (is_parametric) "eta_sq" else "epsilon_sq"

      if (is_parametric) {
        fit <- tryCatch(lm(value ~ donor_status, data = ddf), error = function(e) NULL)
        if (!is.null(fit)) {
          at <- tryCatch(anova(fit), error = function(e) NULL)
          if (!is.null(at)) donor_p <- suppressWarnings(as.numeric(at[["Pr(>F)"]][1]))
          eta_sq <- eta_squared(fit)
        }
      } else {
        kt <- tryCatch(kruskal.test(value ~ donor_status, data = ddf), error = function(e) NULL)
        if (!is.null(kt)) donor_p <- kt$p.value
        if (!is.null(kt) && n_donors > 1) eta_sq <- unname(kt$statistic / (n_donors - 1))
      }

      # ---- Pairwise comparisons with effect sizes (donor-level, conservative) ----
      groups <- levels(ddf$donor_status)
      pairs <- combn(groups, 2, simplify = FALSE)

      pairwise_rows <- lapply(pairs, function(pair) {
        g1_vals <- ddf$value[ddf$donor_status == pair[1]]
        g2_vals <- ddf$value[ddf$donor_status == pair[2]]
        contrast_label <- paste(pair[1], "vs", pair[2])

        if (is_parametric) {
          tt <- tryCatch(t.test(g1_vals, g2_vals), error = function(e) NULL)
          p_val <- if (!is.null(tt)) tt$p.value else NA_real_
        } else {
          wt <- tryCatch(wilcox.test(g1_vals, g2_vals, exact = FALSE), error = function(e) NULL)
          p_val <- if (!is.null(wt)) wt$p.value else NA_real_
        }

        cd <- cohens_d(g1_vals, g2_vals)

        data.frame(
          contrast = contrast_label,
          p_value = p_val,
          cohens_d = cd$d,
          d_ci_lo = cd$ci_lo,
          d_ci_hi = cd$ci_hi,
          # Unstandardised difference + Welch CI, in the feature's own units.
          # This is what the forest plot now shows: it does not depend on which
          # SD you standardise by, which matters because the group SDs here can
          # differ by 50x (floor effects) and the pooled-SD d then lands at an
          # arbitrary point between the two Glass deltas.
          raw_diff = cd$raw_diff,
          raw_ci_lo = cd$raw_ci_lo,
          raw_ci_hi = cd$raw_ci_hi,
          sd1 = cd$sd1, sd2 = cd$sd2,
          glass1 = cd$glass1, glass2 = cd$glass2,
          var_ratio = cd$var_ratio,
          var_ratio_warn = cd$var_ratio_warn,
          n1 = length(g1_vals[is.finite(g1_vals)]),
          n2 = length(g2_vals[is.finite(g2_vals)]),
          stringsAsFactors = FALSE
        )
      })

      pairwise_df <- do.call(rbind, pairwise_rows)
      if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
        pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "BH")
        pairwise_df$sig <- ifelse(!is.na(pairwise_df$p_adj),
                                  pairwise_df$p_adj < alpha_num,
                                  FALSE)
      }

      # ---- Per-group descriptive summary (donor level) ----
      # The spread within each group is what makes a standardised effect size
      # readable, and on this cohort it is wildly unequal (donor-mean SDs of
      # 4.06 / 10.75 / 0.575 for prop_Beta cell). Surfacing it next to the
      # effect sizes is the difference between "d = 6.15" meaning something and
      # meaning nothing.
      group_summary <- tryCatch({
        ddf %>%
          dplyr::group_by(donor_status) %>%
          dplyr::summarise(
            n_donors = dplyr::n(),
            n_islets = sum(n_islets, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            sd = stats::sd(value, na.rm = TRUE),
            median = stats::median(value, na.rm = TRUE),
            .groups = "drop"
          )
      }, error = function(e) NULL)

      # ---- Islet-size companion models ----
      # See size_companion_models() in utils_stats.R for why: the primary model
      # conditions on neither diameter nor a group-varying size slope, yet both
      # matter, and the diameter distribution differs by group.
      size_models <- size_companion_models(rdf)

      # ---- Age-group companion (Steenbuck et al. 2025 covariate form) ----
      # `value ~ donor_status + age_group + (1|donor)`. See
      # age_group_companion_model() in utils_stats.R for why HbA1c and
      # C-peptide are deliberately excluded (mediators, not confounders).
      age_models <- age_group_companion_model(rdf)

      # ---- Per-bin analysis (exploratory, Section 3) ----
      # Bin width is FIXED at STATS_BIN_WIDTH_UM — a tunable width is a
      # multiple-comparison p-hacking surface. Always computed; the Section 3 "Show
      # exploratory analyses" toggle is a client-side display gate only.
      bin_anova <- NULL
      bin_kendall <- NULL
      bin_pairwise <- NULL

      {
        bw <- STATS_BIN_WIDTH_UM

        rdf_binned <- bin_islet_sizes(rdf, "islet_diam_um", bw)

        bin_anova <- per_bin_donor_anova(rdf_binned, "diam_bin", "donor_status", "value")
        bin_kendall <- per_bin_donor_kendall(rdf_binned, "diam_bin", "donor_status", "value")
        bin_pairwise <- per_bin_donor_pairwise(rdf_binned, "diam_bin", "donor_status", "value",
                                                parametric = is_parametric)

        # ---- Per-bin coverage, so the heatmap can state what each cell rests on ----
        # Coverage is strongly group-dependent at the large-diameter end: bin 225
        # is carried by 4 ND islets from 2 donors versus 14 Aab+ islets from 5,
        # and it rendered with exactly the same visual weight as bin 25's 4,282
        # islets from 15 donors. Bins failing the >=2 groups x >=2 donors rule
        # were dropped silently, so the size range actually compared was capped
        # by whichever group ran out of large islets first — invisibly.
        bin_coverage <- tryCatch({
          rdf_binned %>%
            dplyr::filter(is.finite(diam_mid)) %>%
            dplyr::group_by(diam_mid, donor_status) %>%
            dplyr::summarise(n_donors = dplyr::n_distinct(`Case ID`),
                             n_islets = dplyr::n(), .groups = "drop")
        }, error = function(e) NULL)

        # Bins with no omnibus row, and WHY.
        #
        # The reason matches bin_coverage_ok() in utils_stats.R, which all three
        # per-bin tests now gate on: >= 2 groups present AND *every* group present
        # in the bin having >= 2 donors. A bin with ND=1, Aab+=3, T1D=2 donors is
        # dropped on ND's single donor, and — since pairwise is gated on the same
        # rule — contributes no contrasts either. No pairwise post-hoc runs in a
        # bin with no omnibus behind it.
        #
        # The `status` field stays as a DRIFT ASSERTION: it should now always read
        # "not tested at all". If it ever reports surviving pairwise contrasts,
        # the three rules have diverged again.
        bin_dropped <- tryCatch({
          if (is.null(bin_coverage)) NULL else {
            tested <- if (!is.null(bin_anova)) bin_anova$mid else numeric(0)
            pw_mids <- if (!is.null(bin_pairwise) && nrow(bin_pairwise) > 0)
              unique(bin_pairwise$mid) else numeric(0)
            miss <- setdiff(sort(unique(bin_coverage$diam_mid)), tested)
            if (!length(miss)) NULL else {
              do.call(rbind, lapply(miss, function(m) {
                cv <- bin_coverage[bin_coverage$diam_mid == m, , drop = FALSE]
                thin <- cv$donor_status[cv$n_donors < 2]
                reason <- if (nrow(cv) < 2) {
                  sprintf("only %d donor-status group%s present (need 2)",
                          nrow(cv), if (nrow(cv) == 1) "" else "s")
                } else if (length(thin)) {
                  sprintf("%s below 2 donors; the omnibus rule requires every group present in the bin to have >= 2",
                          paste(as.character(thin), collapse = " and "))
                } else "coverage rule not met"
                pairs_here <- if (m %in% pw_mids)
                  bin_pairwise$pair[bin_pairwise$mid == m] else character(0)
                data.frame(
                  mid = m,
                  n_islets = sum(cv$n_islets),
                  detail = paste(sprintf("%s %dd/%di", cv$donor_status,
                                         cv$n_donors, cv$n_islets), collapse = ", "),
                  reason = reason,
                  status = if (length(pairs_here))
                    paste0("BUG: pairwise contrasts survived without an omnibus (",
                           paste(pairs_here, collapse = ", "),
                           ") — the per-bin coverage rules have diverged")
                  else "not tested at all",
                  stringsAsFactors = FALSE)
              }))
            }
          }
        }, error = function(e) NULL)

        # BH correction across bins to control false discovery rate
        if (!is.null(bin_anova) && nrow(bin_anova) > 0) {
          bin_anova$p_adj <- p.adjust(bin_anova$p_anova, method = "BH")
        }
        if (!is.null(bin_kendall) && nrow(bin_kendall) > 0) {
          bin_kendall$p_adj <- p.adjust(bin_kendall$p_kendall, method = "BH")
        }
        if (!is.null(bin_pairwise) && nrow(bin_pairwise) > 0) {
          bin_pairwise$p_adj <- p.adjust(bin_pairwise$p_value, method = "BH")
        }
      }

      # ---- Demographics (NULL if columns absent) ----
      demo_summary <- NULL
      age_corr <- NULL
      age_model_p <- NULL
      gender_strat <- NULL

      has_age <- "age" %in% colnames(ddf)
      has_gender <- "gender" %in% colnames(ddf)

      if (has_age || has_gender) {
        # Demographic summary per donor group
        demo_rows <- lapply(levels(ddf$donor_status), function(g) {
          sub <- ddf[ddf$donor_status == g, , drop = FALSE]
          n_d <- nrow(sub)
          islets_total <- sum(sub$n_islets, na.rm = TRUE)

          age_med <- NA_real_; age_range_str <- ""
          if (has_age) {
            ages <- as.numeric(sub$age)
            ages <- ages[is.finite(ages)]
            if (length(ages) > 0) {
              age_med <- median(ages)
              age_range_str <- paste0(min(ages), "-", max(ages))
            }
          }

          pct_m <- NA_real_; pct_f <- NA_real_
          if (has_gender) {
            genders <- as.character(sub$gender)
            genders <- genders[!is.na(genders) & nzchar(genders)]
            if (length(genders) > 0) {
              pct_m <- round(100 * sum(genders == "M") / length(genders), 1)
              pct_f <- round(100 * sum(genders == "F") / length(genders), 1)
            }
          }

          data.frame(
            donor_status = g, n = n_d, n_islets = islets_total,
            age_median = age_med, age_range = age_range_str,
            pct_male = pct_m, pct_female = pct_f,
            stringsAsFactors = FALSE
          )
        })
        demo_summary <- do.call(rbind, demo_rows)

        # Age-feature correlation across ALL islets (note: correlated within donors)
        if (has_age) {
          a <- as.numeric(rdf$age)
          v <- rdf$value
          keep <- is.finite(a) & is.finite(v)
          age_corr <- if (sum(keep) >= 3) {
            ct <- tryCatch(cor.test(a[keep], v[keep], method = "pearson"), error = function(e) NULL)
            if (!is.null(ct)) {
              data.frame(cor = unname(ct$estimate), p_value = ct$p.value,
                         n = sum(keep), stringsAsFactors = FALSE)
            }
          }
        }

        # Covariate-adjusted model on DONOR-level data (correct level for this)
        if (has_age && has_gender) {
          ddf_cov <- ddf
          ddf_cov$age_num <- as.numeric(ddf_cov$age)
          ddf_cov$gender_f <- factor(ddf_cov$gender)
          cov_fit <- tryCatch(lm(value ~ donor_status + age_num + gender_f, data = ddf_cov), error = function(e) NULL)
          if (!is.null(cov_fit)) {
            at <- tryCatch(anova(cov_fit), error = function(e) NULL)
            if (!is.null(at) && "donor_status" %in% rownames(at)) {
              age_model_p <- as.numeric(at["donor_status", "Pr(>F)"])
            }
          }
        }

        # Sex-stratified tests on donor-level data.
        # Pairwise tests are run UNADJUSTED within each sex, then BH-corrected
        # across all 6 contrasts (3 pairs x 2 sexes) so the family-wise FDR is
        # controlled the same way as the primary pairwise table.
        if (has_gender) {
          genders_present <- unique(na.omit(as.character(ddf$gender)))
          gender_rows <- lapply(genders_present, function(gen) {
            sub <- ddf[as.character(ddf$gender) == gen, , drop = FALSE]
            if (n_distinct(sub$donor_status) < 2) return(NULL)
            n_per_grp <- table(sub$donor_status)
            low_power <- any(n_per_grp[n_per_grp > 0] < 3)
            pw <- if (is_parametric) {
              # Welch t-tests over the SAME contrast list Section 2 uses, not
              # `pairwise.t.test`. Two reasons, both measured:
              #  1. `pairwise.t.test` pools the SD across ALL groups. Donor-level
              #     SDs here differ by up to 50:1 (ND 4.06 / Aab+ 10.75 / T1D
              #     0.575 for prop_Beta cell), so pooling is exactly the wrong
              #     assumption -- and it disagreed with Section 2's Welch result
              #     (Males ND-vs-T1D 0.104 pooled vs 0.088 Welch).
              #  2. Its p-value matrix is indexed [later, earlier], so unpacking
              #     it with as.table() produced reversed labels ("Aab+ vs ND"
              #     where the rest of the tab says "ND vs Aab+").
              # Iterating `pairs` fixes the statistic and the label together.
              tryCatch({
                rows <- lapply(pairs, function(pair) {
                  v1 <- sub$value[sub$donor_status == pair[1]]
                  v2 <- sub$value[sub$donor_status == pair[2]]
                  if (length(v1) < 2 || length(v2) < 2) return(NULL)
                  tt <- tryCatch(t.test(v1, v2), error = function(e) NULL)
                  if (is.null(tt)) return(NULL)
                  data.frame(group1 = pair[1], group2 = pair[2],
                             p_value = unname(tt$p.value),
                             stringsAsFactors = FALSE)
                })
                mat <- do.call(rbind, Filter(Negate(is.null), rows))
                if (is.null(mat)) NULL else mat[!is.na(mat$p_value), , drop = FALSE]
              }, error = function(e) NULL)
            } else {
              # Same [later, earlier] matrix orientation as pairwise.t.test, so
              # iterate `pairs` here too and keep the labels in level order.
              tryCatch({
                rows <- lapply(pairs, function(pair) {
                  v1 <- sub$value[sub$donor_status == pair[1]]
                  v2 <- sub$value[sub$donor_status == pair[2]]
                  if (length(v1) < 1 || length(v2) < 1) return(NULL)
                  wt <- tryCatch(wilcox.test(v1, v2, exact = FALSE),
                                 error = function(e) NULL)
                  if (is.null(wt)) return(NULL)
                  data.frame(group1 = pair[1], group2 = pair[2],
                             p_value = unname(wt$p.value),
                             stringsAsFactors = FALSE)
                })
                mat <- do.call(rbind, Filter(Negate(is.null), rows))
                if (is.null(mat)) NULL else mat[!is.na(mat$p_value), , drop = FALSE]
              }, error = function(e) NULL)
            }
            if (is.null(pw) || nrow(pw) == 0) return(NULL)
            pw$gender <- gen
            pw$contrast <- paste(pw$group1, "vs", pw$group2)
            pw$low_power <- low_power
            pw
          })
          gender_strat <- do.call(rbind, gender_rows)
          if (!is.null(gender_strat) && nrow(gender_strat) > 0) {
            gender_strat$p_adj <- p.adjust(gender_strat$p_value, method = "BH")
          }
        }
      }

      # ---- Autoantibody analysis (Aab+ islets only) ----
      aab_cols <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
      aab_avail <- intersect(aab_cols, colnames(rdf))
      aab_summary <- NULL
      aab_donor_summary <- NULL
      aab_kw <- NULL

      if (length(aab_avail) > 0) {
        aab_sub <- rdf[rdf$donor_status == "Aab+", , drop = FALSE]
        if (nrow(aab_sub) > 0) {
          # Count AAbs per islet (same for all islets from same donor).
          # `sapply` returns a plain vector when only one AAb column survives the
          # intersect above, which rowSums() cannot take — force a matrix.
          aab_sub$n_aab <- rowSums(
            matrix(unlist(lapply(aab_avail, function(col) as.numeric(as.logical(aab_sub[[col]])))),
                   nrow = nrow(aab_sub), ncol = length(aab_avail)),
            na.rm = TRUE)
          aab_summary <- aab_sub

          # ---- AAb-count test, at the DONOR level ----
          # The AAb count is a donor attribute: every islet from a donor carries
          # the same value. Testing it across islets is pure pseudoreplication —
          # on this cohort that produced KW p = 1.6e-15 from 1,666 islets that
          # come from just 5 donors split 1 / 3 / 1 across the count groups,
          # while the donor-level test gives p = 0.67. Aggregate first, and carry
          # the per-group donor counts so the plot can state how thin this is.
          aab_kw <- tryCatch({
            grp_of <- function(n) ifelse(n >= 3, "3+", as.character(n))
            ddf_aab <- aab_sub %>%
              dplyr::mutate(.aab_group = grp_of(n_aab)) %>%
              dplyr::group_by(`Case ID`, .aab_group) %>%
              dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
            n_donor_per_group <- table(ddf_aab$.aab_group)
            p_val <- if (length(n_donor_per_group) >= 2) {
              kt <- kruskal.test(value ~ factor(.aab_group), data = ddf_aab)
              unname(kt$p.value)
            } else NA_real_
            list(p_value = p_val,
                 n_donors = nrow(ddf_aab),
                 n_islets = nrow(aab_sub),
                 donors_per_group = n_donor_per_group,
                 # Any group carried by a single donor makes the comparison a
                 # between-donor difference of one, not a group effect.
                 underpowered = any(n_donor_per_group < 2))
          }, error = function(e) NULL)

          # Per-donor AAb profile table
          aab_donor_rows <- lapply(unique(aab_sub$`Case ID`), function(cid) {
            dsub <- aab_sub[aab_sub$`Case ID` == cid, , drop = FALSE]
            row <- data.frame(
              `Case ID` = cid,
              n_islets = nrow(dsub),
              mean_value = mean(dsub$value, na.rm = TRUE),
              n_aab = dsub$n_aab[1],
              check.names = FALSE, stringsAsFactors = FALSE
            )
            # Add individual AAb flags
            for (ac in aab_avail) {
              short_name <- sub("^AAb_", "", ac)
              row[[short_name]] <- if (as.logical(dsub[[ac]][1])) "+" else "-"
            }
            row
          })
          aab_donor_summary <- do.call(rbind, aab_donor_rows)
        }
      }

      # ---- AUC: trapezoidal integration of the binned mean, per donor group ----
      #
      # Computed from THIS tab's filtered frame (`rdf_binned`, fixed
      # STATS_BIN_WIDTH_UM). It used to integrate the Plot module's summary_df(),
      # which made it the one card on the page that answered a different
      # question than everything above it:
      #   * it ignored every Advanced filter — with min-cells + diameter + zero +
      #     outlier filters cutting 5,212 islets to 131, the AUC was byte-identical
      #     to the unfiltered run, and the "T1D shows X% decrease" sentence
      #     described a dataset the rest of the page had discarded;
      #   * it obeyed the Plot sidebar's `binwidth` slider (hidden on this tab),
      #     swinging ND from 9,639 to 10,483 to 9,124 across 10/50/75 µm, while
      #     Section 3 advertises a fixed 50 µm bin width;
      #   * it obeyed the Plot sidebar's Mean/Median radio, so with "Median + IQR"
      #     selected it integrated medians while the card said "mean" (T1D 465 -> 186).
      auc_by_group <- tryCatch({
        rdf_binned %>%
          dplyr::filter(is.finite(diam_mid), is.finite(value)) %>%
          dplyr::group_by(donor_status, diam_mid) %>%
          dplyr::summarise(y = mean(value, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(donor_status, diam_mid) %>%
          dplyr::group_by(donor_status) %>%
          dplyr::summarise(
            auc = if (dplyr::n() < 2) 0
                  else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
            n_bins = dplyr::n(),
            .groups = "drop"
          )
      }, error = function(e) NULL)

      list(
        n_islets         = n_islets_filtered,
        n_donors         = n_donors,
        donor_df         = ddf,
        feature_name     = feature_name,
        feature_col      = feature_col,
        # Primary test: mixed-effects
        lmer_p           = lmer_p,
        icc              = icc,
        lmer_pairwise    = lmer_pairwise,
        lmer_error       = lmer_error,
        # Sensitivity test: donor-level
        donor_test       = donor_test,
        donor_p          = donor_p,
        eta_sq           = eta_sq,
        eta_sq_label     = eta_sq_label,
        alpha            = alpha_num,
        is_parametric    = is_parametric,
        pairwise_df      = pairwise_df,
        bin_anova        = bin_anova,
        bin_kendall      = bin_kendall,
        bin_pairwise     = bin_pairwise,
        bin_coverage     = bin_coverage,
        bin_dropped      = bin_dropped,
        group_summary    = group_summary,
        size_models      = size_models,
        age_models       = age_models,
        auc_by_group     = auc_by_group,
        demo_summary     = demo_summary,
        age_corr         = age_corr,
        age_model_p      = age_model_p,
        gender_strat     = gender_strat,
        norm_result      = norm_result,
        aab_summary      = aab_summary,
        aab_donor_summary = aab_donor_summary,
        aab_kw           = aab_kw,
        rdf              = rdf,
        outlier_summary  = outlier_summary_pre,
        exclude_zero     = exclude_zero,
        zero_dropped     = zero_dropped,
        # Data-subset knobs, recorded so the CSV header can report the settings
        # that actually produced this run. The UI promises "every setting is
        # recorded in the downloaded CSV header" (mod_statistics_ui.R) — these
        # two were the ones it silently omitted.
        min_cells        = min_cells,
        min_cells_applied = min_cells_applied,
        diam_range       = diam_range_used,
        # df method actually used for the mixed-effects pairwise contrasts
        lmer_pairwise_df_method = lmer_pairwise_df_method,
        # The Plot sidebar's own filters run BEFORE this tab sees the data
        # (`raw_df_prezero`), and its controls are hidden while the Statistics tab
        # is active — so donor groups, min-cells, diameter, age and sex set over
        # there silently change every number here. Capture the description so the
        # banner and the CSV can name the inherited subset instead of leaving the
        # user to infer it.
        plot_subset_desc = tryCatch(get_selection_description(), error = function(e) character(0))
      )
    }, ignoreInit = TRUE)

    # ===========================================================================
    # Card 1: Overview Banner
    # ===========================================================================

    output$overview_banner <- renderUI({
      st <- stats_run()
      if (is.null(st)) {
        # Distinguish "not run yet" from "ran, but the filters left nothing".
        # Both used to render the same grey "click Run Statistics" line, so an
        # over-restrictive filter looked exactly like an idle tab.
        reason <- stats_empty_reason()
        if (!is.null(reason)) {
          return(tags$div(
            style = paste0("padding: 10px 14px; background: #fdecea; ",
                           "border: 1px solid #f5c6cb; border-radius: 5px; color: #721c24;"),
            tags$strong("Statistics ran, but produced no result."), tags$br(),
            tags$span(paste0("Reason: ", reason))))
        }
        return(tags$p(style = "margin: 0",
                      "Select a feature in the sidebar, then click 'Run Statistics'."))
      }

      # Descriptive only — the primary p-value now lives in Section 2's headline.
      zero_badge <- if (isTRUE(st$exclude_zero)) {
        tags$span(style = paste0("background: #fff3e0; color: #6a3a00; padding: 4px 10px; ",
                                 "border: 1px solid #ffcc80; border-radius: 4px;"),
                  sprintf("Zero-valued excluded (%s islet%s dropped)",
                          formatC(st$zero_dropped, format = "d", big.mark = ","),
                          if (st$zero_dropped == 1) "" else "s"))
      } else NULL

      # Warn when a min-cells value was typed but could not be honoured (Excel
      # fallback path, where the per-islet cell counts are unavailable). Silently
      # ignoring it is what made the Microenvironment min-cells bug invisible.
      min_cells_warn <- if (!isTRUE(st$min_cells_applied)) {
        tags$span(style = paste0("background: #fdecea; color: #8a1c12; padding: 4px 10px; ",
                                 "border: 1px solid #f5b7b1; border-radius: 4px;"),
                  sprintf("Min cells/islet = %g NOT applied — per-islet cell counts unavailable for this data source",
                          st$min_cells))
      } else NULL

      # Filters set on the Plot sidebar run upstream of this tab and its controls
      # are hidden here, so name them rather than letting them act invisibly.
      inherited <- st$plot_subset_desc
      inherited <- inherited[!grepl("^Bin Width:", inherited)]
      inherited_block <- if (length(inherited)) {
        tags$details(style = "margin-top: 8px; width: 100%; font-size: 12px",
          tags$summary(style = "cursor: pointer; font-weight: 600",
            "Data subset inherited from the Plot sidebar ",
            tags$span(style = "font-weight: 400",
                      "— applied before this tab's filters")),
          tags$ul(style = "margin: 6px 0 0 0; padding-left: 20px",
                  lapply(inherited, tags$li)))
      } else NULL

      tags$div(
        tags$div(style = "display: flex; gap: 15px; align-items: center; flex-wrap: wrap",
          tags$span(style = "padding: 4px 10px; font-weight: 600",
                    st$feature_name),
          tags$span(style = "padding: 4px 10px",
                    paste0(formatC(st$n_islets, format = "d", big.mark = ","),
                           " islets, ", st$n_donors, " donors")),
          tags$span(style = "font-size: 12px",
                    "Primary result is shown in Section 2 below."),
          zero_badge,
          min_cells_warn
        ),
        inherited_block
      )
    })

    # ===========================================================================
    # Section 2 headline: the single pre-specified primary result
    # ===========================================================================

    output$primary_headline <- renderUI({
      st <- stats_run()
      if (is.null(st)) {
        return(tags$p(style = "margin: 0",
                      "Click 'Run Statistics' to compute the primary test."))
      }

      primary_p <- st$lmer_p
      sig <- !is.na(primary_p) && primary_p < st$alpha
      p_col <- if (is.na(primary_p)) "#999" else if (sig) "#d62728" else "#2e7d32"
      icc_str <- if (is.na(st$icc)) "N/A" else sprintf("%.2f", st$icc)

      # ---- Per-group mean / SD / n, with an unequal-spread flag ----
      gs <- st$group_summary
      group_block <- if (!is.null(gs) && nrow(gs) > 0) {
        vr <- {
          v <- gs$sd[is.finite(gs$sd) & gs$sd > 0]^2
          if (length(v) >= 2) max(v) / min(v) else NA_real_
        }
        spread_note <- if (is.finite(vr) && vr > 4) {
          tags$div(style = paste0("margin-top: 6px; padding: 6px 8px; background: #fff8e1;",
                                  " border-left: 3px solid #f0ad4e; border-radius: 3px;",
                                  " font-size: 11px; color: #6a5200; line-height: 1.4;"),
            sprintf("Group spreads differ %.0f-fold. ", vr),
            "A pooled-SD effect size (Cohen's d) divides by a common scale that does not ",
            "exist here — read the raw difference in the forest plot, and treat d as indicative only.")
        } else NULL
        tagList(
          tags$div(style = "margin-top: 12px; font-size: 12px",
            tags$strong("Group spread (donor means):"),
            tags$table(style = "width: 100%; margin-top: 4px; font-size: 12px; border-collapse: collapse",
              tags$tr(style = "text-align: right",
                tags$th(style = "text-align: left; font-weight: 600", ""),
                tags$th(style = "font-weight: 600", "donors"),
                tags$th(style = "font-weight: 600", "islets"),
                tags$th(style = "font-weight: 600", "mean"),
                tags$th(style = "font-weight: 600", "SD")),
              lapply(seq_len(nrow(gs)), function(i) {
                tags$tr(style = "text-align: right",
                  tags$td(style = "text-align: left; font-weight: 600",
                          as.character(gs$donor_status[i])),
                  tags$td(gs$n_donors[i]),
                  tags$td(formatC(gs$n_islets[i], format = "d", big.mark = ",")),
                  tags$td(sprintf("%.3g", gs$mean[i])),
                  tags$td(if (is.finite(gs$sd[i])) sprintf("%.3g", gs$sd[i]) else "—"))
              }))),
          spread_note)
      } else NULL

      # ---- Islet-size companion models ----
      sm <- st$size_models
      size_block <- if (!is.null(sm) && any(is.finite(c(sm$p_status_adj, sm$p_interaction)))) {
        interaction_hit <- is.finite(sm$p_interaction) && sm$p_interaction < st$alpha
        warn <- if (interaction_hit) {
          tags$div(style = paste0("margin-top: 6px; padding: 6px 8px; background: #fdecea;",
                                  " border-left: 3px solid #d9534f; border-radius: 3px;",
                                  " font-size: 11px; color: #8a1c12; line-height: 1.4;"),
            tags$strong("The group difference depends on islet size. "),
            "The feature-vs-diameter slope differs by donor status, so the single p above is a ",
            "marginal average over a size distribution that itself differs between groups. ",
            "Read the size-stratified analysis in Section 3 alongside it, not as an afterthought.")
        } else NULL
        tagList(
          tags$div(style = "margin-top: 12px; font-size: 12px; line-height: 1.6",
            tags$strong("Islet-size companions:"), tags$br(),
            tags$span("diameter effect: "),
            format_p_value(sm$p_diam), tags$br(),
            tags$span("donor status adjusted for diameter: "),
            format_p_value(sm$p_status_adj),
            tags$span(
                      sprintf(" (vs %s unadjusted)", format_p_value(primary_p))), tags$br(),
            tags$span("status × diameter interaction: "),
            tags$strong(style = if (interaction_hit) "color: #d9534f;" else "",
                        format_p_value(sm$p_interaction))),
          warn)
      } else if (!is.null(sm) && !is.null(sm$error_msg)) {
        tags$div(style = "margin-top: 10px; font-size: 11px",
                 "Islet-size companion models unavailable: ", sm$error_msg)
      } else NULL

      # ---- Age-group companion ----
      # 5 of our 15 donors are under 5 years old, so "is this a disease effect or
      # an age effect?" is a live question rather than a formality.
      am <- st$age_models
      age_block <- if (!is.null(am) && is.finite(am$p_status_adj)) {
        age_drives <- is.finite(am$p_age_group) && am$p_age_group < st$alpha
        tagList(
          tags$div(style = "margin-top: 12px; font-size: 12px; line-height: 1.6",
            tags$strong("Age-group companion:"),
            tags$span(
                      sprintf(" (<%g y: %d donors | \u2265%g y: %d donors)",
                              am$cut, am$n_young, am$cut, am$n_old)), tags$br(),
            tags$span("donor status adjusted for age group: "),
            format_p_value(am$p_status_adj),
            tags$span(
                      sprintf(" (vs %s unadjusted)", format_p_value(primary_p))), tags$br(),
            tags$span("age group alone: "),
            tags$strong(style = if (age_drives) "color: #d9534f;" else "",
                        format_p_value(am$p_age_group))),
          if (age_drives) {
            tags$div(style = paste0("margin-top: 6px; padding: 6px 8px; background: #fdecea;",
                                    " border-left: 3px solid #d9534f; border-radius: 3px;",
                                    " font-size: 11px; color: #8a1c12; line-height: 1.4;"),
              tags$strong("Age group independently predicts this feature. "),
              "With 5 of 15 donors under 5 years old, read the adjusted p above as the ",
              "disease-stage result and treat the unadjusted headline with caution.")
          } else NULL)
      } else if (!is.null(am) && !is.null(am$error_msg)) {
        tags$div(style = "margin-top: 10px; font-size: 11px",
                 "Age-group companion unavailable: ", am$error_msg)
      } else NULL

      tags$div(
        tags$div(style = "font-size: 13px; margin-bottom: 8px; line-height: 1.4",
          "Linear mixed-effects model — donor status with a per-donor random intercept ",
          tags$em("(the pre-specified primary test).")),
        tags$div(style = paste0("display: inline-block; background: ", p_col,
                                "; color: white; padding: 8px 16px; border-radius: 6px;",
                                " font-size: 20px; font-weight: 700;"),
          # format_p_value(): fixed notation down to 1e-4, scientific below.
          # This used to be an unconditional formatC(format = "e"), so p = 0.3476
          # rendered as "3.48e-01" — trivially misread as 0.0348 — while the
          # hypothesis table below printed the same value as "0.3476".
          paste0("p = ", format_p_value(primary_p))),
        tags$span(style = "margin-left: 10px; font-size: 14px",
          paste0(if (is.na(primary_p)) "—" else if (sig) "Significant" else "Not significant",
                 " at α = ", st$alpha, " (fixed)")),
        tags$div(style = "margin-top: 12px; font-size: 13px",
          tags$strong("ICC = "), icc_str,
          tags$span(
                    "  (fraction of variance from donor identity)")),
        tags$div(style = "margin-top: 6px; font-size: 12px; line-height: 1.4",
          "Sensitivity test auto-selected from normality: ", tags$strong(st$donor_test),
          " on donor means. Full secondary tests are in the ",
          tags$em("Sensitivity check"), " panel below."),

        # ---- Per-group spread ----
        # Effect sizes below are unreadable without it: a standardised d divides
        # by a pooled SD, which only means something when the group SDs are
        # comparable. Show them so the reader can judge that themselves.
        group_block,

        # ---- Islet-size companion models ----
        size_block,

        # ---- Age-group companion model ----
        age_block
      )
    })

    # ===========================================================================
    # Normality test results
    # ===========================================================================

    output$normality_result <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$norm_result)) return(NULL)

      nr <- st$norm_result
      items <- lapply(seq_len(nrow(nr)), function(i) {
        g <- nr$donor_status[i]
        w_val <- nr$W[i]
        p_val <- nr$p_value[i]
        n_val <- nr$n[i]
        is_norm <- nr$is_normal[i]
        if (is.na(is_norm)) {
          col <- "#999"
          label <- paste0(g, ": n=", n_val, " (too few)")
        } else if (is_norm) {
          col <- "#2e7d32"
          label <- paste0(g, ": W=", sprintf("%.3f", w_val),
                          ", p=", sprintf("%.3f", p_val), " (normal)")
        } else {
          col <- "#e65100"
          label <- paste0(g, ": W=", sprintf("%.3f", w_val),
                          ", p=", sprintf("%.3f", p_val), " (non-normal)")
        }
        tags$span(style = paste0("color: ", col, "; margin-right: 12px;"), label)
      })

      # Report the SAME decision the engine made (st$is_parametric is the single
      # source of truth) so the readout never contradicts the test actually run.
      auto_parametric <- isTRUE(st$is_parametric)
      all_na <- all(is.na(nr$is_normal))
      reason <- if (auto_parametric) "all groups passed Shapiro-Wilk"
                else if (all_na) "normality indeterminate (too few donors)"
                else "at least one group is non-normal"
      suggestion <- tags$span(
        style = paste0("font-weight: 600; color: ", if (auto_parametric) "#2e7d32" else "#e65100", ";"),
        paste0("→ Sensitivity test auto-selected: ", st$donor_test,
               " on donor means (", reason, ")."))

      tags$div(
        style = "padding: 8px 12px; font-size: 12px; line-height: 1.8",
        tags$strong("Shapiro-Wilk normality (donor means):"),
        tags$br(),
        tagList(items),
        tagList(tags$br(), suggestion),
        tags$br(),
        tags$span(style = "font-size: 11px",
                  "Note: with n=5/group, Shapiro-Wilk has limited power to detect departures from normality. ",
                  "The test type is chosen automatically and is not a user setting.")
      )
    })

    # ===========================================================================
    # Pseudoreplication info banner
    # ===========================================================================

    output$pseudorep_banner <- renderUI({
      st <- stats_run()
      if (is.null(st)) return(NULL)

      icc_str <- if (is.na(st$icc)) "N/A" else sprintf("%.2f", st$icc)
      donor_str <- if (is.na(st$donor_p)) "N/A" else
        paste0("p = ", format_p_value(st$donor_p))

      icc_interp <- if (!is.na(st$icc) && st$icc > 0.3) {
        " — substantial clustering: donor identity explains much of the variance."
      } else if (!is.na(st$icc) && st$icc > 0.1) {
        " — moderate clustering."
      } else if (!is.na(st$icc)) {
        " — low clustering: most variation is between islets within donors."
      } else {
        "."
      }

      # Show warning if mixed-effects model failed
      lmer_warn <- NULL
      if (!is.null(st$lmer_error)) {
        lmer_warn <- tags$div(
          style = "padding: 8px 12px; margin-bottom: 8px; font-size: 12px",
          tags$strong("Mixed-effects model issue: "), st$lmer_error,
          " — donor-level results below remain valid."
        )
      }

      tagList(
        lmer_warn,
        tags$div(
          style = "padding: 10px 14px; margin-bottom: 12px; font-size: 13px; line-height: 1.5",
          tags$strong("Statistical approach: "),
          paste0("The primary test is a mixed-effects model with donor as a random intercept. ",
                 "This uses all islets while accounting for the fact that islets from the same donor are correlated. ",
                 "Islets of different sizes are biologically distinct and contribute meaningful signal."),
          tags$br(),
          paste0("ICC = ", icc_str, icc_interp, " "),
          tags$br(),
          tags$strong("Donor-level sensitivity check: "),
          paste0("As a conservative cross-check, islets are averaged within each donor to yield one value per donor (N=",
                 st$n_donors, "), then tested with ", st$donor_test, ": ", donor_str, ". ",
                 "If both approaches agree, the result is robust.")
        )
      )
    })

    # ===========================================================================
    # Card 2: Hypothesis Testing — table + forest plot
    # ===========================================================================

    output$test_results_table <- renderTable({
      st <- stats_run()
      if (is.null(st)) return(NULL)

      # Shared with the Section 2 headline so the same p never appears in two
      # different notations on one page.
      fmt_p <- function(x) format_p_value(x)
      fmt_d <- function(x) ifelse(is.na(x), "N/A", sprintf("%.2f", x))

      # Effect-size columns: the unstandardised difference leads (it does not
      # depend on which SD you standardise by), Cohen's d follows as the
      # standardised companion, and Glass's delta shows the range d is picking a
      # point inside of when the group SDs differ.
      blank_row <- function(contrast, p, padj, sig) data.frame(
        Contrast = contrast, `p-value` = p, `p (adj)` = padj,
        `Diff (raw)` = "-", `95% CI (raw)` = "-",
        `Cohen's d` = "-", `Glass's d range` = "-", Sig = sig,
        check.names = FALSE, stringsAsFactors = FALSE)

      rows <- data.frame(
        Contrast = character(), `p-value` = character(), `p (adj)` = character(),
        `Diff (raw)` = character(), `95% CI (raw)` = character(),
        `Cohen's d` = character(), `Glass's d range` = character(), Sig = character(),
        check.names = FALSE, stringsAsFactors = FALSE
      )

      # PRIMARY: Mixed-effects global (islet-level, donor as random intercept)
      rows <- rbind(rows, blank_row(
        "Overall (mixed-effects, donor random intercept)",
        fmt_p(st$lmer_p), "-",
        if (!is.na(st$lmer_p) && st$lmer_p < st$alpha) "*" else ""))

      # Size companions, so the primary block states what it does and does not
      # condition on right where the primary p is read.
      sm <- st$size_models
      if (!is.null(sm) && any(is.finite(c(sm$p_status_adj, sm$p_interaction)))) {
        rows <- rbind(rows, blank_row(
          "  \u21b3 donor status adjusted for islet diameter",
          fmt_p(sm$p_status_adj), "-",
          if (is.finite(sm$p_status_adj) && sm$p_status_adj < st$alpha) "*" else ""))
        rows <- rbind(rows, blank_row(
          "  \u21b3 donor status \u00d7 islet diameter (is the effect size-dependent?)",
          fmt_p(sm$p_interaction), "-",
          if (is.finite(sm$p_interaction) && sm$p_interaction < st$alpha) "*" else ""))
      }

      am <- st$age_models
      if (!is.null(am) && is.finite(am$p_status_adj)) {
        rows <- rbind(rows, blank_row(
          "  \u21b3 donor status adjusted for age group (<13 / \u226513 y)",
          fmt_p(am$p_status_adj), "-",
          if (is.finite(am$p_status_adj) && am$p_status_adj < st$alpha) "*" else ""))
        rows <- rbind(rows, blank_row(
          "  \u21b3 age group (does age itself predict the feature?)",
          fmt_p(am$p_age_group), "-",
          if (is.finite(am$p_age_group) && am$p_age_group < st$alpha) "*" else ""))
      }

      # Pairwise from mixed model (estimated marginal means)
      lmer_pw <- st$lmer_pairwise
      if (!is.null(lmer_pw) && nrow(lmer_pw) > 0) {
        # emmeans(pairwise ~ ...) Tukey-adjusts `p.value` by default, so this
        # column belongs under "p (adj)". It was previously printed under
        # "p-value" with "p (adj)" showing "-", which read as an uncorrected
        # result and understated the correction that had in fact been applied.
        lmer_pw_rows <- blank_row(
          paste0(lmer_pw$contrast, " (mixed-effects pairwise, Tukey)"),
          "-", fmt_p(lmer_pw$p.value),
          ifelse(!is.na(lmer_pw$p.value) & lmer_pw$p.value < st$alpha, "*", ""))
        rows <- rbind(rows, lmer_pw_rows)
      }

      # SENSITIVITY: Donor-level means (collapses islets to one value per donor)
      rows <- rbind(rows, blank_row(
        paste0("Donor-level means (", st$donor_test, ", N=", st$n_donors, " donors)"),
        fmt_p(st$donor_p), "-",
        if (!is.na(st$donor_p) && st$donor_p < st$alpha) "*" else ""))

      # Donor-level pairwise with effect sizes
      pw <- st$pairwise_df
      if (!is.null(pw) && nrow(pw) > 0) {
        test_label <- if (st$is_parametric) "t-test" else "Wilcoxon"
        # A dagger marks contrasts whose group variances differ more than 4-fold:
        # there the pooled-SD d is not on a meaningful common scale, and the raw
        # difference is the number to read.
        pw_rows <- data.frame(
          Contrast = paste0(pw$contrast, " (donor-level ", test_label, ")",
                            # NB: isTRUE() collapses a vector to a single FALSE,
                            # so the marker never fired. Compare element-wise.
                            ifelse(pw$var_ratio_warn %in% TRUE, " \u2020", "")),
          `p-value` = fmt_p(pw$p_value),
          `p (adj)` = fmt_p(pw$p_adj),
          `Diff (raw)` = ifelse(is.na(pw$raw_diff), "N/A", sprintf("%.3g", pw$raw_diff)),
          `95% CI (raw)` = ifelse(is.na(pw$raw_ci_lo), "N/A",
                                  paste0("[", sprintf("%.3g", pw$raw_ci_lo), ", ",
                                         sprintf("%.3g", pw$raw_ci_hi), "]")),
          `Cohen's d` = fmt_d(pw$cohens_d),
          `Glass's d range` = ifelse(is.na(pw$glass1) | is.na(pw$glass2), "N/A",
                                     paste0(sprintf("%.2f", pmin(pw$glass1, pw$glass2)),
                                            " to ", sprintf("%.2f", pmax(pw$glass1, pw$glass2)))),
          Sig = ifelse(!is.na(pw$p_adj) & pw$p_adj < st$alpha, "*", ""),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        rows <- rbind(rows, pw_rows)
      }

      rows
    # 8 columns: Contrast, p-value, p (adj), Diff (raw), 95% CI (raw),
    # Cohen's d, Glass's d range, Sig. `align` must have exactly one character
    # per column or renderTable errors.
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lccccccc")

    output$forest_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)
      pw <- st$pairwise_df
      if (is.null(pw) || nrow(pw) == 0 || all(is.na(pw$raw_diff))) return(NULL)

      # Plot the UNSTANDARDISED difference with a Welch CI, in the feature's own
      # units. It used to plot Cohen's d — but d divides by a pooled SD, and when
      # the group SDs differ 50-fold (routine here: floor effects put a whole
      # group at ~0 with almost no spread) the pooled SD is not a common scale.
      # d = 6.15 with CI [3.18, 9.12] looked precise while Glass's delta ranged
      # 4.39 to 31.03 for the same contrast. The raw difference has no such
      # ambiguity; d and the Glass range stay in the table beside it.
      pw <- pw[!is.na(pw$raw_diff), , drop = FALSE]
      pw$contrast <- factor(pw$contrast, levels = rev(pw$contrast))
      pw$sig_label <- ifelse(pw$sig, "Significant", "NS")
      pw$hover <- sprintf(
        "%s<br>difference: %.3g<br>95%% CI: [%.3g, %.3g]<br>Cohen's d: %s<br>SDs: %.3g vs %.3g (ratio %.1f)",
        as.character(pw$contrast), pw$raw_diff, pw$raw_ci_lo, pw$raw_ci_hi,
        ifelse(is.na(pw$cohens_d), "N/A", sprintf("%.2f", pw$cohens_d)),
        pw$sd1, pw$sd2, pw$var_ratio)

      zero_suffix <- if (isTRUE(st$exclude_zero) && st$zero_dropped > 0) {
        sprintf(" — %s zero-valued dropped",
                formatC(st$zero_dropped, format = "d", big.mark = ","))
      } else ""
      warn_suffix <- if (any(pw$var_ratio_warn %in% TRUE)) {
        " — † group SDs differ >4×"
      } else ""
      forest_title <- paste0("Pairwise Differences (donor means)",
                             outlier_title_suffix(st$outlier_summary),
                             zero_suffix, warn_suffix)
      g <- ggplot(pw, aes(x = raw_diff, y = contrast, color = sig_label, text = hover)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#999") +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = raw_ci_lo, xmax = raw_ci_hi), height = 0.2) +
        scale_color_manual(values = c("Significant" = "#d62728", "NS" = "#1f77b4"),
                           name = "") +
        labs(x = paste0("Difference in ", st$feature_name, " (Welch 95% CI)"),
             y = NULL, title = forest_title) +
        theme_islet() +
        theme(legend.position = "none")

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = "text") %>%
        layout(font = fonts$global,
               title = list(font = fonts$title),
               xaxis = list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick,
                            automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick,
                            automargin = TRUE),
               margin = list(b = 80),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25,
                             font = fonts$legend_text,
                             title = list(font = fonts$legend_title))) %>%
        plotly_islet_config("statistics_forest_plot")
    })

    # ===========================================================================
    # Card 3: Per-Bin Significance Heatmap
    # ===========================================================================

    output$bin_heatmap <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)

      ba <- st$bin_anova
      bk <- st$bin_kendall
      bp <- st$bin_pairwise

      if ((is.null(ba) || nrow(ba) == 0) && (is.null(bk) || nrow(bk) == 0)) return(NULL)

      # Collect all bin midpoints
      all_mids <- sort(unique(c(ba$mid, bk$mid,
                                if (!is.null(bp) && nrow(bp) > 0) bp$mid)))
      if (length(all_mids) == 0) return(NULL)

      n_mids <- length(all_mids)

      # Build p-value vectors for global tests
      anova_p <- rep(NA_real_, n_mids)
      kendall_p <- rep(NA_real_, n_mids)

      if (!is.null(ba) && nrow(ba) > 0) {
        idx <- match(ba$mid, all_mids)
        anova_p[idx[!is.na(idx)]] <- ba$p_adj[!is.na(idx)]
      }
      if (!is.null(bk) && nrow(bk) > 0) {
        idx <- match(bk$mid, all_mids)
        kendall_p[idx[!is.na(idx)]] <- bk$p_adj[!is.na(idx)]
      }

      # Build p-value vectors for pairwise tests
      pair_names <- c("ND vs Aab+", "ND vs T1D", "Aab+ vs T1D")
      pair_p_list <- lapply(pair_names, function(pn) {
        pv <- rep(NA_real_, n_mids)
        if (!is.null(bp) && nrow(bp) > 0) {
          sub <- bp[bp$pair == pn, , drop = FALSE]
          if (nrow(sub) > 0) {
            idx <- match(sub$mid, all_mids)
            pv[idx[!is.na(idx)]] <- sub$p_adj[!is.na(idx)]
          }
        }
        pv
      })

      # Convert to -log10(q) and build matrices
      to_z <- function(p) ifelse(!is.na(p) & p > 0, -log10(p), NA_real_)
      star_fn <- function(p) {
        ifelse(is.na(p), "",
               ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*", ""))))
      }

      # Row order: ANOVA, pairwise (3), Kendall
      row_names <- c("ANOVA", pair_names, "Kendall \u03c4")
      all_p <- list(anova_p, pair_p_list[[1]], pair_p_list[[2]], pair_p_list[[3]], kendall_p)

      z_mat <- do.call(rbind, lapply(all_p, to_z))
      text_mat <- do.call(rbind, lapply(all_p, star_fn))

      # ---- Per-CELL coverage, so a cell states what it rests on ----
      # Coverage is group-dependent at the large-diameter end (bin 225 was 4 ND
      # islets from 2 donors vs 14 Aab+ from 5) yet every cell rendered with the
      # same visual weight as bin 25's 4,282 islets from 15 donors. The
      # granularity is genuinely per-cell, not per-column: for a pairwise row
      # only the two groups in that contrast matter.
      cov <- st$bin_coverage
      cov_for <- function(mid, groups) {
        if (is.null(cov)) return("")
        sub <- cov[cov$diam_mid == mid & as.character(cov$donor_status) %in% groups, , drop = FALSE]
        if (!nrow(sub)) return("no donors in this bin")
        paste(sprintf("%s: %d donor%s, %d islet%s", as.character(sub$donor_status),
                      sub$n_donors, ifelse(sub$n_donors == 1, "", "s"),
                      sub$n_islets, ifelse(sub$n_islets == 1, "", "s")),
              collapse = "<br>")
      }
      all_groups <- c("ND", "Aab+", "T1D")
      row_groups <- list(all_groups,                       # ANOVA
                         c("ND", "Aab+"), c("ND", "T1D"), c("Aab+", "T1D"),  # pairwise
                         all_groups)                       # Kendall
      cover_mat <- do.call(rbind, lapply(row_groups, function(gs) {
        vapply(all_mids, function(m) cov_for(m, gs), character(1))
      }))

      # Bin midpoints only. Stacking "15d / 4,424i" under each midpoint was tried
      # and reverted: at this card's ~400 px width the second line of adjacent
      # bins collides into an unreadable run ("15d / 4,4285 / 216i / 168i / 24i")
      # and the taller tick block pushed the x-axis title off the canvas and over
      # the coverage note below. The per-bin counts live in the cell hover and in
      # the coverage strip under the plot, both of which have room for them.
      mid_labels <- as.character(round(all_mids))

      fonts <- plotly_islet_fonts()
      plot_ly(
        x = mid_labels, y = row_names, z = z_mat,
        type = "heatmap",
        colorscale = list(c(0, "#f7fbff"), c(0.5, "#6baed6"), c(1, "#08306b")),
        colorbar = list(title = list(text = "-log10(q)", font = fonts$colorbar_title),
                        tickfont = fonts$axis_tick),
        text = text_mat, texttemplate = "%{text}",
        textfont = list(size = ISLET_FONT_SIZES$dense_min),
        customdata = cover_mat,
        hovertemplate = paste0("Bin: %{x} \u00b5m<br>Test: %{y}<br>-log10(q): %{z:.2f}",
                               "<br><br>%{customdata}<extra></extra>")
      ) %>%
        layout(
          font = fonts$global,
          # Layout notes \u2014 this plot lives in a col-4 (~400 px) card, which is
          # the binding constraint and the reason three things were clipped:
          #   * tickangle 0, NOT -45. The tick labels became two lines when the
          #     donor/islet counts were added; rotating a two-line label made
          #     them 76 px tall, which swallowed the bottom margin and pushed
          #     the x-axis title 23 px below the 520 px canvas. Horizontal,
          #     5-7 short labels fit fine across 400 px.
          #   * explicit `margin`, not automargin alone. automargin does not
          #     measure HTML-styled multi-line ticks reliably, and it never
          #     reserves room for the plot title.
          #   * no in-plot title. "Stratified Significance (BH-corrected,
          #     donor-level)" rendered 453 px wide in a 401 px plot and was
          #     clipped at BOTH ends. The card h5 already names the plot and
          #     the caption already states the BH correction, so the in-plot
          #     title was redundant as well as broken.
          margin = list(t = 20, b = 105, l = 70, r = 20),
          # type = "category" is REQUIRED. The tick labels are numeric-looking
          # strings ("25", "75", ...), so plotly auto-detects a linear axis and
          # draws its own ticks (0, 100, 200) that no longer line up one-to-one
          # with the heatmap columns. Forcing category gives exactly one tick per
          # bin. (This was masked while the labels carried HTML, which is never
          # numeric.)
          xaxis = list(title = list(text = "Diameter bin midpoint (\u00b5m)",
                                    font = fonts$axis_title),
                       type = "category",
                       # -45 is fine now that the labels are ONE line. Five 3-digit
                       # midpoints collide horizontally in this ~400 px card when
                       # laid flat ("125175225"). The earlier clipping came from
                       # rotating a TWO-line label, not from the rotation itself.
                       tickfont = fonts$axis_tick, tickangle = -45, automargin = TRUE),
          yaxis = list(title = list(text = "", font = fonts$axis_title),
                       tickfont = fonts$axis_tick, automargin = TRUE,
                       categoryorder = "array", categoryarray = rev(row_names))
        ) %>%
        # "show", not the default "hide": these are single-glyph significance
        # stars. A hidden star reads as "not significant" \u2014 the one failure mode
        # this plot cannot afford.
        plotly_islet_config("statistics_bin_heatmap", uniform_text_mode = "show")
    })

    # Name the bins that exist in the data but were never tested. They used to
    # just be absent from the heatmap, which reads as "there are no islets that
    # big" rather than "the coverage rule excluded them" — and because ND runs
    # out of large islets first, the size range actually compared was being
    # capped by the narrowest group, silently.
    output$bin_coverage_note <- renderUI({
      st <- tryCatch(stats_run(), error = function(e) NULL)
      if (is.null(st)) return(NULL)
      bd <- st$bin_dropped

      # Per-bin donor / islet counts. These were tried as a second line on each
      # x-axis tick and were unreadable at this card's width (adjacent bins'
      # counts collided into one run and shoved the axis title off the canvas).
      # As a collapsible text strip they get the full column width to lay out in.
      cov <- st$bin_coverage
      cov_strip <- if (!is.null(cov) && nrow(cov) > 0) {
        per_bin <- lapply(sort(unique(cov$diam_mid)), function(m) {
          sub <- cov[cov$diam_mid == m, , drop = FALSE]
          tags$li(tags$strong(sprintf("%g µm: ", m)),
                  paste(sprintf("%s %dd/%di", as.character(sub$donor_status),
                                sub$n_donors, sub$n_islets), collapse = ", "))
        })
        tags$details(style = "margin-top: 10px; font-size: 11px",
          tags$summary(style = "cursor: pointer; font-weight: 600",
                       "Donors / islets per bin"),
          tags$ul(style = "margin: 4px 0 0 0; padding-left: 18px; line-height: 1.5",
                  per_bin))
      } else NULL

      if (is.null(bd) || nrow(bd) == 0) {
        return(tagList(
          cov_strip,
          tags$p(style = "font-size: 11px; margin: 6px 0 0 0",
                 "Every diameter bin present in this subset met the coverage rule ",
                 "(≥ 2 donor-status groups with ≥ 2 donors each).")))
      }
      tagList(cov_strip,
      tags$div(style = paste0("margin-top: 8px; padding: 8px 10px; background: #fbfcfd;",
                              " border: 1px solid #e0e0e0; border-radius: 4px;",
                              " font-size: 11px; color: #666; line-height: 1.5;"),
        tags$strong(sprintf("%d diameter bin%s without an omnibus (ANOVA / Kendall) row:",
                            nrow(bd), if (nrow(bd) == 1) "" else "s")),
        tags$ul(style = "margin: 4px 0 0 0; padding-left: 18px",
          lapply(seq_len(nrow(bd)), function(i) {
            tags$li(
              tags$strong(sprintf("%g µm", bd$mid[i])),
              sprintf(" — %d islet%s (%s). %s. ",
                      bd$n_islets[i], if (bd$n_islets[i] == 1) "" else "s",
                      bd$detail[i], bd$reason[i]),
              tags$em(bd$status[i]))
          })),
        tags$p(style = "margin: 6px 0 0 0; font-style: italic",
          "The compared size range is therefore capped by whichever group runs out of ",
          "islets first — usually ND at the large end.")))
    })

    # ===========================================================================
    # Card 4: Trend Analysis
    # ===========================================================================

    output$trend_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)
      bk <- st$bin_kendall
      if (is.null(bk) || nrow(bk) == 0) return(NULL)

      bk$sig_label <- ifelse(!is.na(bk$p_adj) & bk$p_adj < st$alpha, "Significant", "NS")

      g <- ggplot(bk, aes(x = mid, y = tau, color = sig_label)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#999") +
        geom_line(color = "#333", alpha = 0.4) +
        geom_point(size = 3) +
        scale_color_manual(values = c("Significant" = "#d62728", "NS" = "#aaaaaa"), name = "") +
        labs(x = "Diameter bin midpoint (\u00b5m)", y = "Kendall \u03c4") +
        theme_islet() +
        theme(legend.position = "none")

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(font = fonts$global,
               title = list(font = fonts$title),
               xaxis = list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick),
               yaxis = list(title = list(font = fonts$axis_title), tickfont = fonts$axis_tick),
               margin = list(b = 110),
               legend = list(orientation = "h", x = 0.5, xanchor = "center",
                             y = -0.3, font = fonts$legend_text,
                             title = list(font = fonts$legend_title)),
               showlegend = TRUE) %>%
        plotly_islet_config("statistics_trend_plot")
    })

    # ===========================================================================
    # Card 5: Demographics (conditional — NULL when no age/sex)
    # ===========================================================================

    output$demographics_card <- renderUI({
      st <- stats_run()
      # Show placeholder before first run
      if (is.null(st)) {
        rdf <- raw_df()
        if (is.null(rdf) || !("age" %in% colnames(rdf) || "gender" %in% colnames(rdf))) return(NULL)
        return(div(class = "card", style = "padding: 15px; margin-bottom: 15px",
                   h5("Demographics Analysis"),
                   tags$p("Run statistics to see demographic analysis.")))
      }
      if (is.null(st$demo_summary)) return(NULL)

      # Build AAb section conditionally
      aab_section <- NULL
      if (!is.null(st$aab_summary) && nrow(st$aab_summary) > 0) {
        aab_section <- div(style = "margin-top: 18px",
          h6("Autoantibody Profile (Aab+ only)", style = "margin-bottom: 10px"),
          fluidRow(
            column(6,
              tableOutput(ns("aab_donor_table"))
            ),
            column(6,
              plotlyOutput(ns("aab_plot"), height = "380px")
            )
          )
        )
      }

      div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible",
        h5("Demographics Analysis"),
        div(style = "margin-bottom: 18px",
          h6("Donor Summary", style = "margin-bottom: 10px"),
          tableOutput(ns("demo_summary_table"))
        ),
        div(style = "margin-bottom: 18px",
          fluidRow(
            column(6,
              h6("Age vs Feature (islet-level)", style = "margin-bottom: 10px"),
              plotlyOutput(ns("age_scatter"), height = "420px")
            ),
            column(6,
              h6("Sex vs Feature", style = "margin-bottom: 10px"),
              plotlyOutput(ns("gender_plot"), height = "420px")
            )
          )
        ),
        aab_section,
        div(style = "margin-top: 18px",
          uiOutput(ns("covariate_results"))
        )
      )
    })

    output$demo_summary_table <- renderTable({
      st <- stats_run()
      if (is.null(st) || is.null(st$demo_summary)) return(NULL)
      ds <- st$demo_summary
      display <- data.frame(
        `Donor Status` = ds$donor_status,
        `N Islets` = ifelse(is.na(ds$n_islets), "-", formatC(ds$n_islets, format = "d", big.mark = ",")),
        `N Donors` = ds$n,
        `Age Median` = ifelse(is.na(ds$age_median), "-", sprintf("%.0f", ds$age_median)),
        `Age Range` = ifelse(nzchar(ds$age_range), ds$age_range, "-"),
        `% Male` = ifelse(is.na(ds$pct_male), "-", sprintf("%.0f%%", ds$pct_male)),
        `% Female` = ifelse(is.na(ds$pct_female), "-", sprintf("%.0f%%", ds$pct_female)),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      display
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$age_scatter <- renderPlotly({
      st <- stats_run()
      if (is.null(st) || is.null(st$age_corr)) return(NULL)
      rdf_plot <- st$rdf
      if (!("age" %in% colnames(rdf_plot))) return(NULL)

      rdf_plot$age_num <- as.numeric(rdf_plot$age)
      rdf_plot <- rdf_plot[is.finite(rdf_plot$age_num) & is.finite(rdf_plot$value), , drop = FALSE]
      if (nrow(rdf_plot) == 0) return(NULL)

      rdf_plot$donor_status <- factor(rdf_plot$donor_status, levels = c("ND", "Aab+", "T1D"))

      corr_label <- sprintf("r = %.2f, p = %s (N=%d islets; correlated within donors)",
                            st$age_corr$cor,
                            format_p_value(st$age_corr$p_value),
                            st$age_corr$n)

      g <- ggplot(rdf_plot, aes(x = age_num, y = value, color = donor_status)) +
        geom_point(alpha = 0.4, size = 1.5) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, aes(group = 1), color = "#333") +
        scale_color_manual(values = donor_colors()) +
        labs(x = "Donor Age (years)", y = st$feature_name, color = NULL,
             title = "Age vs Feature (all islets)",
             subtitle = corr_label) +
        theme_islet() +
        theme(legend.position = "bottom")

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(font = fonts$global,
               # ggplotly drops labs(subtitle=), so `corr_label` — the whole
               # point of this plot — never reached the page. Re-attach it here.
               title = list(text = plotly_title_with_sub("Age vs Feature (all islets)",
                                                         corr_label),
                            font = fonts$title),
               margin = list(t = 80),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               legend = list(orientation = "h", x = 0.2, y = -0.2,
                             font = fonts$legend_text,
                             title = list(font = fonts$legend_title))) %>%
        plotly_islet_config("statistics_age_scatter")
    })

    output$gender_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)
      rdf_plot <- st$rdf
      if (!("gender" %in% colnames(rdf_plot))) return(NULL)

      rdf_plot$gender <- as.character(rdf_plot$gender)
      rdf_plot <- rdf_plot[!is.na(rdf_plot$gender) & nzchar(rdf_plot$gender) &
                           is.finite(rdf_plot$value), , drop = FALSE]
      if (nrow(rdf_plot) == 0) return(NULL)

      rdf_plot$donor_status <- factor(rdf_plot$donor_status, levels = c("ND", "Aab+", "T1D"))
      rdf_plot$gender <- factor(rdf_plot$gender, levels = c("M", "F"), labels = c("Male", "Female"))

      # Build subtitle with stratified p-values (if available).
      # p_adj is BH-corrected across the 6 contrasts (3 pairs x 2 sexes); we
      # report the minimum adjusted p per sex so users see whether ANY contrast
      # survived the family-wise correction.
      gs <- st$gender_strat
      sub_label <- ""
      if (!is.null(gs) && nrow(gs) > 0) {
        p_labels <- c()
        for (gen in unique(as.character(gs$gender))) {
          sub_gs <- gs[gs$gender == gen, , drop = FALSE]
          if (nrow(sub_gs) > 0) {
            min_padj_col <- if ("p_adj" %in% colnames(sub_gs)) "p_adj" else "p_value"
            min_p <- min(sub_gs[[min_padj_col]], na.rm = TRUE)
            lp <- sub_gs$low_power[1]
            suffix <- if (isTRUE(lp)) " (low power)" else ""
            p_labels <- c(p_labels, sprintf("%s: min q=%s%s", gen,
                                            format_p_value(min_p), suffix))
          }
        }
        sub_label <- paste0(paste(p_labels, collapse = "  |  "),
                            "   (BH-corrected across 6 contrasts)")
      }

      g <- ggplot(rdf_plot, aes(x = donor_status, y = value, fill = donor_status)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        facet_wrap(~ gender) +
        scale_fill_manual(values = donor_colors()) +
        labs(x = NULL, y = st$feature_name,
             title = "Feature by Donor Status & Sex",
             subtitle = sub_label) +
        theme_islet() +
        theme(legend.position = "none")

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = c("y")) %>%
        # This is the app's only facet. ggplotly turns the "Male"/"Female" strip
        # labels into layout annotations that neither layout(font=) nor
        # theme_islet()'s strip.text reaches, so they rendered at plotly's ~11px
        # annotation default.
        plotly_islet_strip_fonts() %>%
        layout(font = fonts$global,
               title = list(text = plotly_title_with_sub("Feature by Donor Status & Sex",
                                                         sub_label),
                            font = fonts$title),
               margin = list(b = 60, t = 80),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE)) %>%
        plotly_islet_config("statistics_sex_box")
    })

    output$aab_donor_table <- renderTable({
      st <- stats_run()
      if (is.null(st) || is.null(st$aab_donor_summary)) return(NULL)
      ads <- st$aab_donor_summary
      # Reorder: Case ID, individual AAbs, n_aab, n_islets, mean_value
      aab_short <- sub("^AAb_", "", intersect(c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA"),
                                               colnames(st$rdf)))
      base_cols <- c("Case ID", aab_short, "n_aab", "n_islets")
      base_cols <- intersect(base_cols, colnames(ads))
      display <- ads[, base_cols, drop = FALSE]
      display$`Mean Value` <- sprintf("%.3f", ads$mean_value)
      colnames(display)[colnames(display) == "n_aab"] <- "Total AAb"
      colnames(display)[colnames(display) == "n_islets"] <- "N Islets"
      display
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$aab_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st) || is.null(st$aab_summary) || nrow(st$aab_summary) == 0) return(NULL)

      aab_df <- st$aab_summary
      aab_df <- aab_df[is.finite(aab_df$value), , drop = FALSE]
      if (nrow(aab_df) == 0) return(NULL)

      # Group by AAb count (1, 2, 3+)
      aab_df$aab_group <- ifelse(aab_df$n_aab >= 3, "3+", as.character(aab_df$n_aab))
      aab_df$aab_group <- factor(aab_df$aab_group, levels = c("0", "1", "2", "3+"))
      # Drop levels with no data
      aab_df$aab_group <- droplevels(aab_df$aab_group)

      if (length(levels(aab_df$aab_group)) < 2) {
        # Only one AAb count level — not informative
        return(NULL)
      }

      # Kruskal-Wallis on DONOR means (computed in stats_run as st$aab_kw).
      # This box plot shows the islet-level spread — but the AAb count is a donor
      # attribute, so a test across islets is pseudoreplicated: on this cohort it
      # returned p = 1.6e-15 from 1,666 islets belonging to 5 donors, where the
      # donor-level test gives p = 0.67. The subtitle now states both the level
      # the test was run at and how many donors sit behind each box.
      kw <- st$aab_kw
      p_label <- if (!is.null(kw) && is.finite(kw$p_value)) {
        sprintf("Kruskal-Wallis on donor means: p = %s.", format_p_value(kw$p_value))
      } else "Donor-level test not computable (fewer than 2 AAb-count groups)."
      donors_label <- if (!is.null(kw) && length(kw$donors_per_group)) {
        paste0(" Donors per group: ",
               paste(sprintf("%s=%d", names(kw$donors_per_group),
                             as.integer(kw$donors_per_group)), collapse = ", "), ".")
      } else ""
      power_label <- if (!is.null(kw) && isTRUE(kw$underpowered)) {
        " At least one group is a single donor — treat as descriptive only."
      } else ""

      aab_sub <- paste0("Boxes show ", nrow(aab_df), " islets from ",
                        length(unique(aab_df$`Case ID`)), " Aab+ donors. ",
                        p_label, donors_label, power_label)

      g <- ggplot(aab_df, aes(x = aab_group, y = value, fill = aab_group)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        scale_fill_brewer(palette = "Set2") +
        labs(x = "Number of Autoantibodies", y = st$feature_name,
             title = "Feature by AAb Count (Aab+ islets)",
             subtitle = aab_sub) +
        theme_islet() +
        theme(legend.position = "none")

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = c("y")) %>%
        layout(font = fonts$global,
               # ggplotly drops labs(subtitle=) — re-attach the N + KW p-value.
               title = list(text = plotly_title_with_sub("Feature by AAb Count (Aab+ islets)",
                                                         aab_sub),
                            font = fonts$title),
               margin = list(t = 80),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE)) %>%
        plotly_islet_config("statistics_aab_plot")
    })

    output$covariate_results <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$age_model_p)) return(NULL)
      p_cov <- st$age_model_p
      survived <- !is.na(p_cov) && p_cov < st$alpha
      # "Survives adjustment" only means something if there was an effect to
      # survive. When the UNADJUSTED donor-level test was already null, the old
      # wording ("does NOT survive adjustment for age and sex") read as though
      # adjustment had destroyed a real effect -- e.g. Capillary, 0.348
      # unadjusted -> 0.387 adjusted, where nothing was ever there. Branch on
      # the unadjusted result so the three cases read differently.
      p_unadj <- st$donor_p
      had_effect <- !is.na(p_unadj) && p_unadj < st$alpha
      tags$div(style = "margin-top: 10px; padding: 10px",
        tags$strong("Covariate-Adjusted Model (donor-level)"),
        tags$p(style = "margin: 5px 0 0 0",
          paste0("lm(value ~ donor_status + age + sex, data = donor_means): donor_status p = ",
                 format_p_value(p_cov)),
          tags$br(),
          if (!had_effect) {
            tags$span(
                      paste0("No effect to adjust away — the unadjusted donor-level ",
                             st$donor_test, " was already non-significant (p = ",
                             format_p_value(p_unadj), "). This model neither ",
                             "supports nor refutes a disease-stage effect."))
          } else if (survived) {
            tags$span(style = "font-weight: 600",
                      "Donor status effect SURVIVES adjustment for age and sex.")
          } else {
            tags$span(
                      paste0("Donor status effect does NOT survive adjustment for age ",
                             "and sex (unadjusted p = ", format_p_value(p_unadj),
                             "), so age/sex may account for it."))
          }
        )
      )
    })

    # ===========================================================================
    # Card 6: AUC Analysis
    # ===========================================================================

    output$auc_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st) || is.null(st$auc_by_group) || nrow(st$auc_by_group) == 0) return(NULL)

      auc_df <- st$auc_by_group
      auc_df$donor_status <- factor(auc_df$donor_status, levels = c("ND", "Aab+", "T1D"))

      g <- ggplot(auc_df, aes(x = donor_status, y = auc, fill = donor_status)) +
        geom_col(width = 0.7, color = "black", alpha = 0.8) +
        scale_fill_manual(values = donor_colors()) +
        # No in-plot title: at 18 pt "Integrated Area by Donor Group" measured
        # 340 px inside this 401 px col-4 card and was clipped at both ends. The
        # card's h5 ("Area Under Curve by Donor Group") already names it.
        labs(x = "Donor Status", y = "Area Under Curve (AUC)") +
        theme_islet() +
        theme(legend.position = "none")

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = c("x", "y")) %>%
        # ggplotly does not reserve room for the axis titles at this font size in
        # a narrow card: the x-axis title rendered 28 px BELOW the 380 px canvas.
        # Explicit margins, same fix as the bin heatmap beside it.
        layout(font = fonts$global,
               margin = list(t = 20, b = 95, l = 80, r = 20),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE)) %>%
        plotly_islet_config("statistics_auc_plot")
    })

    output$auc_table <- renderTable({
      st <- stats_run()
      if (is.null(st) || is.null(st$auc_by_group) || nrow(st$auc_by_group) == 0) return(NULL)

      auc_df <- st$auc_by_group

      result <- data.frame(
        "Donor Group" = as.character(auc_df$donor_status),
        "AUC" = sprintf("%.2f", auc_df$auc),
        "N Bins" = as.character(auc_df$n_bins),
        check.names = FALSE, stringsAsFactors = FALSE
      )

      # Add T1D vs ND percentage change row (with division-by-zero guard)
      nd_auc <- auc_df$auc[auc_df$donor_status == "ND"]
      t1d_auc <- auc_df$auc[auc_df$donor_status == "T1D"]
      if (length(nd_auc) > 0 && length(t1d_auc) > 0 && is.finite(nd_auc) && nd_auc != 0) {
        pct_change <- ((t1d_auc - nd_auc) / nd_auc) * 100
        result <- rbind(result, data.frame(
          "Donor Group" = "T1D vs ND",
          "AUC" = sprintf("%.1f%%", pct_change),
          "N Bins" = "",
          check.names = FALSE, stringsAsFactors = FALSE
        ))
      }

      result
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$auc_interpretation <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$auc_by_group) || nrow(st$auc_by_group) == 0) return(NULL)

      auc_df <- st$auc_by_group
      nd_auc <- auc_df$auc[auc_df$donor_status == "ND"]
      t1d_auc <- auc_df$auc[auc_df$donor_status == "T1D"]

      if (length(nd_auc) > 0 && length(t1d_auc) > 0 && is.finite(nd_auc) && nd_auc != 0) {
        pct <- ((t1d_auc - nd_auc) / nd_auc) * 100
        direction <- if (pct < 0) "decrease" else "increase"
        tags$p(style = "margin-top: 10px; font-size: 90%",
          sprintf("The AUC integrates the binned mean feature value across islet diameters (trapezoidal rule). T1D shows a %.1f%% %s relative to ND.",
                  abs(pct), direction))
      } else {
        NULL
      }
    })

    # ===========================================================================
    # Card 7: Methods & Interpretation
    # ===========================================================================

    output$stats_explanation <- renderUI({
      st <- stats_run()
      if (is.null(st)) {
        return(tags$div(style = "padding: 10px",
          tags$p("Click 'Run Statistics' to perform analysis. ",
                 "The sidebar controls which feature and subset of data are analyzed.")))
      }

      alpha_str <- as.character(st$alpha)
      test_desc <- if (st$is_parametric) {
        paste0("Auto-selected parametric (all donor-mean groups passed Shapiro-Wilk): ",
               "one-way ANOVA (global), pairwise t-tests (BH-corrected).")
      } else {
        paste0("Auto-selected non-parametric (at least one group failed Shapiro-Wilk, ",
               "or normality could not be established): Kruskal-Wallis (global), ",
               "pairwise Wilcoxon rank-sum (BH-corrected).")
      }

      outlier_para <- {
        os <- st$outlier_summary
        local_thr <- {
          v <- suppressWarnings(as.numeric(input$stats_outlier_z))
          if (is.finite(v) && v >= 0.5) v else 3
        }
        if (is.null(os)) {
          tagList(
            tags$strong("Outlier handling: "),
            sprintf("For each donor status (ND / Aab+ / T1D), islets whose robust |z| (median/MAD modified z-score) exceeds %g are excluded from all donor-level aggregation, mixed-effects modelling, and stratified tests. Detection is two-sided — both high and low extremes are flagged. When MAD = 0 (zero-inflated groups) the rule falls back to classical mean/SD. ", local_thr),
            "The control is local to the Statistics tab (Section 1 ", tags$em("Exclude outliers"),
            " checkbox + threshold slider) and operates independently of the Plot / Trajectory / Spatial synced toggle. ",
            "On the Plot and Trajectory scatters, outliers still render as faint grey × markers when their own filter is on.",
            tags$br()
          )
        } else if (isTRUE(os$enabled)) {
          tagList(
            tags$strong("Outlier handling: "),
            sprintf("This run excluded %d islet%s (%s, robust median/MAD; mean/SD fallback when MAD = 0) before donor-level aggregation, mixed-effects modelling, and stratified tests. ",
                    os$n_excluded, if (os$n_excluded == 1) "" else "s", os$label),
            "The control is local to the Statistics tab (Section 1 ", tags$em("Exclude outliers"),
            " checkbox + threshold slider); the Plot / Trajectory / Spatial tabs share their own synced controls. ",
            "Spatial Card A uses Tukey 1.5×IQR per donor status, since enrichment z is already a standardised quantity.",
            tags$br()
          )
        } else {
          tagList(
            tags$strong("Outlier handling: "),
            sprintf("Filter is OFF for this run — all islets contribute to the tests below, including extreme values. Re-enable the "),
            tags$em("Exclude outliers"), sprintf(" checkbox in Section 1 to apply the |z| > %g (per donor status, median/MAD) rule.", local_thr),
            tags$br()
          )
        }
      }

      zero_para <- if (isTRUE(st$exclude_zero)) {
        tagList(
          tags$strong("Zero-value exclusion: "),
          sprintf("This run dropped %s islet%s with value == 0 before donor-level aggregation, mixed-effects modelling, demographics, and per-bin tests. ",
                  formatC(st$zero_dropped, format = "d", big.mark = ","),
                  if (st$zero_dropped == 1) "" else "s"),
          "Useful for zero-inflated features — rare phenotypes whose core or peri-islet proportion is exactly 0 in most islets. ",
          tags$strong("Note: "), "AUC is recomputed from this same filtered subset, so it moves with this filter too.",
          tags$br()
        )
      } else {
        tagList(
          tags$strong("Zero-value exclusion: "),
          "Filter is OFF for this run — all finite values contribute to the tests, including zero-valued islets. ",
          "Toggle ", tags$em("Exclude zero-valued islets"),
          " in Section 1 to re-run on the non-zero subset (recommended for zero-inflated features).",
          tags$br()
        )
      }

        ieCore::ie_explain(
        tags$p(style = paste0("background: #eef6ee; border-left: 3px solid #2e7d32; ",
                              "padding: 8px 12px; border-radius: 4px;"),
          tags$strong("Design safeguards against p-hacking. "),
          "This tab reports a single pre-specified primary test — the mixed-effects model. ",
          "The significance threshold is fixed at α = 0.05, the parametric vs. non-parametric ",
          "form of the sensitivity test is auto-selected from normality (not a user toggle), and ",
          "the size-stratified analysis uses a fixed ", STATS_BIN_WIDTH_UM, " µm bin width. ",
          "Optional data-subset filters (outliers, zero-values, min-cells, diameter range) are ",
          "confined to an Advanced panel and recorded in the exported CSV, so any robustness ",
          "probing is explicit and auditable rather than a hidden researcher degree of freedom."),
        tags$p(outlier_para),
        tags$p(zero_para),
        tags$p(
          tags$strong("Primary test (mixed-effects model): "),
          paste0("A linear mixed-effects model tests whether donor status (ND / Aab+ / T1D) predicts the feature value across all ",
                 formatC(st$n_islets, big.mark = ","),
                 " islets, with donor included as a random intercept. ",
                 "This accounts for the fact that islets from the same donor share biology, tissue processing, and imaging conditions, ",
                 "while still leveraging the full islet-level variation — important because islets of different sizes are ",
                 "structurally and functionally distinct. ",
                 "Pairwise group comparisons use estimated marginal means from the mixed model, ",
                 "with Satterthwaite degrees of freedom \u2014 the same approximation lmerTest uses for ",
                 "the omnibus p above. (emmeans defaults to asymptotic z-tests once the model exceeds ",
                 "3,000 observations, which is anti-conservative for a 15-donor design; that limit is ",
                 "raised explicitly so the method does not change with dataset size.) ",
                 "This run used: ", st$lmer_pairwise_df_method %||% "N/A", "."),
          tags$br(), tags$br(),
          tags$strong("ICC (intra-class correlation): "),
          "Measures what fraction of total variation is explained by donor identity. ",
          "Low ICC means most variation is between islets within the same donor, supporting islet-level analysis. ",
          "High ICC means donor identity dominates the signal.",
          tags$br(), tags$br(),
          tags$strong("Donor-level sensitivity check: "), test_desc,
          paste0(" As a cross-check, all islets within each donor are averaged to one value per donor (N = ", st$n_donors,
                 "), then tested directly. This is deliberately conservative — it collapses size-dependent variation ",
                 "to a single mean per donor. If both the mixed-effects model and this donor-level test agree, ",
                 "the result is robust."),
          tags$br(),
          tags$strong("Significance threshold: "), paste0("\u03b1 = ", alpha_str, " (fixed)"),
          tags$br(),
          tags$strong("Effect sizes: "),
          "The forest plot reports the ", tags$strong("unstandardised difference"),
          " between donor-level group means, in the feature's own units, with a Welch 95% CI ",
          "(no equal-variance assumption). Cohen's d and the Glass's \u0394 range are reported ",
          "alongside in the Sensitivity check table. d divides by a ", tags$em("pooled"),
          " SD, which is only a meaningful common scale when the groups have comparable ",
          "spread \u2014 on this cohort they frequently do not (floor effects can leave one group ",
          "at ~0 with almost no variance, giving 50-fold variance ratios). Where that happens ",
          "the contrast is marked \u2020 and d should be read as indicative only; the raw ",
          "difference is the quantity that does not depend on which SD you standardise by. ",
          "Per-group means and SDs are shown in the primary card so this is checkable.",
          tags$br(), tags$br(),
          tags$strong("Islet size (companion models): "),
          "The primary model conditions on neither islet diameter nor a group-varying size ",
          "slope, but both matter: diameter is a strong predictor of most features, and the ",
          "diameter distribution itself differs between donor groups. Two companions are ",
          "therefore fitted on the same data and reported next to the primary p \u2014 ",
          tags$code("value ~ donor_status + diameter + (1|donor)"), " (is the group difference ",
          "still there at a common islet size?) and ",
          tags$code("value ~ donor_status \u00d7 diameter + (1|donor)"),
          " (does the group difference depend on size?). Diameter is mean-centred, so the ",
          "adjusted donor-status term is the difference at the mean islet diameter. ",
          tags$strong("When the interaction is significant, the single omnibus p is a marginal ",
          "average over a size distribution that differs between groups"),
          " \u2014 the size-stratified analysis in Section 3 is then part of the primary reading, ",
          "not an optional extra. Neither companion replaces the pre-specified primary test; ",
          "they are reported to make its scope explicit.",
          tags$br(), tags$br(),
          tags$strong("Age (companion model): "),
          "5 of the 15 donors are under 5 years old, so a companion model fits ",
          tags$code("value ~ donor_status + age_group + (1|donor)"),
          " following Steenbuck et al. 2025, who used the same form to test whether ",
          "disease-stage differences survive adjustment for age-associated endotypes. ",
          "Their <13 / \u226513 year cut is data-supported on this cohort rather than ",
          "arbitrary: our ages are bimodal with a gap from 4.22 to 14.56 years, and 13 ",
          "falls inside it. Age comes from the H5AD, which is canonical for modelling ",
          "(it disagrees with the nPOD record for one donor, both under the cut). ",
          tags$strong("HbA1c and C-peptide are deliberately NOT adjusted for: "),
          "both are downstream consequences of beta-cell loss rather than confounders ",
          "of it, so conditioning on them is mediator adjustment and would bias the ",
          "disease effect toward the null. They appear in the demographics table and ",
          "are useful as external-validity checks, not as covariates. Disease duration ",
          "is likewise a cohort descriptor only \u2014 all T1D donors here are ",
          "recent-onset (<1 year), which leaves it no variance to model.",
          tags$br(),
          tags$strong("Normality: "), "Shapiro-Wilk test on donor-level means per group. ",
          "With n = 5, this test has limited power, so failure to reject normality does not guarantee it.",
          tags$br(),
          tags$strong("Stratified analysis (exploratory): "), "ANOVA, pairwise tests, and Kendall \u03c4 on donor-level means within each fixed 50 \u00b5m diameter bin. ",
          "Kendall \u03c4 correlates disease stage (ND=0, Aab+=1, T1D=2) with the feature value: ",
          "\u03c4 > 0 means the feature increases with disease progression, \u03c4 < 0 means it decreases. ",
          "P-values are BH-corrected across bins. A bin is tested only if at least 2 donor-status groups ",
          "each have \u2265 2 donors with finite values (same rule for the heatmap and the trend plot). ",
          "The same rule gates the per-bin ANOVA, Kendall \u03c4 AND the pairwise contrasts: a pairwise ",
          "post-hoc is never reported for a bin that has no omnibus test behind it. ",
          "Coverage is group-dependent \u2014 the large-diameter bins can rest on a handful of islets ",
          "from two donors while the smallest bin carries thousands from every donor \u2014 so each ",
          "bin's donor / islet count is on the axis label, the per-group coverage behind each ",
          "individual cell is in its hover, and bins that failed the rule are listed under the ",
          "plot with the reason rather than silently omitted.",
          tags$br(),
          tags$strong("Multiple comparisons: "),
          "Benjamini-Hochberg (FDR) correction is applied ", tags$em("within"),
          " each test family — donor-level pairwise contrasts (3 pairs), per-bin ANOVA / Kendall / ",
          "pairwise (across bins), and sex-stratified pairwise (6 contrasts). It is ",
          tags$strong("not"), " applied across families, nor across the many features you might test. ",
          "The Section 3–4 exploratory views are therefore hypothesis-generating, not confirmatory; ",
          "confirm any finding with the Section 2 primary test on a single, pre-specified feature.",
          tags$br(),
          tags$strong("AUC (exploratory): "), paste0("Area under curve computed by trapezoidal integration of the per-status mean feature value across fixed ", STATS_BIN_WIDTH_UM, " \u00b5m islet-diameter bins, on the same filtered subset as every test above. Descriptive summary \u2014 no inferential test.")
        ),
        if (!is.null(st$demo_summary)) {
          tags$p(style = "margin-top: 8px",
            tags$strong("Demographic covariates: "),
            "Age-feature scatter shows all islets (note: correlated within donors). ",
            "Sex analysis uses box plots by donor status, faceted by sex. ",
            "Sex-stratified pairwise tests are BH-corrected across all 6 contrasts (3 pairs × 2 sexes); ",
            "the subtitle reports the minimum adjusted q per sex (low power with ~2-3 donors per sex per group). ",
            if (!is.null(st$age_model_p)) {
              "A linear model on donor means adjusting for age and sex tests whether donor_status retains significance."
            } else {
              ""
            },
            if (!is.null(st$aab_summary)) {
              paste0(" Autoantibody analysis shows the islet-level feature distribution by AAb count (number of positive AAbs) for Aab+ donors only; the accompanying Kruskal-Wallis p is computed on DONOR means, because AAb count is a donor attribute and testing it across islets is pseudoreplication.")
            } else {
              ""
            }
          )
        },
        tags$p(style = "font-size: 90%; margin-top: 8px",
          tags$strong("Interpretation: "),
          "A significant result (p < \u03b1) suggests the feature differs between donor groups. ",
          # The donor counts come from the analysis that just ran, not from a literal: this text
          # used to say "N = 5 per group", which is a fact about the CODEX cohort and is wrong
          # for any other dataset the app serves.
          sprintf("With %s, power is limited; non-significant results should not be over-interpreted. ",
                  local({
                    gs <- st$group_summary
                    if (is.null(gs) || !nrow(gs) || !"n_donors" %in% names(gs)) "these group sizes"
                    else paste0("n = ", paste(sprintf("%d %s", as.integer(gs$n_donors),
                                                      as.character(gs$donor_status)),
                                              collapse = ", "), " donors")
                  })),
          "Effect sizes (Cohen's d) provide magnitude information independent of sample size."
        ),
          label = "How these statistics are computed, and the safeguards"
      )
    })

    # ===========================================================================
    # Download handler
    # ===========================================================================

    output$download_stats <- downloadHandler(
      filename = function() {
        paste0("statistics_", gsub("[^0-9A-Za-z]+", "_", Sys.time()), ".csv")
      },
      content = function(file) {
        # `stats_run` is an eventReactive(..., ignoreInit = TRUE): before the Run
        # button is pressed it raises rather than returning NULL, which crashed
        # the download instead of producing the "no results" fallback below.
        st <- tryCatch(stats_run(), error = function(e) NULL)

        # get_selection_description() belongs to the PLOT module, and since the
        # cohort panel landed it carries two things: the Plot tab's SELECTION
        # (layer, feature, metric) and the panel's COHORT SCOPE, generated by
        # cohort_provenance() from the spec plus the audit. Both filter the data
        # before this tab sees it.
        #
        # Both blocks used to emit a bare "Min cells per islet (core+peri): N"
        # line, so a file could carry two contradictory values for the same named
        # quantity (Plot 50 / Statistics 1) with nothing to say which one bound.
        # Both are explicitly attributed now, and the panel's line cannot claim a
        # filter that was skipped -- cohort_provenance() marks those
        # "*** NOT APPLIED ***". "Bin Width" is dropped: it controls only the
        # Plot scatter and no longer reaches any number here.
        descriptions <- if (!is.null(st)) st$plot_subset_desc else
          tryCatch(get_selection_description(), error = function(e) character(0))
        descriptions <- descriptions[!grepl("^Bin Width:", descriptions)]

        lines <- c("# COHORT SCOPE AND PLOT SELECTION (both applied upstream of this tab)",
                   paste("#  ", descriptions))
        lines <- c(lines, "")

        if (!is.null(st)) {
          # Methodology notes
          lines <- c(lines, "# METHODOLOGY")
          # Record BOTH forms: the readable label matches what is on screen, the
          # raw column is what an analyst needs to reproduce the selection.
          lines <- c(lines, paste0("# Feature: ", st$feature_name,
                                   "  [column: ", st$feature_col %||% "unknown", "]"))
          lines <- c(lines, paste0("# Primary: mixed-effects lmer (islet-level with donor random intercept)"))
          lines <- c(lines, paste0("# Alpha (fixed): ", st$alpha))
          lines <- c(lines, paste0("# Sensitivity test (auto-selected from normality): ", st$donor_test,
                                   " (", if (st$is_parametric) "parametric" else "non-parametric", ")"))
          lines <- c(lines, paste0("# Per-bin width (fixed): ", STATS_BIN_WIDTH_UM, " um"))
          os <- st$outlier_summary
          if (!is.null(os)) {
            lines <- c(lines, paste0("# Outlier rule: ",
                                     if (isTRUE(os$enabled)) os$label else "OFF (filter disabled)",
                                     "; excluded = ", os$n_excluded,
                                     " of ", os$n_total, " input rows."))
          }
          lines <- c(lines, paste0("# Zero-value filter: ",
                                   if (isTRUE(st$exclude_zero)) "ON" else "OFF",
                                   "; N dropped = ", st$zero_dropped,
                                   ". (This tab owns its own zero filter; the Plot",
                                   " sidebar's checkbox does not affect these results.)"))
          lines <- c(lines, paste0("# [Statistics tab] Min cells per islet (core+peri): ",
                                   st$min_cells %||% 1,
                                   if (!isTRUE(st$min_cells_applied))
                                     "  *** NOT APPLIED: per-islet cell counts unavailable for this data source ***"
                                   else ""))
          lines <- c(lines, paste0("# [Statistics tab] Diameter range: ",
                                   if (is.null(st$diam_range)) "no filter"
                                   else paste0(st$diam_range[1], " - ", st$diam_range[2],
                                               " um (intersected with the Plot sidebar range above)")))
          lines <- c(lines, paste0("# Total islets: ",
                                   formatC(st$n_islets, big.mark = ","),
                                   ", Donors: ", st$n_donors))
          lines <- c(lines, paste0("# ICC (intra-class correlation): ",
                                   if (is.na(st$icc)) "N/A" else sprintf("%.4f", st$icc)))
          lines <- c(lines, paste0("# Mixed-effects p-value (primary): ",
                                   if (is.na(st$lmer_p)) "N/A" else formatC(st$lmer_p, format = "e", digits = 4)))
          lines <- c(lines, paste0("# Mixed-effects pairwise d.f. method: ",
                                   st$lmer_pairwise_df_method %||% "N/A"))
          sm <- st$size_models
          if (!is.null(sm)) {
            lines <- c(lines, paste0("# Size companion - diameter effect p: ",
                                     format_p_value(sm$p_diam)))
            lines <- c(lines, paste0("# Size companion - donor_status adjusted for diameter p: ",
                                     format_p_value(sm$p_status_adj)))
            lines <- c(lines, paste0("# Size companion - donor_status x diameter interaction p: ",
                                     format_p_value(sm$p_interaction),
                                     if (is.finite(sm$p_interaction) && sm$p_interaction < st$alpha)
                                       "  *** effect is size-dependent; the omnibus p is a marginal average ***"
                                     else ""))
            lines <- c(lines, paste0("# Size companion - diameter mean-centred at: ",
                                     if (is.finite(sm$mean_diam)) sprintf("%.2f um", sm$mean_diam) else "N/A"))
          }
          am <- st$age_models
          if (!is.null(am)) {
            lines <- c(lines, paste0("# Age companion - donor_status adjusted for age_group p: ",
                                     format_p_value(am$p_status_adj)))
            lines <- c(lines, paste0("# Age companion - age_group (<", am$cut, " / >=", am$cut,
                                     " y) p: ", format_p_value(am$p_age_group),
                                     "   donors: ", am$n_young, " young / ", am$n_old, " old"))
            lines <- c(lines, paste0("# Age companion - NOT adjusted for HbA1c or C-peptide:",
                                     " both are consequences of beta-cell loss, so adjusting",
                                     " for them is mediator adjustment and biases toward null."))
            if (!is.null(am$error_msg)) {
              lines <- c(lines, paste0("# Age companion - note: ", am$error_msg))
            }
          }
          lines <- c(lines, "")

          # Sensitivity test
          lines <- c(lines, "# Sensitivity: Donor-Level Test (conservative)")
          # `format_p_value()` is the single p-formatting rule (see the Jul 2026
          # correctness pass): fixed to 4 dp down to 1e-4, scientific below. A
          # bare formatC(format="e") here rendered p = 0.3476 as "3.48e-01".
          # The effect size is named by its branch -- ANOVA gives eta^2,
          # Kruskal-Wallis gives epsilon^2; they are not the same statistic.
          lines <- c(lines, paste0("# ", st$donor_test, ", p = ",
                                   format_p_value(st$donor_p),
                                   ", ", st$eta_sq_label %||% "eta_sq", " = ",
                                   round(st$eta_sq, 4)))
          lines <- c(lines, "")

          # Write header lines
          writeLines(lines, file)

          # Normality tests
          if (!is.null(st$norm_result)) {
            cat("# Shapiro-Wilk Normality Tests (donor-level)\n", file = file, append = TRUE)
            write.table(st$norm_result, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # `gender` is the DATA column name (H5AD obs) and stays that way
          # in memory; on export it becomes `sex` to match every user-facing
          # label in the app.
          as_sex <- function(df) {
            if (!is.null(df) && "gender" %in% colnames(df)) {
              colnames(df)[colnames(df) == "gender"] <- "sex"
            }
            df
          }

          # Per-group descriptive summary — the spread needed to read the effect sizes
          if (!is.null(st$group_summary)) {
            cat("\n# Group Summary (donor level; SD is across donor means)\n",
                file = file, append = TRUE)
            write.table(st$group_summary, file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }

          # Per-bin coverage + the bins that failed the coverage rule
          if (!is.null(st$bin_coverage)) {
            cat("\n# Per-Bin Coverage (donors and islets per group in each diameter bin)\n",
                file = file, append = TRUE)
            write.table(st$bin_coverage, file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }
          if (!is.null(st$bin_dropped) && nrow(st$bin_dropped) > 0) {
            cat(paste0("\n# Bins Not Tested (coverage rule: >= 2 donor-status groups present AND every\n",
                       "#   group present in the bin having >= 2 donors). ANOVA, Kendall and pairwise all\n",
                       "#   gate on this rule, so these bins contribute nothing to the stratified analysis.\n"),
                file = file, append = TRUE)
            write.table(st$bin_dropped, file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }

          # Donor-level data
          if (!is.null(st$donor_df)) {
            cat("\n# Donor-Level Means\n", file = file, append = TRUE)
            write.table(as_sex(st$donor_df), file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }

          # ---- POST-HOC: mixed-effects pairwise (the PRIMARY test's post-hoc) --
          # emmeans(pairwise ~ donor_status) on the lmer fit. `pairwise ~` applies
          # a Tukey adjustment to p.value by default, so that column IS adjusted.
          # These contrasts are shown on screen but were never exported: only the
          # omnibus p and ICC survived, as header comments — which made the
          # per-contrast primary results unrecoverable from the CSV.
          lpw <- st$lmer_pairwise
          if (!is.null(lpw) && nrow(lpw) > 0) {
            lpw_out <- lpw
            # `p.value` from emmeans pairwise is Tukey-adjusted; name it so.
            if ("p.value" %in% colnames(lpw_out)) {
              colnames(lpw_out)[colnames(lpw_out) == "p.value"] <- "p_tukey"
              lpw_out$sig <- !is.na(lpw_out$p_tukey) & lpw_out$p_tukey < st$alpha
            }
            cat("\n# Post-hoc: Mixed-Effects Pairwise (emmeans, Tukey-adjusted) -- post-hoc for the PRIMARY test\n",
                file = file, append = TRUE)
            write.table(lpw_out, file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }

          # Pairwise
          cat("\n# Post-hoc: Pairwise Comparisons (donor-level, BH-corrected across 3 contrasts)\n",
              file = file, append = TRUE)
          pw <- st$pairwise_df
          if (!is.null(pw) && nrow(pw) > 0) {
            write.table(pw, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Covariate-adjusted model (donor-level). Rendered on screen but never
          # exported.
          if (!is.null(st$age_model_p)) {
            cat("\n# Covariate-Adjusted Model: lm(value ~ donor_status + age + sex), donor-level\n",
                file = file, append = TRUE)
            write.table(
              data.frame(term = "donor_status",
                         p_value = st$age_model_p,
                         sig = !is.na(st$age_model_p) && st$age_model_p < st$alpha,
                         stringsAsFactors = FALSE),
              file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Per-bin ANOVA
          if (!is.null(st$bin_anova) && nrow(st$bin_anova) > 0) {
            cat("\n# Per-Bin ANOVA (donor-level)\n", file = file, append = TRUE)
            write.table(st$bin_anova, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Per-bin Kendall
          if (!is.null(st$bin_kendall) && nrow(st$bin_kendall) > 0) {
            cat("\n# Per-Bin Kendall (donor-level)\n", file = file, append = TRUE)
            write.table(st$bin_kendall, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # ---- POST-HOC: per-bin pairwise ------------------------------------
          # Rows 2-4 of the on-screen bin_heatmap (rows 1 and 5 are the ANOVA and
          # Kendall sections above). Computed, BH-corrected and plotted, but
          # previously absent from the export — so three fifths of that figure
          # could not be reproduced from the CSV.
          if (!is.null(st$bin_pairwise) && nrow(st$bin_pairwise) > 0) {
            cat("\n# Post-hoc: Per-Bin Pairwise (donor-level, BH-corrected across bins)\n",
                file = file, append = TRUE)
            write.table(st$bin_pairwise, file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }

          # AUC
          if (!is.null(st$auc_by_group) && nrow(st$auc_by_group) > 0) {
            cat("\n# AUC by Donor Group\n", file = file, append = TRUE)
            write.table(st$auc_by_group, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Demographics
          if (!is.null(st$demo_summary)) {
            cat("\n# Demographics Summary\n", file = file, append = TRUE)
            write.table(st$demo_summary, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          if (!is.null(st$age_corr)) {
            cat("\n# Age-Feature Correlation (all islets)\n", file = file, append = TRUE)
            write.table(st$age_corr, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          if (!is.null(st$gender_strat) && nrow(st$gender_strat) > 0) {
            cat("\n# Post-hoc: Sex-Stratified Pairwise (BH-corrected across 6 contrasts: 3 pairs x 2 sexes)\n",
                file = file, append = TRUE)
            write.table(as_sex(st$gender_strat), file, append = TRUE, sep = ",",
                        row.names = FALSE, quote = TRUE)
          }

          if (!is.null(st$aab_donor_summary)) {
            cat("\n# Autoantibody Profile (Aab+ donors)\n", file = file, append = TRUE)
            write.table(st$aab_donor_summary, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # AAb-count test — donor-level, matching the on-screen subtitle. The
          # plot showed a p-value that was never exported, so the figure could
          # not be reproduced from the CSV.
          kw <- st$aab_kw
          if (!is.null(kw)) {
            cat("\n# AAb Count Test (Kruskal-Wallis on DONOR means; the box plot shows islets)\n",
                file = file, append = TRUE)
            write.table(
              data.frame(p_value = kw$p_value, n_donors = kw$n_donors,
                         n_islets = kw$n_islets,
                         donors_per_group = paste(sprintf("%s=%d", names(kw$donors_per_group),
                                                          as.integer(kw$donors_per_group)),
                                                  collapse = "; "),
                         underpowered = kw$underpowered,
                         stringsAsFactors = FALSE),
              file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }
        } else {
          writeLines(c(lines, "# No results computed"), file)
        }
      }
    )
  })
}
