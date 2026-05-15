# ---------- Statistics module server ----------
# Exports: statistics_server(id, raw_df, summary_df, get_selection_description)
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
#   cohens_d, eta_squared, pairwise_wilcox           -- utils_stats.R
#   bin_islet_sizes                                  -- data_loading.R

statistics_server <- function(id, raw_df, summary_df, get_selection_description,
                              donor_colors = reactive(DONOR_COLORS),
                              remove_outliers = reactive(TRUE)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Helper: section heading (mirrors the one in UI but available server-side)
    # `tip` (NULL or character): if supplied, renders a small "(?)" icon with a
    # native browser tooltip describing the section's contents.
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

    # Dynamic section numbers: shift down by 1 when binning is skipped
    section_offset <- reactive({
      if (isTRUE(input$stats_no_binning)) -1L else 0L
    })

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

    output$section_confounders_heading <- renderUI({
      step <- 4L + section_offset()
      section_heading(as.character(step), "Confounders & Deeper Analysis",
                      "Check whether demographics explain the effect, and compare integrated area under the curve across groups.",
                      tip = paste0("Demographic confounders (age, sex, AAb count) and a covariate-adjusted ",
                                   "donor-level linear model. Sex-stratified pairwise tests are BH-corrected ",
                                   "across the 6 contrasts (3 pairs x 2 sexes)."))
    })

    output$section_methods_heading <- renderUI({
      step <- 5L + section_offset()
      section_heading(as.character(step), "Methods Reference",
                      "Technical details about tests, corrections, and guidelines.",
                      tip = paste0("Full methods text. Updates dynamically with the current test type, ",
                                   "outlier rule, zero-exclusion state, and binning settings."))
    })

    # ===========================================================================
    # Core computation: eventReactive triggered by Run button
    # ===========================================================================

    stats_run <- eventReactive(input$run_tests, {
      rdf <- raw_df()
      if (is.null(rdf) || !nrow(rdf)) return(NULL)

      alpha_num <- as.numeric(input$alpha)
      is_parametric <- identical(input$test_type, "Parametric")

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
      outlier_summary_pre <- summarize_outlier_filter(rdf)
      if (stats_filter_on) {
        rdf <- rdf[!rdf$is_outlier | is.na(rdf$is_outlier), , drop = FALSE]
      }

      # ---- Zero-value exclusion (sensitivity for zero-inflated features) ----
      # Removes rows whose feature value == 0 before any donor-level aggregation,
      # mixed-effects modelling, demographics, and per-bin tests. AUC is read
      # from the Plot tab's summary_df and is NOT affected.
      exclude_zero <- isTRUE(input$stats_exclude_zero)
      zero_dropped <- 0L
      if (exclude_zero) {
        n_before <- nrow(rdf)
        rdf <- rdf[is.finite(rdf$value) & rdf$value != 0, , drop = FALSE]
        zero_dropped <- as.integer(n_before - nrow(rdf))
      }

      # ---- Min-cells filter (like spatial tab) ----
      min_cells <- input$stats_min_cells %||% 1
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
      if (!is.null(input$stats_diam_range) && length(input$stats_diam_range) == 2) {
        dmin <- as.numeric(input$stats_diam_range[1])
        dmax <- as.numeric(input$stats_diam_range[2])
        if (is.finite(dmin) && is.finite(dmax)) {
          rdf <- rdf[is.finite(rdf$islet_diam_um) &
                      rdf$islet_diam_um >= dmin &
                      rdf$islet_diam_um <= dmax, , drop = FALSE]
        }
      }

      if (nrow(rdf) == 0) return(NULL)

      rdf$donor_status <- factor(rdf$donor_status, levels = c("ND", "Aab+", "T1D"))
      n_islets_filtered <- nrow(rdf)

      # Detect feature name from raw_df attributes
      feature_name <- attr(raw_df(), "selection_used") %||% "Unknown"

      # ---- Donor-level aggregation (conservative sensitivity) ----
      ddf <- aggregate_to_donor(rdf, agg_fn = "mean")
      ddf$donor_status <- factor(ddf$donor_status, levels = c("ND", "Aab+", "T1D"))
      n_donors <- nrow(ddf)

      # ---- Normality test on donor-level data ----
      norm_result <- normality_tests(ddf)

      # ---- PRIMARY: Mixed-effects model (preserves islet variation, accounts for donor clustering) ----
      lmer_result <- tryCatch(lmer_test_donor(rdf), error = function(e) {
        list(p_value = NA_real_, icc = NA_real_, fit = NULL,
             error_msg = paste("Unexpected error:", e$message))
      })
      lmer_p <- if (!is.null(lmer_result)) lmer_result$p_value else NA_real_
      icc <- if (!is.null(lmer_result)) lmer_result$icc else NA_real_
      lmer_error <- lmer_result$error_msg  # NULL if no error

      # Pairwise from mixed model via emmeans
      lmer_pairwise <- NULL
      if (!is.null(lmer_result) && !is.null(lmer_result$fit)) {
        lmer_pairwise <- tryCatch({
          emm <- emmeans::emmeans(lmer_result$fit, pairwise ~ donor_status)
          ct <- as.data.frame(emm$contrasts)
          # ct has columns: contrast, estimate, SE, df, t.ratio, p.value
          ct
        }, error = function(e) NULL)
      }

      # ---- SENSITIVITY: Donor-level test (conservative — collapses size variation) ----
      donor_test <- if (is_parametric) "ANOVA" else "Kruskal-Wallis"
      donor_p <- NA_real_
      eta_sq <- NA_real_

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
        if (!is.null(kt) && n_donors > 1) eta_sq <- kt$statistic / (n_donors - 1)
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

      # ---- Per-bin analysis (conditional on no-binning checkbox) ----
      bin_anova <- NULL
      bin_kendall <- NULL
      bin_pairwise <- NULL

      if (!isTRUE(input$stats_no_binning)) {
        bw <- suppressWarnings(as.numeric(input$bin_width))
        if (!is.finite(bw) || bw <= 0) bw <- 50

        rdf_binned <- bin_islet_sizes(rdf, "islet_diam_um", bw)

        bin_anova <- per_bin_donor_anova(rdf_binned, "diam_bin", "donor_status", "value")
        bin_kendall <- per_bin_donor_kendall(rdf_binned, "diam_bin", "donor_status", "value")
        bin_pairwise <- per_bin_donor_pairwise(rdf_binned, "diam_bin", "donor_status", "value",
                                                parametric = is_parametric)

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
              tryCatch({
                pt <- pairwise.t.test(sub$value, sub$donor_status, p.adjust.method = "none")
                mat <- as.data.frame(as.table(pt$p.value))
                colnames(mat) <- c("group1", "group2", "p_value")
                mat[!is.na(mat$p_value), , drop = FALSE]
              }, error = function(e) NULL)
            } else {
              tryCatch({
                pw_raw <- pairwise.wilcox.test(sub$value, sub$donor_status,
                                               p.adjust.method = "none", exact = FALSE)
                mat <- as.data.frame(as.table(pw_raw$p.value))
                colnames(mat) <- c("group1", "group2", "p_value")
                mat <- mat[!is.na(mat$p_value), , drop = FALSE]
                mat$group1 <- as.character(mat$group1)
                mat$group2 <- as.character(mat$group2)
                mat
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

      if (length(aab_avail) > 0) {
        aab_sub <- rdf[rdf$donor_status == "Aab+", , drop = FALSE]
        if (nrow(aab_sub) > 0) {
          # Count AAbs per islet (same for all islets from same donor)
          aab_sub$n_aab <- rowSums(sapply(aab_avail, function(col) {
            as.numeric(as.logical(aab_sub[[col]]))
          }), na.rm = TRUE)
          aab_summary <- aab_sub

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

      # ---- AUC: trapezoidal integration of binned summary curve per donor group ----
      sdf <- summary_df()
      auc_by_group <- tryCatch({
        if (is.null(sdf) || nrow(sdf) == 0) NULL
        else {
          sdf %>%
            dplyr::group_by(donor_status) %>%
            dplyr::arrange(diam_mid) %>%
            dplyr::summarise(
              auc = if (n() < 2) 0 else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
              n_bins = n(),
              .groups = "drop"
            )
        }
      }, error = function(e) NULL)

      list(
        n_islets         = n_islets_filtered,
        n_donors         = n_donors,
        donor_df         = ddf,
        feature_name     = feature_name,
        # Primary test: mixed-effects
        lmer_p           = lmer_p,
        icc              = icc,
        lmer_pairwise    = lmer_pairwise,
        lmer_error       = lmer_error,
        # Sensitivity test: donor-level
        donor_test       = donor_test,
        donor_p          = donor_p,
        eta_sq           = eta_sq,
        alpha            = alpha_num,
        is_parametric    = is_parametric,
        pairwise_df      = pairwise_df,
        bin_anova        = bin_anova,
        bin_kendall      = bin_kendall,
        bin_pairwise     = bin_pairwise,
        auc_by_group     = auc_by_group,
        demo_summary     = demo_summary,
        age_corr         = age_corr,
        age_model_p      = age_model_p,
        gender_strat     = gender_strat,
        norm_result      = norm_result,
        aab_summary      = aab_summary,
        aab_donor_summary = aab_donor_summary,
        rdf              = rdf,
        outlier_summary  = outlier_summary_pre,
        exclude_zero     = exclude_zero,
        zero_dropped     = zero_dropped
      )
    }, ignoreInit = TRUE)

    # ===========================================================================
    # Card 1: Overview Banner
    # ===========================================================================

    output$overview_banner <- renderUI({
      st <- stats_run()
      if (is.null(st)) {
        return(tags$p(style = "color: #666; margin: 0;",
                      "Select a feature in the sidebar, then click 'Run Statistics'."))
      }

      # Primary p-value from mixed-effects model
      primary_p <- st$lmer_p
      p_col <- if (is.na(primary_p)) "#999" else if (primary_p < st$alpha) "#d62728" else "#2e7d32"
      icc_label <- if (is.na(st$icc)) "" else sprintf("  ICC = %.2f", st$icc)

      zero_badge <- if (isTRUE(st$exclude_zero)) {
        tags$span(style = paste0("background: #fff3e0; color: #6a3a00; padding: 4px 10px; ",
                                 "border: 1px solid #ffcc80; border-radius: 4px;"),
                  sprintf("Zero-valued excluded (%s islet%s dropped)",
                          formatC(st$zero_dropped, format = "d", big.mark = ","),
                          if (st$zero_dropped == 1) "" else "s"))
      } else NULL

      tags$div(style = "display: flex; gap: 15px; align-items: center; flex-wrap: wrap;",
        tags$span(style = "background: #e3f2fd; padding: 4px 10px; border-radius: 4px; font-weight: 600;",
                  st$feature_name),
        tags$span(style = "padding: 4px 10px;",
                  paste0(formatC(st$n_islets, format = "d", big.mark = ","),
                         " islets, ", st$n_donors, " donors", icc_label)),
        tags$span(style = paste0("background: ", p_col, "; color: white; padding: 4px 10px; border-radius: 4px;"),
                  paste0("Mixed-effects p = ",
                         if (is.na(primary_p)) "N/A" else formatC(primary_p, format = "e", digits = 2))),
        zero_badge
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

      all_normal <- all(nr$is_normal, na.rm = TRUE)
      any_non_normal <- any(!nr$is_normal, na.rm = TRUE)
      suggestion <- if (any_non_normal) {
        tags$span(style = "color: #e65100; font-weight: 600;",
                  "Consider non-parametric tests.")
      } else if (all_normal) {
        tags$span(style = "color: #2e7d32;",
                  "Parametric tests may be appropriate.")
      } else {
        NULL
      }

      tags$div(
        style = "background: #f5f5f5; padding: 8px 12px; border-radius: 5px; font-size: 12px; line-height: 1.8;",
        tags$strong("Shapiro-Wilk normality (donor means):"),
        tags$br(),
        tagList(items),
        if (!is.null(suggestion)) tagList(tags$br(), suggestion),
        tags$br(),
        tags$span(style = "color: #999; font-size: 11px;",
                  "Note: with n=5/group, Shapiro-Wilk has limited power to detect departures from normality.")
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
        paste0("p = ", formatC(st$donor_p, format = "e", digits = 2))

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
          style = "background: #fff3cd; padding: 8px 12px; border-radius: 4px; margin-bottom: 8px; font-size: 12px; color: #856404;",
          tags$strong("Mixed-effects model issue: "), st$lmer_error,
          " — donor-level results below remain valid."
        )
      }

      tagList(
        lmer_warn,
        tags$div(
          style = "background: #e8f4fd; padding: 10px 14px; border-radius: 5px; margin-bottom: 12px; font-size: 13px; line-height: 1.5;",
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

      fmt_p <- function(x) {
        ifelse(is.na(x), "N/A",
               ifelse(x < 1e-5, formatC(x, format = "e", digits = 2),
                      sprintf("%.4f", x)))
      }
      fmt_d <- function(x) ifelse(is.na(x), "N/A", sprintf("%.2f", x))

      rows <- data.frame(
        Contrast = character(), `p-value` = character(), `p (adj)` = character(),
        `Cohen's d` = character(), `95% CI` = character(), Sig = character(),
        check.names = FALSE, stringsAsFactors = FALSE
      )

      # PRIMARY: Mixed-effects global (islet-level, donor as random intercept)
      rows <- rbind(rows, data.frame(
        Contrast = "Overall (mixed-effects, donor random intercept)",
        `p-value` = fmt_p(st$lmer_p),
        `p (adj)` = "-",
        `Cohen's d` = "-",
        `95% CI` = "-",
        Sig = if (!is.na(st$lmer_p) && st$lmer_p < st$alpha) "*" else "",
        check.names = FALSE, stringsAsFactors = FALSE
      ))

      # Pairwise from mixed model (estimated marginal means)
      lmer_pw <- st$lmer_pairwise
      if (!is.null(lmer_pw) && nrow(lmer_pw) > 0) {
        lmer_pw_rows <- data.frame(
          Contrast = paste0(lmer_pw$contrast, " (mixed-effects pairwise)"),
          `p-value` = fmt_p(lmer_pw$p.value),
          `p (adj)` = "-",
          `Cohen's d` = "-",
          `95% CI` = "-",
          Sig = ifelse(!is.na(lmer_pw$p.value) & lmer_pw$p.value < st$alpha, "*", ""),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        rows <- rbind(rows, lmer_pw_rows)
      }

      # SENSITIVITY: Donor-level means (collapses islets to one value per donor)
      rows <- rbind(rows, data.frame(
        Contrast = paste0("Donor-level means (", st$donor_test, ", N=", st$n_donors, " donors)"),
        `p-value` = fmt_p(st$donor_p),
        `p (adj)` = "-",
        `Cohen's d` = "-",
        `95% CI` = "-",
        Sig = if (!is.na(st$donor_p) && st$donor_p < st$alpha) "*" else "",
        check.names = FALSE, stringsAsFactors = FALSE
      ))

      # Donor-level pairwise with effect sizes
      pw <- st$pairwise_df
      if (!is.null(pw) && nrow(pw) > 0) {
        test_label <- if (st$is_parametric) "t-test" else "Wilcoxon"
        pw_rows <- data.frame(
          Contrast = paste0(pw$contrast, " (donor-level ", test_label, ")"),
          `p-value` = fmt_p(pw$p_value),
          `p (adj)` = fmt_p(pw$p_adj),
          `Cohen's d` = fmt_d(pw$cohens_d),
          `95% CI` = ifelse(is.na(pw$d_ci_lo), "N/A",
                            paste0("[", sprintf("%.2f", pw$d_ci_lo), ", ", sprintf("%.2f", pw$d_ci_hi), "]")),
          Sig = ifelse(!is.na(pw$p_adj) & pw$p_adj < st$alpha, "*", ""),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        rows <- rbind(rows, pw_rows)
      }

      rows
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lccccc")

    output$forest_plot <- renderPlotly({
      st <- stats_run()
      if (is.null(st)) return(NULL)
      pw <- st$pairwise_df
      if (is.null(pw) || nrow(pw) == 0 || all(is.na(pw$cohens_d))) return(NULL)

      pw <- pw[!is.na(pw$cohens_d), , drop = FALSE]
      pw$contrast <- factor(pw$contrast, levels = rev(pw$contrast))
      pw$sig_label <- ifelse(pw$sig, "Significant", "NS")

      zero_suffix <- if (isTRUE(st$exclude_zero) && st$zero_dropped > 0) {
        sprintf(" — %s zero-valued dropped",
                formatC(st$zero_dropped, format = "d", big.mark = ","))
      } else ""
      forest_title <- paste0("Pairwise Effect Sizes (conservative)",
                             outlier_title_suffix(st$outlier_summary),
                             zero_suffix)
      g <- ggplot(pw, aes(x = cohens_d, y = contrast, color = sig_label)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "#999") +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = d_ci_lo, xmax = d_ci_hi), height = 0.2) +
        scale_color_manual(values = c("Significant" = "#d62728", "NS" = "#1f77b4"),
                           name = "") +
        labs(x = "Cohen's d (donor means)", y = NULL, title = forest_title) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(margin = list(b = 60),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))
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

      mid_labels <- as.character(round(all_mids))

      plot_ly(
        x = mid_labels, y = row_names, z = z_mat,
        type = "heatmap",
        colorscale = list(c(0, "#f7fbff"), c(0.5, "#6baed6"), c(1, "#08306b")),
        colorbar = list(title = "-log10(q)"),
        text = text_mat, texttemplate = "%{text}",
        hovertemplate = "Bin: %{x} \u00b5m<br>Test: %{y}<br>-log10(q): %{z:.2f}<extra></extra>"
      ) %>%
        layout(
          xaxis = list(title = "Diameter bin midpoint (\u00b5m)", tickangle = -45),
          yaxis = list(title = "", categoryorder = "array", categoryarray = rev(row_names)),
          title = "Stratified Significance (BH-corrected, donor-level)"
        )
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
        theme_minimal(base_size = 13) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(margin = list(b = 100),
               legend = list(orientation = "h", x = 0.5, xanchor = "center",
                             y = -0.3, font = list(size = 12)),
               showlegend = TRUE)
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
        return(div(class = "card", style = "padding: 15px; margin-bottom: 15px;",
                   h5("Demographics Analysis"),
                   tags$p(style = "color: #666;", "Run statistics to see demographic analysis.")))
      }
      if (is.null(st$demo_summary)) return(NULL)

      # Build AAb section conditionally
      aab_section <- NULL
      if (!is.null(st$aab_summary) && nrow(st$aab_summary) > 0) {
        aab_section <- div(style = "margin-top: 18px;",
          h6("Autoantibody Profile (Aab+ only)", style = "color: #555; margin-bottom: 10px;"),
          fluidRow(
            column(6,
              tableOutput(ns("aab_donor_table"))
            ),
            column(6,
              plotlyOutput(ns("aab_plot"), height = "260px")
            )
          )
        )
      }

      div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible;",
        h5("Demographics Analysis"),
        div(style = "margin-bottom: 18px;",
          h6("Donor Summary", style = "color: #555; margin-bottom: 10px;"),
          tableOutput(ns("demo_summary_table"))
        ),
        div(style = "margin-bottom: 18px;",
          fluidRow(
            column(6,
              h6("Age vs Feature (islet-level)", style = "color: #555; margin-bottom: 10px;"),
              plotlyOutput(ns("age_scatter"), height = "300px")
            ),
            column(6,
              h6("Sex vs Feature", style = "color: #555; margin-bottom: 10px;"),
              plotlyOutput(ns("gender_plot"), height = "300px")
            )
          )
        ),
        aab_section,
        div(style = "margin-top: 18px;",
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
                            formatC(st$age_corr$p_value, format = "e", digits = 2),
                            st$age_corr$n)

      g <- ggplot(rdf_plot, aes(x = age_num, y = value, color = donor_status)) +
        geom_point(alpha = 0.4, size = 1.5) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, aes(group = 1), color = "#333") +
        scale_color_manual(values = donor_colors()) +
        labs(x = "Donor Age (years)", y = st$feature_name, color = NULL,
             title = "Age vs Feature (all islets)",
             subtitle = corr_label) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")

      ggplotly(g, tooltip = c("x", "y", "colour")) %>%
        layout(legend = list(orientation = "h", x = 0.2, y = -0.2))
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
                                            formatC(min_p, format = "e", digits = 1), suffix))
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
        theme_minimal(base_size = 12) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("y")) %>%
        layout(margin = list(b = 50))
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

      # Kruskal-Wallis test for trend
      kt <- tryCatch(kruskal.test(value ~ aab_group, data = aab_df), error = function(e) NULL)
      p_label <- if (!is.null(kt)) sprintf("KW p = %s", formatC(kt$p.value, format = "e", digits = 1)) else ""

      g <- ggplot(aab_df, aes(x = aab_group, y = value, fill = aab_group)) +
        geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
        scale_fill_brewer(palette = "Set2") +
        labs(x = "Number of Autoantibodies", y = st$feature_name,
             title = "Feature by AAb Count (Aab+ islets)",
             subtitle = paste0("N = ", nrow(aab_df), " islets from ",
                               length(unique(aab_df$`Case ID`)), " Aab+ donors. ", p_label)) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("y"))
    })

    output$covariate_results <- renderUI({
      st <- stats_run()
      if (is.null(st) || is.null(st$age_model_p)) return(NULL)
      p_cov <- st$age_model_p
      survived <- !is.na(p_cov) && p_cov < st$alpha
      tags$div(style = "margin-top: 10px; padding: 10px; background: #f8f9fa; border-radius: 5px;",
        tags$strong("Covariate-Adjusted Model (donor-level)"),
        tags$p(style = "margin: 5px 0 0 0;",
          paste0("lm(value ~ donor_status + age + sex, data = donor_means): donor_status p = ",
                 if (is.na(p_cov)) "N/A" else formatC(p_cov, format = "e", digits = 2)),
          tags$br(),
          if (survived) {
            tags$span(style = "color: #d62728; font-weight: 600;",
                      "Donor status effect SURVIVES adjustment for age and sex.")
          } else {
            tags$span(style = "color: #666;",
                      "Donor status effect does NOT survive adjustment for age and sex.")
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
        labs(x = "Donor Status", y = "Area Under Curve (AUC)",
             title = "Integrated Area by Donor Group") +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")

      ggplotly(g, tooltip = c("x", "y"))
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
        tags$p(style = "margin-top: 10px; font-size: 90%; color: #555;",
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
        return(tags$div(style = "padding: 10px; color: #666;",
          tags$p("Click 'Run Statistics' to perform analysis. ",
                 "The sidebar controls which feature and subset of data are analyzed.")))
      }

      alpha_str <- as.character(st$alpha)
      test_desc <- if (st$is_parametric) {
        "Parametric tests: one-way ANOVA (global), pairwise t-tests (BH-corrected)."
      } else {
        "Non-parametric tests: Kruskal-Wallis (global), pairwise Wilcoxon rank-sum (BH-corrected)."
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
            "Spatial Card A uses Tukey 1.5×IQR per donor status to tolerate zero-inflated phenotypes.",
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
          "Useful for zero-inflated features (rare phenotypes or peri-islet enrichment z-scores where many islets register 0). ",
          tags$strong("Note: "), "AUC is computed from the Plot-tab summary curves and is not affected by this filter.",
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

      tags$div(style = "padding: 10px;",
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
                 "Pairwise group comparisons use estimated marginal means from the mixed model."),
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
          tags$strong("Significance threshold: "), paste0("\u03b1 = ", alpha_str),
          tags$br(),
          tags$strong("Effect sizes: "), "Cohen's d with 95% CI computed on donor-level means. ",
          "These are conservative estimates.",
          tags$br(),
          tags$strong("Normality: "), "Shapiro-Wilk test on donor-level means per group. ",
          "With n = 5, this test has limited power, so failure to reject normality does not guarantee it.",
          tags$br(),
          tags$strong("Stratified analysis: "), "ANOVA, pairwise tests, and Kendall \u03c4 on donor-level means within each diameter bin. ",
          "Kendall \u03c4 correlates disease stage (ND=0, Aab+=1, T1D=2) with the feature value: ",
          "\u03c4 > 0 means the feature increases with disease progression, \u03c4 < 0 means it decreases. ",
          "P-values are BH-corrected across bins. A bin is tested only if at least 2 donor-status groups ",
          "each have \u2265 2 donors with finite values (same rule for the heatmap and the trend plot).",
          tags$br(),
          tags$strong("Multiple comparisons: "), "Benjamini-Hochberg (FDR) correction applied to pairwise p-values.",
          tags$br(),
          tags$strong("AUC: "), "Area under curve computed by trapezoidal integration of binned mean values across islet diameters, per donor group."
        ),
        if (!is.null(st$demo_summary)) {
          tags$p(style = "margin-top: 8px;",
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
              " Autoantibody analysis shows feature distribution by AAb count (number of positive AAbs) for Aab+ donors only."
            } else {
              ""
            }
          )
        },
        tags$p(style = "font-size: 90%; color: #666; margin-top: 8px;",
          tags$strong("Interpretation: "),
          "A significant result (p < \u03b1) suggests the feature differs between donor groups. ",
          "With N = 5 per group, power is limited; non-significant results should not be over-interpreted. ",
          "Effect sizes (Cohen's d) provide magnitude information independent of sample size."
        )
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
        st <- stats_run()
        descriptions <- get_selection_description()

        lines <- c(paste("#", descriptions))
        lines <- c(lines, "")

        if (!is.null(st)) {
          # Methodology notes
          lines <- c(lines, "# METHODOLOGY")
          lines <- c(lines, paste0("# Primary: mixed-effects lmer (islet-level with donor random intercept)"))
          os <- st$outlier_summary
          if (!is.null(os)) {
            lines <- c(lines, paste0("# Outlier rule: ",
                                     if (isTRUE(os$enabled)) os$label else "OFF (filter disabled)",
                                     "; excluded = ", os$n_excluded,
                                     " of ", os$n_total, " input rows."))
          }
          lines <- c(lines, paste0("# Zero-value filter: ",
                                   if (isTRUE(st$exclude_zero)) "ON" else "OFF",
                                   "; N dropped = ", st$zero_dropped))
          lines <- c(lines, paste0("# Total islets: ",
                                   formatC(st$n_islets, big.mark = ","),
                                   ", Donors: ", st$n_donors))
          lines <- c(lines, paste0("# ICC (intra-class correlation): ",
                                   if (is.na(st$icc)) "N/A" else sprintf("%.4f", st$icc)))
          lines <- c(lines, paste0("# Mixed-effects p-value (primary): ",
                                   if (is.na(st$lmer_p)) "N/A" else formatC(st$lmer_p, format = "e", digits = 4)))
          lines <- c(lines, "")

          # Sensitivity test
          lines <- c(lines, "# Sensitivity: Donor-Level Test (conservative)")
          lines <- c(lines, paste0("# ", st$donor_test, ", p = ",
                                   formatC(st$donor_p, format = "e", digits = 4),
                                   ", eta_sq = ", round(st$eta_sq, 4)))
          lines <- c(lines, "")

          # Write header lines
          writeLines(lines, file)

          # Normality tests
          if (!is.null(st$norm_result)) {
            cat("# Shapiro-Wilk Normality Tests (donor-level)\n", file = file, append = TRUE)
            write.table(st$norm_result, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Donor-level data
          if (!is.null(st$donor_df)) {
            cat("\n# Donor-Level Means\n", file = file, append = TRUE)
            write.table(st$donor_df, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          # Pairwise
          cat("\n# Pairwise Comparisons (donor-level)\n", file = file, append = TRUE)
          pw <- st$pairwise_df
          if (!is.null(pw) && nrow(pw) > 0) {
            write.table(pw, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
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
            cat("\n# Sex-Stratified Pairwise (BH-corrected across 6 contrasts: 3 pairs x 2 sexes)\n",
                file = file, append = TRUE)
            write.table(st$gender_strat, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }

          if (!is.null(st$aab_donor_summary)) {
            cat("\n# Autoantibody Profile (Aab+ donors)\n", file = file, append = TRUE)
            write.table(st$aab_donor_summary, file, append = TRUE, sep = ",", row.names = FALSE, quote = TRUE)
          }
        } else {
          writeLines(c(lines, "# No results computed"), file)
        }
      }
    )
  })
}
