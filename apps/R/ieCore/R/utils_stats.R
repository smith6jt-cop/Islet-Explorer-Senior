# ---------- Unified outlier flagging (used by Plot, Trajectory, Statistics, Spatial) ----------
#
# Single source of truth for outlier identification across the app. Each tab
# previously had its own inline loop with slightly different rules (z-score on
# raw values vs Tukey IQR vs 99th-percentile clip). The unified rule is:
#   For each donor_status group (ND, Aab+, T1D), exclude rows whose value is
#   more than `threshold` SDs from that group's mean. Outliers are flagged on
#   the returned data.frame; downstream code decides whether to render them as
#   ghost points or drop them. The IQR fallback is kept for cards plotting an
#   already-standardised, roughly symmetric quantity (Spatial Card A on
#   enrichment z-scores).
#
# `exclude_zeros = TRUE` removes zero-valued rows from the scoring population
# entirely \u2014 they neither contribute to the centre/scale nor receive a score
# (their outlier_score stays NA). Use it for genuinely zero-inflated features
# where a spike at 0 would collapse the MAD; do NOT use it on quantities that
# merely take negative values (e.g. enrichment z-scores contain no exact zeros,
# so the flag is a no-op there).
#
# Returns the input df with two appended columns:
#   is_outlier    -- logical, TRUE for rows beyond the threshold
#   outlier_score -- numeric, z-score (method="zscore") or signed distance to
#                    the nearer IQR fence in IQR units (method="iqr"); NA for
#                    rows we couldn't score (single-row groups, etc.)
flag_outliers <- function(df,
                          value_col = "value",
                          group_col = "donor_status",
                          method = c("zscore", "iqr"),
                          threshold = NULL,
                          transform = c("identity", "log1p"),
                          exclude_zeros = FALSE,
                          enabled = TRUE) {
  method <- match.arg(method)
  transform <- match.arg(transform)
  if (is.null(threshold)) {
    threshold <- if (method == "zscore") 3 else 1.5
  }

  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    if (is.data.frame(df)) {
      df$is_outlier <- logical(0)
      df$outlier_score <- numeric(0)
    }
    attr(df, "outlier_method") <- method
    attr(df, "outlier_threshold") <- threshold
    attr(df, "outlier_enabled") <- isTRUE(enabled)
    return(df)
  }

  df$is_outlier <- FALSE
  df$outlier_score <- NA_real_

  if (!isTRUE(enabled)) {
    attr(df, "outlier_method") <- method
    attr(df, "outlier_threshold") <- threshold
    attr(df, "outlier_enabled") <- FALSE
    return(df)
  }

  if (!value_col %in% colnames(df)) {
    attr(df, "outlier_method") <- method
    attr(df, "outlier_threshold") <- threshold
    attr(df, "outlier_enabled") <- TRUE
    return(df)
  }

  vals_raw <- suppressWarnings(as.numeric(df[[value_col]]))
  vals <- if (transform == "log1p") {
    ifelse(is.finite(vals_raw) & vals_raw >= 0, log1p(vals_raw), NA_real_)
  } else {
    vals_raw
  }

  if (group_col %in% colnames(df)) {
    groups <- as.character(df[[group_col]])
  } else {
    groups <- rep("__all__", nrow(df))
  }

  # Track which scaling rule(s) actually fired across groups. `stat_seen`
  # accumulates the names; the attribute set below names "modified_z",
  # "classical_z", or "modified_z (some groups: classical_z fallback)".
  stat_seen <- character(0)

  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(vals))
    if (length(idx) < 3) next
    grp_vals <- vals[idx]

    # Helper: compute a robust signed z-like score for method="zscore". Uses
    # Iglewicz-Hoaglin modified z = 0.6745 * (x - median) / MAD, which is
    # robust to right-skew (mean/SD pulled by long tails leaves the lower
    # bound mathematically unreachable on non-negative data; median/MAD does
    # not). Falls back to classical mean/SD when MAD = 0 (zero-inflated /
    # constant groups). 0.6745 rescales MAD to a normal-distribution-equivalent
    # SD so the slider threshold is comparable to a classical z.
    robust_zscore <- function(ref_vals, score_vals) {
      med <- stats::median(ref_vals, na.rm = TRUE)
      mad <- stats::median(abs(ref_vals - med), na.rm = TRUE)
      if (is.finite(mad) && mad > 0) {
        list(z = 0.6745 * (score_vals - med) / mad, statistic = "modified_z")
      } else {
        m <- mean(ref_vals, na.rm = TRUE)
        s <- sd(ref_vals, na.rm = TRUE)
        if (!is.finite(s) || s == 0) return(NULL)
        list(z = (score_vals - m) / s, statistic = "classical_z")
      }
    }

    if (exclude_zeros) {
      # Zeros are removed from the SCORING POPULATION entirely: they neither
      # contribute to the centre/scale nor receive a score. Scoring a zero
      # against a distribution it was excluded from produces a large spurious
      # negative (e.g. 0.6745 * (0 - median_nz) / MAD_nz) that then surfaces in
      # the outlier tables and tooltips as if it meant something. NA is the
      # honest answer: "not evaluated under this rule".
      is_zero <- vals[idx] == 0
      nz_idx <- idx[!is_zero]
      nz_vals <- vals[nz_idx]
      # Same minimum-group rule as the default path below (>= 3 scorable values).
      if (length(nz_vals) < 3) next
      if (method == "zscore") {
        rz <- robust_zscore(nz_vals, vals[nz_idx])
        if (is.null(rz)) next
        z <- rz$z
        stat_seen <- unique(c(stat_seen, rz$statistic))
        df$outlier_score[nz_idx] <- z
        df$is_outlier[nz_idx] <- !is.na(z) & abs(z) > threshold
      } else {
        q <- stats::quantile(nz_vals, c(0.25, 0.75), na.rm = TRUE)
        iqr <- q[2] - q[1]
        if (!is.finite(iqr) || iqr == 0) next
        lo <- q[1] - threshold * iqr
        hi <- q[2] + threshold * iqr
        fence_dist <- ifelse(nz_vals < q[1],
                             (nz_vals - q[1]) / iqr,
                             ifelse(nz_vals > q[2], (nz_vals - q[2]) / iqr, 0))
        df$outlier_score[nz_idx] <- fence_dist
        df$is_outlier[nz_idx] <- !is.na(nz_vals) & (nz_vals < lo | nz_vals > hi)
      }
      # Zero rows keep the initialised is_outlier = FALSE / outlier_score = NA.
    } else {
      if (method == "zscore") {
        rz <- robust_zscore(grp_vals, grp_vals)
        if (is.null(rz)) next
        z <- rz$z
        stat_seen <- unique(c(stat_seen, rz$statistic))
        df$outlier_score[idx] <- z
        df$is_outlier[idx] <- !is.na(z) & abs(z) > threshold
      } else {
        q <- stats::quantile(grp_vals, c(0.25, 0.75), na.rm = TRUE)
        iqr <- q[2] - q[1]
        if (!is.finite(iqr) || iqr == 0) next
        lo <- q[1] - threshold * iqr
        hi <- q[2] + threshold * iqr
        fence_dist <- ifelse(grp_vals < q[1],
                             (grp_vals - q[1]) / iqr,
                             ifelse(grp_vals > q[2], (grp_vals - q[2]) / iqr, 0))
        df$outlier_score[idx] <- fence_dist
        df$is_outlier[idx] <- !is.na(grp_vals) & (grp_vals < lo | grp_vals > hi)
      }
    }
  }

  attr(df, "outlier_method") <- method
  attr(df, "outlier_threshold") <- threshold
  attr(df, "outlier_enabled") <- TRUE
  attr(df, "outlier_statistic") <- if (method != "zscore") {
    "iqr"
  } else if (length(stat_seen) == 0) {
    "modified_z"
  } else if (length(stat_seen) == 1) {
    stat_seen
  } else {
    "modified_z (classical_z fallback used for some groups)"
  }
  df
}

# Build a one-line summary of an outlier-flagged data frame for badges/tooltips.
summarize_outlier_filter <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !("is_outlier" %in% colnames(df))) {
    return(list(n_total = 0L, n_excluded = 0L, n_kept = 0L,
                method = NA_character_, threshold = NA_real_,
                enabled = FALSE, label = ""))
  }
  n_total <- nrow(df)
  n_excluded <- sum(isTRUE(attr(df, "outlier_enabled")) & df$is_outlier, na.rm = TRUE)
  if (!isTRUE(attr(df, "outlier_enabled"))) n_excluded <- 0L
  method <- attr(df, "outlier_method") %||% "zscore"
  threshold <- attr(df, "outlier_threshold") %||% NA_real_
  enabled <- isTRUE(attr(df, "outlier_enabled"))
  label <- if (!enabled) {
    "outlier filter OFF"
  } else if (method == "zscore") {
    sprintf("|z| > %g per donor status", threshold)
  } else {
    sprintf("Tukey %g\u00d7IQR per donor status", threshold)
  }
  list(n_total = n_total, n_excluded = as.integer(n_excluded),
       n_kept = n_total - as.integer(n_excluded),
       method = method, threshold = threshold,
       enabled = enabled, label = label)
}

# Format a one-line suffix for plot titles. Empty string when no exclusion.
outlier_title_suffix <- function(summary) {
  if (is.null(summary) || !isTRUE(summary$enabled) || summary$n_excluded == 0) {
    return("")
  }
  sprintf(" \u2014 %d outlier%s excluded (%s)",
          summary$n_excluded,
          if (summary$n_excluded == 1) "" else "s",
          summary$label)
}

# Convenience: shared tooltip text for the global outlier control.
GLOBAL_OUTLIER_TOOLTIP <- paste0(
  "For each donor status (ND / Aab+ / T1D), islets whose robust |z| from the ",
  "group's median and MAD (modified z-score) exceeds the slider threshold are ",
  "excluded from trend lines, donor-level means, and statistical tests. ",
  "Detection is two-sided \u2014 both high and low extremes are flagged. When MAD = 0 ",
  "(zero-inflated groups), the rule falls back to classical mean/SD. Excluded ",
  "points still appear on scatter and box plots as faint grey \u00d7 markers; counts ",
  "are shown next to each plot title. Spatial Card A (infiltration) uses Tukey ",
  "1.5\u00d7IQR instead, since enrichment z-scores are already a standardised, ",
  "roughly symmetric quantity; Spatial Card C distance plots apply the z-score ",
  "rule on log(1 + distance)."
)

# ---------- Plotly hover-text scrubber ----------
#
# ggplot geoms with a hard-coded `color = "black"` parameter (e.g. shape=21
# summary points, "overall" LOESS lines) cause ggplotly to surface the literal
# string "colour: black" in the hover tooltip. The same conversion also emits
# both the `fill` and inherited `colour` aesthetics for shape-21 points, so the
# donor_status label can appear twice. This helper post-processes a ggplotly
# object's traces and removes those artefacts. Safe to apply unconditionally.
strip_plotly_color_leak <- function(gg) {
  if (is.null(gg) || is.null(gg$x) || is.null(gg$x$data)) return(gg)

  scrub <- function(s) {
    if (is.null(s) || length(s) == 0) return(s)
    s <- as.character(s)
    # Drop any "colour: black" or "color: black" line (case-insensitive),
    # collapsing the surrounding <br> separator so the tooltip stays compact.
    s <- gsub("(?i)<br\\s*/?\\s*>\\s*colou?r:\\s*black\\s*(?=<br|$)", "", s, perl = TRUE)
    s <- gsub("(?i)^\\s*colou?r:\\s*black\\s*<br\\s*/?\\s*>", "", s, perl = TRUE)
    s <- gsub("(?i)^\\s*colou?r:\\s*black\\s*$", "", s, perl = TRUE)
    # Collapse adjacent duplicate donor_status lines (handles the ggplotly
    # case where both `fill` and `colour` aesthetics emit the same value).
    s <- gsub("(donor_status:\\s*([^<]+))(<br\\s*/?\\s*>donor_status:\\s*\\2)+",
              "\\1", s, perl = TRUE)
    s
  }

  for (i in seq_along(gg$x$data)) {
    tr <- gg$x$data[[i]]
    if (!is.null(tr$text))          gg$x$data[[i]]$text          <- scrub(tr$text)
    if (!is.null(tr$hovertemplate)) gg$x$data[[i]]$hovertemplate <- scrub(tr$hovertemplate)
    if (!is.null(tr$hovertext))     gg$x$data[[i]]$hovertext     <- scrub(tr$hovertext)
  }
  gg
}

summary_stats <- function(df, group_cols, value_col, stat = c("mean_se","mean_sd","median_iqr")) {
  stat <- match.arg(stat)
  df %>% group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      se = ifelse(n > 0, sd / sqrt(pmax(n, 1)), NA_real_),
      median = median(.data[[value_col]], na.rm = TRUE),
      q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE, type = 7),
      q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE, type = 7),
      .groups = "drop"
    ) %>%
    mutate(
      y = dplyr::case_when(
        stat == "mean_se" ~ mean,
        stat == "mean_sd" ~ mean,
        stat == "median_iqr" ~ median,
        TRUE ~ mean
      ),
      ymin = dplyr::case_when(
        stat == "mean_se" ~ mean - se,
        stat == "mean_sd" ~ mean - sd,
        stat == "median_iqr" ~ q1,
        TRUE ~ mean - se
      ),
      ymax = dplyr::case_when(
        stat == "mean_se" ~ mean + se,
        stat == "mean_sd" ~ mean + sd,
        stat == "median_iqr" ~ q3,
        TRUE ~ mean + se
      )
    )
}

# Map configured group levels to an ordinal trend code used by Kendall tau.
# Values outside the declared vocabulary stay unavailable (NA), never zero.
code_group_status <- function(x, levels = NULL) {
  if (is.null(levels)) {
    levels <- if (exists("domain_group_levels", mode = "function")) {
      domain_group_levels()
    } else {
      c("ND", "Aab+", "T1D")
    }
  }
  code <- match(as.character(x), as.character(levels)) - 1L
  code[is.na(match(as.character(x), as.character(levels)))] <- NA_integer_
  as.numeric(code)
}

# Compatibility name retained for the existing module API.
code_donor_status <- function(x) code_group_status(x)

per_bin_anova <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    if (n_distinct(na.omit(sub[[group_col]])) < 2 || sum(is.finite(sub[[value_col]])) < 2) return(NULL)
    fit <- tryCatch(aov(reformulate(group_col, response = value_col), data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pval <- tryCatch({
      at <- anova(fit)
      as.numeric(at[["Pr(>F)"]][1])
    }, error = function(e) NA_real_)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_anova = pval)
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  arrange(out, mid)
}

per_bin_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    x <- code_donor_status(sub[[group_col]])
    y <- suppressWarnings(as.numeric(sub[[value_col]]))
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    if (length(unique(x)) < 2 || length(y) < 3) return(NULL)
    ct <- tryCatch(cor.test(x, y, method = "kendall", exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_kendall = unname(ct$p.value), tau = unname(ct$estimate))
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(out)
  arrange(out, mid)
}

# ---------- Donor-level aggregation (pseudoreplication fix) ----------

#' Aggregate islet-level data to donor-level means
#' Returns one row per donor with mean value, n_islets, and preserved metadata
aggregate_to_donor <- function(rdf, agg_fn = c("mean", "median")) {
  agg_fn <- match.arg(agg_fn)
  fn <- if (agg_fn == "mean") mean else median
  if (!"Case ID" %in% colnames(rdf)) return(rdf)
  has_age <- "age" %in% colnames(rdf)
  has_gender <- "gender" %in% colnames(rdf)
  out <- rdf %>%
    dplyr::group_by(`Case ID`, donor_status) %>%
    dplyr::summarise(
      value = fn(value, na.rm = TRUE),
      n_islets = dplyr::n(),
      .groups = "drop"
    )
  if (has_age) {
    age_lookup <- rdf %>%
      dplyr::group_by(`Case ID`) %>%
      dplyr::summarise(age = dplyr::first(na.omit(age)), .groups = "drop")
    out <- dplyr::left_join(out, age_lookup, by = "Case ID")
  }
  if (has_gender) {
    gender_lookup <- rdf %>%
      dplyr::group_by(`Case ID`) %>%
      dplyr::summarise(gender = dplyr::first(na.omit(as.character(gender))), .groups = "drop")
    out <- dplyr::left_join(out, gender_lookup, by = "Case ID")
  }
  out
}

#' Mixed-effects test: value ~ donor_status + (1 | Case ID)
#' Returns list(p_value, icc, fit) or NULL on failure
lmer_test_donor <- function(rdf) {
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    return(list(p_value = NA_real_, icc = NA_real_, fit = NULL,
                error_msg = "lmerTest package not installed"))
  }
  if (!"Case ID" %in% colnames(rdf)) {
    return(list(p_value = NA_real_, icc = NA_real_, fit = NULL,
                error_msg = "No 'Case ID' column in data"))
  }
  rdf$case_id_f <- factor(rdf$`Case ID`)

  # Fit mixed-effects model
  fit <- tryCatch(
    lmerTest::lmer(value ~ donor_status + (1 | case_id_f), data = rdf),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(p_value = NA_real_, icc = NA_real_, fit = NULL,
                error_msg = paste("Model failed to fit:", fit$message)))
  }

  # Extract p-value via lmerTest's anova (Satterthwaite F-test)
  # With a single fixed effect, Type I = Type III, so no need for car::Anova
  p_val <- tryCatch({
    at <- anova(fit)  # lmerTest anova method
    if ("donor_status" %in% rownames(at)) {
      as.numeric(at["donor_status", "Pr(>F)"])
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  # ICC: donor variance / total variance
  vc <- tryCatch(as.data.frame(lme4::VarCorr(fit)), error = function(e) NULL)
  icc <- NA_real_
  if (!is.null(vc)) {
    donor_var <- vc$vcov[vc$grp == "case_id_f"]
    resid_var <- vc$vcov[vc$grp == "Residual"]
    if (length(donor_var) > 0 && length(resid_var) > 0 && (donor_var + resid_var) > 0) {
      icc <- donor_var / (donor_var + resid_var)
    }
  }

  error_msg <- NULL
  if (is.na(p_val)) error_msg <- "Could not extract p-value from model ANOVA table"
  if (lme4::isSingular(fit)) error_msg <- paste0(error_msg %||% "", "Model is singular (zero donor variance)")

  list(p_value = p_val, icc = icc, fit = fit, error_msg = error_msg)
}

#' Islet-size companion models to the primary mixed-effects test.
#'
#' The primary model is `value ~ donor_status + (1|donor)` \u2014 it conditions on
#' neither islet diameter nor a group-varying size slope. Both matter on this
#' cohort:
#'   * diameter is a strong predictor (main effect p < 2e-16 for most features);
#'   * the diameter DISTRIBUTION differs by group (ND vs T1D Kolmogorov-Smirnov
#'     D = 0.091, p = 2.9e-06);
#'   * and for some features the feature-vs-size SLOPE differs by disease status
#'     (donor_status x diameter p < 2e-16 for prop_Beta cell and Ins_any).
#'
#' The third point is the important one. When the interaction is significant the
#' group difference is not a single number \u2014 it depends on islet size \u2014 and the
#' primary omnibus p is a marginal average over a size mix that is itself
#' group-dependent. That is an estimand issue, not bias, but it has to be visible
#' rather than left for the user to infer from the Section 3 per-bin views.
#'
#' Diameter is MEAN-CENTRED so that, under lmerTest's Type III anova, the
#' `donor_status` term in the additive model is the group difference at the mean
#' islet diameter rather than at the unobservable diameter of zero.
#'
#' Returns list(p_status_adj, p_diam, p_interaction, mean_diam, error_msg) with
#' NA fields on failure \u2014 never NULL, so callers can render unconditionally.
size_companion_models <- function(rdf, diam_col = "islet_diam_um") {
  out <- list(p_status_adj = NA_real_, p_diam = NA_real_, p_interaction = NA_real_,
              mean_diam = NA_real_, error_msg = NULL)
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    out$error_msg <- "lmerTest not installed"; return(out)
  }
  if (!all(c(diam_col, "Case ID", "value", "donor_status") %in% colnames(rdf))) {
    out$error_msg <- "required columns missing"; return(out)
  }
  d <- rdf
  d$.diam <- suppressWarnings(as.numeric(d[[diam_col]]))
  d <- d[is.finite(d$.diam) & is.finite(d$value), , drop = FALSE]
  # Need at least 2 donor-status levels present and some size spread.
  if (nrow(d) < 10 || dplyr::n_distinct(d$donor_status) < 2 ||
      !is.finite(stats::sd(d$.diam)) || stats::sd(d$.diam) == 0) {
    out$error_msg <- "insufficient data or no diameter variation"; return(out)
  }
  out$mean_diam <- mean(d$.diam)
  d$.diam_c <- d$.diam - out$mean_diam
  d$.cf <- factor(d$`Case ID`)
  d$donor_status <- droplevels(factor(d$donor_status))

  p_of <- function(fit, term) {
    tryCatch({
      at <- anova(fit)
      if (term %in% rownames(at)) as.numeric(at[term, "Pr(>F)"]) else NA_real_
    }, error = function(e) NA_real_)
  }

  fit_add <- tryCatch(
    lmerTest::lmer(value ~ donor_status + .diam_c + (1 | .cf), data = d),
    error = function(e) e)
  if (inherits(fit_add, "error")) {
    out$error_msg <- paste("size-adjusted model failed:", fit_add$message)
  } else {
    out$p_status_adj <- p_of(fit_add, "donor_status")
    out$p_diam <- p_of(fit_add, ".diam_c")
  }

  fit_int <- tryCatch(
    lmerTest::lmer(value ~ donor_status * .diam_c + (1 | .cf), data = d),
    error = function(e) e)
  if (inherits(fit_int, "error")) {
    out$error_msg <- paste(c(out$error_msg,
                             paste("interaction model failed:", fit_int$message)),
                           collapse = "; ")
  } else {
    out$p_interaction <- p_of(fit_int, "donor_status:.diam_c")
  }
  out
}

#' Age-group companion to the primary mixed-effects test.
#'
#' Fits `value ~ donor_status + age_group + (1 | donor)`, the covariate form used
#' by Steenbuck et al. 2025 (Bodenmiller lab, bioRxiv 2025.03.05.641526) to test
#' whether disease-stage differences survive adjustment for age-associated
#' endotypes. They pooled donors at 13 years; that cut is data-supported on this
#' cohort rather than arbitrary, because our ages are bimodal with a gap from
#' 4.22 to 14.56 years and 13 falls inside it.
#'
#' Age is read from the H5AD `age` column, NOT from the nPOD metadata join. The
#' two disagree for exactly one donor (112: H5AD 3.00 vs nPOD 3.72) and the H5AD
#' value is canonical for modelling so existing results do not shift. Both fall
#' under 13, so the binning is identical either way \u2014 but the rule matters if
#' another donor is ever added.
#'
#' Deliberately NOT offered here: HbA1c and C-peptide. Both are downstream
#' consequences of beta-cell loss rather than confounders of it, so adjusting for
#' them is mediator adjustment and would bias the disease effect toward the null.
#' They belong in the demographics table and as external-validity checks (does
#' clinical C-peptide track our imaging-derived beta-cell proportion?), not in
#' this model. Age, sex and BMI are the legitimate confounders.
#'
#' Returns list(p_status_adj, p_age_group, n_young, n_old, cell_counts,
#' error_msg) with NA fields on failure \u2014 never NULL, so callers can render
#' unconditionally.
age_group_companion_model <- function(rdf, cut = 13) {
  out <- list(p_status_adj = NA_real_, p_age_group = NA_real_,
              n_young = NA_integer_, n_old = NA_integer_,
              cell_counts = NULL, cut = cut, error_msg = NULL)
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    out$error_msg <- "lmerTest not installed"; return(out)
  }
  if (!all(c("age", "Case ID", "value", "donor_status") %in% colnames(rdf))) {
    out$error_msg <- "age / Case ID column missing"; return(out)
  }
  d <- rdf
  d$.age <- suppressWarnings(as.numeric(d$age))
  d <- d[is.finite(d$.age) & is.finite(d$value), , drop = FALSE]
  if (nrow(d) < 10) { out$error_msg <- "too few rows with a finite age"; return(out) }

  # ">=13" first so it is the reference level: the age_group coefficient then
  # reads as the effect of being in the young group.
  d$age_group <- factor(ifelse(d$.age < cut, "<13", ">=13"), levels = c(">=13", "<13"))
  d$donor_status <- droplevels(factor(d$donor_status))
  d$.cf <- factor(d$`Case ID`)

  donors <- unique(d[, c("Case ID", "donor_status", "age_group")])
  out$cell_counts <- table(donors$age_group, donors$donor_status)
  out$n_young <- sum(donors$age_group == "<13")
  out$n_old <- sum(donors$age_group == ">=13")

  if (dplyr::n_distinct(d$age_group) < 2) {
    out$error_msg <- sprintf("all donors fall on one side of the %g-year cut", cut)
    return(out)
  }
  if (dplyr::n_distinct(d$donor_status) < 2) {
    out$error_msg <- "fewer than 2 donor-status groups"; return(out)
  }

  fit <- tryCatch(
    lmerTest::lmer(value ~ donor_status + age_group + (1 | .cf), data = d),
    error = function(e) e)
  if (inherits(fit, "error")) {
    out$error_msg <- paste("age-adjusted model failed:", fit$message); return(out)
  }
  p_of <- function(term) tryCatch({
    at <- anova(fit)
    if (term %in% rownames(at)) as.numeric(at[term, "Pr(>F)"]) else NA_real_
  }, error = function(e) NA_real_)
  out$p_status_adj <- p_of("donor_status")
  out$p_age_group <- p_of("age_group")
  out
}

#' Shared per-bin coverage rule \u2014 the single definition of "this bin is analysable".
#'
#' A bin qualifies only if, after aggregating to donor means within the bin, at
#' least 2 donor_status groups are present AND every group present in the bin has
#' >= 2 donors.
#'
#' This lives in one function because it previously did not. `per_bin_donor_anova`
#' required every present group to clear 2 donors; `per_bin_donor_kendall` applied
#' an equivalent rule on the status codes; but `per_bin_donor_pairwise` had no bin
#' rule at all and only checked the two groups in the contrast. The 275 \u00b5m bin
#' (ND = 1 donor, Aab+ = 3, T1D = 2) therefore produced an `Aab+ vs T1D` star on
#' the stratified heatmap with no ANOVA or Kendall row above it \u2014 a pairwise
#' post-hoc with no omnibus behind it.
#'
#' All three now gate on this function, so a bin either supports the whole
#' stratified analysis or contributes nothing to it.
#'
#' @param group_values donor_status labels, one per DONOR (already aggregated)
#' @return TRUE if the bin can be analysed
bin_coverage_ok <- function(group_values) {
  g <- as.character(group_values)
  g <- g[!is.na(g) & nzchar(g)]
  if (!length(g)) return(FALSE)
  n_per <- table(g)
  length(n_per) >= 2 && all(n_per >= 2)
}

#' Per-bin ANOVA on donor-level means within each bin
#
# Coverage: bin_coverage_ok() \u2014 see there. Bins that fail are omitted, so the
# heatmap (ANOVA + pairwise rows) and the trend plot (Kendall tau) cover exactly
# the same bins.
per_bin_donor_anova <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  if (!"Case ID" %in% colnames(df)) return(per_bin_anova(df, bin_col, group_col, value_col, mid_col))
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    # Aggregate to donor means within this bin
    ddf <- sub %>%
      dplyr::group_by(`Case ID`, .data[[group_col]]) %>%
      dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    colnames(ddf)[colnames(ddf) == group_col] <- "grp"
    if (!bin_coverage_ok(ddf$grp)) return(NULL)
    fit <- tryCatch(aov(value ~ grp, data = ddf), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pval <- tryCatch({
      at <- anova(fit)
      as.numeric(at[["Pr(>F)"]][1])
    }, error = function(e) NA_real_)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_anova = pval)
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  arrange(out, mid)
}

#' Per-bin Kendall tau on donor-level means within each bin
#
# Coverage rule matches per_bin_donor_anova: a bin is tested only if at least
# 2 donor_status groups each have >= 2 donors with finite values. This keeps
# the trend plot's bins aligned with the heatmap.
per_bin_donor_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  if (!"Case ID" %in% colnames(df)) return(per_bin_kendall(df, bin_col, group_col, value_col, mid_col))
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    # Aggregate to donor means within this bin
    ddf <- sub %>%
      dplyr::group_by(`Case ID`, .data[[group_col]]) %>%
      dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    grp <- as.character(ddf[[group_col]])
    x <- code_donor_status(ddf[[group_col]])
    y <- ddf$value
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]; grp <- grp[keep]
    if (length(y) < 3) return(NULL)
    if (!bin_coverage_ok(grp)) return(NULL)
    ct <- tryCatch(cor.test(x, y, method = "kendall", exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_kendall = unname(ct$p.value), tau = unname(ct$estimate))
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(out)
  arrange(out, mid)
}

#' Per-bin pairwise tests on donor-level means within each bin
#' Returns data.frame with bin, mid, pair, p_value
per_bin_donor_pairwise <- function(df, bin_col, group_col, value_col,
                                    mid_col = "diam_mid", parametric = TRUE) {
  if (nrow(df) == 0 || !"Case ID" %in% colnames(df)) {
    return(tibble(bin = character(), mid = numeric(), pair = character(), p_value = numeric()))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  pairs_list <- list(c("ND", "Aab+"), c("ND", "T1D"), c("Aab+", "T1D"))
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    ddf <- sub %>%
      dplyr::group_by(`Case ID`, .data[[group_col]]) %>%
      dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    # Gate the whole bin on the omnibus coverage rule: a pairwise post-hoc does
    # not run in a bin that has no ANOVA behind it. Without this the heatmap drew
    # pairwise stars in bins whose omnibus row was blank.
    if (!bin_coverage_ok(ddf[[group_col]])) return(NULL)
    pair_rows <- lapply(pairs_list, function(pair) {
      g1 <- ddf$value[ddf[[group_col]] == pair[1]]
      g2 <- ddf$value[ddf[[group_col]] == pair[2]]
      if (length(g1) < 2 || length(g2) < 2) return(NULL)
      p_val <- if (parametric) {
        tt <- tryCatch(t.test(g1, g2), error = function(e) NULL)
        if (!is.null(tt)) tt$p.value else NA_real_
      } else {
        wt <- tryCatch(wilcox.test(g1, g2, exact = FALSE), error = function(e) NULL)
        if (!is.null(wt)) wt$p.value else NA_real_
      }
      tibble(bin = as.character(b), mid = bmeta$mid[i],
             pair = paste(pair[1], "vs", pair[2]), p_value = p_val)
    })
    bind_rows(pair_rows)
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), pair = character(), p_value = numeric()))
  arrange(out, mid, pair)
}

#' Shapiro-Wilk normality test per donor group on donor-level values
#' Returns data.frame with donor_status, W, p_value, is_normal
normality_tests <- function(ddf) {
  groups <- levels(ddf$donor_status)
  if (is.null(groups)) groups <- unique(ddf$donor_status)
  rows <- lapply(groups, function(g) {
    vals <- ddf$value[ddf$donor_status == g]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 3) {
      return(data.frame(donor_status = g, W = NA_real_, p_value = NA_real_,
                        is_normal = NA, n = length(vals), stringsAsFactors = FALSE))
    }
    sw <- tryCatch(shapiro.test(vals), error = function(e) NULL)
    if (is.null(sw)) {
      return(data.frame(donor_status = g, W = NA_real_, p_value = NA_real_,
                        is_normal = NA, n = length(vals), stringsAsFactors = FALSE))
    }
    data.frame(donor_status = g, W = unname(sw$statistic), p_value = sw$p.value,
               is_normal = sw$p.value > 0.05, n = length(vals), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---------- p-value formatting ----------

#' Format a p-value for display \u2014 one rule, used everywhere.
#'
#' Scientific notation is only readable when the exponent is large. The
#' Statistics headline used to call formatC(format = "e") unconditionally, so
#' p = 0.3476 rendered as "3.48e-01" \u2014 one keystroke away from being read as
#' 0.0348 \u2014 while the hypothesis table three lines below printed the same number
#' as "0.3476". Fixed notation down to `sci_below`, scientific under it.
#'
#' @param x numeric p-value (NA-safe, vectorised)
#' @param digits significant digits for the scientific branch
#' @param sci_below switch to scientific notation below this value
#' @param na label for NA / non-finite input
format_p_value <- function(x, digits = 2, sci_below = 1e-4, na = "N/A") {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.finite(x), na,
         ifelse(x < sci_below,
                formatC(x, format = "e", digits = digits),
                formatC(x, format = "f", digits = 4)))
}

# ---------- Effect size utilities (Statistics tab) ----------

#' Effect size for a two-sample comparison.
#'
#' Returns Cohen's d (pooled SD, 95% CI by normal approximation) AND the
#' unstandardised difference with a Welch CI, plus the per-group SDs and their
#' variance ratio.
#'
#' Why the extra fields: Cohen's d divides by a POOLED SD, which is only
#' meaningful when the two groups have comparable spread. On this cohort they
#' routinely do not \u2014 for prop_Beta cell the donor-level SDs are ND 4.06 vs
#' T1D 0.575, a 50:1 variance ratio (a floor effect: T1D beta-cell proportions
#' sit at ~0). The pooled SD lands at 2.90, giving d = 6.15 with CI [3.18, 9.12]
#' \u2014 but Glass's delta is 4.39 against the ND SD and 31.03 against the T1D SD.
#' d = 6.15 is an arbitrary point in that range, presented with a CI that implies
#' precision it does not have. `raw_diff` (17.85 percentage points here) is the
#' quantity that does not depend on which SD you standardise by.
#'
#' `var_ratio_warn` flags var ratio > 4, the point past which the pooled SD is
#' not a reasonable common scale.
#'
#' Backwards compatible: `$d`, `$ci_lo`, `$ci_hi` are unchanged.
cohens_d <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  empty <- list(d = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                raw_diff = NA_real_, raw_ci_lo = NA_real_, raw_ci_hi = NA_real_,
                mean1 = NA_real_, mean2 = NA_real_, sd1 = NA_real_, sd2 = NA_real_,
                glass1 = NA_real_, glass2 = NA_real_,
                var_ratio = NA_real_, var_ratio_warn = FALSE)
  if (n1 < 2 || n2 < 2) return(empty)

  m1 <- mean(x); m2 <- mean(y)
  s1 <- sd(x);   s2 <- sd(y)
  raw_diff <- m1 - m2

  # Welch CI on the raw difference \u2014 does not assume equal variances.
  wc <- tryCatch(t.test(x, y)$conf.int, error = function(e) c(NA_real_, NA_real_))

  v1 <- var(x); v2 <- var(y)
  var_ratio <- if (is.finite(v1) && is.finite(v2) && min(v1, v2) > 0) {
    max(v1, v2) / min(v1, v2)
  } else NA_real_

  # Glass's delta against each group's SD in turn \u2014 the two ends of the range
  # the pooled-SD d silently picks a point inside of.
  glass1 <- if (is.finite(s1) && s1 > 0) raw_diff / s1 else NA_real_
  glass2 <- if (is.finite(s2) && s2 > 0) raw_diff / s2 else NA_real_

  pooled_sd <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    out <- empty
    out$raw_diff <- raw_diff; out$raw_ci_lo <- wc[1]; out$raw_ci_hi <- wc[2]
    out$mean1 <- m1; out$mean2 <- m2; out$sd1 <- s1; out$sd2 <- s2
    out$glass1 <- glass1; out$glass2 <- glass2
    out$var_ratio <- var_ratio
    out$var_ratio_warn <- isTRUE(var_ratio > 4)
    return(out)
  }
  d <- raw_diff / pooled_sd
  se_d <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))

  list(d = d, ci_lo = d - 1.96 * se_d, ci_hi = d + 1.96 * se_d,
       raw_diff = raw_diff, raw_ci_lo = wc[1], raw_ci_hi = wc[2],
       mean1 = m1, mean2 = m2, sd1 = s1, sd2 = s2,
       glass1 = glass1, glass2 = glass2,
       var_ratio = var_ratio, var_ratio_warn = isTRUE(var_ratio > 4))
}

#' Eta-squared from anova() output: SS_effect / SS_total for the first factor
eta_squared <- function(fit) {
  at <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(at)) return(NA_real_)
  ss <- at[["Sum Sq"]]
  if (length(ss) < 2) return(NA_real_)
  ss[1] / sum(ss)
}

# ---------------------------------------------------------------------------
# Carried over from the Senior lane during the app merge (2026-08-16). Junior's
# copy of this file never had it. It has no caller today -- mod_statistics_server.R
# only names it in a header comment -- but it is a working, self-contained helper
# and deleting it silently during a merge would be a loss, not a simplification.
# ---------------------------------------------------------------------------
pairwise_wilcox <- function(df, group_col, value_col) {
  pw <- tryCatch(
    pairwise.wilcox.test(df[[value_col]], df[[group_col]], p.adjust.method = "BH", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(pw)) return(data.frame(group1 = character(), group2 = character(), p_value = numeric(), stringsAsFactors = FALSE))
  mat <- as.data.frame(as.table(pw$p.value))
  colnames(mat) <- c("group1", "group2", "p_value")
  mat <- mat[!is.na(mat$p_value), , drop = FALSE]
  mat$group1 <- as.character(mat$group1)
  mat$group2 <- as.character(mat$group2)
  mat
}
