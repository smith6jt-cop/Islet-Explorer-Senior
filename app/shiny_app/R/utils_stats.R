# ---------- Unified outlier flagging (used by Plot, Trajectory, Statistics, Spatial) ----------
#
# Single source of truth for outlier identification across the app. Each tab
# previously had its own inline loop with slightly different rules (z-score on
# raw values vs Tukey IQR vs 99th-percentile clip). The unified rule is:
#   For each donor_status group (ND, Aab+, T1D), exclude rows whose value is
#   more than `threshold` SDs from that group's mean. Outliers are flagged on
#   the returned data.frame; downstream code decides whether to render them as
#   ghost points or drop them. The IQR fallback is kept for cards that need to
#   tolerate zero-inflated distributions (Spatial Card A on phenotypes that are
#   absent from most islets).
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

  for (grp in unique(groups)) {
    idx <- which(groups == grp & is.finite(vals))
    if (length(idx) < 3) next
    grp_vals <- vals[idx]

    if (exclude_zeros) {
      nz_idx <- idx[grp_vals != 0]
      nz_vals <- vals[nz_idx]
      if (length(nz_vals) < 4) next
      if (method == "zscore") {
        m <- mean(nz_vals, na.rm = TRUE)
        s <- sd(nz_vals, na.rm = TRUE)
        if (!is.finite(s) || s == 0) next
        z <- (vals[idx] - m) / s
        df$outlier_score[idx] <- z
        df$is_outlier[idx] <- !is.na(z) & (vals[idx] != 0) & abs(z) > threshold
      } else {
        q <- stats::quantile(nz_vals, c(0.25, 0.75), na.rm = TRUE)
        iqr <- q[2] - q[1]
        if (!is.finite(iqr) || iqr == 0) next
        lo <- q[1] - threshold * iqr
        hi <- q[2] + threshold * iqr
        fence_dist <- ifelse(vals[idx] < q[1],
                             (vals[idx] - q[1]) / iqr,
                             ifelse(vals[idx] > q[2], (vals[idx] - q[2]) / iqr, 0))
        df$outlier_score[idx] <- fence_dist
        df$is_outlier[idx] <- !is.na(vals[idx]) & (vals[idx] != 0) &
                              (vals[idx] < lo | vals[idx] > hi)
      }
    } else {
      if (method == "zscore") {
        m <- mean(grp_vals, na.rm = TRUE)
        s <- sd(grp_vals, na.rm = TRUE)
        if (!is.finite(s) || s == 0) next
        z <- (grp_vals - m) / s
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
    sprintf("z > %g per donor status", threshold)
  } else {
    sprintf("Tukey %g×IQR per donor status", threshold)
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
  sprintf(" — %d outlier%s excluded (%s)",
          summary$n_excluded,
          if (summary$n_excluded == 1) "" else "s",
          summary$label)
}

# Convenience: shared tooltip text for the global outlier control.
GLOBAL_OUTLIER_TOOLTIP <- paste0(
  "For each donor status (ND / Aab+ / T1D), islets whose value is more than ",
  "3 standard deviations from that group's mean are excluded from trend lines, ",
  "donor-level means, and statistical tests. Excluded points still appear on ",
  "scatter and box plots as faint grey × markers. Counts are shown next to ",
  "each plot title. Spatial Card A (infiltration) uses Tukey 1.5×IQR instead, ",
  "to tolerate zero-inflated phenotypes; Spatial Card C distance plots apply ",
  "the z-score rule on log(1 + distance)."
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

# Map donor_status to ordinal trend code used by Kendall tau.
# ND -> 0, Aab+ -> 1, T1D -> 2. Other labels return NA.
# Used by per_bin_kendall() and per_bin_donor_kendall() so the encoding lives
# in one place. If donor cohort labels ever change, update only here.
code_donor_status <- function(x) {
  x <- as.character(x)
  ifelse(x == "ND", 0,
    ifelse(x == "Aab+", 1,
      ifelse(x == "T1D", 2, NA_real_)))
}

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

#' Per-bin ANOVA on donor-level means within each bin
#
# Shared coverage rule (used by per_bin_donor_anova and per_bin_donor_kendall):
#   A bin is tested only if, after aggregating to donor means within the bin,
#   at least 2 donor_status groups each have >= 2 donors. Bins that fail this
#   rule are silently omitted so the heatmap (ANOVA + pairwise rows) and the
#   trend plot (Kendall tau) include the same bins.
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
    # Need >= 2 groups with >= 2 donors each
    grp_n <- ddf %>% dplyr::group_by(grp) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop")
    if (nrow(grp_n) < 2 || any(grp_n$n < 2)) return(NULL)
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
    x <- code_donor_status(ddf[[group_col]])
    y <- ddf$value
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    if (length(y) < 3) return(NULL)
    # Match per_bin_donor_anova coverage: at least 2 distinct groups present,
    # each present group with >= 2 donors. ANOVA cannot estimate within-group
    # variance otherwise; we apply the same rule here so the heatmap and trend
    # plot include the same bins.
    code_counts <- table(x)
    if (length(code_counts) < 2 || any(code_counts < 2)) return(NULL)
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

# ---------- Effect size utilities (Statistics tab) ----------

#' Cohen's d for two-sample comparison with 95% CI (normal approximation)
cohens_d <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(d = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  pooled_sd <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(list(d = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_))
  d <- (mean(x) - mean(y)) / pooled_sd
  se_d <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  list(d = d, ci_lo = d - 1.96 * se_d, ci_hi = d + 1.96 * se_d)
}

#' Eta-squared from anova() output: SS_effect / SS_total for the first factor
eta_squared <- function(fit) {
  at <- tryCatch(anova(fit), error = function(e) NULL)
  if (is.null(at)) return(NA_real_)
  ss <- at[["Sum Sq"]]
  if (length(ss) < 2) return(NA_real_)
  ss[1] / sum(ss)
}

#' Non-parametric pairwise Wilcoxon tests with BH correction
#' Returns data.frame with group1, group2, p_value columns
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
