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
  code_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    x <- code_group(sub[[group_col]])
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
per_bin_donor_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  if (!"Case ID" %in% colnames(df)) return(per_bin_kendall(df, bin_col, group_col, value_col, mid_col))
  code_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    # Aggregate to donor means within this bin
    ddf <- sub %>%
      dplyr::group_by(`Case ID`, .data[[group_col]]) %>%
      dplyr::summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
    x <- code_group(ddf[[group_col]])
    y <- ddf$value
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
