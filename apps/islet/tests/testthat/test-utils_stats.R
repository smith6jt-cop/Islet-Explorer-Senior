# Unit tests for the statistical utilities in R/utils_stats.R.
# These functions back the Statistics tab (donor-level aggregation, mixed-effects
# test, per-bin analysis, effect sizes, outlier flagging) and are pure R — no app
# data, no Python bridge — so they are the natural CI test surface.

test_that("code_donor_status maps disease stage to ordinal trend codes", {
  expect_equal(code_donor_status(c("ND", "Aab+", "T1D")), c(0, 1, 2))
  expect_true(is.na(code_donor_status("Other")))
  # accepts factors (it coerces with as.character first)
  expect_equal(code_donor_status(factor("T1D", levels = c("ND", "Aab+", "T1D"))), 2)
})

test_that("cohens_d returns d with ci_lo/ci_hi and guards small samples", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(2, 3, 4, 5, 6)
  res <- cohens_d(x, y)
  # Contract: names are d / ci_lo / ci_hi (NOT ci_lower / ci_upper).
  # The Jul 2026 group-comparability pass WIDENED this return value with the
  # unstandardised difference + Welch CI, per-group means/SDs, both Glass's
  # deltas and the variance ratio -- but did not update this test, so the suite
  # went red at that commit and stayed red. Assert the required fields are
  # present rather than pinning the exact set, so adding a field is not a
  # breaking change while a RENAME still fails loudly.
  expect_true(all(c("d", "ci_lo", "ci_hi") %in% names(res)))
  expect_true(all(c("raw_diff", "raw_ci_lo", "raw_ci_hi", "mean1", "mean2",
                    "sd1", "sd2", "glass1", "glass2", "var_ratio",
                    "var_ratio_warn") %in% names(res)))
  # pooled SD = sqrt(2.5); mean diff = -1  =>  d = -1 / sqrt(2.5)
  expect_equal(res$d, -1 / sqrt(2.5), tolerance = 1e-6)
  expect_lt(res$ci_lo, res$d)
  expect_gt(res$ci_hi, res$d)
  # the unstandardised difference is the plain mean difference
  expect_equal(res$raw_diff, -1, tolerance = 1e-6)
  # fewer than 2 obs per group -> all NA
  na_res <- cohens_d(c(1), c(2, 3))
  expect_true(is.na(na_res$d))
  expect_true(is.na(na_res$ci_lo))
})

test_that("eta_squared is bounded in [0,1] and large for separated groups", {
  set.seed(1)
  g <- factor(rep(c("a", "b", "c"), each = 10))
  val <- c(rnorm(10, 0), rnorm(10, 10), rnorm(10, 20))  # strong separation
  es <- eta_squared(lm(val ~ g))
  expect_gte(es, 0)
  expect_lte(es, 1)
  expect_gt(es, 0.8)
})

test_that("aggregate_to_donor collapses islets to one row per donor", {
  rdf <- tibble(
    `Case ID`    = c("A", "A", "A", "B", "B"),
    donor_status = c("ND", "ND", "ND", "T1D", "T1D"),
    value        = c(1, 3, 5, 10, 20),
    age          = c(30, 30, 30, 50, 50),
    gender       = c("M", "M", "M", "F", "F")
  )
  out <- aggregate_to_donor(rdf, "mean")
  expect_equal(nrow(out), 2)
  a <- out[out$`Case ID` == "A", ]
  expect_equal(a$value, 3)      # mean(1, 3, 5)
  expect_equal(a$n_islets, 3)
  expect_equal(a$age, 30)       # metadata preserved
  expect_equal(a$gender, "M")
  # no 'Case ID' column -> returned unchanged
  bare <- tibble(donor_status = "ND", value = 1)
  expect_identical(aggregate_to_donor(bare), bare)
})

test_that("normality_tests reports per-group W/p and flags undersized groups", {
  set.seed(2)
  ddf <- data.frame(
    donor_status = factor(c(rep("ND", 5), rep("Aab+", 5), rep("T1D", 2)),
                          levels = c("ND", "Aab+", "T1D")),
    value        = c(rnorm(5, 5), rnorm(5, 6), c(1, 2))
  )
  nr <- normality_tests(ddf)
  expect_true(all(c("donor_status", "W", "p_value", "is_normal", "n") %in% names(nr)))
  expect_setequal(as.character(nr$donor_status), c("ND", "Aab+", "T1D"))
  # T1D has only 2 finite values (< 3) -> untestable -> NA flag, n = 2
  t1d <- nr[nr$donor_status == "T1D", ]
  expect_true(is.na(t1d$is_normal))
  expect_equal(t1d$n, 2)
  # ND has 5 -> a real W statistic
  nd <- nr[nr$donor_status == "ND", ]
  expect_equal(nd$n, 5)
  expect_false(is.na(nd$W))
})

test_that("flag_outliers flags extreme values per group and honours enabled=FALSE", {
  # 1..9 are inliers, 100 is a clear outlier (modified z >> 3 with MAD > 0)
  df <- data.frame(donor_status = rep("ND", 10), value = c(1:9, 100))
  fl <- flag_outliers(df, "value", "donor_status", method = "zscore",
                      threshold = 3, enabled = TRUE)
  expect_true(all(c("is_outlier", "outlier_score") %in% names(fl)))
  expect_true(fl$is_outlier[10])
  expect_false(any(fl$is_outlier[1:9]))
  # disabled -> nothing flagged
  off <- flag_outliers(df, "value", "donor_status", enabled = FALSE)
  expect_false(any(off$is_outlier))
  # empty df -> zero-length flag columns, no error
  empty <- flag_outliers(data.frame(donor_status = character(), value = numeric()))
  expect_length(empty$is_outlier, 0)
})

test_that("flag_outliers(exclude_zeros=TRUE) removes zeros from the scoring population", {
  # 12 zeros + the inlier run 1..9 + one clear outlier. With zeros included the
  # median/MAD collapse toward 0; excluding them must derive the scale from the
  # non-zero values only.
  df <- data.frame(donor_status = rep("ND", 22),
                   value = c(rep(0, 12), 1:9, 100))
  fl <- flag_outliers(df, "value", "donor_status", method = "zscore",
                      threshold = 3, exclude_zeros = TRUE, enabled = TRUE)

  # Zeros are not evaluated at all: no score, never flagged.
  expect_true(all(is.na(fl$outlier_score[1:12])))
  expect_false(any(fl$is_outlier[1:12]))

  # Non-zero rows are scored, and 100 is still the outlier.
  expect_true(all(is.finite(fl$outlier_score[13:22])))
  expect_true(fl$is_outlier[22])
  expect_false(any(fl$is_outlier[13:21]))

  # Centre/scale come from the non-zero subset, not the full vector.
  nz <- c(1:9, 100)
  med <- median(nz); mad_nz <- median(abs(nz - med))
  expect_equal(fl$outlier_score[13], 0.6745 * (1 - med) / mad_nz)

  # Contrast: with zeros included the MAD collapses to 0, so the helper falls
  # back to classical mean/SD — a different statistic on a different scale.
  incl <- flag_outliers(df, "value", "donor_status", method = "zscore",
                        threshold = 3, exclude_zeros = FALSE, enabled = TRUE)
  expect_equal(attr(fl, "outlier_statistic"), "modified_z")
  expect_equal(attr(incl, "outlier_statistic"), "classical_z")
  expect_true(all(is.finite(incl$outlier_score)))          # zeros ARE scored here
  # Same rows, different scores — this is why the filter order matters.
  expect_false(isTRUE(all.equal(incl$outlier_score[13:22], fl$outlier_score[13:22])))
})

test_that("flag_outliers(exclude_zeros=TRUE, method='iqr') fences on the non-zero subset", {
  df <- data.frame(donor_status = rep("ND", 20),
                   value = c(rep(0, 10), 10, 11, 12, 13, 14, 15, 16, 17, 18, 500))
  fl <- flag_outliers(df, "value", "donor_status", method = "iqr",
                      threshold = 1.5, exclude_zeros = TRUE, enabled = TRUE)
  expect_true(all(is.na(fl$outlier_score[1:10])))
  expect_false(any(fl$is_outlier[1:10]))
  expect_true(fl$is_outlier[20])          # 500 is beyond the upper fence
  expect_false(any(fl$is_outlier[11:19]))
})

test_that("flag_outliers(exclude_zeros=TRUE) skips groups with <3 non-zero values", {
  df <- data.frame(donor_status = rep("ND", 8), value = c(rep(0, 6), 1, 99))
  fl <- flag_outliers(df, "value", "donor_status", method = "zscore",
                      threshold = 3, exclude_zeros = TRUE, enabled = TRUE)
  expect_true(all(is.na(fl$outlier_score)))
  expect_false(any(fl$is_outlier))
})

test_that("summarize_outlier_filter + outlier_title_suffix report exclusion counts", {
  df <- data.frame(donor_status = rep("ND", 10), value = c(1:9, 100))
  fl <- flag_outliers(df, "value", "donor_status", method = "zscore",
                      threshold = 3, enabled = TRUE)
  s <- summarize_outlier_filter(fl)
  expect_equal(s$n_total, 10)
  expect_equal(s$n_excluded, 1)
  expect_equal(s$n_kept, 9)
  expect_true(s$enabled)
  expect_match(outlier_title_suffix(s), "1 outlier excluded")
  # filter OFF -> empty suffix
  off <- flag_outliers(df, enabled = FALSE)
  expect_equal(outlier_title_suffix(summarize_outlier_filter(off)), "")
})

test_that("summary_stats produces mean/se and median/iqr summaries", {
  df <- data.frame(g = rep(c("a", "b"), each = 4),
                   v = c(1, 2, 3, 4, 10, 20, 30, 40))
  ms <- summary_stats(df, "g", "v", stat = "mean_se")
  expect_true(all(c("n", "mean", "sd", "se", "median", "q1", "q3",
                    "y", "ymin", "ymax") %in% names(ms)))
  a <- ms[ms$g == "a", ]
  expect_equal(a$mean, 2.5)
  expect_equal(a$y, 2.5)                    # mean_se -> y = mean
  # unname(): q1/q3 inherit a "25%"/"75%" name from quantile() that propagates
  # through case_when() into ymin; we only care about the numeric value here.
  expect_equal(unname(a$ymin), unname(a$mean - a$se))
  mi <- summary_stats(df, "g", "v", stat = "median_iqr")
  b <- mi[mi$g == "a", ]
  expect_equal(unname(b$y), unname(b$median))   # median_iqr -> y = median
  expect_equal(unname(b$ymin), unname(b$q1))
})

# ---- Per-bin donor-level tests (Statistics Section 3) ----

# Helper: two islets for one donor in a bin, centred on `base`.
.donor_islets <- function(bin, mid, cid, status, base) {
  data.frame(diam_bin = bin, diam_mid = mid, `Case ID` = cid,
             donor_status = status, value = c(base, base + 0.1),
             check.names = FALSE, stringsAsFactors = FALSE)
}

.make_bin_df <- function() {
  b25 <- rbind(
    .donor_islets("b25", 25, "A", "ND",   1.0),
    .donor_islets("b25", 25, "B", "ND",   1.2),
    .donor_islets("b25", 25, "C", "Aab+", 5.0),
    .donor_islets("b25", 25, "D", "Aab+", 5.2),
    .donor_islets("b25", 25, "E", "T1D",  10.0),
    .donor_islets("b25", 25, "F", "T1D",  10.2)
  )
  # bin b75 has only ND donors -> fails the >=2 groups coverage rule -> skipped
  b75 <- rbind(
    .donor_islets("b75", 75, "A", "ND", 1.0),
    .donor_islets("b75", 75, "B", "ND", 1.1)
  )
  rbind(b25, b75)
}

test_that("per_bin_donor_anova tests only bins with adequate donor coverage", {
  df <- .make_bin_df()
  ba <- per_bin_donor_anova(df, "diam_bin", "donor_status", "value")
  expect_true(all(c("bin", "mid", "p_anova") %in% names(ba)))
  expect_equal(nrow(ba), 1)          # only b25 is testable
  expect_equal(ba$bin, "b25")
  expect_lt(ba$p_anova, 0.05)        # strong ND/Aab+/T1D separation
})

test_that("per_bin_donor_kendall returns tau and matches the coverage rule", {
  df <- .make_bin_df()
  bk <- per_bin_donor_kendall(df, "diam_bin", "donor_status", "value")
  expect_true(all(c("bin", "mid", "p_kendall", "tau") %in% names(bk)))
  expect_equal(nrow(bk), 1)
  expect_gt(bk$tau, 0)               # value increases ND -> Aab+ -> T1D
})

test_that("per_bin_donor_pairwise emits the three named contrasts", {
  df <- .make_bin_df()
  bp <- per_bin_donor_pairwise(df, "diam_bin", "donor_status", "value",
                               parametric = TRUE)
  expect_true(all(c("bin", "mid", "pair", "p_value") %in% names(bp)))
  expect_setequal(unique(bp$pair), c("ND vs Aab+", "ND vs T1D", "Aab+ vs T1D"))
})

test_that("lmer_test_donor returns the expected result structure", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  set.seed(4)
  rdf <- do.call(rbind, lapply(1:6, function(i) {
    status <- c("ND", "ND", "Aab+", "Aab+", "T1D", "T1D")[i]
    base   <- c(1, 1.2, 5, 5.2, 10, 10.1)[i]
    data.frame(`Case ID` = paste0("D", i), donor_status = status,
               value = base + rnorm(8, 0, 0.3), check.names = FALSE)
  }))
  res <- lmer_test_donor(rdf)
  expect_true(all(c("p_value", "icc", "fit", "error_msg") %in% names(res)))
  expect_true(is.na(res$p_value) || (res$p_value >= 0 && res$p_value <= 1))
  expect_true(is.na(res$icc) || (res$icc >= 0 && res$icc <= 1))
  # missing 'Case ID' -> graceful error_msg, no crash
  bad <- lmer_test_donor(data.frame(donor_status = "ND", value = 1))
  expect_true(is.na(bad$p_value))
  expect_false(is.null(bad$error_msg))
})
