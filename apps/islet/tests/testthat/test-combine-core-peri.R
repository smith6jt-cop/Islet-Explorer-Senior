# The Plot tab's Core+Peri measure: one proportion made from two, each over its own cell total.
#
# THE DEFECT THIS PINS (found 2026-09-21). mod_plot_server.R's combined_prop_ branch read
#
#   core_n <- ifelse(is.finite(core_prop) & is.finite(core_total), core_prop * core_total, 0)
#   ctot   <- ifelse(is.finite(core_total), core_total, 0)
#   value  <- ifelse(combined_total > 0, 100 * combined_n / combined_total, NA_real_)
#
# so a NaN proportion on a side that HAS cells became 0 cells over a kept denominator, and the
# islet rendered at exactly 0. build_galectin_layer.py writes NaN for a blank count or a
# classifier disagreement (719 of 15,557 islet-regions on the 2026-09-14 layer), and
# data_loading.R pairs every peri_prop_* with its prop_*, so the galectin classes reach that
# branch. A missing measurement is unavailable -- not zero.
#
# THE BOUNDARY (2026-09-22). Each column is read once, in ie_combine_core_peri(), and refused
# there: a total is a count (finite, non-negative) or null; a proportion is in [0, 1] or null.
# A side is "no cells there" only when BOTH its total and its proportion are null (total 0 is
# a count of no cells and reads the same way) -- "an islet with no peri zone still contributes
# its core, and vice versa", the comment that stood on the fold, and it stays. A null total
# beside a finite proportion is a missing measurement: the proportion was measured and the
# denominator it needs was not.

.combined_pct <- function(x) ifelse(x$total > 0, 100 * x$n / x$total, NA_real_)

test_that("both sides finite gives the cell-weighted combination, as before", {
  x <- ie_combine_core_peri(core_prop = 0.2, peri_prop = 0.5, core_total = 10, peri_total = 20)
  expect_s3_class(x, "data.frame")
  expect_identical(nrow(x), 1L)
  expect_equal(x$n, 0.2 * 10 + 0.5 * 20)     # 12 cells
  expect_equal(x$total, 30)
  expect_equal(.combined_pct(x), 40)         # 100 * 12 / 30, the pre-2026-09-21 arithmetic
  expect_true(is.na(x$reason))
})

test_that("a non-finite core proportion over a positive core total is unavailable, not 0", {
  x <- ie_combine_core_peri(core_prop = NaN, peri_prop = 0.5, core_total = 5, peri_total = 20)
  expect_true(is.na(x$n))
  expect_true(is.na(.combined_pct(x)))
  # The denominator is a fact about the islet and is left as measured; the absence is in n.
  expect_equal(x$total, 25)
  expect_match(x$reason, "core", fixed = TRUE)
  expect_match(x$reason, "5", fixed = TRUE)
  # NA reads the same as NaN here: the proportion is not a number either way.
  y <- ie_combine_core_peri(NA_real_, 0.5, 5, 20)
  expect_true(is.na(y$n))
})

test_that("a non-finite peri proportion over a positive peri total is unavailable too", {
  x <- ie_combine_core_peri(core_prop = 0.2, peri_prop = NaN, core_total = 10, peri_total = 3)
  expect_true(is.na(x$n))
  expect_match(x$reason, "peri-islet", fixed = TRUE)
  both <- ie_combine_core_peri(NaN, NaN, 10, 3)
  expect_true(is.na(both$n))
  expect_match(both$reason, "core", fixed = TRUE)
  expect_match(both$reason, "peri-islet", fixed = TRUE)
})

test_that("a core total of 0 with a finite peri side gives the peri-only value, as before", {
  x <- ie_combine_core_peri(core_prop = 0.2, peri_prop = 0.5, core_total = 0, peri_total = 20)
  expect_equal(x$n, 10)
  expect_equal(x$total, 20)
  expect_equal(.combined_pct(x), 50)
  expect_true(is.na(x$reason))
})

# ITEM 1 (decided 2026-09-22). Two cases that the previous test in this slot read as one.

test_that("an NA total beside a finite proportion is a missing measurement, not an empty side", {
  # The peri proportion 0.5 was measured; the peri total it is a proportion OF is blank. The
  # previous test here expected the core-only value (20%) -- that folded a lost denominator
  # into "no cells there" and drew the islet at its core value alone.
  peri_lost <- ie_combine_core_peri(0.2, 0.5, 10, NA_real_)
  expect_true(is.na(peri_lost$n))
  expect_true(is.na(.combined_pct(peri_lost)))
  expect_equal(peri_lost$total, 10)          # the core's cells are still a fact about the islet
  expect_match(peri_lost$reason, "peri-islet expansion total", fixed = TRUE)
  expect_match(peri_lost$reason, "0.5", fixed = TRUE)
  core_lost <- ie_combine_core_peri(0.2, 0.5, NA_real_, 20)
  expect_true(is.na(core_lost$n))
  expect_equal(core_lost$total, 20)
  expect_match(core_lost$reason, "core total", fixed = TRUE)
  expect_match(core_lost$reason, "0.2", fixed = TRUE)
  # A blank in the source column arrives as "" or a non-number and reads the same as NA.
  expect_true(is.na(ie_combine_core_peri(0.2, 0.5, 10, "")$n))
  expect_true(is.na(ie_combine_core_peri(0.2, 0.5, 10, "n/a")$n))
})

test_that("an NA total beside an NA proportion is no cells there, as the branch always has", {
  # "an islet with no peri zone still contributes its core, and vice versa"
  core_only <- ie_combine_core_peri(0.2, NA_real_, 10, NA_real_)
  expect_equal(core_only$n, 2)
  expect_equal(core_only$total, 10)
  expect_equal(.combined_pct(core_only), 20)
  peri_only <- ie_combine_core_peri(NaN, 0.5, NA_real_, 20)
  expect_equal(peri_only$n, 10)
  expect_equal(peri_only$total, 20)
  expect_equal(.combined_pct(peri_only), 50)
  expect_true(is.na(core_only$reason) && is.na(peri_only$reason))
})

test_that("a NaN proportion on a side with NO cells is a measured absence, not unavailable", {
  # build_galectin_layer.py writes NaN exactly where n_cells_qupath is 0; that side has nothing
  # to measure, and the other side stands alone.
  x <- ie_combine_core_peri(core_prop = NaN, peri_prop = 0.5, core_total = 0, peri_total = 20)
  expect_equal(x$n, 10)
  expect_equal(.combined_pct(x), 50)
  expect_true(is.na(x$reason))
  none <- ie_combine_core_peri(NaN, NaN, 0, 0)
  expect_equal(none$n, 0)
  expect_equal(none$total, 0)
  expect_true(is.na(.combined_pct(none)))   # no denominator, as before
  expect_true(is.na(none$reason))
})

# ITEM 2 (decided 2026-09-22). Values that cannot be what the column holds are refused where
# the column is read, once, with a reason; nothing downstream is hardened site by site.

test_that("a negative or infinite total is not a count and reads as NA, then the total rule applies", {
  # -5 beside a finite proportion: refused, then a null total beside a measured proportion.
  neg <- ie_combine_core_peri(0.2, 0.5, -5, 20)
  expect_true(is.na(neg$n))
  expect_equal(neg$total, 20)                # the refused side counts 0 in the denominator
  expect_match(neg$reason, "core total -5 is not a count", fixed = TRUE)
  inf <- ie_combine_core_peri(0.2, 0.5, 10, Inf)
  expect_true(is.na(inf$n))
  expect_equal(inf$total, 10)
  expect_match(inf$reason, "peri-islet expansion total Inf is not a count", fixed = TRUE)
  expect_true(is.na(ie_combine_core_peri(0.2, 0.5, 10, -Inf)$n))
  # The composition the decision spells out: a refused total is a null total, and a null total
  # beside a null proportion is no cells there. The refusal is then not visible in `reason`,
  # which records why n is NA and nothing else.
  neg_and_nan <- ie_combine_core_peri(0.2, NaN, 10, -5)
  expect_equal(neg_and_nan$n, 2)
  expect_equal(neg_and_nan$total, 10)
  expect_true(is.na(neg_and_nan$reason))
})

test_that("a proportion outside [0, 1] is not a proportion and the islet is unavailable", {
  over <- ie_combine_core_peri(1.5, 0.5, 10, 20)
  expect_true(is.na(over$n))
  expect_equal(over$total, 30)               # the denominator is left as measured
  expect_match(over$reason, "core proportion 1.5 is not a proportion", fixed = TRUE)
  under <- ie_combine_core_peri(0.2, -0.2, 10, 20)
  expect_true(is.na(under$n))
  expect_match(under$reason, "peri-islet expansion proportion -0.2 is not a proportion", fixed = TRUE)
  # Unconditional: the refusal is about the value read, not about whether cells were counted
  # beside it -- 1.5 over 0 cells is still not a proportion.
  expect_true(is.na(ie_combine_core_peri(1.5, 0.5, 0, 20)$n))
  expect_true(is.na(ie_combine_core_peri(1.5, 0.5, NA_real_, 20)$n))
  # The bounds themselves are proportions.
  edge <- ie_combine_core_peri(0, 1, 10, 20)
  expect_equal(edge$n, 20)
  expect_true(is.na(edge$reason))
  # One reason per refused column, joined when both sides have one.
  two <- ie_combine_core_peri(1.5, 0.5, 10, -5)
  expect_true(is.na(two$n))
  expect_match(two$reason, "core proportion 1.5 is not a proportion", fixed = TRUE)
  expect_match(two$reason, "peri-islet expansion total -5 is not a count", fixed = TRUE)
})

test_that("it is vectorised over islets and recycles a scalar side", {
  x <- ie_combine_core_peri(core_prop = c(0.2, NaN, 0.4, NaN, 0.2, 1.5),
                            peri_prop = c(0.5, 0.5, NaN, 0.5, 0.5, 0.5),
                            core_total = c(10, 5, 10, 0, 10, 10),
                            peri_total = c(20, 20, 0, 20, NA, 20))
  expect_identical(nrow(x), 6L)
  expect_equal(x$n, c(12, NA, 4, 10, NA, NA))
  expect_equal(x$total, c(30, 25, 10, 20, 10, 30))
  expect_identical(is.na(x$reason), c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  # The module passes a scalar NA_real_ when a side's column is absent from the frame; every
  # islet with cells on that side is then unavailable, not 0.
  s <- ie_combine_core_peri(NA_real_, c(0.5, 0.5), c(5, 0), c(20, 20))
  expect_equal(s$n, c(NA, 10))
  # A missing total column beside a measured proportion column: every islet is unavailable.
  t <- ie_combine_core_peri(0.2, c(0.5, 0.5), c(5, 0), NA_real_)
  expect_true(all(is.na(t$n)))
  expect_identical(nrow(ie_combine_core_peri(numeric(0), numeric(0), numeric(0), numeric(0))), 0L)
})

test_that("the Plot module's combined branch uses the helper and no longer folds NaN to 0", {
  src <- readLines(.find_source("R/mod_plot_server.R"), warn = FALSE)
  expect_true(any(grepl("ie_combine_core_peri(", src, fixed = TRUE)))
  expect_false(any(grepl("core_prop * core_total, 0)", src, fixed = TRUE)))
  expect_false(any(grepl("peri_prop * peri_total, 0)", src, fixed = TRUE)))
})
