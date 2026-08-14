test_that("compute_subpanel_coverage returns one row per subpanel + ancillary", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  cov <- compute_subpanel_coverage(panels, xen, min_detection_pct = 0)

  # Every subpanel + custom + xenium5k row
  expect_equal(nrow(cov),
               length(panels$subpanels) + 2L)
  expect_named(cov, c("subpanel", "n_genes", "n_present", "n_passing",
                      "frac_present", "frac_passing",
                      "mean_detection_pct", "median_detection_pct",
                      "mean_expr", "median_expr"))

  # Synthetic genes are drawn from subpanels ⊂ xenium5k, so xenium5k_in_audit
  # row should have n_present == n_data ≤ 4992.
  xen5k <- cov[cov$subpanel == "xenium5k_in_audit", ]
  expect_equal(xen5k$n_genes, 4992L)
  expect_lte(xen5k$n_present, 4992L)
  expect_gt(xen5k$n_present, 0L)

  # frac_present is in [0,1]
  expect_true(all(cov$frac_present >= 0 & cov$frac_present <= 1, na.rm = TRUE))
  expect_true(all(cov$frac_passing >= 0 & cov$frac_passing <= 1, na.rm = TRUE))
})

test_that("min_detection_pct cutoff monotonically reduces n_passing", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  c0  <- compute_subpanel_coverage(panels, xen, min_detection_pct = 0)
  c50 <- compute_subpanel_coverage(panels, xen, min_detection_pct = 50)

  expect_equal(c0$subpanel, c50$subpanel)
  expect_true(all(c50$n_passing <= c0$n_passing))
  expect_true(all(c50$n_present == c0$n_present))   # presence is cutoff-free
})

test_that("top_n_by_subpanel returns a long data.frame respecting N", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  df <- top_n_by_subpanel(panels, xen, n = 5L, min_detection_pct = 0)
  expect_named(df, c("subpanel", "gene", "detection_pct", "mean_expr"))

  # Each subpanel contributes at most N rows
  per <- table(df$subpanel)
  expect_true(all(per <= 5L))
})

test_that("per_gene_detection_pct is a named numeric vector summing reasonably", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  d <- per_gene_detection_pct(xen)
  expect_type(d, "double")
  expect_named(d)
  expect_equal(length(d), nrow(xen))
  expect_true(all(d >= 0 & d <= 100))
})
