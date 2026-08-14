test_that("compute_markers returns the standard presto column set", {
  skip_if_not_installed("presto")
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out    <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10, resolutions = 0.5
  ))

  m <- compute_markers(out, "seurat_clusters_res_0.5")
  expect_s3_class(m, "data.frame")
  expect_setequal(names(m),
                  c("feature", "group", "avgExpr", "logFC",
                    "statistic", "auc", "pval", "padj",
                    "pct_in", "pct_out"))
  # One row per (cluster × gene)
  expect_equal(nrow(m),
               nrow(out) * length(unique(out@meta.data[["seurat_clusters_res_0.5"]])))

  # AUC is in [0, 1]
  expect_true(all(m$auc >= 0 & m$auc <= 1))
})

test_that("top_markers respects the per-group N and the filter cutoffs", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out    <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10, resolutions = 0.5
  ))
  m <- compute_markers(out, "seurat_clusters_res_0.5")

  top <- top_markers(m, n = 5L, min_pct_in = 0, max_padj = 1)
  per <- table(top$group)
  expect_true(all(per <= 5L))

  # Tighter cutoffs return ≤ as many rows
  loose <- top_markers(m, n = 30L, min_pct_in = 0,    max_padj = 1)
  tight <- top_markers(m, n = 30L, min_pct_in = 0.10, max_padj = 0.05)
  expect_lte(nrow(tight), nrow(loose))
})

test_that("compute_markers errors when fewer than two groups exist", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  xen$one_group <- "a"
  expect_error(compute_markers(xen, "one_group"),
               "at least 2 groups")
})
