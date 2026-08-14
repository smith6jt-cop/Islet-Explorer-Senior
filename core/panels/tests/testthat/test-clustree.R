test_that("clustree() returns a ggplot from a fixture with res cols", {
  skip_if_not_installed("clustree")
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  out <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10,
    resolutions = c(0.2, 0.5, 0.8)
  ))

  g <- clustree::clustree(out, prefix = "seurat_clusters_res_")
  expect_s3_class(g, "ggplot")
})

test_that("compute_stability_summary has the expected shape and ranges", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10,
    resolutions = c(0.2, 0.5, 0.8)
  ))

  st <- compute_stability_summary(out)
  expect_named(st, c("resolution", "n_clusters",
                     "ARI_vs_prev", "NMI_vs_prev",
                     "frac_changed_vs_prev"))
  expect_equal(nrow(st), 3L)

  # First row has no "previous" → all *_vs_prev are NA
  expect_true(is.na(st$ARI_vs_prev[1]))
  expect_true(is.na(st$NMI_vs_prev[1]) || requireNamespace("aricode",
                                                            quietly = TRUE))
  # n_clusters strictly positive
  expect_true(all(st$n_clusters > 0))
  # frac_changed_vs_prev in [0,1] for non-NA rows
  fch <- st$frac_changed_vs_prev[!is.na(st$frac_changed_vs_prev)]
  expect_true(all(fch >= 0 & fch <= 1))
})

test_that("clustree_edge_table parses edges between adjacent resolutions", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10,
    resolutions = c(0.2, 0.5, 0.8)
  ))

  edges <- clustree_edge_table(out)
  expect_named(edges, c("from_resolution", "from",
                        "to_resolution", "to",
                        "Freq", "in_prop"))
  expect_gt(nrow(edges), 0L)
  expect_true(all(edges$Freq > 0))
  expect_true(all(edges$in_prop > 0 & edges$in_prop <= 1 + 1e-9))
})

test_that(".adjusted_rand returns 1 for identical labellings, 0 for random-equivalent", {
  expect_equal(.adjusted_rand(c(1,1,2,2), c(1,1,2,2)), 1)
  expect_lt(.adjusted_rand(c(1,2,3,4,5,6), c(1,1,2,2,3,3)), 1)
})
