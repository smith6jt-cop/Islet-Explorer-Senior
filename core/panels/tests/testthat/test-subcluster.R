test_that("subcluster pipeline produces a Seurat ⊂ parent with res cols", {
  skip_if_not_installed("Seurat")
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  parent <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10,
    resolutions = c(0.2, 0.5, 0.8)
  ))

  cluster_col <- "seurat_clusters_res_0.5"
  cl_levels   <- as.character(unique(parent@meta.data[[cluster_col]]))
  pick_cluster <- cl_levels[1]
  keep <- which(as.character(parent@meta.data[[cluster_col]]) == pick_cluster)
  skip_if(length(keep) < 30L,
          "First cluster has < 30 cells; skip test on this synthetic seed")

  sub_in <- parent[, keep]
  for (c0 in cluster_resolution_columns(sub_in)) {
    sub_in@meta.data[[c0]] <- NULL
  }

  child <- run_cluster_pipeline(sub_in, panels, list(
    use_all_features = TRUE, npcs = 5,
    resolutions = c(0.2, 0.5)
  ))

  expect_s4_class(child, "Seurat")
  # Cells of child ⊂ cells of parent
  expect_true(all(colnames(child) %in% colnames(parent)))
  expect_lt(ncol(child), ncol(parent))
  # Res cols start fresh
  expect_setequal(cluster_resolution_columns(child),
                  c("seurat_clusters_res_0.2", "seurat_clusters_res_0.5"))
  # Pipeline history includes one new run on the child
  expect_true("pipeline_history" %in% names(child@misc))
  expect_gte(length(child@misc$pipeline_history), 1L)
})
