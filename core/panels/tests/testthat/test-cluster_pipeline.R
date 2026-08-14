test_that("run_cluster_pipeline produces ≥2 cluster columns and a UMAP", {
  skip_if_not_installed("Seurat")
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  out <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE,
    npcs        = 15,
    resolutions = c(0.2, 0.5, 0.8),
    seed        = 42L
  ))

  expect_s4_class(out, "Seurat")
  expect_setequal(SeuratObject::Reductions(out), c("pca", "umap"))

  res_cols <- cluster_resolution_columns(out)
  expect_gte(length(res_cols), 2L)
  expect_setequal(res_cols, c("seurat_clusters_res_0.2",
                              "seurat_clusters_res_0.5",
                              "seurat_clusters_res_0.8"))
  for (col in res_cols) {
    expect_gt(length(unique(out@meta.data[[col]])), 0L)
  }

  # The intermediate column is dropped
  expect_false("seurat_clusters" %in% names(out@meta.data))

  # Run log is recorded
  expect_true("pipeline_history" %in% names(out@misc))
  expect_equal(length(out@misc$pipeline_history), 1L)
  expect_true(nchar(out@misc$last_run_id) > 0)
})

test_that("output object survives qs2::qs_save / qs_read round trip", {
  skip_if_not_installed("qs2")
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  out <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE,
    npcs = 10, resolutions = c(0.3, 0.7)
  ))
  tmp <- tempfile(fileext = ".qs2")
  qs2::qs_save(out, tmp)
  back <- qs2::qs_read(tmp)
  expect_s4_class(back, "Seurat")
  expect_setequal(cluster_resolution_columns(back),
                  cluster_resolution_columns(out))
  expect_equal(SeuratObject::Embeddings(back, "umap"),
               SeuratObject::Embeddings(out,  "umap"))
})

test_that("subpanel restriction selects only matching genes", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)

  out <- run_cluster_pipeline(xen, panels, list(
    subpanels   = "01_pancreas_endocrine",
    npcs        = 5,
    resolutions = 0.5
  ))
  # Variable features must be ⊂ (subpanel ∩ data)
  vf <- SeuratObject::VariableFeatures(out)
  sp_genes <- panels$subpanels[["01_pancreas_endocrine"]]$gene
  expect_true(all(vf %in% sp_genes))
})

test_that("nCount filter trims cells", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  hi <- max(xen$nCount_Xenium)
  cutoff <- as.integer(stats::median(xen$nCount_Xenium))

  out <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE,
    nCount_min = cutoff,
    npcs = 5, resolutions = 0.5
  ))
  expect_lt(ncol(out), ncol(xen))
  expect_true(all(out$nCount_Xenium >= cutoff))
})

test_that("Leiden falls back to Louvain when leidenAlg is missing", {
  skip_if(requireNamespace("leidenAlg", quietly = TRUE),
          "leidenAlg installed; fallback path not exercised")
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  expect_warning(
    run_cluster_pipeline(xen, panels, list(
      use_all_features  = TRUE,
      cluster_algorithm = "Leiden",
      npcs = 5, resolutions = 0.5
    )),
    "leidenAlg not installed"
  )
})
