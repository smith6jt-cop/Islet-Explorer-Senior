test_that("build_cluster_csv returns one row per cell with res cols", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out    <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10, resolutions = c(0.3, 0.7)
  ))

  df <- build_cluster_csv(out)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), ncol(out))
  expect_true("cell" %in% names(df))
  expect_true(all(c("seurat_clusters_res_0.3",
                    "seurat_clusters_res_0.7") %in% names(df)))
  # round-trip via tempfile to confirm CSV-clean
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(df, tmp, row.names = FALSE)
  back <- utils::read.csv(tmp, stringsAsFactors = FALSE)
  expect_equal(nrow(back), nrow(df))
})

test_that("build_markers_csv concatenates with source / group_col cols", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out    <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10, resolutions = 0.5
  ))
  m <- compute_markers(out, "seurat_clusters_res_0.5")
  cache <- list(`run_a::seurat_clusters_res_0.5` = m,
                `run_b::seurat_clusters_res_0.5` = m)
  df <- build_markers_csv(cache)
  expect_equal(nrow(df), 2 * nrow(m))
  expect_true(all(c("source_run_id", "group_col") %in% names(df)))
  expect_setequal(unique(df$source_run_id), c("run_a", "run_b"))
})

test_that("build_cluster_csv errors without a clustered object", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  expect_error(build_cluster_csv(xen), "No seurat_clusters_res_")
  expect_error(build_cluster_csv(NULL), "No clustered Seurat")
})

test_that("render_session_report writes a non-empty self-contained HTML", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out    <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10, resolutions = 0.5
  ))

  app_state <- list(
    xen              = xen,
    xen_path         = "/tmp/synthetic.rds",
    xen_clustered    = out,
    cluster_stack    = list(list(obj = out, label = "root")),
    markers_cache    = list(),
    markers_last_key = NULL
  )

  tmp <- tempfile(fileext = ".html")
  render_session_report(tmp, panels, app_state)
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 1000)

  body <- readLines(tmp, warn = FALSE)
  body <- paste(body, collapse = "\n")
  expect_match(body, "Audit overview")
  expect_match(body, "Loaded dataset")
  expect_match(body, "Cluster pipeline")
  expect_match(body, "Session")    # session info card
})

test_that("session report includes markers when the cache has entries", {
  panels <- test_load_panels()
  xen    <- make_test_seurat(panels = panels)
  out    <- run_cluster_pipeline(xen, panels, list(
    use_all_features = TRUE, npcs = 10, resolutions = 0.5
  ))
  m <- compute_markers(out, "seurat_clusters_res_0.5")
  app_state <- list(
    xen              = xen,
    xen_path         = "synthetic",
    xen_clustered    = out,
    cluster_stack    = list(list(obj = out, label = "root")),
    markers_cache    = list(`syn::seurat_clusters_res_0.5` = m),
    markers_last_key = "syn::seurat_clusters_res_0.5"
  )
  tmp <- tempfile(fileext = ".html")
  render_session_report(tmp, panels, app_state)
  body <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  expect_match(body, "Top-10 markers per cluster")
})
