test_that(".rds ingest returns a valid Seurat object and caches", {
  skip_if_not_installed("Seurat")
  cache <- tempfile("cache_rds_")
  rds   <- make_test_rds()

  x1 <- load_xenium(rds, cache_dir = cache)
  expect_s4_class(x1, "Seurat")
  expect_equal(ncol(x1), 800L)
  expect_equal(nrow(x1), 500L)
  cache_info <- attr(x1, "load_xenium_cache")
  expect_false(isTRUE(cache_info$hit))
  expect_true(file.exists(cache_info$file))

  # Second call must be a cache hit
  x2 <- load_xenium(rds, cache_dir = cache)
  expect_true(isTRUE(attr(x2, "load_xenium_cache")$hit))
  expect_equal(ncol(x2), ncol(x1))
})

test_that("Xenium-bundle ingest reads h5 + parquet into a Seurat object", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("arrow")
  bdir  <- make_test_bundle()
  cache <- tempfile("cache_bundle_")

  x <- load_xenium(bdir, cache_dir = cache)
  expect_s4_class(x, "Seurat")
  expect_equal(ncol(x), 800L)
  expect_equal(nrow(x), 500L)
  expect_true("x_centroid" %in% names(x@meta.data))
  expect_true("y_centroid" %in% names(x@meta.data))
  expect_equal(x@misc$xenium_bundle_dir, bdir)

  # Cache hit on a second call
  x2 <- load_xenium(bdir, cache_dir = cache)
  expect_true(isTRUE(attr(x2, "load_xenium_cache")$hit))
})

test_that("panel_validate returns the expected report shape", {
  panels <- test_load_panels()
  obj    <- make_test_seurat(panels = panels)
  rep    <- panel_validate(rownames(obj), panels)

  expect_named(rep, c("n_data", "n_reference", "intersection",
                      "missing_from_data", "extra_in_data",
                      "n_intersection", "n_missing_from_data",
                      "n_extra_in_data", "pct_reference_covered",
                      "partial", "has_custom"))
  expect_true(rep$has_custom)
  expect_equal(rep$n_data, 500L)
  expect_equal(rep$n_reference, 5092L)
  expect_equal(rep$n_extra_in_data, 0L)         # synthetic genes ⊂ subpanels ⊂ 5K
  expect_lt(rep$pct_reference_covered, 1)        # 500 << 5092
  expect_true(rep$partial)
  expect_match(panel_validate_summary(rep), "Loaded dataset:")
})

test_that("h5ad ingest dispatches to .load_xenium_h5ad", {
  skip_if_not_installed("zellkonverter")
  tmp <- tempfile(fileext = ".h5ad")
  file.create(tmp)
  # Empty file => zellkonverter/rhdf5 can't open it. We only care that
  # the dispatcher routes to the h5ad path and surfaces an error.
  expect_error(load_xenium(tmp))
})

test_that("unknown extensions are rejected", {
  tmp <- tempfile(fileext = ".bin")
  file.create(tmp)
  expect_error(load_xenium(tmp), "unsupported file extension")
})
