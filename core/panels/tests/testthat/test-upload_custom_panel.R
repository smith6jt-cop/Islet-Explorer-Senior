tissue_root_present <- function() test_tissue_present("pancreas")

test_that("custom_panel_canonical_cols composes columns from reference_runs", {
  cols <- custom_panel_canonical_cols(c("0041323", "0041326"))
  expect_length(cols, 10L)
  expect_true("gene" %in% cols)
  expect_true("detection_pct_0041323" %in% cols)
  expect_true("log2_detection_ratio"  %in% cols)
  expect_false("log2_detection_ratio_326_over_323" %in% cols)

  # No reference_runs -> only the biology columns.
  expect_equal(custom_panel_canonical_cols(character()),
               c("gene", "category", "rationale", "single_gene_panel_note",
                 "exclude_recommended"))
})

test_that("upload_custom_panel enriches a gene-only CSV from references", {
  skip_if_not(tissue_root_present(), "tissues/pancreas/ missing")
  panels <- test_load_panels("pancreas")

  # One gene known to be in panels$custom (T1D-GWAS), one in 5K only,
  # one made-up unmatched gene.
  g_custom <- panels$custom$gene[1]
  g_5k     <- setdiff(panels$xenium5k$gene, panels$custom$gene)[1]
  g_none   <- "DEFINITELY_NOT_A_REAL_GENE_XYZZY"

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(
    data.frame(gene = c(g_custom, g_5k, g_none),
               stringsAsFactors = FALSE),
    csv, row.names = FALSE
  )

  res <- upload_custom_panel(csv, panels)
  expect_equal(res$n_genes, 3L)
  expect_equal(res$source_name, basename(csv))
  expect_equal(res$n_enriched_from_custom, 1L)
  expect_equal(res$n_enriched_from_5k, 1L)
  expect_equal(res$n_unmatched, 1L)

  # Canonical columns present, in canonical order at the front.
  expect_equal(head(names(res$df), 10L),
               custom_panel_canonical_cols(panels$reference_runs))

  row_custom <- res$df[res$df$gene == g_custom, ]
  row_5k     <- res$df[res$df$gene == g_5k, ]
  row_none   <- res$df[res$df$gene == g_none, ]
  expect_false(is.na(row_custom$detection_pct_0041323))
  expect_false(is.na(row_5k$detection_pct_0041323))
  expect_true(is.na(row_none$detection_pct_0041323))
})

test_that("upload_custom_panel preserves user-supplied non-NA values", {
  skip_if_not(tissue_root_present(), "tissues/pancreas/ missing")
  panels <- test_load_panels("pancreas")

  g <- panels$custom$gene[1]
  user_category <- "USER_OVERRIDE_CATEGORY"

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(
    data.frame(gene = g, category = user_category,
               stringsAsFactors = FALSE),
    csv, row.names = FALSE
  )

  res <- upload_custom_panel(csv, panels)
  expect_equal(res$df$category[1], user_category)
  expect_false(is.na(res$df$detection_pct_0041323[1]))
})

test_that("upload_custom_panel accepts case-insensitive gene column", {
  skip_if_not(tissue_root_present(), "tissues/pancreas/ missing")
  panels <- test_load_panels("pancreas")

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(
    data.frame(Gene = panels$custom$gene[1:3],
               stringsAsFactors = FALSE),
    csv, row.names = FALSE
  )

  res <- upload_custom_panel(csv, panels)
  expect_true("gene" %in% names(res$df))
  expect_equal(res$n_genes, 3L)
})

test_that("upload_custom_panel deduplicates and drops empty gene rows", {
  skip_if_not(tissue_root_present(), "tissues/pancreas/ missing")
  panels <- test_load_panels("pancreas")

  g <- panels$custom$gene[1]
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  writeLines(
    c("gene,extra",
      paste0(g, ",1"),
      paste0(g, ",2"),
      ",3",
      " ,4"),
    csv
  )

  res <- suppressWarnings(upload_custom_panel(csv, panels))
  expect_equal(res$n_genes, 1L)
  expect_equal(res$df$gene[1], g)
})

test_that("upload_custom_panel errors when `gene` column is missing", {
  skip_if_not(tissue_root_present(), "tissues/pancreas/ missing")
  panels <- test_load_panels("pancreas")

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(
    data.frame(symbol = c("AAA", "BBB"), stringsAsFactors = FALSE),
    csv, row.names = FALSE
  )

  expect_error(upload_custom_panel(csv, panels),
               regexp = "gene", ignore.case = TRUE)
})

test_that("upload_custom_panel errors on a file with no valid gene rows", {
  skip_if_not(tissue_root_present(), "tissues/pancreas/ missing")
  panels <- test_load_panels("pancreas")

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  writeLines(c("gene,extra", "NA,1", "NA,2"), csv)

  expect_error(
    suppressWarnings(upload_custom_panel(csv, panels)),
    regexp = "empty|valid|gene", ignore.case = TRUE
  )
})

test_that("upload_custom_panel works for tissues without a default custom panel", {
  td <- make_test_tissue(has_custom = FALSE)
  on.exit(unlink(td$root, recursive = TRUE), add = TRUE)

  panels <- load_panels(
    tissue_id              = td$tissue,
    reference_5k_path      = file.path(td$root, "reference_5k"),
    subpanels_shared_path  = file.path(td$root, "subpanels_shared"),
    tissues_root           = file.path(td$root, "tissues")
  )
  expect_null(panels$custom)

  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(
    data.frame(gene = c("AAA", "BBB", "NOT_IN_5K"),
               stringsAsFactors = FALSE),
    csv, row.names = FALSE
  )

  res <- upload_custom_panel(csv, panels)
  expect_equal(res$n_genes, 3L)
  # No default panel means n_enriched_from_custom is zero by construction.
  expect_equal(res$n_enriched_from_custom, 0L)
  # The 5K audit still provides enrichment for AAA / BBB; NOT_IN_5K is unmatched.
  expect_equal(res$n_enriched_from_5k, 2L)
  expect_equal(res$n_unmatched, 1L)
})

test_that("custom_panel_label returns stem from upload, manifest, or fallback", {
  # No upload, no panels arg → generic fallback.
  expect_equal(custom_panel_label(NULL), "custom_panel")
  expect_equal(custom_panel_label(list()), "custom_panel")

  # Uploaded source overrides everything.
  st <- list(source_name = "my_t1d_v2.csv")
  expect_equal(custom_panel_label(st), "my_t1d_v2")
  expect_equal(custom_panel_label(list(source_name = "panel_alt")),
               "panel_alt")

  # Manifest display_name wins when no upload is active.
  fake_panels <- list(tissue = list(manifest = list(
    custom_panel = list(display_name = "T1D-GWAS custom 100"))))
  expect_equal(custom_panel_label(NULL, fake_panels),
               "T1D-GWAS custom 100")
  expect_equal(custom_panel_label(st, fake_panels), "my_t1d_v2")
})

test_that("panel_display_label rewrites only the custom slot key", {
  st <- list(source_name = "my_t1d_v2.csv")
  expect_equal(panel_display_label("custom_panel", st),
               "my_t1d_v2")
  expect_equal(panel_display_label("custom_panel", NULL),
               "custom_panel")
  expect_equal(panel_display_label("01_pancreas_endocrine", st),
               "01_pancreas_endocrine")
  expect_equal(panel_display_label("xenium5k_in_audit", st),
               "xenium5k_in_audit")
})
