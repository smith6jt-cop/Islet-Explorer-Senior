test_that("load_panels composes a tissue with the legacy panels shape", {
  skip_if_not(test_tissue_present("pancreas"),
              "no tissue with tissue_id 'pancreas' on disk")

  panels <- test_load_panels("pancreas")

  # Tissue metadata
  expect_equal(panels$tissue$id, "pancreas")
  expect_equal(panels$tissue$display_name, "Human Pancreas")
  expect_setequal(panels$reference_runs, c("0041323", "0041326"))

  # 5K reference is constant (gene + biology cols, no detection_pct).
  expect_equal(nrow(panels$reference_5k), 4992L)
  expect_true("gene"      %in% names(panels$reference_5k))
  expect_true("gene_id"   %in% names(panels$reference_5k))
  expect_true("full_name" %in% names(panels$reference_5k))
  expect_false(any(grepl("^detection_pct_", names(panels$reference_5k))))

  # Subpanels: tissue-local (manifest sets use_shared_subpanels: false); each
  # already carries audit + biology columns on disk, re-asserted by the join.
  expect_gt(length(panels$subpanels), 50)
  for (key in names(panels$subpanels)) {
    df <- panels$subpanels[[key]]
    expect_true("gene" %in% names(df), info = key)
    expect_gt(nrow(df), 0L)
    # Audit columns appear on every subpanel after the load-time join.
    expect_true("detection_pct_0041323" %in% names(df), info = key)
    expect_true("log2_detection_ratio"  %in% names(df), info = key)
  }

  # All pancreas subpanels are materialized under tissues/pancreas/subpanels/
  # (cell-type, pathway, and residual), so both pancreas-specific and the
  # general biology subpanels resolve there.
  expect_true("01_pancreas_endocrine" %in% names(panels$subpanels))
  expect_true("03_immune_T_cell" %in% names(panels$subpanels))

  # Custom panel
  expect_equal(nrow(panels$custom), 100L)
  expect_true("detection_pct_0041323" %in% names(panels$custom))

  # xenium5k (joined audit) preserves the legacy reference-table use.
  expect_gte(nrow(panels$xenium5k), 4992L)
  expect_true("gene" %in% names(panels$xenium5k))

  # Excluded list
  expect_equal(nrow(panels$excluded), 106L)

  # hIO supplementary tables
  expect_setequal(
    names(panels$hIO),
    c("hIO_genes_vs_5K_status",
      "hIO_genes_unique_vs_5K",
      "hIO_genes_gained_vs_current_panel",
      "hIO_vs_5K_subpanel_coverage",
      "hIO_vs_5K_immune_lost_high_specificity",
      "hIO_vs_5K_immune_losses_detailed")
  )
})

test_that("available_tissues returns every directory with a manifest", {
  tids <- available_tissues(test_tissues_root())
  expect_true("pancreas" %in% tids)
  # thymus stub is present in this checkout but exposes no subpanels yet.
  expect_true("thymus" %in% tids)
})

test_that("load_panels rejects unknown tissues", {
  expect_error(test_load_panels("not_a_real_tissue"),
               regexp = "not found",
               ignore.case = TRUE)
})

test_that("load_panels handles a tissue with no custom panel", {
  td <- make_test_tissue(has_custom = FALSE)
  on.exit(unlink(td$root, recursive = TRUE), add = TRUE)

  panels <- load_panels(
    tissue_id              = td$tissue,
    reference_5k_path      = file.path(td$root, "reference_5k"),
    subpanels_shared_path  = file.path(td$root, "subpanels_shared"),
    tissues_root           = file.path(td$root, "tissues")
  )
  expect_null(panels$custom)
  expect_length(panels$reference_runs, 2L)
  expect_true("01_shared_demo" %in% names(panels$subpanels))
  expect_true("02_tissue_demo" %in% names(panels$subpanels))
})

test_that("load_panels surfaces tissue-specific overrides over shared ones", {
  # Build a tissue whose subpanels/ contains a file with the same key as
  # the shared layer; tissue wins on collision.
  td <- make_test_tissue()
  override <- data.frame(gene = c("DDD"), category = "tissue_override",
                         rationale = "test",
                         single_gene_panel_note = NA_character_,
                         stringsAsFactors = FALSE)
  data.table::fwrite(override,
                    file.path(td$tissue_dir, "subpanels",
                              "01_shared_demo.csv"), na = "")
  on.exit(unlink(td$root, recursive = TRUE), add = TRUE)

  panels <- load_panels(
    tissue_id              = td$tissue,
    reference_5k_path      = file.path(td$root, "reference_5k"),
    subpanels_shared_path  = file.path(td$root, "subpanels_shared"),
    tissues_root           = file.path(td$root, "tissues")
  )
  expect_equal(panels$subpanels[["01_shared_demo"]]$category[1],
               "tissue_override")
})

test_that("resolve_subpanel_key maps long summary names to file stems", {
  skip_if_not(test_tissue_present("pancreas"))
  panels <- test_load_panels("pancreas")

  expect_equal(
    resolve_subpanel_key("99d_truly_unannotated_subset_of_99c", panels),
    "99d_truly_unannotated"
  )
  expect_equal(
    resolve_subpanel_key("01_pancreas_endocrine", panels),
    "01_pancreas_endocrine"
  )
  expect_true(is.na(resolve_subpanel_key("not_a_real_subpanel", panels)))
})

test_that("manifest tissue_id is authoritative — folder name is arbitrary", {
  # Build a tissue at <root>/tissues/some_arbitrary_folder/ whose
  # manifest declares tissue_id: real_id. Loader should look it up by
  # tissue_id, not folder name.
  td <- make_test_tissue(tissue = "real_id")
  on.exit(unlink(td$root, recursive = TRUE), add = TRUE)

  weird_path <- file.path(td$root, "tissues", "some_arbitrary_folder")
  file.rename(td$tissue_dir, weird_path)

  ids <- available_tissues(file.path(td$root, "tissues"))
  expect_equal(ids, "real_id")

  panels <- load_panels(
    tissue_id              = "real_id",
    reference_5k_path      = file.path(td$root, "reference_5k"),
    subpanels_shared_path  = file.path(td$root, "subpanels_shared"),
    tissues_root           = file.path(td$root, "tissues")
  )
  expect_equal(panels$tissue$id, "real_id")
  expect_equal(basename(panels$tissue$dir), "some_arbitrary_folder")
})

test_that("tissues_index warns on duplicate tissue_ids and keeps the first", {
  td1 <- make_test_tissue(tissue = "dup")
  td2 <- make_test_tissue(tissue = "dup",
                          root = tempfile("xen_tissue_dup_"))
  on.exit(unlink(c(td1$root, td2$root), recursive = TRUE), add = TRUE)

  # Move td2's tissue dir alongside td1's so they share the same root.
  alt <- file.path(td1$root, "tissues", "another_dup")
  file.rename(td2$tissue_dir, alt)

  expect_warning(idx <- tissues_index(file.path(td1$root, "tissues")),
                 "duplicate tissue_id")
  expect_length(idx, 1L)
  expect_equal(names(idx), "dup")
})

test_that("audit_detection_cols returns the column-names contract", {
  expect_equal(audit_detection_cols(c("A", "B")),
               c("detection_pct_A", "detection_pct_B",
                 "n_cells_A", "n_cells_B",
                 "log2_detection_ratio"))
  expect_equal(audit_detection_cols(character()), character())
})
