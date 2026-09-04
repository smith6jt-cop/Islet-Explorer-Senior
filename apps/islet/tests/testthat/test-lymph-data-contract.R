test_that("the lymph lane is configuration plus synthetic data, not an app fork", {
  lymph_root <- dirname(dirname(.find_repo_file("apps/lymph/config/lymphoid_default.yml")))
  data_root <- file.path(lymph_root, "data", "app_data")

  expect_false(dir.exists(file.path(lymph_root, "R")))
  expect_false(file.exists(file.path(lymph_root, "app.R")))
  expect_true(file.exists(file.path(data_root, "master_results.xlsx")))
  expect_true(file.exists(file.path(data_root, "follicle_spatial_lookup.csv")))
  expect_true(file.exists(file.path(data_root, "phenotype_rules.csv")))
  expect_length(list.files(file.path(data_root, "cells"), pattern = "[.]csv$"), 24L)
  expect_length(list.files(file.path(data_root, "json"), pattern = "[.]geojson$"), 4L)
})

test_that("the committed lymph example satisfies tab-facing schemas", {
  data_root <- dirname(.find_repo_file("apps/lymph/data/app_data/phenotype_rules.csv"))
  lookup <- read.csv(file.path(data_root, "follicle_spatial_lookup.csv"), check.names = FALSE)
  expect_true(all(c(
    "case_id", "follicle_key", "centroid_x_um", "centroid_y_um", "area_um2"
  ) %in% names(lookup)))
  expect_false(anyDuplicated(paste(lookup$case_id, lookup$follicle_key, sep = "|")) > 0L)

  cell_files <- list.files(file.path(data_root, "cells"), pattern = "[.]csv$", full.names = TRUE)
  cell_columns <- names(read.csv(cell_files[[1L]], nrows = 1L, check.names = FALSE))
  expect_true(all(c(
    "X_centroid", "Y_centroid", "phenotype", "cell_region", "Cell Area", "Nucleus Area"
  ) %in% cell_columns))

  skip_if_not_installed("readxl")
  workbook <- file.path(data_root, "master_results.xlsx")
  expect_identical(
    readxl::excel_sheets(workbook),
    c("Follicle_Targets", "Follicle_Markers", "Follicle_Composition", "LGALS3")
  )
  targets <- readxl::read_excel(workbook, sheet = "Follicle_Targets", n_max = 1L)
  expect_true(all(c(
    "Case ID", "Donor Status", "region", "type", "class", "area_um2", "region_um2"
  ) %in% names(targets)))
})
