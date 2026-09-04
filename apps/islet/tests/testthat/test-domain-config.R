test_that("islet DOMAIN preserves the established display contract", {
  domain <- ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml"))

  expect_identical(domain_app_name(domain), "Islet Explorer Senior")
  expect_equal(domain_pixel_size(domain), 0.5078)
  expect_identical(domain_group_levels(domain), c("ND", "Aab+", "T1D"))
  expect_identical(domain_key_col(domain), "islet_key")
  expect_identical(domain_id_regex(domain), "Islet_\\d+")
  expect_length(domain_phenotype_order(domain), 34L)
  expect_setequal(names(domain_phenotype_colors(domain)), domain_phenotype_order(domain))
  expect_true(all(vapply(
    domain_group_palettes(domain),
    function(palette) all(domain_group_levels(domain) %in% names(palette)),
    logical(1)
  )))
})

test_that("lymph DOMAIN is independently loadable", {
  domain <- ie_load_domain(.find_repo_file("apps/lymph/config/lymphoid_default.yml"))

  expect_identical(domain_app_name(domain), "LymphVizkit")
  expect_equal(domain_pixel_size(domain), 0.5078)
  expect_identical(domain_group_levels(domain), c("Resting", "Reactive", "Involuted"))
  expect_identical(domain_key_col(domain), "follicle_key")
  expect_identical(domain_id_regex(domain), "Follicle_\\d+")
  expect_identical(domain_pseudotime_root(domain), "BCL6")
  expect_length(domain_phenotype_order(domain), 19L)
  expect_setequal(names(domain_phenotype_colors(domain)), domain_phenotype_order(domain))
})

test_that("every selectable islet palette covers the DOMAIN vocabulary", {
  env <- new.env(parent = globalenv())
  sys.source(.find_repo_file("apps/islet/R/drilldown_helpers.R"), envir = env)
  expected <- domain_phenotype_order()
  expect_true(all(vapply(
    env$PHENOTYPE_PALETTES,
    function(palette) all(expected %in% names(palette)),
    logical(1)
  )))
})

test_that("invalid or ambiguous DOMAIN values fail closed", {
  domain <- ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml"))
  domain$grouping$levels <- c("ND", "ND", "T1D")
  expect_error(ie_domain_validate(domain), "non-empty and unique", fixed = TRUE)

  domain <- ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml"))
  domain$phenotypes$colors[[domain$phenotypes$display_order[[1L]]]] <- NULL
  expect_error(ie_domain_validate(domain), "does not cover", fixed = TRUE)
})

test_that("unit-key extraction and ordinal coding follow active DOMAIN", {
  islet <- ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml"))
  lymph <- ie_load_domain(.find_repo_file("apps/lymph/config/lymphoid_default.yml"))
  on.exit(ie_set_domain(islet), add = TRUE)

  ie_set_domain(islet)
  got <- add_unit_key(data.frame(region = c("Islet_7_core", NA_character_)))
  expect_identical(got$islet_key, c("Islet_7", NA_character_))
  expect_equal(code_group_status(c("ND", "Aab+", "T1D", "missing")), c(0, 1, 2, NA))

  ie_set_domain(lymph)
  got <- add_unit_key(data.frame(name = c("Follicle_12", "region 9")))
  expect_identical(got$follicle_key, c("Follicle_12", "Follicle_9"))
  expect_equal(
    code_group_status(c("Resting", "Reactive", "Involuted", "missing")),
    c(0, 1, 2, NA)
  )
})
