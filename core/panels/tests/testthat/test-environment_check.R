test_that("check_environment returns the expected shape on a healthy host", {
  env <- check_environment()
  expect_named(env, c("errors", "warnings"))
  expect_type(env$errors,   "character")
  expect_type(env$warnings, "character")
})

test_that("check_environment flags a missing required package", {
  env <- check_environment(required = c("shiny", "definitely_not_a_package"))
  expect_match(paste(env$warnings, collapse = " "),
               "definitely_not_a_package")
})

test_that("check_environment flags an unrealistically high R minimum", {
  env <- check_environment(required = "shiny", min_r_version = "99.0.0")
  expect_match(paste(env$warnings, collapse = " "), "R \\d")
})
