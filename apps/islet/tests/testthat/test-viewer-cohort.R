# Like scripts/test_viewer_helpers.R, load function definitions without running
# startup code, registering resource paths, or starting a TissUUmaps service.
.viewer_helpers_path <- .find_source("R/viewer_helpers.R")
.viewer_server_path <- .find_source("R/mod_viewer_server.R")

.viewer_helpers <- function() {
  env <- new.env(parent = globalenv())
  for (expr in parse(.viewer_helpers_path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, envir = env)
    }
  }
  env$project_root <- NULL
  env
}

.viewer_cohort <- function(ids = c("112", "7001", "7002")) {
  list(donor_ids = ids, availability = ie_avail("measured", n = length(ids)))
}

.viewer_files <- function(root) {
  dirs <- list(www = file.path(root, "www", "local_images"),
               local = file.path(root, "local"), projects = file.path(root, "projects"))
  for (dir in dirs) dir.create(dir, recursive = TRUE)
  file.create(file.path(dirs$www, c("0112.ome.tiff", "7002.OME.TIF", "7003.ome.tiff",
                                   "8112.ome.tiff", "sample_112.ome.tiff", "notes.txt")))
  file.create(file.path(dirs$local, c("0112.ome.tiff", "7001.ome.tif", "7003.ome.tiff")))
  file.create(file.path(dirs$projects, c("0112.tmap", "7001.TMAP", "7003.tmap", "notes.txt")))
  dirs
}

.viewer_server_env <- function(backend, dirs, cohort = .viewer_cohort()) {
  env <- .viewer_helpers()
  for (name in c("moduleServer", "reactive", "reactiveVal", "observe", "observeEvent",
                 "renderUI", "isolate", "tagList", "tags", "fluidRow", "column",
                 "selectInput", "HTML")) {
    env[[name]] <- getExportedValue("shiny", name)
  }
  env$VIEWER_BACKEND <- backend
  env$VIEWER_DEBUG_ENABLED <- FALSE
  env$channel_names_vec <- NULL
  env$default_image_url <- NULL
  env$local_images_root <- dirs$local
  env$read_viewer_cohort <- function() cohort
  env$resolve_avivator_base <- function() "avivator/index.html"
  env$resolve_tissuumaps_project_dir <- function() dirs$projects
  env$ensure_tissuumaps <- function() TRUE
  sys.source(.viewer_server_path, envir = env)
  env
}

test_that("cohort filtering preserves basenames and full paths and matches exact stems", {
  env <- .viewer_helpers()
  images <- c("0112.ome.tiff", "7001.OME.TIF", "7002.tif", "7003.ome.tiff",
              "8112.ome.tiff", "sample_112.ome.tiff", "112_extra.ome.tiff", NA_character_)
  for (prefix in c("", "/images/7003/")) {
    paths <- if (nzchar(prefix)) paste0(prefix, images) else images
    expect_identical(env$filter_viewer_images(paths, .viewer_cohort()), paths[1:3])
  }
  projects <- c("0112", "7001", "7003", "sample_112", "8112")
  expect_identical(env$filter_viewer_images(projects, .viewer_cohort()), projects[1:2])
  expect_identical(env$filter_viewer_images(paste0(projects, ".TMAP"), .viewer_cohort()),
                   paste0(projects[1:2], ".TMAP"))
  expect_identical(env$normalize_viewer_donor_id(c(" 0112 ", "000", "7001", "", NA)),
                   c("112", "0", "7001", NA_character_, NA_character_))
  expect_identical(env$filter_viewer_images(character(0), .viewer_cohort()), character(0))
})

test_that("case lookup filters both roots and does not use substring or regex matches", {
  root <- withr::local_tempdir()
  dirs <- .viewer_files(root)
  env <- .viewer_helpers()
  env$read_viewer_cohort <- .viewer_cohort
  expect_identical(env$find_image_for_case("112", dirs$www, dirs$local),
                   file.path(dirs$www, "0112.ome.tiff"))
  expect_identical(env$find_image_for_case("07001", dirs$www, dirs$local),
                   file.path(dirs$local, "7001.ome.tif"))
  expect_identical(env$find_image_for_case("112", NULL, dirs$local),
                   file.path(dirs$local, "0112.ome.tiff"))
  expect_null(env$find_image_for_case("7003", dirs$www, dirs$local))
  expect_null(env$find_image_for_case("7003", NULL, dirs$local))
  expect_null(env$find_image_for_case(".*", dirs$www, dirs$local))
  expect_null(env$find_image_for_case("12", dirs$www, dirs$local))
})

test_that("TissUUmaps filters project files and the OME-TIFF fallback", {
  root <- withr::local_tempdir()
  dirs <- .viewer_files(root)
  env <- .viewer_helpers()
  env$read_viewer_cohort <- .viewer_cohort
  env$resolve_tissuumaps_project_dir <- function() dirs$projects
  withr::local_dir(root)
  expect_identical(env$list_available_tmap_projects(), c("0112", "7001"))
  env$read_viewer_cohort <- function() .viewer_cohort("7002")
  expect_identical(env$list_available_tmap_projects(), character(0))
  inventory <- env$tmap_project_inventory()
  expect_identical(inventory$n_found, 3L)
  expect_identical(inventory$n_excluded, 3L)
  expect_s3_class(inventory$note, "ie_guard_note")
  expect_match(format(inventory$note), "[viewer_cohort_assets_excluded]", fixed = TRUE)
  expect_match(format(inventory$note), "Assets found: 3; excluded by cohort membership: 3", fixed = TRUE)
  env$read_viewer_cohort <- .viewer_cohort
  env$resolve_tissuumaps_project_dir <- function() NULL
  expect_identical(env$list_available_tmap_projects(), c("0112", "7002"))
  env$read_viewer_cohort <- function() .viewer_cohort("9000")
  expect_identical(env$list_available_tmap_projects(), character(0))
  inventory <- env$tmap_project_inventory()
  expect_identical(inventory$n_found, 5L)
  expect_identical(inventory$n_excluded, 5L)
  expect_s3_class(inventory$note, "ie_guard_note")
  expect_match(format(inventory$note), "Assets found: 5; excluded by cohort membership: 5", fixed = TRUE)
})

test_that("both viewer backends reject excluded dropdown and forced selections", {
  skip_if_not_installed("shiny")
  root <- withr::local_tempdir()
  dirs <- .viewer_files(root)
  withr::local_envvar(c(LOCAL_IMAGE_ROOT = dirs$local, TISSUUMAPS_PROJECT_PATH = ""))

  for (backend in c("avivator", "tissuumaps")) {
    env <- .viewer_server_env(backend, dirs)
    expected <- if (backend == "avivator") {
      c("0112.ome.tiff", "7002.OME.TIF", "7001.ome.tif")
    } else c("0112", "7001")
    excluded <- if (backend == "avivator") "7003.ome.tiff" else "7003"
    forced <- shiny::reactiveVal(excluded)
    withr::local_envvar(c(ISLET_VIEWER_BACKEND = backend))
    withr::local_dir(root)
    suppressMessages(shiny::testServer(env$viewer_server,
      args = list(forced_image = forced, current_tab = shiny::reactive("Viewer")), {
        session$flushReact()
        expect_identical(available_images_r(), expected)
        expect_null(cohort_note_r())
        expect_identical(selected_basename_r(), expected[[1]])
        expect_identical(viewer_info()$selection, expected[[1]])
        expect_false(grepl("7003", output$image_selector$html, fixed = TRUE))
        expect_false(grepl("7003", output$viewer_frame$html, fixed = TRUE))

        forced(NULL)
        session$setInputs(selected_image = expected[[2]])
        expect_identical(selected_basename_r(), expected[[2]])
        expect_identical(viewer_info()$selection, expected[[2]])
        forced(excluded)
        session$flushReact()
        expect_identical(selected_basename_r(), expected[[2]])
        session$setInputs(selected_image = excluded)
        expect_identical(selected_basename_r(), expected[[1]])
        forced(expected[[2]])
        session$flushReact()
        expect_identical(selected_basename_r(), expected[[2]])
      }))
  }
})

test_that("assets excluded by a measured cohort render counted notes on both backends", {
  skip_if_not_installed("shiny")
  cases <- list(
    list(backend = "avivator", roots = "www", n = 1L),
    list(backend = "avivator", roots = "local", n = 1L),
    list(backend = "avivator", roots = c("www", "local"), n = 2L),
    list(backend = "tissuumaps", roots = c("www", "projects"), n = 1L),
    list(backend = "tissuumaps", roots = "www", n = 1L)
  )
  for (case in cases) {
    root <- withr::local_tempdir()
    dirs <- list(www = file.path(root, "www", "local_images"),
                 local = file.path(root, "local"), projects = file.path(root, "projects"))
    for (dir in dirs) dir.create(dir, recursive = TRUE)
    for (site in case$roots) {
      filename <- if (site == "projects") "6539.tmap" else "6539.ome.tiff"
      file.create(file.path(dirs[[site]], filename))
    }
    # Existing project files stay authoritative even with eligible TIFFs nearby.
    if ("projects" %in% case$roots) file.create(file.path(dirs$www, "0112.ome.tiff"))
    env <- .viewer_server_env(case$backend, dirs, .viewer_cohort(c("112", "6450")))
    withr::local_envvar(c(LOCAL_IMAGE_ROOT = dirs$local, ISLET_VIEWER_BACKEND = case$backend))
    withr::local_dir(root)
    forced <- shiny::reactive(if (case$backend == "avivator") "6539.ome.tiff" else "6539")
    suppressMessages(shiny::testServer(env$viewer_server,
      args = list(forced_image = forced, current_tab = shiny::reactive("Viewer")), {
        session$flushReact()
        expect_identical(viewer_cohort_r()$availability$state, "measured")
        expect_identical(available_images_r(), character(0))
        expect_null(selected_basename_r())
        expect_identical(image_inventory_r()$n_found, case$n)
        expect_identical(image_inventory_r()$n_excluded, case$n)
        expect_s3_class(cohort_note_r(), "ie_guard_note")
        expect_null(output$image_selector)
        html <- output$viewer_frame$html
        expect_match(html, "[viewer_cohort_assets_excluded]", fixed = TRUE)
        expect_match(html, sprintf("Assets found: %d; excluded by cohort membership: %d", case$n, case$n),
                     fixed = TRUE)
        expect_match(html, "eligible: 0", fixed = TRUE)
        expect_false(grepl("<iframe|Build one per donor|No images found|No TissUUmaps projects found", html))
      }))
  }
})

test_that("empty asset directories keep their missing-asset message on both backends", {
  skip_if_not_installed("shiny")
  root <- withr::local_tempdir()
  dirs <- list(www = file.path(root, "www", "local_images"),
               local = file.path(root, "local"), projects = file.path(root, "projects"))
  for (dir in dirs) dir.create(dir, recursive = TRUE)
  withr::local_envvar(c(LOCAL_IMAGE_ROOT = dirs$local))
  withr::local_dir(root)
  for (backend in c("avivator", "tissuumaps")) {
    env <- .viewer_server_env(backend, dirs)
    withr::local_envvar(c(ISLET_VIEWER_BACKEND = backend))
    shiny::testServer(env$viewer_server,
      args = list(forced_image = shiny::reactive(NULL), current_tab = shiny::reactive("Viewer")), {
        session$flushReact()
        expect_identical(available_images_r(), character(0))
        expect_identical(image_inventory_r()$n_found, 0L)
        expect_identical(image_inventory_r()$n_excluded, 0L)
        expect_null(cohort_note_r())
        message <- if (backend == "avivator") "No images found." else "No TissUUmaps projects found."
        expect_match(output$image_selector$html, message, fixed = TRUE)
        expect_false(grepl("viewer_cohort_assets_excluded|<iframe", output$viewer_frame$html))
      })
  }
})

test_that("missing cohort membership renders a named note on both backends", {
  skip_if_not_installed("shiny")
  root <- withr::local_tempdir()
  dirs <- .viewer_files(root)
  unavailable <- list(donor_ids = NULL, availability = ie_avail(
    "unavailable", "viewer_cohort_source_missing", source = "missing-build"))
  withr::local_envvar(c(LOCAL_IMAGE_ROOT = dirs$local))
  for (backend in c("avivator", "tissuumaps")) {
    env <- .viewer_server_env(backend, dirs, unavailable)
    withr::local_envvar(c(ISLET_VIEWER_BACKEND = backend))
    withr::local_dir(root)
    shiny::testServer(env$viewer_server,
      args = list(forced_image = shiny::reactive(NULL), current_tab = shiny::reactive("Viewer")), {
        session$flushReact()
        expect_length(available_images_r(), 0L)
        expect_null(selected_basename_r())
        expect_match(output$viewer_frame$html, "viewer_cohort_source_missing", fixed = TRUE)
        expect_match(output$viewer_frame$html, "[viewer_cohort]", fixed = TRUE)
        expect_false(grepl("<iframe", output$viewer_frame$html, fixed = TRUE))
      })
  }
  expect_error(.viewer_helpers()$filter_viewer_images("0112.ome.tif", unavailable),
               "viewer_cohort_source_missing")
})

.viewer_parquet <- function(root, ids) {
  dir.create(file.path(root, "parquet"), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbWriteTable(con, "cohort", data.frame(imageid = ids))
  path <- file.path(root, "parquet", "islets.parquet")
  DBI::dbExecute(con, paste("COPY cohort TO", DBI::dbQuoteString(con, path), "(FORMAT PARQUET)"))
  path
}

test_that("cohort membership comes from the active build, not the broader donor table", {
  skip_if_not_installed("duckdb")
  root <- withr::local_tempdir()
  path <- .viewer_parquet(root, c("0112", "112", "7001"))
  write.csv(data.frame(case_id = c("112", "7001", "7003")),
            file.path(root, "donor_metadata.csv"), row.names = FALSE)
  env <- .viewer_helpers()
  env$domain_data_root <- function() root
  cohort <- env$read_viewer_cohort()
  expect_setequal(cohort$donor_ids, c("112", "7001"))
  expect_identical(cohort$availability$state, "measured")
  expect_identical(cohort$availability$n, 2L)
  expect_identical(cohort$availability$source, path)
  expect_identical(env$filter_viewer_images(c("0112.tmap", "7003.tmap")), "0112.tmap")

  other_root <- withr::local_tempdir()
  .viewer_parquet(other_root, 7002L)
  env$domain_data_root <- function() other_root
  expect_identical(env$read_viewer_cohort()$donor_ids, "7002")
  expect_identical(env$filter_viewer_images(c("0112.tmap", "7002.tmap")), "7002.tmap")
})

test_that("the domain workbook supplies membership when Parquet is absent", {
  skip_if_not_installed("writexl")
  root <- withr::local_tempdir()
  path <- file.path(root, "cohort.xlsx")
  writexl::write_xlsx(list(Islet_Composition = data.frame(
    `Case ID` = c("0112", "112", "7002"), check.names = FALSE)), path)
  env <- .viewer_helpers()
  env$domain_data_root <- function() root
  env$domain_data_file <- function(key, default) path
  cohort <- env$read_viewer_cohort()
  expect_setequal(cohort$donor_ids, c("112", "7002"))
  expect_identical(cohort$availability$state, "measured")
  expect_identical(cohort$availability$source, path)
})

test_that("missing, empty, invalid and unreadable cohort sources never allow images", {
  skip_if_not_installed("duckdb")
  root <- withr::local_tempdir()
  env <- .viewer_helpers()
  env$domain_data_root <- function() root
  env$domain_data_file <- function(key, default) file.path(root, "missing.xlsx")
  cohort <- env$read_viewer_cohort()
  expect_identical(cohort$availability$state, "unavailable")
  expect_identical(cohort$availability$reason, "viewer_cohort_source_missing")
  expect_null(cohort$donor_ids)
  expect_error(env$filter_viewer_images("0112.ome.tiff"), "viewer_cohort_source_missing")

  path <- .viewer_parquet(root, character(0))
  expect_identical(env$read_viewer_cohort()$availability$reason, "viewer_cohort_has_no_donors")
  unlink(path)
  .viewer_parquet(root, c("112", NA_character_))
  expect_identical(env$read_viewer_cohort()$availability$reason, "viewer_cohort_invalid_donor_ids")
  writeLines("not a parquet file", path)
  cohort <- env$read_viewer_cohort()
  expect_identical(cohort$availability$state, "failed")
  expect_match(cohort$availability$reason, "viewer_cohort_read_failed", fixed = TRUE)
  expect_null(cohort$donor_ids)
  expect_error(env$filter_viewer_images("0112.ome.tiff"), "viewer_cohort_read_failed")
})
