.function_is_declared <- function(path, function_name) {
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  pattern <- paste0("(?m)^", function_name, "\\s*<-\\s*function\\s*\\(")
  grepl(pattern, text, perl = TRUE)
}

test_that("every established tab retains one UI and one server implementation", {
  contract <- data.frame(
    tab = c("Plot", "Trajectory", "Viewer", "Statistics", "Spatial", "Inflammation", "AI"),
    ui_file = c(
      "mod_plot_ui.R", "mod_trajectory_ui.R", "mod_viewer_ui.R",
      "mod_statistics_ui.R", "mod_spatial_ui.R", "mod_inflammation_ui.R",
      "mod_ai_assistant_ui.R"
    ),
    ui_function = c(
      "plot_main_ui", "trajectory_ui", "viewer_ui", "statistics_ui",
      "spatial_ui", "inflammation_ui", "ai_assistant_ui"
    ),
    server_file = c(
      "mod_plot_server.R", "mod_trajectory_server.R", "mod_viewer_server.R",
      "mod_statistics_server.R", "mod_spatial_server.R", "mod_inflammation_server.R",
      "mod_ai_assistant_server.R"
    ),
    server_function = c(
      "plot_server", "trajectory_server", "viewer_server", "statistics_server",
      "spatial_server", "inflammation_server", "ai_assistant_server"
    )
  )

  app_text <- paste(
    readLines(.find_repo_file("apps/islet/app.R"), warn = FALSE),
    collapse = "\n"
  )
  for (i in seq_len(nrow(contract))) {
    ui_path <- .find_repo_file(file.path("apps/islet/R", contract$ui_file[[i]]))
    server_path <- .find_repo_file(file.path("apps/islet/R", contract$server_file[[i]]))
    expect_true(.function_is_declared(ui_path, contract$ui_function[[i]]),
                info = paste(contract$tab[[i]], "UI declaration"))
    expect_true(.function_is_declared(server_path, contract$server_function[[i]]),
                info = paste(contract$tab[[i]], "server declaration"))
    expect_match(app_text, paste0(contract$ui_function[[i]], "("), fixed = TRUE,
                 info = paste(contract$tab[[i]], "UI wiring"))
    expect_match(app_text, paste0(contract$server_function[[i]], "("), fixed = TRUE,
                 info = paste(contract$tab[[i]], "server wiring"))
  }
})

test_that("shared helper implementations live only in ieCore", {
  shared <- c(
    "safe_left_join", "add_unit_key", "add_islet_key", "compute_diameter_um",
    "flag_outliers", "summary_stats", "aggregate_to_donor", "format_p_value",
    "theme_islet", "plotly_islet_fonts", "plotly_islet_config"
  )
  roots <- c(
    dirname(.find_repo_file("apps/R/ieCore/R/domain_config.R")),
    dirname(.find_repo_file("apps/islet/R/00_globals.R"))
  )
  files <- unlist(lapply(roots, list.files, pattern = "[.]R$", full.names = TRUE))
  definitions <- vapply(shared, function(name) {
    pattern <- paste0("(?m)^", name, "\\s*<-\\s*function\\s*\\(")
    sum(vapply(files, function(path) {
      grepl(pattern, paste(readLines(path, warn = FALSE), collapse = "\n"), perl = TRUE)
    }, logical(1)))
  }, integer(1))

  expect_identical(unname(definitions), rep(1L, length(shared)),
                   info = paste(names(definitions), definitions, collapse = "; "))
  expect_false(file.exists(.find_repo_file("apps/islet/R/00_globals.R") |> dirname() |>
                             file.path("utils_stats.R")))
})
