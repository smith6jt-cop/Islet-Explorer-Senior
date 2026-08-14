test_that("end-to-end smoke: app starts, navbar lists every tab", {
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  # chromote auto-detects the Chrome binary; if it fails, skip.
  bin <- tryCatch(chromote::find_chrome(), error = function(e) NULL)
  skip_if(is.null(bin) || !nzchar(bin), "no Chrome binary detected")

  app <- shinytest2::AppDriver$new(
    app_dir = rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
    name    = "xenium_panel_app",
    seed    = 1L,
    timeout = 60 * 1000,
    load_timeout = 60 * 1000,
    options = list(shiny.testmode = TRUE)
  )
  on.exit(app$stop(), add = TRUE)

  expected <- c("Overview", "Panel Browser", "Load Xenium",
                "Panel-vs-Data", "Cluster", "Subcluster",
                "Markers", "Clustree", "Export")
  html <- app$get_html("body")
  for (tab in expected) {
    expect_match(html, tab, fixed = TRUE,
                 info = sprintf("expected nav label '%s'", tab))
  }
})
