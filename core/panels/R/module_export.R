#' Export module — Seurat object, cluster assignments, marker tables,
#' and a self-contained HTML session report.
#'
#' Reads from `app_state` only — never mutates state. Each download is a
#' separate `downloadHandler` so the user picks what they want.

export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::h5("Source"),
      shiny::radioButtons(ns("seurat_source"),
                          "Seurat object source",
                          inline = TRUE,
                          choices = c("Root (Cluster tab)" = "root",
                                      "Stack top (Subcluster)" = "top")),
      shiny::hr(),
      shiny::h5("Downloads"),
      shiny::p(shiny::strong("Seurat object")),
      shiny::downloadButton(ns("dl_seurat"), "Download .qs2",
                            class = "btn-primary"),
      shiny::br(), shiny::br(),
      shiny::p(shiny::strong("Cluster assignments")),
      shiny::downloadButton(ns("dl_clusters"), "Download CSV"),
      shiny::br(), shiny::br(),
      shiny::p(shiny::strong("Marker tables (combined)")),
      shiny::downloadButton(ns("dl_markers"), "Download CSV"),
      shiny::br(), shiny::br(),
      shiny::p(shiny::strong("Session report")),
      shiny::downloadButton(ns("dl_report"), "Download HTML"),
      shiny::hr(),
      shiny::p(shiny::em("Cluster + marker downloads pull from the source ",
                         "selected above. The session report includes the ",
                         "loaded dataset, validation, the cluster pipeline ",
                         "run-log, the subcluster stack, and the most ",
                         "recently computed marker table."))
    ),
    bslib::card(
      bslib::card_header("Status"),
      shiny::uiOutput(ns("status"))
    )
  )
}

export_server <- function(id, panels, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    pick_seurat <- function() {
      if (identical(input$seurat_source, "top")) {
        stk <- app_state$cluster_stack %||% list()
        if (length(stk)) stk[[length(stk)]]$obj else NULL
      } else {
        app_state$xen_clustered
      }
    }

    output$status <- shiny::renderUI({
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      x  <- pick_seurat()
      cache <- app_state$markers_cache %||% list()
      shiny::tagList(
        shiny::p(shiny::strong("Now: "), ts),
        shiny::p(shiny::strong("Loaded dataset: "),
                 if (!is.null(app_state$xen))
                   sprintf("%d cells × %d genes", ncol(app_state$xen),
                           nrow(app_state$xen))
                 else "(none)"),
        shiny::p(shiny::strong("Selected source: "),
                 if (is.null(x)) "(no clustered run)"
                 else sprintf("%d cells, %d resolutions",
                              ncol(x),
                              length(cluster_resolution_columns(x)))),
        shiny::p(shiny::strong("Marker cache entries: "),
                 length(cache)),
        shiny::p(shiny::strong("Subcluster stack depth: "),
                 length(app_state$cluster_stack %||% list()))
      )
    })

    output$dl_seurat <- shiny::downloadHandler(
      filename = function() sprintf("xenium_%s.qs2",
                                     format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        x <- pick_seurat()
        if (is.null(x)) {
          stop("No clustered Seurat object available — run the Cluster ",
               "or Subcluster tab first.")
        }
        qs2::qs_save(x, file)
      })

    output$dl_clusters <- shiny::downloadHandler(
      filename = function() sprintf("clusters_%s.csv",
                                     format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        utils::write.csv(build_cluster_csv(pick_seurat()),
                         file, row.names = FALSE)
      })

    output$dl_markers <- shiny::downloadHandler(
      filename = function() sprintf("markers_%s.csv",
                                     format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        cache <- app_state$markers_cache %||% list()
        utils::write.csv(build_markers_csv(cache),
                         file, row.names = FALSE)
      })

    output$dl_report <- shiny::downloadHandler(
      filename = function() sprintf("session_report_%s.html",
                                     format(Sys.time(), "%Y%m%dT%H%M%S")),
      content  = function(file) {
        render_session_report(
          file, panels(), app_state,
          custom_label = custom_panel_label(app_state$custom_panel_status,
                                            panels()))
      })
  })
}
