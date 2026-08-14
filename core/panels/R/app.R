#' Launch the Xenium Panel Explorer Shiny app.
#'
#' Returns a `shiny.appobj` that callers pass to `shiny::runApp` (or print).
#' The navbar wires every tab module declared in §5 of CLAUDE.md; M2 lights
#' up Overview and Panel Browser, the rest remain placeholder cards.
#'
#' @export
xenium_panel_app <- function() {
  env_check <- check_environment()
  print_environment_check(env_check)

  ui <- bslib::page_navbar(
    id    = "main_nav",
    title = "Xenium Panel Explorer",
    theme = bslib::bs_theme(version = 5, preset = "shiny"),
    bslib::nav_panel("Overview",            overview_ui("overview")),
    bslib::nav_panel("Panel Browser",       panel_browser_ui("panel_browser")),
    bslib::nav_panel("Load Xenium",         load_xenium_ui("load_xenium")),
    bslib::nav_panel("Panel-vs-Data",       panel_compare_ui("panel_compare")),
    bslib::nav_panel("Cluster",             cluster_ui("cluster")),
    bslib::nav_panel("Subcluster",          subcluster_ui("subcluster")),
    bslib::nav_panel("Markers",             marker_ui("marker")),
    bslib::nav_panel("Clustree",            clustree_ui("clustree")),
    bslib::nav_panel("Export",              export_ui("export"))
  )

  server <- function(input, output, session) {
    # Default tissue: the first manifest under data/tissues/ alphabetically.
    initial_tissues <- available_tissues()
    initial_tissue  <- if (length(initial_tissues)) initial_tissues[1] else NULL

    # App-wide reactive state. `selected_subpanel` is the canonical key
    # (filename stem), e.g. "01_pancreas_endocrine". `nav_target` is the
    # `value` of the nav_panel to switch to.
    app_state <- shiny::reactiveValues(
      selected_tissue       = initial_tissue,
      selected_subpanel     = NULL,
      nav_target            = NULL,
      xen                   = NULL,
      xen_path              = NULL,
      xen_clustered         = NULL,
      cluster_error         = NULL,
      cluster_jump_res      = NULL,
      compare_min_det       = 0,
      compare_topn          = 10,
      cluster_stack         = list(),
      subcluster_error      = NULL,
      markers_error         = NULL,
      markers_cache         = list(),
      markers_last_key      = NULL,
      env_warnings          = env_check$warnings,
      env_errors            = env_check$errors,
      # Session-scoped override of `panels$custom`. When non-NULL, spliced
      # in by the `panels` reactive below. Cleared on tissue switch.
      custom_panel_override = NULL,
      custom_panel_status   = NULL
    )

    # Whenever a fresh root run lands in xen_clustered, reset the
    # subcluster stack so it begins again at the new root.
    shiny::observeEvent(app_state$xen_clustered, {
      x <- app_state$xen_clustered
      if (is.null(x)) {
        app_state$cluster_stack <- list()
        return()
      }
      app_state$cluster_stack <- list(list(
        obj            = x,
        label          = "root",
        parent_res     = NA_character_,
        parent_cluster = NA_character_,
        opts           = NULL
      ))
    }, ignoreNULL = FALSE)

    # Re-runs whenever the selected tissue changes. Hits the disk again
    # so a user editing data/tissues/<t>/ at runtime is picked up.
    panels_default <- shiny::reactive({
      shiny::req(app_state$selected_tissue)
      load_panels(app_state$selected_tissue)
    })
    # `panels` splices any user-uploaded custom-panel override into
    # `panels$custom` so every downstream module sees it.
    panels <- shiny::reactive({
      base <- panels_default()
      ov   <- app_state$custom_panel_override
      if (!is.null(ov)) base$custom <- ov
      base
    })

    # Clearing the custom-panel override on tissue switch is required:
    # an upload enriched against (say) pancreas T1D-GWAS does not make
    # sense as the "custom panel" for thymus.
    shiny::observeEvent(app_state$selected_tissue, {
      if (!is.null(app_state$custom_panel_override) ||
          !is.null(app_state$custom_panel_status)) {
        app_state$custom_panel_override <- NULL
        app_state$custom_panel_status   <- NULL
      }
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    overview_server("overview",          panels, app_state)
    panel_browser_server("panel_browser", panels, panels_default, app_state)
    load_xenium_server("load_xenium",    panels, app_state)
    panel_compare_server("panel_compare", panels, app_state)
    cluster_server("cluster", panels, app_state)
    subcluster_server("subcluster", panels, app_state)
    marker_server("marker", panels, app_state)
    clustree_server("clustree", panels, app_state)
    export_server("export", panels, app_state)

    # Cross-tab routing: when a module sets app_state$nav_target, switch
    # the navbar to that panel.
    shiny::observeEvent(app_state$nav_target, {
      target <- app_state$nav_target
      if (is.null(target) || !nzchar(target)) return()
      bslib::nav_select(id = "main_nav", selected = target)
      app_state$nav_target <- NULL
    }, ignoreNULL = TRUE)

    # Sweet-alert toasts. The status alert cards inside each module
    # remain as the persistent record; this is a transient pop-up
    # confirmation when an error transitions from NULL to a string.
    .toast_error <- function(slot_name, title) {
      shiny::observeEvent(app_state[[slot_name]], {
        msg <- app_state[[slot_name]]
        if (is.null(msg) || !nzchar(msg)) return()
        shinyWidgets::sendSweetAlert(session = session,
                                     title   = title,
                                     text    = msg,
                                     type    = "error",
                                     btn_labels = "OK")
      }, ignoreNULL = TRUE, ignoreInit = TRUE)
    }
    .toast_error("cluster_error",    "Cluster pipeline failed")
    .toast_error("subcluster_error", "Subcluster pipeline failed")
    .toast_error("markers_error",    "Marker computation failed")
  }

  shiny::shinyApp(ui, server)
}
