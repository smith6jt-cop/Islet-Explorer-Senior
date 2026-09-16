# Post-Xenium IF overlay. Only prebuilt layer artifacts are read in the app.
# Shared restain readers also serve card_x3_identity.R, sourced after this file.

ie_restain_root <- function() {
  file.path(Sys.getenv("ISLET_LAYERS_DIR", unset = "/data/islet_layers"), "restain")
}

ie_restain_normalize_donor <- function(x) {
  x <- sub("^0+", "", trimws(as.character(x)))
  x[x == ""] <- "0"
  x
}

ie_restain_failure <- function(reason, state = "unavailable") {
  list(availability = ie_avail(state, reason = reason), data = NULL)
}

ie_restain_index <- function(root = ie_restain_root()) {
  path <- file.path(root, "availability.parquet")
  if (!file.exists(path)) return(ie_restain_failure("restain_artifact_not_built"))
  tryCatch({
    d <- as.data.frame(arrow::read_parquet(path))
    needed <- c("donor_id", "state", "reason", "n", "frame", "pixel_size_um")
    if (!all(needed %in% names(d)) || !nrow(d)) {
      stop("restain availability schema is missing or empty")
    }
    d$donor_id <- ie_restain_normalize_donor(d$donor_id)
    if (anyDuplicated(d$donor_id) || any(!d$state %in% IE_STATES) ||
        any(is.na(d$state)) ||
        any(d$state != "measured" & (is.na(d$reason) | !nzchar(d$reason)))) {
      stop("invalid restain donor availability")
    }
    list(availability = ie_avail("measured", n = nrow(d), source = path), data = d)
  }, error = function(e) ie_restain_failure(
    paste0("invalid_restain_artifact: ", conditionMessage(e)), "failed"))
}

ie_restain_ledger <- function(donor, root = ie_restain_root()) {
  path <- file.path(dirname(root), "availability")
  if (!dir.exists(path)) return(ie_avail("unavailable", "availability_artifact_not_built"))
  tryCatch({
    d <- as.data.frame(arrow::open_dataset(path, partitioning = NULL))
    if (!all(c("donor_id", "capability", "state", "reason", "n") %in% names(d))) {
      stop("invalid donor capability ledger")
    }
    d <- d[ie_restain_normalize_donor(d$donor_id) == donor &
             d$capability == "restain_if", , drop = FALSE]
    if (nrow(d) != 1L) return(ie_avail("unavailable", "no_restain_availability_record"))
    ie_avail(d$state[[1]], reason = d$reason[[1]], n = d$n[[1]], source = path)
  }, error = function(e) ie_avail("failed", paste0("invalid_availability_artifact: ",
                                                 conditionMessage(e))))
}

ie_restain_selection <- function(index, donor, root = ie_restain_root()) {
  if (!ie_is_measured(index$availability)) return(index)
  if (is.null(donor) || length(donor) != 1L || !nzchar(donor)) {
    return(ie_restain_failure("no_donor_selected"))
  }
  donor <- ie_restain_normalize_donor(donor)
  row <- index$data[index$data$donor_id == donor, , drop = FALSE]
  if (nrow(row) != 1L) return(ie_restain_failure("no_restain_pairing"))
  available <- ie_avail(row$state[[1]], reason = row$reason[[1]], n = row$n[[1]])
  if (ie_is_measured(available)) available <- ie_restain_ledger(donor, root)
  list(availability = available, data = row)
}

ie_restain_note <- function(availability, guard) {
  note <- ie_require_measured(availability, guard = guard,
    artifact = file.path(ie_restain_root(), "availability.parquet"),
    builder = "scripts/senior/build_restain_layer.py")
  if (is.null(note)) return(NULL)
  shiny::tags$p(class = "ie-card-note", format(note))
}

ie_restain_cells <- function(row, root = ie_restain_root()) {
  donor <- as.character(row$donor_id[[1]])
  # The path is derived from a validated identifier, never from an input path.
  if (!grepl("^[0-9]+$", donor)) stop("invalid restain donor identifier")
  d <- as.data.frame(arrow::read_parquet(file.path(root, "cells", paste0(donor, ".parquet"))))
  expected <- ie_frame("if_um", donor, "panc")
  needed <- c("donor_id", "cell_id", "object_id", "x_um", "y_um", "frame",
              "object_id_state", "object_id_reason", "feat__if_INS", "feat__if_GCG",
              "feat__if_INS_state", "feat__if_INS_reason",
              "feat__if_GCG_state", "feat__if_GCG_reason")
  if (!all(needed %in% names(d)) || !nrow(d) || anyNA(d$frame) ||
      any(d$frame != expected) || anyNA(d$donor_id) || any(d$donor_id != donor) ||
      anyNA(d$cell_id) || any(!nzchar(d$cell_id)) || anyDuplicated(d$cell_id) ||
      any(!is.na(d$object_id)) || anyNA(d$object_id_state) || anyNA(d$object_id_reason) ||
      any(d$object_id_state != "unavailable") || any(d$object_id_reason != "no_qupath_object_id") ||
      any(!is.finite(d$x_um)) || any(!is.finite(d$y_um)) || nrow(d) != row$n[[1]]) {
    stop("invalid IF-frame restain cells or count")
  }
  for (marker in c("INS", "GCG")) {
    field <- paste0("feat__if_", marker)
    state <- d[[paste0(field, "_state")]]
    reason <- d[[paste0(field, "_reason")]]
    if (anyNA(state) || any(!state %in% IE_STATES) ||
        any(state != "measured" & (is.na(reason) | !nzchar(reason))) ||
        any((state == "measured") != is.finite(d[[field]]))) {
      stop("invalid IF intensity availability")
    }
  }
  ie_attach_frame(d, expected)
}

ie_restain_reactives <- function(input, output, session, id, cohort, active_tab) {
  index <- shiny::reactive({
    shiny::req(identical(active_tab(), "CrossModal"))
    ie_restain_index()
  })
  choices <- shiny::reactive({
    d <- index()$data
    if (is.null(d)) return(character())
    ids <- as.character(d$donor_id)
    if (!is.null(cohort) && is.function(cohort$donor_ids)) {
      ids <- intersect(ids, ie_restain_normalize_donor(cohort$donor_ids()))
    }
    ids
  })
  shiny::observeEvent(choices(), {
    ids <- choices()
    selected <- shiny::isolate(input[[paste0(id, "_donor")]])
    if (is.null(selected) || !selected %in% ids) selected <- if (length(ids)) ids[[1]] else ""
    shiny::updateSelectInput(session, paste0(id, "_donor"), choices = ids, selected = selected)
  })
  shiny::reactive({
    selected <- input[[paste0(id, "_donor")]]
    current <- index()
    if (!ie_is_measured(current$availability)) return(current)
    if (is.null(selected) || !selected %in% choices()) {
      return(ie_restain_failure("no_donor_selected_in_cohort"))
    }
    ie_restain_selection(current, selected)
  })
}

card_x2_overlay_ui <- function(ns) {
  shiny::tags$div(class = "ie-card", style = "overflow: visible;",
    shiny::tags$h3(class = "ie-card-title", "Post-Xenium IF overlay"),
    shiny::selectInput(ns("x2_overlay_donor"), "Donor", choices = character()),
    shiny::selectInput(ns("x2_overlay_marker"), "IF intensity", choices = c("INS", "GCG")),
    shiny::uiOutput(ns("x2_overlay_content")))
}

card_x2_overlay_server <- function(input, output, session, prepared = NULL, cohort = NULL,
                                   active_tab = function() "CrossModal") {
  selection <- ie_restain_reactives(input, output, session, "x2_overlay", cohort, active_tab)
  cells <- shiny::reactive({
    selected <- selection()
    if (!ie_is_measured(selected$availability)) return(selected)
    tryCatch(list(availability = selected$availability, data = ie_restain_cells(selected$data)),
      error = function(e) ie_restain_failure(paste0("invalid_restain_cells: ",
                                                   conditionMessage(e)), "failed"))
  })
  drawing <- shiny::reactive({
    loaded <- cells()
    if (!ie_is_measured(loaded$availability)) return(loaded)
    marker <- input$x2_overlay_marker
    if (is.null(marker) || !marker %in% c("INS", "GCG")) marker <- "INS"
    all <- loaded$data
    field <- paste0("feat__if_", marker)
    valid <- all[[paste0(field, "_state")]] == "measured"
    if (!any(valid)) return(ie_restain_failure(paste0("no_finite_IF_", marker)))
    # Deterministic display sampling only; all counts come from the complete artifact.
    take <- unique(round(seq(1, nrow(all), length.out = min(nrow(all), 50000L))))
    background <- all[take, , drop = FALSE]
    foreground <- background[background[[paste0(field, "_state")]] == "measured", , drop = FALSE]
    foreground$intensity <- foreground[[field]]
    ie_assert_same_frame(background, foreground, context = "post-Xenium IF overlay")
    caption <- sprintf(
      "IF frame; Xenium centroids warped into IF. Showing %s of %s nucleated cells. %s: %s measured; %s unavailable. CSV pixel size: %.4f µm/px. Hormone calls are reported separately.",
      format(nrow(background), big.mark = ","), format(nrow(all), big.mark = ","), marker,
      format(sum(valid), big.mark = ","), format(sum(!valid), big.mark = ","),
      selection()$data$pixel_size_um[[1]])
    plot <- ggplot2::ggplot(background, ggplot2::aes(x = x_um, y = y_um)) +
      ggplot2::geom_point(colour = "#9AA0A6", size = 0.2) +
      ggplot2::geom_point(data = foreground, ggplot2::aes(colour = intensity), size = 0.2) +
      ggplot2::scale_colour_viridis_c(name = paste(marker, "IF intensity"),
                                    trans = scales::pseudo_log_trans(base = 10)) +
      ggplot2::scale_y_reverse() + ggplot2::coord_fixed() +
      ggplot2::labs(x = "IF x (µm)", y = "IF y (µm)") + theme_islet()
    list(availability = loaded$availability, plot = plot, caption = caption)
  })
  output$x2_overlay_content <- shiny::renderUI({
    d <- drawing()
    if (!ie_is_measured(d$availability)) return(ie_restain_note(d$availability, "x2_overlay"))
    shiny::tagList(shiny::plotOutput(session$ns("x2_overlay_plot"), height = "500px"),
      shiny::tags$p(class = "ie-card-note", d$caption),
      shiny::downloadButton(session$ns("x2_overlay_download"), "Download plot"))
  })
  output$x2_overlay_plot <- shiny::renderPlot({
    d <- drawing()
    shiny::validate(shiny::need(ie_is_measured(d$availability), "IF overlay unavailable [x2_overlay]"))
    d$plot
  })
  output$x2_overlay_download <- shiny::downloadHandler(
    filename = function() paste0("restain_", input$x2_overlay_donor, ".png"),
    content = function(file) {
      d <- drawing(); shiny::req(ie_is_measured(d$availability))
      ggplot2::ggsave(file, d$plot + ggplot2::labs(caption = d$caption),
                      width = 9, height = 8, dpi = 180)
    })
  invisible(list(selection = selection, cells = cells, drawing = drawing))
}
