# Small layer artifacts only. Both cards use these readers; no cell-tree scan
# or tissue-coordinate transform belongs in a Shiny reactive.

.ie_distance_norm <- function(x) {
  value <- sub("^0+", "", trimws(as.character(x)))
  value[!is.na(value) & value == ""] <- "0"
  value
}

.ie_distance_note <- function(value, guard, artifact = "distance_profiles/data.parquet",
                               builder = "scripts/senior/build_distance_profiles.py") {
  note <- ie_require_measured(
    value, guard = guard, artifact = artifact, builder = builder
  )
  if (is.null(note)) return(NULL)
  shiny::tags$p(class = "ie-card-note", role = "status",
                paste0(value$state, "(", value$reason, "): ", format(note)))
}

.ie_distance_read <- function() {
  root <- Sys.getenv("ISLET_LAYERS_DIR", unset = "/data/islet_layers")
  path <- file.path(root, "distance_profiles", "data.parquet")
  absent <- function(reason, state = "unavailable") {
    list(value = ie_avail(state, reason, source = path))
  }
  if (!file.exists(path)) return(absent("distance_profiles_export_absent"))
  if (!requireNamespace("arrow", quietly = TRUE)) return(absent("arrow_reader_unavailable"))
  tryCatch({
    table <- arrow::read_parquet(path, as_data_frame = FALSE)
    # The external layers directory is explicit in the artifact contract. Its
    # producer must still identify the active domain, so CODEX cannot consume
    # PhenoCycler measurements simply because some donor ids overlap.
    source_root <- table$schema$metadata[["source_domain_root"]]
    if (is.null(source_root) || !exists("domain_data_root", mode = "function")) {
      return(absent("distance_profile_domain_not_recorded"))
    }
    if (!identical(normalizePath(source_root, mustWork = FALSE),
                   normalizePath(domain_data_root(), mustWork = FALSE))) {
      return(absent("distance_profile_domain_mismatch"))
    }
    data <- as.data.frame(table)
    required <- c("donor_id", "donor_status", "distance_bin_um", "n_cells", "n_immune",
                  "immune_fraction", "immune_subtype", "resolved_by", "frame",
                  "compartment", "state", "reason")
    if (!all(required %in% names(data))) return(absent("distance_profile_schema_invalid", "failed"))
    if (!nrow(data)) return(absent("distance_profile_export_empty"))
    if (anyNA(data$state) || any(!data$state %in% IE_STATES) ||
        anyNA(data$donor_id) || any(!nzchar(data$donor_id))) {
      return(absent("distance_profile_state_or_donor_invalid", "failed"))
    }
    data$donor_id <- .ie_distance_norm(data$donor_id)
    for (column in c("distance_bin_um", "n_cells", "n_immune", "immune_fraction")) {
      data[[column]] <- as.numeric(data[[column]])
    }
    measured <- data[data$state == "measured", , drop = FALSE]
    if (any(!is.finite(measured$n_cells) | measured$n_cells <= 0 |
            !is.finite(measured$n_immune) | measured$n_immune < 0 |
            measured$n_immune > measured$n_cells |
            !is.finite(measured$distance_bin_um) |
            !is.finite(measured$immune_fraction)) ||
        any(abs(measured$immune_fraction - measured$n_immune / measured$n_cells) > 1e-12) ||
        anyNA(measured$donor_status) || anyNA(measured$compartment)) {
      return(absent("distance_profile_measurement_invalid", "failed"))
    }
    unmeasured <- data[data$state != "measured", , drop = FALSE]
    if (anyNA(unmeasured$reason) || any(!nzchar(unmeasured$reason)) ||
        any(!is.na(unmeasured$n_cells) | !is.na(unmeasured$n_immune) |
            !is.na(unmeasured$immune_fraction))) {
      return(absent("distance_profile_absence_invalid", "failed"))
    }
    paths <- list.files(file.path(root, "availability"), pattern = "[.]parquet$",
                        recursive = TRUE, full.names = TRUE)
    ledger <- if (!length(paths)) {
      ie_avail("unavailable", "availability_export_absent")
    } else {
      as.data.frame(arrow::open_dataset(paths, format = "parquet", partitioning = NULL))
    }
    if (is.data.frame(ledger)) {
      if (!all(c("donor_id", "capability", "state", "reason", "n", "source") %in% names(ledger))) {
        return(absent("availability_schema_invalid", "failed"))
      }
      ledger$donor_id <- .ie_distance_norm(ledger$donor_id)
    }
    list(value = ie_avail("measured", n = nrow(data), source = path), data = data,
         ledger = ledger, bin_width = table$schema$metadata[["bin_width_um"]])
  }, error = function(e) absent(paste0("distance_profile_read_failed: ", conditionMessage(e)), "failed"))
}

.ie_distance_donors <- function(layer, cohort) {
  if (!ie_is_measured(layer$value)) return(character())
  donors <- unique(layer$data$donor_id)
  if (!is.null(cohort) && is.function(cohort$donor_ids)) {
    donors <- intersect(donors, .ie_distance_norm(cohort$donor_ids()))
  }
  sort(donors)
}

.ie_distance_pair <- function(data) paste(data$immune_subtype, data$resolved_by, sep = " | ")

.ie_distance_select <- function(layer, donor, pair = "Immune | all") {
  if (!ie_is_measured(layer$value)) return(list(value = layer$value))
  unavailable <- function(reason) list(value = ie_avail("unavailable", reason))
  if (is.null(donor) || !length(donor)) return(unavailable("no_donor_in_current_selection"))
  data <- layer$data[layer$data$donor_id == donor, , drop = FALSE]
  if (!nrow(data)) return(unavailable("distance_profile_donor_absent"))
  if (any(data$state != "measured")) {
    row <- data[which(data$state != "measured")[[1L]], , drop = FALSE]
    return(list(value = ie_avail(row$state, row$reason)))
  }
  if (inherits(layer$ledger, "ie_avail")) return(list(value = layer$ledger))
  for (capability in c("phenocycler_panc", "islet_outline")) {
    entry <- layer$ledger[layer$ledger$donor_id == donor &
                            layer$ledger$capability == capability, , drop = FALSE]
    if (nrow(entry) != 1L) return(unavailable(paste0("availability_row_missing_or_ambiguous:", capability)))
    value <- ie_avail(entry$state, entry$reason, as.numeric(entry$n), entry$source)
    if (!ie_is_measured(value)) return(list(value = value))
  }
  provenance <- sort(unique(data$resolved_by[data$resolved_by != "all" & data$n_immune > 0]))
  data <- data[.ie_distance_pair(data) == pair, , drop = FALSE]
  if (!nrow(data)) return(unavailable("immune_subtype_not_recorded"))
  list(value = ie_avail("measured"), data = data, provenance = provenance)
}

.ie_distance_framed <- function(data) {
  frames <- unique(data$frame)
  if (length(frames) != 1L || is.na(frames) || !nzchar(frames)) {
    stop("Distance profile contains mixed or missing frames", call. = FALSE)
  }
  ie_attach_frame(data, frames[[1L]])
}

.ie_distance_controls <- function(ns, prefix, layer, donors, donor, pair) {
  if (!length(donors)) return(NULL)
  if (is.null(donor) || !donor %in% donors) donor <- donors[[1L]]
  data <- layer$data[layer$data$donor_id == donor & layer$data$state == "measured", , drop = FALSE]
  pairs <- unique(c("Immune | all", .ie_distance_pair(data)))
  if (is.null(pair) || !pair %in% pairs) pair <- pairs[[1L]]
  shiny::fluidRow(
    shiny::column(4, shiny::selectInput(ns(paste0(prefix, "_donor")), "Donor", donors, donor)),
    shiny::column(8, shiny::selectInput(ns(paste0(prefix, "_pair")),
                                       "Immune call | resolved by", pairs, pair))
  )
}

.ie_distance_theme <- function() {
  if (exists("theme_islet", mode = "function")) theme_islet() else ggplot2::theme_minimal()
}

card_x1_gradient_ui <- function(ns) {
  shiny::div(class = "ie-card", style = "overflow: visible;",
             shiny::tags$h3(class = "ie-card-title", "Peri-islet immune gradient"),
             shiny::uiOutput(ns("x1_gradient_controls")),
             shiny::uiOutput(ns("x1_gradient_body")))
}

card_x1_gradient_server <- function(input, output, session, prepared = NULL,
                                    cohort = NULL, active_tab = NULL, tab_key = NULL, ...) {
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                             identical(active_tab(), tab_key))
  layer <- shiny::reactive({
    shiny::req(active())
    .ie_distance_read()
  })
  donors <- shiny::reactive(.ie_distance_donors(layer(), cohort))
  selected <- shiny::reactive({
    donor <- input$x1_gradient_donor
    if (is.null(donor) || !donor %in% donors()) donor <- head(donors(), 1L)
    pair <- input$x1_gradient_pair
    if (is.null(pair)) pair <- "Immune | all"
    .ie_distance_select(layer(), donor, pair)
  })
  output$x1_gradient_controls <- shiny::renderUI({
    .ie_distance_controls(session$ns, "x1_gradient", layer(), donors(),
                          input$x1_gradient_donor, input$x1_gradient_pair)
  })
  output$x1_gradient_body <- shiny::renderUI({
    value <- selected()
    if (!ie_is_measured(value$value)) return(.ie_distance_note(value$value, "guard_x1_gradient"))
    shiny::tagList(
      shiny::plotOutput(session$ns("x1_gradient_plot")),
      shiny::tags$p("Existing immune calls / all QC-eligible, non-excluded cells in each bin. ",
                    "Signed distance is negative inside the islet and positive outside; zero is retained. ",
                    "Resolved by: ", if (length(value$provenance)) paste(value$provenance, collapse = ", ")
                    else "no immune calls observed", "."),
      shiny::downloadButton(session$ns("x1_gradient_download"), "Download plot"),
      shiny::tableOutput(session$ns("x1_gradient_counts"))
    )
  })
  plot_data <- shiny::reactive({
    value <- selected()
    shiny::req(ie_is_measured(value$value))
    source <- .ie_distance_framed(value$data)
    reference <- ie_attach_frame(data.frame(), ie_frame("signed_distance_um"))
    ie_assert_same_frame(source, reference, "peri-islet immune gradient")
    data <- stats::aggregate(cbind(n_cells, n_immune) ~ donor_id + donor_status + distance_bin_um,
                             data = source, FUN = sum)
    data$immune_fraction <- data$n_immune / data$n_cells
    ie_attach_frame(data[order(data$distance_bin_um), ], ie_frame_of(source))
  })
  plot <- shiny::reactive({
    data <- plot_data()
    reference <- ie_attach_frame(data.frame(distance_bin_um = 0), ie_frame("signed_distance_um"))
    ie_assert_same_frame(data, reference, "peri-islet immune gradient plot")
    ggplot2::ggplot(data, ggplot2::aes(distance_bin_um, immune_fraction)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = "Signed distance to islet edge (µm)", y = "Immune fraction of eligible cells",
                    title = paste("Donor", data$donor_id[[1L]], "—", data$donor_status[[1L]]),
                    subtitle = .ie_distance_pair(selected()$data)[[1L]],
                    caption = paste("Bin labels are outer edges; width", layer()$bin_width,
                                    "µm. Zero has its own bin. Existing typing, no refitted thresholds.")) +
      .ie_distance_theme()
  })
  output$x1_gradient_plot <- shiny::renderPlot(plot())
  output$x1_gradient_counts <- shiny::renderTable(plot_data(), digits = 6)
  output$x1_gradient_download <- shiny::downloadHandler(
    filename = function() "peri-islet-immune-gradient.pdf",
    content = function(file) ggplot2::ggsave(file, plot = plot(), width = 8, height = 5)
  )
  invisible(list(active = active, selected = selected, plot = plot, data = plot_data))
}
