# S2 card interface for the layer registry: card_s2_trajectory_ui(ns) and
# card_s2_trajectory_server(input, output, session, ...). Reads layer artifacts only.
# The builder preserves the source fit and ranks core medians within each donor.

.s2_trajectory_note <- function(state, reason, guard = "s2_trajectory", artifact = "trajectory/") {
  ie_require_measured(ie_avail(state, reason = reason), guard = guard,
                      artifact = artifact, builder = "scripts/senior/build_trajectory_layer.py")
}

.s2_trajectory_note_ui <- function(note) {
  shiny::tags$p(class = "ie-card-note", role = "status", format(note))
}

.s2_trajectory_norm_donor <- function(x) {
  x <- sub("^0+", "", trimws(as.character(x)))
  x[!is.na(x) & !nzchar(x)] <- "0"
  x
}

.s2_trajectory_read <- function(root, name, required) {
  path <- file.path(root, name)
  if (!dir.exists(path)) {
    return(list(note = .s2_trajectory_note("unavailable", "layer_artifact_not_found",
                                          artifact = path)))
  }
  tryCatch({
    # Layer roots are the spec's explicit exception to the app-data path convention.
    # trajectory/domain_data_root is checked against the active domain before display.
    data <- as.data.frame(arrow::open_dataset(path, format = "parquet", partitioning = NULL))
    missing <- setdiff(required, names(data))
    if (length(missing)) stop("missing columns: ", paste(missing, collapse = ", "))
    if (!nrow(data)) {
      return(list(note = .s2_trajectory_note("unavailable", "layer_artifact_has_no_rows",
                                            artifact = path)))
    }
    data$donor_id <- .s2_trajectory_norm_donor(data$donor_id)
    list(data = data, note = NULL)
  }, error = function(e) {
    list(note = .s2_trajectory_note("failed", conditionMessage(e), artifact = path))
  })
}

.s2_trajectory_domain_note <- function(data, active_root = domain_data_root()) {
  recorded <- unique(as.character(data$domain_data_root))
  if (length(recorded) != 1L || is.na(recorded) || !nzchar(recorded) ||
      !identical(normalizePath(recorded, mustWork = FALSE),
                 normalizePath(active_root, mustWork = FALSE))) {
    return(.s2_trajectory_note("unavailable", "trajectory_belongs_to_another_domain"))
  }
  NULL
}

.s2_trajectory_availability_note <- function(ledger, donor) {
  if (!is.null(ledger$note)) return(ledger$note)
  row <- ledger$data[!is.na(ledger$data$donor_id) & ledger$data$donor_id == donor &
                       ledger$data$capability == "trajectory", , drop = FALSE]
  if (nrow(row) != 1L) {
    return(.s2_trajectory_note("unavailable", "no_unique_trajectory_availability_for_donor",
                               artifact = "availability/"))
  }
  if (is.na(row$state) || !row$state %in% IE_STATES ||
      (row$state != "measured" && (is.na(row$reason) || !nzchar(row$reason)))) {
    return(.s2_trajectory_note("failed", "invalid_trajectory_availability_state",
                               artifact = "availability/"))
  }
  ie_require_measured(ie_avail(row$state, reason = row$reason),
                      guard = "s2_trajectory_availability", artifact = "availability/")
}

.s2_trajectory_caption <- function(data) {
  row <- data[1L, , drop = FALSE]
  number <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  sampling <- if (identical(row$sampling_state, "measured")) {
    sprintf("Ordering was fit on %s of %s cells.",
            number(row$n_ordered_cells), number(row$n_population_cells))
  } else {
    sprintf("Sampling %s: %s.", row$sampling_state, row$sampling_reason)
  }
  batch <- if (identical(row$batch_eta_sq_state, "measured")) {
    sprintf(paste0("Within-status batch η² (stain group; %s): %.6g; ",
                   "%s sampled cells, %s stain groups."),
            row$donor_status, row$batch_eta_sq_within_status,
            number(row$batch_n_cells), number(row$batch_n_groups))
  } else {
    sprintf("Within-status batch η² %s: %s.",
            row$batch_eta_sq_state, row$batch_eta_sq_reason)
  }
  paste(sampling, batch,
        "Islets are ordered by core median within the selected donor; ties share a rank.",
        "Ranks summarize the existing sampled-cell fit; the fit was not recomputed per donor.",
        "Pseudotime is an ordering, not elapsed time or evidence of disease progression.",
        sprintf("%s of %s islets have a core ordering; %s have a measured peri-minus-core difference.",
                number(sum(data$state == "measured")), number(sum(!is.na(data$islet_key))),
                number(sum(data$delta_state == "measured"))))
}

.s2_trajectory_plot_data <- function(data) {
  donors <- unique(data$donor_id)
  if (length(donors) != 1L || is.na(donors)) stop("trajectory plot requires exactly one donor")
  frames <- unique(data$frame)
  if (length(frames) != 1L || is.na(frames) || frames != "none") {
    stop("trajectory ordering requires the declared scale-free frame 'none'")
  }
  keep <- data$state == "measured" & data$delta_state == "measured"
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) stop("no measured peri-minus-core pairs")
  if (any(!is.finite(data$within_donor_rank)) ||
      any(!is.finite(data$sc_pseudotime_peri_minus_core))) {
    stop("measured trajectory row has no finite rank or paired difference")
  }
  ie_attach_frame(data, "none")
}

card_s2_trajectory_ui <- function(ns) {
  shiny::tags$div(class = "ie-card", style = "overflow: visible;",
    shiny::tags$h3(class = "ie-card-title", "Within-donor islet trajectory"),
    shiny::uiOutput(ns("s2_trajectory_controls")),
    shiny::uiOutput(ns("s2_trajectory_note")),
    shiny::uiOutput(ns("s2_trajectory_caption")),
    shiny::uiOutput(ns("s2_trajectory_plot_ui"))
  )
}

card_s2_trajectory_server <- function(input, output, session, prepared = NULL, cohort = NULL,
                                      active_tab = NULL, active_tabs = "Spatial",
                                      layers_dir = Sys.getenv("ISLET_LAYERS_DIR", "/data/islet_layers")) {
  active <- shiny::reactive(is.null(active_tab) || active_tab() %in% active_tabs)
  layer <- shiny::reactive({
    shiny::req(active())
    result <- .s2_trajectory_read(layers_dir, "trajectory", c(
      "donor_id", "donor_status", "islet_key", "within_donor_rank", "state", "reason",
      "sc_pseudotime_peri_minus_core", "delta_state", "delta_reason", "frame", "domain_data_root",
      "batch_eta_sq_within_status", "batch_eta_sq_state", "batch_eta_sq_reason",
      "batch_n_cells", "batch_n_groups", "n_ordered_cells", "n_population_cells",
      "sampling_state", "sampling_reason"
    ))
    if (is.null(result$note)) result$note <- .s2_trajectory_domain_note(result$data)
    if (is.null(result$note) && all(is.na(result$data$donor_id))) {
      result$note <- .s2_trajectory_note(result$data$state[[1L]], result$data$reason[[1L]])
    }
    if (!is.null(result$note)) return(list(note = result$note))
    result
  })
  ledger <- shiny::reactive({
    shiny::req(active())
    .s2_trajectory_read(layers_dir, "availability", c("donor_id", "capability", "state", "reason"))
  })
  output$s2_trajectory_controls <- shiny::renderUI({
    if (!active()) return(NULL)
    loaded <- layer()
    if (!is.null(loaded$note)) return(NULL)
    donors <- sort(unique(stats::na.omit(loaded$data$donor_id)))
    current <- shiny::isolate(input$s2_trajectory_donor)
    selected <- if (length(current) == 1L && current %in% donors) current else donors[[1L]]
    shiny::selectInput(session$ns("s2_trajectory_donor"), "Donor", donors, selected = selected)
  })
  selected <- shiny::reactive({
    loaded <- layer()
    if (!is.null(loaded$note)) return(loaded)
    donor <- input$s2_trajectory_donor
    if (is.null(donor)) donor <- sort(unique(stats::na.omit(loaded$data$donor_id)))[[1L]]
    data <- loaded$data[!is.na(loaded$data$donor_id) & loaded$data$donor_id == donor, , drop = FALSE]
    if (!nrow(data)) return(list(note = .s2_trajectory_note("unavailable", "donor_has_no_trajectory_rows")))
    note <- .s2_trajectory_availability_note(ledger(), donor)
    if (!is.null(note)) return(list(note = note))
    if (!any(data$state == "measured" & data$delta_state == "measured")) {
      reasons <- unique(c(data$reason, data$delta_reason))
      reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
      note <- .s2_trajectory_note("unavailable", paste(c("no_measured_peri_minus_core_pairs", reasons),
                                                       collapse = "; "))
    }
    list(data = data, note = note)
  })
  output$s2_trajectory_note <- shiny::renderUI({
    if (!active()) return(NULL)
    note <- selected()$note
    if (!is.null(note)) .s2_trajectory_note_ui(note)
  })
  output$s2_trajectory_caption <- shiny::renderUI({
    if (!active()) return(NULL)
    data <- selected()$data
    if (is.null(data)) return(NULL)
    shiny::tags$p(class = "ie-card-caption", .s2_trajectory_caption(data))
  })
  plot <- shiny::reactive({
    shiny::req(active())
    chosen <- selected()
    if (!is.null(chosen$note)) return(NULL)
    data <- .s2_trajectory_plot_data(chosen$data)
    ie_assert_same_frame(data, data, context = "within-donor islet trajectory")
    ggplot2::ggplot(data, ggplot2::aes(x = within_donor_rank, y = sc_pseudotime_peri_minus_core)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
      ggplot2::geom_point(alpha = 0.65) +
      ggplot2::labs(x = "Core pseudotime rank within donor", y = "Peri minus core pseudotime",
                    title = paste("Donor", unique(data$donor_id)),
                    caption = paste(strwrap(.s2_trajectory_caption(chosen$data), width = 100),
                                    collapse = "\n")) +
      theme_islet()
  })
  output$s2_trajectory_plot_ui <- shiny::renderUI({
    if (!active() || !is.null(selected()$note)) return(NULL)
    shiny::tagList(shiny::plotOutput(session$ns("s2_trajectory_plot"), height = 360),
                   shiny::downloadButton(session$ns("s2_trajectory_download"), "Download PNG"))
  })
  output$s2_trajectory_plot <- shiny::renderPlot({ plot() })
  output$s2_trajectory_download <- shiny::downloadHandler(
    filename = function() paste0("islet_trajectory_", input$s2_trajectory_donor, ".png"),
    content = function(file) {
      p <- plot()
      shiny::req(!is.null(p))
      ggplot2::ggsave(file, p, width = 10, height = 7, dpi = 300, bg = "white")
    }
  )
  invisible(list(data = selected, plot = plot))
}
