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

# The standing methods text for this card, in ONE place.
#
# It was five sentences printed above the figure on every render: the sampling, the ranking rule,
# the claim that "the fit was not recomputed per donor", the ordering-is-not-elapsed-time caveat
# and the counts. The Lymph Node column was measured at ~470 words of the same shape and trimmed
# to a single line with a modal behind it; this follows that. Nothing is retracted -- every
# statement still holds and is one click away. What changes is that STANDING text is stated once
# and MEASURED state stays inline, because measured state is a fact about what the reader is
# looking at and standing text is not.
#
# One of those five sentences is also simply gone, because it became false on 2026-09-17: the fit
# IS recomputed per donor now. That is the whole point of the rebuild.
.s2_trajectory_methods_modal <- function() {
  shiny::modalDialog(
    title = "How this ordering is measured", easyClose = TRUE, size = "l",
    footer = shiny::modalButton("Close"),
    shiny::tags$p(shiny::tags$strong("One fit per donor."), " Each donor's cells were ordered on ",
      "that donor's own diffusion map — destiny → slingshot, one fit per donor, no pooled ",
      "embedding. The stored value is that donor's own quantile, so two donors' pseudotimes are ",
      "not on one axis and nothing here sorts or compares donors."),
    shiny::tags$p(shiny::tags$strong("Why within donor."), " Acquisition batch is aliased with ",
      "disease status at the coarse level in this cohort, so a pooled ordering measures ",
      "acquisition and is then read as progression. Fitting inside each donor removes the ",
      "question rather than correcting for it."),
    shiny::tags$p(shiny::tags$strong("Peri minus core."), " The quantity plotted is a difference ",
      "taken inside ONE islet — where its peri-islet expansion sits on the ordering relative ",
      "to its core. A within-islet difference needs no shared axis to mean something, which is ",
      "why it is the quantity carried across donors. The peri-islet expansion is cells at ",
      "0 < distance ≤ 20 µm, derived from the signed distance and not the drawn annotation."),
    shiny::tags$p(shiny::tags$strong("A sample, and it is stated."), " The ordering was fit on a ",
      "per-donor subsample, not on every cell. An islet's value summarises the cells that were ",
      "ordered in it, and the count travels with it. A region with fewer than five ordered cells ",
      "has no value — a median of one or two cells is not a summary of anything."),
    shiny::tags$p(shiny::tags$strong("The sign is a convention."), " Each donor's axis is ",
      "oriented so its core end is low. The cohort-level SIGN of peri-minus-core follows from ",
      "that convention; the SPREAD across islets does not."),
    shiny::tags$p(shiny::tags$strong("Not elapsed time."), " Pseudotime is an ordering. It is not ",
      "time, and it is not evidence of disease progression."),
    shiny::tags$p(shiny::tags$strong("No batch η² is printed here, deliberately."), " ",
      "Every donor is rank-uniform on its own axis, so a pooled batch η² over these ",
      "values is arithmetic about the construction rather than a measurement of batch — and ",
      "a zero printed beside the word \"batch\" would read as \"batch was removed\". The ",
      "cross-donor question is asked separately, on the relation between the twenty orderings, ",
      "against a permutation null over donors.")
  )
}

#' ONE line, measured, with the standing method a click away.
.s2_trajectory_claim <- function(data, modal_link) {
  row <- data[1L, , drop = FALSE]
  number <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  sampling <- if (identical(row$sampling_state, "measured")) {
    sprintf("fit within donor on %s of %s cells", number(row$n_ordered_cells),
            number(row$n_population_cells))
  } else {
    sprintf("fit within donor; sampling %s (%s)", row$sampling_state, row$sampling_reason)
  }
  shiny::tags$p(
    class = "ie-card-note", style = "font-size: 12px",
    sprintf("%s · %s of %s islets ordered · %s with a paired peri-minus-core difference",
            sampling,
            number(sum(data$state == "measured")),
            number(nrow(data)),
            number(sum(data$delta_state == "measured"))),
    shiny::tags$span(" · "),
    shiny::actionLink(modal_link, "how this is measured →", style = "font-size: 12px")
  )
}

#' The islets this donor has that could NOT be ordered, by the reason they carry.
#'
#' They used to leave the figure silently: the plot filters to measured pairs and nothing said how
#' many rows that dropped or why. An islet that cannot be ordered is not a blank and not a zero --
#' it is a named state, and the artifact already carries the name.
.s2_trajectory_unavailable_note <- function(data) {
  missing <- data[data$delta_state != "measured", , drop = FALSE]
  if (!nrow(missing)) return(NULL)
  reasons <- missing$delta_reason
  reasons[is.na(reasons) | !nzchar(reasons)] <- "reason_not_recorded"
  # The artifact names the shortfall per region and per cell count. Rolled up to the region here,
  # because "n_sampled=3" and "n_sampled=4" are the same fact about the same floor.
  rolled <- ifelse(grepl("no_sampled_core_cells", reasons) &
                     grepl("peri_pseudotime_not_reported|no_sampled_peri_cells", reasons),
                   "no ordered cell in either region",
            ifelse(grepl("no_sampled_core_cells", reasons), "no ordered cell in the islet core",
            ifelse(grepl("core_pseudotime_not_reported", reasons) &
                     grepl("peri_pseudotime_not_reported|no_sampled_peri_cells", reasons),
                   "below the 5-ordered-cell floor in both regions",
            ifelse(grepl("core_pseudotime_not_reported", reasons),
                   "below the 5-ordered-cell floor in the core",
            ifelse(grepl("peri_pseudotime_not_reported|no_sampled_peri_cells", reasons),
                   "below the 5-ordered-cell floor in the peri-islet expansion",
                   reasons)))))
  counts <- sort(table(rolled), decreasing = TRUE)
  shiny::tags$p(
    class = "ie-card-note", role = "status", style = "font-size: 12px",
    sprintf("unavailable(%s of %s islets have no paired difference): %s. They are not zero and not absent data.",
            format(nrow(missing), big.mark = ","), format(nrow(data), big.mark = ","),
            paste(sprintf("%s — %s",
                          format(as.integer(counts), big.mark = ",", trim = TRUE),
                          names(counts)), collapse = "; "))
  )
}

#' The fuller caption, for the downloaded PNG only. A figure that travels alone carries its method
#' with it; the on-screen card does not, because the modal is one click away there and the PNG has
#' no modal.
.s2_trajectory_plot_caption <- function(data) {
  row <- data[1L, , drop = FALSE]
  number <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  sampling <- if (identical(row$sampling_state, "measured")) {
    sprintf("Ordering was fit within donor — one destiny/slingshot fit per donor, no pooled embedding — on %s of %s cells.",
            number(row$n_ordered_cells), number(row$n_population_cells))
  } else {
    sprintf("Ordering was fit within donor. Sampling %s: %s.", row$sampling_state, row$sampling_reason)
  }
  paste(sampling,
        "Pseudotime is this donor's own quantile; it is not comparable to another donor's.",
        "Islets are ordered by core median within the selected donor; ties share a rank.",
        "Peri minus core is a difference taken inside one islet, so it needs no shared axis.",
        "Pseudotime is an ordering, not elapsed time or evidence of disease progression.",
        sprintf("%s of %s islets have a core ordering; %s have a measured peri-minus-core difference.",
                number(sum(data$state == "measured")), number(nrow(data)),
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
                                      layers_dir = Sys.getenv("ISLET_LAYERS_DIR", "/data/islet_layers"),
                                      active_tab = NULL, tab_key = NULL, ...) {
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                             identical(active_tab(), tab_key))
  # The standing methods text lives in one modal rather than above every render. Registered here
  # because the shared phenotype-rules loop is keyed by id and does not know about this card.
  shiny::observeEvent(input$s2_show_trajectory_methods, {
    shiny::showModal(.s2_trajectory_methods_modal())
  }, ignoreInit = TRUE)
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
    shiny::tagList(
      .s2_trajectory_claim(data, session$ns("s2_show_trajectory_methods")),
      .s2_trajectory_unavailable_note(data)
    )
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
                    caption = paste(strwrap(.s2_trajectory_plot_caption(chosen$data), width = 100),
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
  invisible(list(active = active, data = selected, plot = plot))
}
