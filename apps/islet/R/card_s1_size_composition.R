# Size must be visible before composition, and each view belongs to one donor.
# Builder contract: scripts/senior/build_islet_summary.py. All reads below are
# small layer artifacts; no typing source tree is opened by the card.

.s1_note <- function(reason, state = "unavailable", artifact = NULL) {
  ie_require_measured(
    ie_avail(state, reason = reason), guard = "s1_size_composition",
    artifact = artifact, builder = "scripts/senior/build_islet_summary.py"
  )
}

.s1_note_ui <- function(note) {
  shiny::tags$p(class = "ie-card-note", role = "status", format(note))
}

card_s1_size_composition_read <- function(
    layers_dir = ie_layers_dir(),
    active_domain_root = domain_data_root()) {
  root <- file.path(layers_dir, "islet_summary")
  receipt_path <- file.path(root, "availability.parquet")
  unavailable <- function(reason, state = "unavailable") {
    list(note = .s1_note(reason, state, root))
  }
  if (!file.exists(receipt_path)) return(unavailable("islet_summary_not_built"))
  tryCatch({
    receipt <- as.data.frame(arrow::read_parquet(receipt_path))
    required <- c("capability", "state", "reason", "source_domain_root")
    if (!all(required %in% names(receipt))) stop("invalid availability artifact")
    source_roots <- unique(receipt$source_domain_root)
    if (length(source_roots) != 1L || is.na(source_roots) || !nzchar(source_roots)) {
      stop("missing or ambiguous source_domain_root")
    }
    # The legacy and run domains share identifiers, so even a donor match cannot
    # authorize displaying an artifact built from the other domain.
    canonical <- function(path) normalizePath(path, mustWork = FALSE)
    if (!identical(canonical(source_roots), canonical(active_domain_root))) {
      return(unavailable("artifact_belongs_to_a_different_domain"))
    }
    size <- receipt[receipt$capability == "islet_size", , drop = FALSE]
    if (nrow(size) != 1L) stop("expected one islet_size availability row")
    if (!identical(size$state[[1]], "measured")) {
      return(unavailable(size$reason[[1]], size$state[[1]]))
    }
    islets_path <- file.path(root, "islets.parquet")
    if (!file.exists(islets_path)) return(unavailable("islets_artifact_missing"))
    islets <- as.data.frame(arrow::read_parquet(islets_path))
    required <- c("donor_id", "islet_key", "n_cells", "size_stratum", "state", "reason",
                  "composition_state", "composition_reason", "outline_state", "outline_reason")
    if (!all(required %in% names(islets))) stop("invalid islets artifact")
    if (!nrow(islets)) return(unavailable("islets_artifact_has_no_rows"))
    if (anyNA(islets$donor_id) || any(!nzchar(islets$donor_id)) ||
        any(duplicated(islets[c("donor_id", "islet_key")]))) {
      stop("invalid or duplicate islet identifiers")
    }
    for (column in c("state", "composition_state", "outline_state")) {
      if (anyNA(islets[[column]]) || any(!islets[[column]] %in% IE_STATES)) {
        stop("invalid availability state in ", column)
      }
    }
    composition_path <- file.path(root, "composition.parquet")
    composition <- NULL
    if (file.exists(composition_path)) {
      composition <- as.data.frame(arrow::read_parquet(composition_path))
      required <- c("donor_id", "islet_key", "phenotype", "phenotype_n_cells", "n_cells",
                    "region_total", "prop", "composition_kind", "state", "frame")
      if (!all(required %in% names(composition))) stop("invalid composition artifact")
      if (anyNA(composition$composition_kind) ||
          any(composition$composition_kind != "phenotype_counts")) {
        stop("composition is not descriptive phenotype counts")
      }
      if (anyNA(composition$frame) || any(composition$frame != "none")) {
        stop("composition counts must carry the scale-free frame 'none'")
      }
    }
    list(note = NULL, islets = islets, composition = composition, root = root)
  }, error = function(e) unavailable(paste0("artifact_read_failed: ", conditionMessage(e)), "failed"))
}

card_s1_size_composition_ui <- function(ns) {
  shiny::div(
    class = "ie-card", id = ns("s1_size_composition"),
    shiny::tags$h3(class = "ie-card-title", "Islet size and phenotype composition"),
    shiny::uiOutput(ns("s1_body"))
  )
}

card_s1_size_composition_server <- function(
    input, output, session, prepared = NULL, cohort = NULL,
    layers_dir = ie_layers_dir(),
    active_domain_root = NULL, active_tab = NULL, tab_key = NULL, ...) {
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                             identical(active_tab(), tab_key))
  ns <- session$ns
  layer <- shiny::reactive({
    shiny::req(active())
    domain_root <- if (is.null(active_domain_root)) domain_data_root() else active_domain_root
    card_s1_size_composition_read(layers_dir, domain_root)
  })
  donor <- shiny::reactive({
    data <- layer()
    shiny::req(is.null(data$note))
    choices <- sort(unique(data$islets$donor_id))
    selected <- input$s1_donor
    if (is.null(selected) || !selected %in% choices) choices[[1L]] else selected
  })
  all_sizes <- shiny::reactive({
    data <- layer()
    shiny::req(is.null(data$note))
    data$islets[data$islets$donor_id == donor(), , drop = FALSE]
  })
  sizes <- shiny::reactive({
    data <- all_sizes()
    # Select only size measurements. Tissue coordinates never enter these plots.
    data <- data[data$state == "measured" & is.finite(data$n_cells),
                 c("donor_id", "islet_key", "n_cells", "size_stratum"), drop = FALSE]
    ie_attach_frame(data, ie_frame("none"))
  })
  composition <- shiny::reactive({
    data <- layer()$composition
    if (is.null(data)) return(NULL)
    data <- data[data$donor_id == donor() & data$state == "measured" &
                   is.finite(data$n_cells), , drop = FALSE]
    stratum <- input$s1_stratum
    if (is.null(stratum)) stratum <- "large"
    if (identical(stratum, "large")) data <- data[data$n_cells > 10, , drop = FALSE]
    if (identical(stratum, "small")) data <- data[data$n_cells <= 10, , drop = FALSE]
    ie_attach_frame(data, ie_frame("none"))
  })
  output$s1_body <- shiny::renderUI({
    data <- layer()
    if (!is.null(data$note)) return(.s1_note_ui(data$note))
    shiny::tagList(
      shiny::selectInput(ns("s1_donor"), "Donor (shown separately)",
                         choices = sort(unique(data$islets$donor_id)),
                         selected = shiny::isolate(donor())),
      shiny::tags$h4("Size distribution — all exported islets for this donor"),
      shiny::uiOutput(ns("s1_size_view")),
      shiny::textOutput(ns("s1_size_caption")),
      shiny::tableOutput(ns("s1_size_counts")),
      shiny::tags$h4("Composition within a size stratum"),
      shiny::radioButtons(
        ns("s1_stratum"), "Cells per islet core",
        choices = c("Above 10 cells" = "large", "0–10 cells" = "small", "All sizes" = "all"),
        selected = "large", inline = TRUE
      ),
      shiny::tags$p(class = "ie-card-note", paste(
        "Small objects can be classifier fragments; inspect the size distribution first.",
        "Endocrine shape comparisons require attention to islets above about 10 cells.",
        "These are descriptive phenotype counts. Thresholded hormone positivity is unavailable",
        "from this source: the argmax hormone table has equal any/single counts and zero",
        "multi-positive and triple-negative counts by construction."
      )),
      shiny::uiOutput(ns("s1_composition_view")),
      shiny::textOutput(ns("s1_composition_caption")),
      shiny::uiOutput(ns("s1_coverage"))
    )
  })
  output$s1_size_view <- shiny::renderUI({
    if (!nrow(sizes())) return(.s1_note_ui(.s1_note("no_measured_core_n_cells_for_donor")))
    shiny::tagList(
      shiny::plotOutput(ns("s1_size_plot"), height = "260px"),
      shiny::downloadButton(ns("s1_size_download"), "Download size plot")
    )
  })
  output$s1_size_caption <- shiny::renderText({
    data <- all_sizes()
    sprintf("Donor %s: %s exported islets; %s measured cell counts; %s unavailable. Size selection below leaves this distribution intact.",
            donor(), nrow(data), sum(data$state == "measured"), sum(data$state != "measured"))
  })
  output$s1_size_counts <- shiny::renderTable({
    data <- sizes()
    shiny::validate(shiny::need(nrow(data) > 0,
                               format(.s1_note("no_measured_core_n_cells_for_donor"))))
    levels <- c("0 cells", "1 cell", "2–10 cells", ">10 cells")
    data.frame(
      Donor = donor(), `Cells per core` = levels,
      Islets = as.integer(table(factor(data$size_stratum, levels = levels))), check.names = FALSE
    )
  }, rownames = FALSE)
  output$s1_composition_view <- shiny::renderUI({
    data <- composition()
    if (is.null(data) || !nrow(data)) {
      reason <- if (is.null(data)) "composition_artifact_unavailable" else "no_composition_in_selected_size_stratum"
      return(.s1_note_ui(.s1_note(reason)))
    }
    shiny::tagList(
      shiny::plotOutput(ns("s1_composition_plot"), height = "480px"),
      shiny::downloadButton(ns("s1_composition_download"), "Download composition plot")
    )
  })
  output$s1_composition_caption <- shiny::renderText({
    data <- composition()
    if (is.null(data) || !nrow(data)) return("No measured composition in this selection.")
    unique_islets <- data[!duplicated(data$islet_key), , drop = FALSE]
    sprintf(paste0(
      "Donor %s: %s islets, %s exported phenotype rows. Each fraction is phenotype_n_cells / region_total. ",
      "The composition denominator differs from the geometry core n_cells in %s selected islets. ",
      "Only exported phenotype rows are shown; no absent measurement is filled with zero."),
      donor(), nrow(unique_islets), nrow(data), sum(unique_islets$region_total != unique_islets$n_cells))
  })
  output$s1_coverage <- shiny::renderUI({
    data <- all_sizes()
    notes <- list()
    for (kind in c("composition", "outline")) {
      state <- data[[paste0(kind, "_state")]]
      reason <- data[[paste0(kind, "_reason")]]
      for (value in unique(reason[state != "measured"])) {
        count <- sum(reason == value & state != "measured")
        notes[[length(notes) + 1L]] <- .s1_note_ui(.s1_note(
          sprintf("%s: %s islets — %s", kind, count, value)
        ))
      }
    }
    shiny::tagList(notes)
  })
  plot_theme <- function() {
    if (exists("theme_islet", mode = "function")) theme_islet() else ggplot2::theme_minimal()
  }
  size_plot <- shiny::reactive({
    data <- sizes()
    shiny::validate(shiny::need(nrow(data) > 0, format(.s1_note("no_measured_core_n_cells_for_donor"))))
    ie_assert_same_frame(data, data, context = "S1 size distribution")
    # Stems remain visible at zero and one on a transformed axis; bar widths
    # derived from the smallest transformed interval would make them disappear.
    frequencies <- stats::aggregate(list(n_islets = rep(1L, nrow(data))),
                                    list(n_cells = data$n_cells), sum)
    ggplot2::ggplot(frequencies, ggplot2::aes(x = n_cells, y = n_islets)) +
      ggplot2::geom_segment(ggplot2::aes(xend = n_cells, yend = 0), colour = "#4477AA") +
      ggplot2::geom_point(size = 1.2, colour = "#4477AA") +
      ggplot2::geom_vline(xintercept = 10, linetype = "dashed", colour = "#555555") +
      ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                                 breaks = c(0, 1, 10, 100, 1000, 10000)) +
      ggplot2::labs(title = paste("Donor", donor()), x = "Cells per islet core (n_cells)",
                    y = "Islets (count)", caption = "Dashed line: 10 cells; zero and one remain visible.") +
      plot_theme()
  })
  composition_plot <- shiny::reactive({
    data <- composition()
    shiny::validate(shiny::need(!is.null(data) && nrow(data) > 0,
                               format(.s1_note("no_composition_in_selected_size_stratum"))))
    ie_assert_same_frame(data, data, context = "S1 size and composition")
    ggplot2::ggplot(data, ggplot2::aes(x = n_cells, y = prop)) +
      ggplot2::geom_point(alpha = 0.45, size = 1.4, colour = "#4477AA") +
      ggplot2::facet_wrap(~phenotype) +
      ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                                 breaks = c(0, 1, 10, 100, 1000, 10000)) +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(title = paste("Donor", donor()), x = "Cells per islet core (n_cells)",
                    y = "Phenotype cells / composition region_total") + plot_theme()
  })
  output$s1_size_plot <- shiny::renderPlot(size_plot())
  output$s1_composition_plot <- shiny::renderPlot(composition_plot())
  output$s1_size_download <- shiny::downloadHandler(
    filename = function() paste0("islet_size_", donor(), ".png"),
    content = function(file) ggplot2::ggsave(file, size_plot(), width = 9, height = 4, dpi = 160)
  )
  output$s1_composition_download <- shiny::downloadHandler(
    filename = function() paste0("islet_composition_", donor(), ".png"),
    content = function(file) ggplot2::ggsave(file, composition_plot(), width = 10, height = 7, dpi = 160)
  )
  invisible(list(active = active, size_plot = size_plot, composition_plot = composition_plot,
                 sizes = sizes, composition = composition, layer = layer))
}
