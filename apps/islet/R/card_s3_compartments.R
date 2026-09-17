# Islet core, the 20 µm peri-islet expansion, and other pancreas retain all
# eligible cells from the same small distance artifact. pLN is a separate
# tissue and remains unavailable until increment 9 exports its measurements.
# No spatial relation crosses tissues.

card_s3_compartments_ui <- function(ns) {
  shiny::div(class = "ie-card", style = "overflow: visible;",
             shiny::tags$h3(class = "ie-card-title", "Immune cells: pancreas compartments and pLN"),
             shiny::uiOutput(ns("s3_compartments_controls")),
             shiny::fluidRow(
               shiny::column(8, shiny::uiOutput(ns("s3_compartments_body"))),
               shiny::column(4, shiny::tags$h4("Pancreatic lymph node (pLN)"),
                             shiny::uiOutput(ns("s3_compartments_pln")))
             ))
}

card_s3_compartments_server <- function(input, output, session, prepared = NULL,
                                       cohort = NULL, active_tab = NULL, tab_key = NULL, ...) {
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                             identical(active_tab(), tab_key))
  layer <- shiny::reactive({
    shiny::req(active())
    .ie_distance_read()
  })
  donors <- shiny::reactive(.ie_distance_donors(layer(), cohort))
  selected <- shiny::reactive({
    donor <- input$s3_compartments_donor
    if (is.null(donor) || !donor %in% donors()) donor <- head(donors(), 1L)
    pair <- input$s3_compartments_pair
    if (is.null(pair)) pair <- "Immune | all"
    .ie_distance_select(layer(), donor, pair)
  })
  output$s3_compartments_controls <- shiny::renderUI({
    .ie_distance_controls(session$ns, "s3_compartments", layer(), donors(),
                          input$s3_compartments_donor, input$s3_compartments_pair)
  })
  output$s3_compartments_pln <- shiny::renderUI({
    shiny::req(active())
    .ie_distance_note(ie_avail("unavailable", "guard_pln_export_absent"),
                      "guard_pln_export_absent", artifact = "pln/",
                      builder = "scripts/senior/build_pln_layer.py (increment 9)")
  })
  compartment_data <- shiny::reactive({
    value <- selected()
    shiny::req(ie_is_measured(value$value))
    source <- .ie_distance_framed(value$data)
    reference <- ie_attach_frame(data.frame(), ie_frame("signed_distance_um"))
    ie_assert_same_frame(source, reference, "pancreas compartment comparison")
    if (!nrow(source)) return(ie_avail("unavailable", "no_eligible_pancreas_cells"))
    # Sum counts, never bin fractions. Each selected subtype/provenance pair
    # repeats the same eligible denominator, so selecting one avoids duplication.
    data <- stats::aggregate(cbind(n_cells, n_immune) ~ donor_id + donor_status + compartment,
                             data = source, FUN = sum)
    data$immune_fraction <- data$n_immune / data$n_cells
    data$compartment <- factor(data$compartment,
                                levels = c("islet_core", "peri_islet", "other_pancreas"),
                                labels = c("Islet core", "Peri-islet expansion", "Other pancreas"))
    ie_attach_frame(data, ie_frame("none"))
  })
  output$s3_compartments_body <- shiny::renderUI({
    value <- selected()
    if (!ie_is_measured(value$value)) return(.ie_distance_note(value$value, "guard_s3_compartments"))
    data <- compartment_data()
    if (inherits(data, "ie_avail")) return(.ie_distance_note(data, "guard_s3_compartments"))
    missing <- setdiff(c("Islet core", "Peri-islet expansion", "Other pancreas"),
                       as.character(data$compartment))
    notes <- lapply(missing, function(name) {
      .ie_distance_note(ie_avail("unavailable", paste0("no_eligible_cells_in_", name)),
                        "guard_s3_compartment_empty")
    })
    shiny::tagList(
      notes, shiny::plotOutput(session$ns("s3_compartments_plot")),
      shiny::tags$p("Core membership follows the exported islet outline. ",
                    "The peri-islet expansion contains non-core cells at 0 < distance ≤ 20 µm. ",
                    "Other pancreas contains all remaining eligible pancreas cells. ",
                    "Fractions use immune counts divided by all eligible, non-excluded cells; ",
                    "bin fractions are never averaged. Resolved by: ",
                    if (length(value$provenance)) paste(value$provenance, collapse = ", ")
                    else "no immune calls observed", "."),
      shiny::tableOutput(session$ns("s3_compartments_counts")),
      shiny::downloadButton(session$ns("s3_compartments_download"), "Download plot")
    )
  })
  plot <- shiny::reactive({
    data <- compartment_data()
    shiny::req(is.data.frame(data))
    reference <- ie_attach_frame(data.frame(), ie_frame("none"))
    ie_assert_same_frame(data, reference, "pancreas compartment fraction plot")
    data$count_label <- paste0(format(data$n_immune, big.mark = ",", trim = TRUE), " / ",
                               format(data$n_cells, big.mark = ",", trim = TRUE))
    ggplot2::ggplot(data, ggplot2::aes(compartment, immune_fraction)) +
      ggplot2::geom_col(fill = "#3A5A8C", width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = count_label), vjust = -0.4) +
      ggplot2::scale_y_continuous(labels = scales::label_percent(),
                                 expand = ggplot2::expansion(mult = c(0, 0.15))) +
      ggplot2::labs(x = NULL, y = "Immune fraction of eligible cells",
                    title = paste("Donor", data$donor_id[[1L]], "—", data$donor_status[[1L]]),
                    subtitle = .ie_distance_pair(selected()$data)[[1L]],
                    caption = "Counts: immune / eligible cells. pLN: unavailable(guard_pln_export_absent).") +
      .ie_distance_theme()
  })
  output$s3_compartments_plot <- shiny::renderPlot(plot())
  output$s3_compartments_counts <- shiny::renderTable({
    data <- compartment_data()
    shiny::req(is.data.frame(data))
    data.frame(Compartment = data$compartment, `Eligible cells` = data$n_cells,
               `Immune cells` = data$n_immune, check.names = FALSE)
  }, digits = 0, rownames = FALSE)
  output$s3_compartments_download <- shiny::downloadHandler(
    filename = function() "pancreas-compartment-immune-fractions.pdf",
    content = function(file) ggplot2::ggsave(file, plot = plot(), width = 8, height = 5)
  )
  invisible(list(active = active, selected = selected, plot = plot, data = compartment_data))
}
