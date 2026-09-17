# ---------------------------------------------------------------------------
# G1 -- the four-cell grid: both tissues, both modalities, per donor.
#
# The app reads ONE of the four cells. The other three have calls and no reader,
# so this card reads the one layer that presents all four at donor level:
# <ISLET_LAYERS_DIR>/grid_composition, built by
# scripts/senior/build_grid_composition.py.
#
# WHAT THIS CARD REFUSES TO DO.
#
# * It never draws a blank panel and never draws a zero for an absence. Every
#   cell of the grid that is not measured for the selected donor renders a NAMED
#   NOTE carrying one of the five availability states and the producer's own
#   reason. A measured zero is a finding; a hole is a hole; they must not print
#   the same.
# * It never names a donor, a cohort size, a status level or a typing method in
#   a string. Every one of those is read out of the layer. A cohort fact written
#   as a literal is a domain leak the test suite cannot see -- that is how the
#   app once printed "(N=15 total)" over a 20-donor dataset with 238 assertions
#   passing.
# * It never joins the two tissues on anything but donor and disease status. A
#   lymph node is a different block with no islet in it, so the comparison here
#   is COMPOSITIONAL. No distance, adjacency, neighbourhood or KDE appears.
# * It never presents one modality as the correction of the other. They are
#   different sections of different blocks: agreement is informative and so is
#   disagreement, and the card shows the offset as it falls.
#
# The artifact records the domain it was built from. The legacy CODEX build and
# the run cohort share five donors and no vocabulary, so a matching donor id is
# not a licence to read: a mismatch is refused by name rather than rendered.

.g1_dir <- function() {
  file.path(Sys.getenv("ISLET_LAYERS_DIR", unset = "/data/islet_layers"), "grid_composition")
}

.g1_builder <- "scripts/senior/build_grid_composition.py"

.g1_note <- function(reason, state = "unavailable", artifact = NULL) {
  ie_require_measured(ie_avail(state, reason = reason), guard = "g1_grid",
                      artifact = if (is.null(artifact) || is.na(artifact)) .g1_dir() else artifact,
                      builder = .g1_builder)
}

.g1_note_ui <- function(note) shiny::tags$p(class = "ie-card-note", role = "status", format(note))

# donor_id is read as text wherever it appears: 115 is a donor name, not a
# number, and a numeric read would silently rename it.
.g1_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    header <- names(utils::read.csv(path, nrows = 1L))
    classes <- if ("donor_id" %in% header) c(donor_id = "character") else NA
    utils::read.csv(path, colClasses = classes, check.names = FALSE)
  }, error = function(e) NULL)
}

card_g1_grid_read <- function(dir = .g1_dir(), active_domain_root = NULL) {
  absent <- function(reason, state = "unavailable") list(note = .g1_note(reason, state, dir))
  manifest_path <- file.path(dir, "manifest.json")
  if (!file.exists(manifest_path)) return(absent("grid_composition_layer_not_built"))
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(absent("json_reader_unavailable"))
  manifest <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = TRUE),
                       error = function(e) NULL)
  if (is.null(manifest)) return(absent("grid_composition_manifest_unreadable", "failed"))
  recorded <- manifest[["source_domain_root"]]
  root <- if (is.null(active_domain_root)) domain_data_root() else active_domain_root
  if (is.null(recorded) || !nzchar(recorded)) return(absent("grid_composition_domain_not_recorded"))
  canonical <- function(path) normalizePath(path, mustWork = FALSE)
  if (!identical(canonical(recorded), canonical(root))) {
    return(absent("grid_composition_belongs_to_a_different_domain"))
  }
  presence <- .g1_read_csv(file.path(dir, "presence.csv"))
  composition <- .g1_read_csv(file.path(dir, "composition.csv"))
  if (is.null(presence) || is.null(composition)) return(absent("grid_composition_tables_unreadable", "failed"))
  required_presence <- c("tissue", "modality", "donor_id", "state", "reason",
                         "n_cells_assessed", "n_no_nucleus_excluded", "method", "source",
                         "citation_conflict", "disease_status")
  required_composition <- c("tissue", "modality", "donor_id", "lineage_common", "bucket_kind",
                            "n_cells", "denom_assessed", "frac_of_assessed", "state", "reason",
                            "pheno_resolution", "xenium_resolution", "is_immune",
                            "resolution_conflict", "disease_status")
  if (!all(required_presence %in% names(presence))) return(absent("presence_schema_invalid", "failed"))
  if (!all(required_composition %in% names(composition))) return(absent("composition_schema_invalid", "failed"))
  if (!nrow(presence) || !nrow(composition)) return(absent("grid_composition_layer_is_empty"))
  for (column in list(presence$state, composition$state)) {
    if (anyNA(column) || any(!column %in% IE_STATES)) return(absent("availability_state_invalid", "failed"))
  }
  # An absence without a recorded reason is indistinguishable from an oversight.
  unresolved <- presence[presence$state != "measured", , drop = FALSE]
  if (nrow(unresolved) && (anyNA(unresolved$reason) || any(!nzchar(unresolved$reason)))) {
    return(absent("presence_reason_missing", "failed"))
  }
  # A missing measurement is unavailable, never zero: a non-measured composition
  # row carrying a count would mean the rule had already been broken upstream.
  holes <- composition[composition$state != "measured", , drop = FALSE]
  if (nrow(holes) && any(!is.na(holes$n_cells))) return(absent("absence_carries_a_count", "failed"))
  # A spatial column would mean a distance had been computed across a tissue
  # boundary somewhere upstream. Refuse the artifact rather than display it.
  spatial <- c("signed_dist_islet_um", "distance_bin_um", "compartment", "neighbourhood")
  if (any(spatial %in% names(composition))) return(absent("composition_carries_spatial_columns", "failed"))
  list(
    note = NULL, dir = dir, manifest = manifest,
    presence = presence, composition = composition,
    crossmodal = .g1_read_csv(file.path(dir, "crossmodal_tnk.csv")),
    eta = .g1_read_csv(file.path(dir, "batch_eta_squared.csv"))
  )
}

# How a lineage's two resolution columns say the comparison should be READ.
#
# The pair does not say whether to compare a lineage -- nothing is dropped for
# disagreeing. It says what kind of comparison it is. A lineage the RNA side
# reaches only through a surrogate panel is ORTHOGONAL evidence against a protein
# gate, not a check of it; a lineage that exists only as a fine subtype on one
# side compares broad to fine; a lineage assigned by fallback on one side is being
# compared against a residual.
.g1_reading <- function(pheno_resolution, xenium_resolution) {
  pair <- c(pheno_resolution, xenium_resolution)
  if (any(pair == "surrogate")) return("orthogonal: one side has no 1:1 marker and uses a surrogate panel")
  if (any(pair == "fine_only")) return("broad to fine: one side reaches this only as a fine subtype")
  if (any(pair == "default_sink")) return("against a residual: one side assigns this by fallback")
  if (any(pair == "subclass")) return("subclass: one side recovers this from per-marker gates, not the broad label")
  if (all(pair == "yes")) return("like for like: both sides resolve this with a positive broad call")
  paste0("mixed: ", paste(pair, collapse = " / "))
}

card_g1_grid_ui <- function(ns) {
  shiny::div(
    class = "ie-card", id = ns("g1_grid"), style = "overflow: visible;",
    shiny::tags$h3(class = "ie-card-title", "Both tissues, both modalities — per donor"),
    shiny::uiOutput(ns("g1_body"))
  )
}

card_g1_grid_server <- function(
    input, output, session, prepared = NULL, cohort = NULL,
    layers_dir = Sys.getenv("ISLET_LAYERS_DIR", "/data/islet_layers"),
    active_domain_root = NULL, active_tab = NULL, tab_key = NULL, ...) {
  ns <- session$ns
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                              identical(active_tab(), tab_key))

  layer <- shiny::reactive({
    shiny::req(active())
    card_g1_grid_read(file.path(layers_dir, "grid_composition"), active_domain_root)
  })

  # The grid's axes come from the layer, never from a list typed here: a fifth
  # cell added upstream must appear, not be silently dropped by a hardcoded pair.
  tissues <- shiny::reactive(sort(unique(layer()$presence$tissue)))
  modalities <- shiny::reactive(sort(unique(layer()$presence$modality)))

  # The layer spans four call trees whose donor sets differ, so it holds donors
  # this domain's cohort roster does not. They are measured and they exist; what
  # they are not is selectable here. They are named rather than dropped in
  # silence -- a donor that vanishes with no note reads as a donor that was
  # never acquired.
  roster <- shiny::reactive({
    if (is.null(cohort)) return(NULL)
    value <- tryCatch(as.character(cohort), error = function(e) character())
    value <- unique(value[grepl("^[0-9]+$", value)])
    if (length(value)) value else NULL
  })

  donors <- shiny::reactive({
    data <- layer()
    shiny::req(is.null(data$note))
    all_donors <- sort(unique(data$presence$donor_id))
    allowed <- roster()
    if (is.null(allowed)) return(all_donors)
    keep <- all_donors[all_donors %in% allowed]
    if (length(keep)) keep else all_donors
  })

  outside_cohort <- shiny::reactive({
    data <- layer()
    shiny::req(is.null(data$note))
    setdiff(sort(unique(data$presence$donor_id)), donors())
  })

  donor <- shiny::reactive({
    choices <- donors()
    selected <- input$g1_donor
    if (is.null(selected) || !selected %in% choices) choices[[1L]] else selected
  })

  donor_presence <- shiny::reactive({
    data <- layer()
    data$presence[data$presence$donor_id == donor(), , drop = FALSE]
  })

  donor_composition <- shiny::reactive({
    data <- layer()
    rows <- data$composition[data$composition$donor_id == donor(), , drop = FALSE]
    rows[rows$state == "measured" & is.finite(rows$frac_of_assessed), , drop = FALSE]
  })

  output$g1_body <- shiny::renderUI({
    data <- layer()
    if (!is.null(data$note)) return(.g1_note_ui(data$note))
    shiny::tagList(
      shiny::selectInput(ns("g1_donor"), "Donor (each cell of the grid shown separately)",
                         choices = donors(), selected = shiny::isolate(donor())),
      shiny::uiOutput(ns("g1_donor_caption")),
      shiny::uiOutput(ns("g1_outside_cohort")),
      shiny::uiOutput(ns("g1_cells")),
      shiny::tags$h4("Which cells of the grid this donor has"),
      shiny::tableOutput(ns("g1_presence_table")),
      shiny::tags$h4("Cross-modality, within a tissue"),
      shiny::uiOutput(ns("g1_crossmodal")),
      shiny::tableOutput(ns("g1_lineage_crossmodal")),
      shiny::uiOutput(ns("g1_crossmodal_reading")),
      shiny::tags$h4("Acquisition batch"),
      shiny::uiOutput(ns("g1_batch"))
    )
  })

  output$g1_donor_caption <- shiny::renderUI({
    rows <- donor_presence()
    status <- unique(rows$disease_status[!is.na(rows$disease_status)])
    measured <- rows[rows$state == "measured", , drop = FALSE]
    shiny::tags$p(
      class = "ie-card-note",
      sprintf(paste0("Donor %s (%s): %d of %d cells of the grid measured, %s cells assessed. ",
                     "The two tissues are different blocks; they are joined here on donor and ",
                     "disease status alone, and nothing spatial crosses between them."),
              donor(),
              if (length(status)) paste(status, collapse = "/") else "status unavailable",
              nrow(measured), nrow(rows),
              format(sum(measured$n_cells_assessed, na.rm = TRUE), big.mark = ","))
    )
  })

  output$g1_outside_cohort <- shiny::renderUI({
    donors_out <- outside_cohort()
    if (!length(donors_out)) return(NULL)
    data <- layer()
    rows <- data$presence[data$presence$donor_id %in% donors_out &
                            data$presence$state == "measured", , drop = FALSE]
    cells <- unique(paste(rows$tissue, rows$modality, sep = "/"))
    .g1_note_ui(.g1_note(sprintf(paste0(
      "%d donors are measured in this layer and are not in this domain's cohort roster, so ",
      "they are not selectable here: %s. Between them they carry %s cells in %s. They are ",
      "named rather than dropped in silence."),
      length(donors_out), paste(donors_out, collapse = ", "),
      format(sum(rows$n_cells_assessed, na.rm = TRUE), big.mark = ","),
      paste(sort(cells), collapse = " and ")), "withheld"))
  })

  # One panel per cell of the grid: a plot when the cell is measured for this
  # donor, a named note when it is not. Never an empty panel, never a zero.
  output$g1_cells <- shiny::renderUI({
    rows <- donor_presence()
    panels <- lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, , drop = FALSE]
      id <- paste0("g1_cell_", row$tissue, "_", row$modality)
      body <- if (identical(row$state, "measured")) {
        shiny::plotOutput(ns(id), height = "300px")
      } else {
        .g1_note_ui(.g1_note(row$reason, row$state, row$source))
      }
      conflict <- row$citation_conflict
      shiny::div(
        class = "ie-card", style = "overflow: visible;",
        shiny::tags$h5(sprintf("%s · %s", row$tissue, row$modality)),
        shiny::tags$p(class = "ie-card-note", sprintf("%s — %s", row$method,
                                                      ie_state_label(row$state, row$n_cells_assessed))),
        body,
        if (!is.null(conflict) && !is.na(conflict) && nzchar(conflict))
          .g1_note_ui(.g1_note(conflict, "unavailable")) else NULL
      )
    })
    shiny::div(class = "ie-card-stack", panels)
  })

  cell_plot <- function(tissue, modality) {
    rows <- donor_composition()
    rows <- rows[rows$tissue == tissue & rows$modality == modality, , drop = FALSE]
    shiny::validate(shiny::need(nrow(rows) > 0,
                                format(.g1_note("no_measured_bucket_for_this_grid_cell"))))
    rows <- stats::aggregate(list(n_cells = rows$n_cells, frac = rows$frac_of_assessed),
                             list(lineage_common = rows$lineage_common,
                                  bucket_kind = rows$bucket_kind), sum)
    rows <- rows[order(rows$frac), , drop = FALSE]
    rows$lineage_common <- factor(rows$lineage_common, levels = rows$lineage_common)
    ggplot2::ggplot(rows, ggplot2::aes(x = lineage_common, y = frac, fill = bucket_kind)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(
        x = NULL, y = "Cells / cells this lane assessed", fill = "Bucket",
        title = sprintf("%s · %s · donor %s", tissue, modality, donor()),
        caption = "Every bucket the lane produced, including the ones that are not a lineage; they sum to the cells assessed."
      ) +
      (if (exists("theme_islet", mode = "function")) theme_islet() else ggplot2::theme_minimal())
  }

  # Renderers are registered for every cell the layer declares, so a cell added
  # upstream draws without this file changing.
  shiny::observeEvent(layer(), {
    data <- layer()
    if (!is.null(data$note)) return(NULL)
    cells <- unique(data$presence[c("tissue", "modality")])
    for (i in seq_len(nrow(cells))) {
      local({
        tissue <- cells$tissue[[i]]
        modality <- cells$modality[[i]]
        id <- paste0("g1_cell_", tissue, "_", modality)
        output[[id]] <- shiny::renderPlot(cell_plot(tissue, modality))
      })
    }
  }, ignoreNULL = TRUE)

  output$g1_presence_table <- shiny::renderTable({
    rows <- donor_presence()
    data.frame(
      Tissue = rows$tissue, Modality = rows$modality,
      State = vapply(seq_len(nrow(rows)), function(i)
        ie_state_label(rows$state[[i]], rows$n_cells_assessed[[i]]), character(1)),
      `Cells with no nucleus, excluded` = ifelse(is.na(rows$n_no_nucleus_excluded), "unavailable",
                                                 format(rows$n_no_nucleus_excluded, big.mark = ",")),
      Reason = ifelse(nzchar(rows$reason), rows$reason, "—"),
      check.names = FALSE
    )
  }, rownames = FALSE)

  output$g1_crossmodal <- shiny::renderUI({
    data <- layer()
    tab <- data$crossmodal
    if (is.null(tab) || !nrow(tab)) return(.g1_note_ui(.g1_note("crossmodal_table_not_exported")))
    notes <- lapply(sort(unique(tab$tissue)), function(tissue) {
      rows <- tab[tab$tissue == tissue & tab$both_modalities %in% c(TRUE, "True", "true"), , drop = FALSE]
      if (!nrow(rows)) {
        return(.g1_note_ui(.g1_note(sprintf(
          "%s: no donor is measured in both modalities, so no cross-modality comparison exists here",
          tissue))))
      }
      pooled <- function(frac, n) sum(frac * n, na.rm = TRUE) / sum(n, na.rm = TRUE)
      protein <- pooled(rows$frac_T_NK_protein, rows$n_resolved_immune_protein)
      rna <- pooled(rows$frac_T_NK_RNA, rows$n_resolved_immune_RNA)
      disagree <- sum(rows$rank_disagrees %in% c(TRUE, "True", "true"))
      shiny::tags$p(class = "ie-card-note", sprintf(paste0(
        "%s, %d donors measured in both: T_NK is %.1f%% of the immune cells the protein lane ",
        "resolved and %.1f%% of those the RNA lane resolved; the per-donor median offset is ",
        "%+.1f pp (RNA minus protein) and the two modalities rank the donors differently in %d ",
        "of %d. These are different sections of different blocks. The disagreement is shown as ",
        "it falls: neither side is corrected toward the other, and neither is the target."),
        tissue, nrow(rows), 100 * protein, 100 * rna,
        stats::median(rows$offset_pp_RNA_minus_protein, na.rm = TRUE),
        disagree, nrow(rows)))
    })
    shiny::tagList(notes)
  })

  # Every lineage both modalities measured in a tissue, not only the immune ones.
  # The comparison is compositional and the donors are those measured in both.
  lineage_crossmodal <- shiny::reactive({
    data <- layer()
    rows <- data$composition
    rows <- rows[rows$state == "measured" & rows$bucket_kind == "lineage" &
                   is.finite(rows$frac_of_assessed), , drop = FALSE]
    out <- list()
    for (tissue in sort(unique(rows$tissue))) {
      here <- rows[rows$tissue == tissue, , drop = FALSE]
      for (lineage in sort(unique(here$lineage_common))) {
        sub <- here[here$lineage_common == lineage, , drop = FALSE]
        protein <- sub[sub$modality == "protein", , drop = FALSE]
        rna <- sub[sub$modality == "RNA", , drop = FALSE]
        shared <- intersect(protein$donor_id, rna$donor_id)
        if (!length(shared)) next
        protein <- protein[match(shared, protein$donor_id), , drop = FALSE]
        rna <- rna[match(shared, rna$donor_id), , drop = FALSE]
        out[[length(out) + 1L]] <- data.frame(
          Tissue = tissue, Lineage = lineage, `Donors in both` = length(shared),
          `Protein %` = 100 * sum(protein$n_cells) / sum(protein$denom_assessed),
          `RNA %` = 100 * sum(rna$n_cells) / sum(rna$denom_assessed),
          `Median offset (pp)` = stats::median(100 * (rna$frac_of_assessed - protein$frac_of_assessed)),
          `How to read it` = .g1_reading(protein$pheno_resolution[[1L]], rna$xenium_resolution[[1L]]),
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }
    }
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  })

  output$g1_lineage_crossmodal <- shiny::renderTable({
    table <- lineage_crossmodal()
    shiny::validate(shiny::need(!is.null(table),
                                format(.g1_note("no_lineage_is_measured_in_both_modalities"))))
    table
  }, rownames = FALSE, digits = 2)

  output$g1_crossmodal_reading <- shiny::renderUI({
    table <- lineage_crossmodal()
    if (is.null(table)) return(NULL)
    kinds <- sort(unique(table$`How to read it`))
    shiny::tags$p(class = "ie-card-note", sprintf(paste0(
      "The two modalities are different sections of different blocks, so a difference here is ",
      "not an error on either side. The crosswalk's resolution pair says what kind of comparison ",
      "each row is, and none is dropped for disagreeing: %s. Percentages are of the cells each ",
      "lane assessed, over the %d donors measured in both."),
      paste(kinds, collapse = "; "),
      max(table$`Donors in both`)))
  })

  output$g1_batch <- shiny::renderUI({
    eta <- layer()$eta
    if (is.null(eta) || !nrow(eta)) return(.g1_note_ui(.g1_note("batch_eta_squared_not_exported")))
    measured <- eta[eta$state == "measured" & is.finite(eta$eta_squared), , drop = FALSE]
    if (!nrow(measured)) {
      reasons <- unique(eta$reason[nzchar(eta$reason)])
      return(shiny::tagList(lapply(reasons, function(r) .g1_note_ui(.g1_note(r, "blank")))))
    }
    variables <- sort(unique(measured$variable))
    lines <- lapply(variables, function(variable) {
      rows <- measured[measured$variable == variable, , drop = FALSE]
      shiny::tags$p(class = "ie-card-note", sprintf(paste0(
        "%s: median eta-squared %.3f over %d cross-donor quantities, against a permutation ",
        "null whose median is %.3f; %d of %d reach a permutation p below 0.05. Coarse batch is ",
        "aliased with disease status, so a large eta-squared beside a large null median is what ",
        "no effect looks like at this donor count. Nothing here is corrected on it."),
        variable, stats::median(rows$eta_squared), nrow(rows),
        stats::median(rows$null_median, na.rm = TRUE),
        sum(rows$p_permutation < 0.05, na.rm = TRUE), nrow(rows)))
    })
    blanks <- unique(eta$reason[eta$state != "measured" & nzchar(eta$reason)])
    shiny::tagList(lines, lapply(blanks, function(r) .g1_note_ui(.g1_note(r, "blank"))))
  })

  invisible(list(active = active, layer = layer, donor = donor,
                 presence = donor_presence, composition = donor_composition,
                 cell_plot = cell_plot))
}
