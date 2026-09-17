# ---------------------------------------------------------------------------
# The compartment gradient artifact: islet core, the 20 um peri-islet expansion,
# the remaining exocrine tissue, and the pLN.
#
# Built by scripts/senior/build_compartment_gradient.py into
# <ISLET_LAYERS_DIR>/compartments. The three pancreas compartments are spatial
# and come from one signed distance on one section. The pLN compartment is
# COMPOSITIONAL ONLY: a lymph node is a different block with no islet in it, so
# a signed distance to one is undefined there rather than zero, and this column
# joins to the pancreas panel on donor and disease status alone.
#
# The artifact records the domain it was built from. Five donors appear in both
# domains under different vocabularies, so a matching donor id is not a licence
# to read: a mismatch is refused by name rather than rendered.
IE_COMPARTMENT_ORDER <- c("core", "peri", "exocrine")
IE_COMPARTMENT_LABEL <- c(core = "Islet core", peri = "Peri-islet expansion",
                          exocrine = "Exocrine tissue")
IE_COMPARTMENT_LINEAGES <- c("T_NK", "B_Plasma", "Myeloid", "Granulocyte")

.ie_compartment_dir <- function() {
  file.path(Sys.getenv("ISLET_LAYERS_DIR", unset = "/data/islet_layers"), "compartments")
}

.ie_compartment_note <- function(value, guard, artifact = "compartments/compartment_composition.csv") {
  .ie_distance_note(value, guard, artifact = artifact,
                    builder = "scripts/senior/build_compartment_gradient.py")
}

# Reads the artifact or says, by name, why it did not. Never a blank panel and
# never a silent empty frame: an absent export and a refused one are different
# statements and both are made out loud.
.ie_compartment_read <- function(dir = .ie_compartment_dir()) {
  path <- file.path(dir, "compartment_composition.csv")
  absent <- function(reason, state = "unavailable") {
    list(value = ie_avail(state, reason, source = path))
  }
  manifest_path <- file.path(dir, "manifest.json")
  if (!file.exists(path)) return(absent("compartment_gradient_export_absent"))
  if (!file.exists(manifest_path)) return(absent("compartment_gradient_manifest_absent"))
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(absent("json_reader_unavailable"))
  manifest <- tryCatch(jsonlite::fromJSON(manifest_path), error = function(e) NULL)
  if (is.null(manifest)) return(absent("compartment_gradient_manifest_unreadable", "failed"))
  recorded <- manifest[["source_domain_root"]]
  if (is.null(recorded) || !exists("domain_data_root", mode = "function")) {
    return(absent("compartment_gradient_domain_not_recorded"))
  }
  if (!identical(normalizePath(recorded, mustWork = FALSE),
                 normalizePath(domain_data_root(), mustWork = FALSE))) {
    return(absent("compartment_gradient_domain_mismatch"))
  }
  data <- tryCatch(utils::read.csv(path, colClasses = c(donor_id = "character")),
                   error = function(e) NULL)
  if (is.null(data)) return(absent("compartment_gradient_unreadable", "failed"))
  required <- c("donor_id", "tissue", "modality", "compartment", "stratum", "lineage",
                "n_cells", "n_lineage", "n_immune", "fraction_of_compartment",
                "fraction_of_immune", "state", "reason", "donor_status")
  if (!all(required %in% names(data))) return(absent("compartment_gradient_schema_invalid", "failed"))
  if (!nrow(data)) return(absent("compartment_gradient_export_empty"))
  if (anyNA(data$state) || any(!data$state %in% IE_STATES)) {
    return(absent("compartment_gradient_state_invalid", "failed"))
  }
  # A pLN row that carried a distance or a compartment of the pancreas would mean
  # a meaningless number had already been computed upstream. Refuse the artifact
  # rather than display it.
  pln <- data[data$tissue != "pancreas", , drop = FALSE]
  if (nrow(pln) && (any(pln$compartment != "pln") ||
                    any(pln$stratum != "not_applicable_no_islets"))) {
    return(absent("compartment_gradient_pln_carries_pancreas_geometry", "failed"))
  }
  if (any(c("distance_bin_um", "signed_dist_islet_um") %in% names(data))) {
    return(absent("compartment_gradient_pln_carries_distance", "failed"))
  }
  measured <- data[data$state == "measured", , drop = FALSE]
  if (nrow(measured) && any(!is.finite(measured$n_cells) | measured$n_cells <= 0 |
                            !is.finite(measured$n_lineage) | measured$n_lineage < 0 |
                            measured$n_lineage > measured$n_cells)) {
    return(absent("compartment_gradient_counts_invalid", "failed"))
  }
  unresolved <- data[data$state != "measured", , drop = FALSE]
  if (nrow(unresolved) && (anyNA(unresolved$reason) || any(!nzchar(unresolved$reason)))) {
    return(absent("compartment_gradient_reason_missing", "failed"))
  }
  strata <- .ie_compartment_read_csv(file.path(dir, "islet_size_summary.csv"))
  eta <- .ie_compartment_read_csv(file.path(dir, "batch_eta_squared.csv"))
  list(value = ie_avail("measured", "compartment_gradient_export_present", source = path),
       data = data, strata = strata, eta = eta, dir = dir)
}

# donor_id is read as text wherever it exists: 115 is a donor name, not a
# number, and a numeric read would silently rename it.
.ie_compartment_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch({
    header <- names(utils::read.csv(path, nrows = 1L))
    classes <- if ("donor_id" %in% header) c(donor_id = "character") else NA
    utils::read.csv(path, colClasses = classes)
  }, error = function(e) NULL)
}

.ie_compartment_donors <- function(layer, cohort = NULL) {
  if (!ie_is_measured(layer$value)) return(character())
  donors <- unique(layer$data$donor_id[layer$data$tissue == "pancreas"])
  donors <- .ie_distance_norm(donors)
  if (!is.null(cohort)) {
    allowed <- .ie_distance_norm(as.character(cohort))
    donors <- donors[donors %in% allowed]
  }
  sort(donors)
}

# The measured fragmentation of one donor's islets, said plainly.
#
# 6539's islets are mostly one-cell classifier fragments; it is not the only
# donor whose median islet is below the size at which an islet-level statement
# is defensible. The flag is measured, not curated, and the donor is shown
# rather than dropped -- dropping it would be a number chosen for being tidy.
.ie_compartment_stratum_note <- function(layer, donor, min_cells = 10L) {
  strata <- layer$strata
  if (is.null(strata) || !nrow(strata)) return(NULL)
  row <- strata[.ie_distance_norm(strata$donor_id) == donor, , drop = FALSE]
  if (!nrow(row)) return(NULL)
  ge <- grep("^n_islets_ge_", names(row), value = TRUE)[[1L]]
  flagged <- nzchar(as.character(row$flag[[1L]]))
  shiny::tags$p(
    class = "ie-card-note", role = "status",
    sprintf(paste0("Islet-size stratification: %s islets, median %s core cells, ",
                   "%s of at least %d core cells, %.1f%% of islets a single cell.%s"),
            format(row$n_islets[[1L]], big.mark = ","),
            format(row$median_core_cells_per_islet[[1L]]),
            format(row[[ge]][[1L]], big.mark = ","), min_cells,
            100 * row$fraction_islets_one_cell[[1L]],
            if (flagged) paste0(" flagged(", row$flag[[1L]], "): this donor's islets are ",
                                "mostly classifier fragments, so no per-islet number ",
                                "from it may be pooled with the cohort untouched.")
            else ""))
}

# Batch eta-squared next to its own permutation null, in the panel and not only
# in a CSV. With 20 donors over 13 acquisition runs a large eta-squared is the
# expected value under no effect at all, which is exactly what the null median
# says -- so both are shown, and neither gates anything.
.ie_compartment_batch_note <- function(layer, compartment, stratum = "all_islets",
                                       lineage = "Immune_total") {
  eta <- layer$eta
  if (is.null(eta) || !nrow(eta)) {
    return(.ie_compartment_note(ie_avail("unavailable", "batch_eta_squared_not_exported"),
                                "guard_s3_batch", artifact = "compartments/batch_eta_squared.csv"))
  }
  rows <- eta[eta$compartment == compartment & eta$stratum == stratum &
                eta$lineage == lineage & eta$grouping %in% c("stain_group", "acq_run"), ,
              drop = FALSE]
  if (!nrow(rows)) {
    return(.ie_compartment_note(ie_avail("unavailable",
                                         paste0("batch_eta_squared_absent_for_", compartment)),
                                "guard_s3_batch", artifact = "compartments/batch_eta_squared.csv"))
  }
  text <- paste(vapply(seq_len(nrow(rows)), function(i) {
    sprintf("%s eta-squared %.3f (null median %.3f, p = %.3f, %d permutations, %d donors)",
            rows$grouping[[i]], rows$eta_squared[[i]], rows$null_median[[i]],
            rows$p_permutation[[i]], rows$n_perm[[i]], rows$n_donors_used[[i]])
  }, character(1)), collapse = "; ")
  shiny::tags$p(class = "ie-card-note",
                paste0("Batch, on the finer variables only — the coarse batch is aliased ",
                       "with disease status and is not corrected on: ", text,
                       ". Nothing on this card is gated on these numbers."))
}

.ie_compartment_rows <- function(layer, donor, tissue, compartment, modality = NULL,
                                 stratum = NULL) {
  data <- layer$data
  rows <- data[.ie_distance_norm(data$donor_id) == donor & data$tissue == tissue &
                 data$compartment == compartment, , drop = FALSE]
  if (!is.null(modality)) rows <- rows[rows$modality == modality, , drop = FALSE]
  if (!is.null(stratum)) rows <- rows[rows$stratum == stratum, , drop = FALSE]
  rows
}

.ie_compartment_composition_table <- function(rows, caption) {
  totals <- rows[rows$lineage == "Immune_total", , drop = FALSE]
  lineages <- rows[rows$lineage %in% IE_COMPARTMENT_LINEAGES, , drop = FALSE]
  shiny::tagList(
    shiny::tags$table(
      class = "ie-pln-composition",
      shiny::tags$tr(lapply(c("Lineage", "Cells", "Compartment cells",
                              "% of compartment", "% of immune"), shiny::tags$th)),
      lapply(seq_len(nrow(lineages)), function(i) {
        row <- lineages[i, , drop = FALSE]
        shiny::tags$tr(
          shiny::tags$td(row$lineage),
          shiny::tags$td(format(row$n_lineage, big.mark = ",", trim = TRUE)),
          shiny::tags$td(format(row$n_cells, big.mark = ",", trim = TRUE)),
          shiny::tags$td(sprintf("%.3f%%", 100 * row$fraction_of_compartment)),
          shiny::tags$td(if (is.na(row$fraction_of_immune)) "unavailable(no immune cell in compartment)"
                         else sprintf("%.1f%%", 100 * row$fraction_of_immune)))
      }),
      if (nrow(totals)) shiny::tags$tr(
        shiny::tags$td(shiny::tags$strong("Immune total")),
        shiny::tags$td(format(totals$n_immune[[1L]], big.mark = ",", trim = TRUE)),
        shiny::tags$td(format(totals$n_cells[[1L]], big.mark = ",", trim = TRUE)),
        shiny::tags$td(sprintf("%.3f%%", 100 * totals$fraction_of_compartment[[1L]])),
        shiny::tags$td("100%")) else NULL
    ),
    shiny::tags$p(class = "ie-card-note", caption))
}

# The pLN column of this card, drawn from the compartment gradient artifact.
# It renders the compositional panel when calls are present and
# unavailable(<measured reason>) when they are not -- never a blank panel, and
# never a distance axis.
#' The pLN boundary claim is emitted ONCE by the caller that composes this column, not here --
#' this panel and `.ie_distance_pln_panel` are stacked together and each prepending it rendered the
#' identical paragraph twice, which reads as a glitch and undermines the guarantee it is making.
.ie_compartment_pln_panel <- function(layer, donor) {
  claim <- NULL
  if (!ie_is_measured(layer$value)) {
    return(shiny::tagList(claim, .ie_compartment_note(layer$value, "guard_s3_pln_composition")))
  }
  rows <- .ie_compartment_rows(layer, donor, "pancreatic_lymph_node", "pln")
  if (!nrow(rows)) {
    return(shiny::tagList(claim, .ie_compartment_note(
      ie_avail("unavailable", "pln_section_not_typed_for_donor"),
      "guard_s3_pln_composition")))
  }
  modalities <- unique(rows$modality)
  panels <- lapply(modalities, function(modality) {
    sub <- rows[rows$modality == modality, , drop = FALSE]
    absent <- sub[sub$state != "measured", , drop = FALSE]
    if (!nrow(sub[sub$state == "measured", , drop = FALSE])) {
      reason <- if (nrow(absent)) absent$reason[[1L]] else "pln_composition_not_measured"
      state <- if (nrow(absent)) absent$state[[1L]] else "unavailable"
      return(shiny::tagList(
        shiny::tags$h5(modality),
        .ie_compartment_note(ie_avail(state, reason), "guard_s3_pln_composition")))
    }
    shiny::tagList(
      shiny::tags$h5(modality),
      .ie_compartment_composition_table(
        sub[sub$state == "measured", , drop = FALSE],
        "Whole-node composition"),
      if (nrow(absent)) lapply(seq_len(nrow(absent)), function(i) {
        .ie_compartment_note(ie_avail(absent$state[[i]], absent$reason[[i]]),
                             "guard_s3_pln_composition")
      }) else NULL)
  })
  shiny::tagList(claim, panels,
                 .ie_compartment_batch_note(layer, "pln", "not_applicable_no_islets"))
}

# The pancreas half of the card: the three compartments and their immune
# composition, each shown twice -- over every islet, and over the islets large
# enough for an islet-level statement.
.ie_compartment_gradient_panel <- function(layer, donor, min_cells = 10L) {
  if (!ie_is_measured(layer$value)) {
    return(.ie_compartment_note(layer$value, "guard_s3_compartment_gradient"))
  }
  stratum_levels <- unique(layer$data$stratum[layer$data$tissue == "pancreas"])
  restricted <- setdiff(stratum_levels, "all_islets")
  panels <- lapply(IE_COMPARTMENT_ORDER, function(compartment) {
    all_rows <- .ie_compartment_rows(layer, donor, "pancreas", compartment,
                                     stratum = "all_islets")
    if (!nrow(all_rows)) {
      return(shiny::tagList(
        shiny::tags$h5(IE_COMPARTMENT_LABEL[[compartment]]),
        .ie_compartment_note(
          ie_avail("unavailable", paste0("no_rows_for_", compartment, "_for_this_donor")),
          "guard_s3_compartment_gradient")))
    }
    measured <- all_rows[all_rows$state == "measured", , drop = FALSE]
    if (!nrow(measured)) {
      return(shiny::tagList(
        shiny::tags$h5(IE_COMPARTMENT_LABEL[[compartment]]),
        .ie_compartment_note(ie_avail(all_rows$state[[1L]], all_rows$reason[[1L]]),
                             "guard_s3_compartment_gradient")))
    }
    strat_rows <- if (length(restricted))
      .ie_compartment_rows(layer, donor, "pancreas", compartment,
                           stratum = restricted[[1L]]) else all_rows[0, , drop = FALSE]
    shiny::tagList(
      shiny::tags$h5(IE_COMPARTMENT_LABEL[[compartment]]),
      .ie_compartment_composition_table(
        measured,
        paste0("Every islet. Compartment cut from the signed distance to the islet ",
               "outline, never from the drawn annotation and never from QuPath's ",
               "Parent column: one donor's drawn expansion measures 23.1 µm against ",
               "18.8–19.3 µm for the rest, and every one of them is named _exp20um.")),
      if (nrow(strat_rows[strat_rows$state == "measured", , drop = FALSE]))
        .ie_compartment_composition_table(
          strat_rows[strat_rows$state == "measured", , drop = FALSE],
          sprintf(paste0("Restricted to islets of at least %d core cells. The endocrine ",
                         "shape signal is real only above roughly this size, so the ",
                         "pooled row above and this one are different claims and are ",
                         "shown as two."), min_cells))
      else if (nrow(strat_rows))
        .ie_compartment_note(ie_avail(strat_rows$state[[1L]], strat_rows$reason[[1L]]),
                             "guard_s3_compartment_stratum")
      else NULL,
      .ie_compartment_batch_note(layer, compartment))
  })
  shiny::tagList(
    shiny::tags$p("Three pancreas compartments cut from one signed distance on one ",
                  "section: the islet core, the 20 µm peri-islet expansion outside it, ",
                  "and the remaining exocrine tissue. Fractions are immune counts over ",
                  "all eligible, non-excluded cells of the compartment; a cell with no ",
                  "nucleus never entered the denominator."),
    .ie_compartment_stratum_note(layer, donor, min_cells),
    panels)
}

# Islet core, the 20 µm peri-islet expansion, and other pancreas retain all
# eligible cells from the same small distance artifact.
#
# The pLN column reads the roi != "panc" rows of that same artifact and is
# COMPOSITIONAL, never spatial. The pLN section on a slide is cut from a
# different block than the pancreas beside it: the two share a pixel grid and
# share nothing anatomically, and a lymph node has no islet to measure a signed
# distance from. So the only join this column makes to the pancreas panel is on
# donor_id and donor_status. It renders either an explicit compositional panel
# or a named unavailable(<reason>) taken from the artifact — never a blank
# panel, and never a distance axis.

#' The standing methods text for this column, in ONE place.
#'
#' It was inline, as a paragraph above every panel on every render: the boundary claim, the
#' counts-not-a-gradient note, the compare-within-one-modality note, the identity-core-does-not-sum
#' note and the membership definitions. Measured at ~470 words of prose on one tab, identical every
#' time, which buries the numbers the card exists to show. None of it is retracted -- every
#' guarantee still holds and is still one click away. What changes is that standing text is stated
#' once and MEASURED state stays inline, because measured state is a fact about what the reader is
#' looking at and standing text is not.
.ie_pln_methods_modal <- function() {
  shiny::modalDialog(
    title = "How this column is measured", easyClose = TRUE, size = "l",
    footer = shiny::modalButton("Close"),
    shiny::tags$p(shiny::tags$strong("Compositional only."), " The pLN section comes from a ",
      "different block than the pancreas on the same slide, so this column joins to the pancreas ",
      "panel on donor and disease status alone. No distance, adjacency, neighbourhood or KDE ",
      "crosses the tissue boundary, and no signed distance to an islet exists here \u2014 there are ",
      "no islets in a lymph node."),
    shiny::tags$p(shiny::tags$strong("Counts, not a gradient."), " Whole-node composition. There ",
      "is no axis to place these cells on."),
    shiny::tags$p(shiny::tags$strong("One modality at a time."), " Each modality is typed ",
      "separately from a different section. Compare within one, never across."),
    shiny::tags$p(shiny::tags$strong("Identity-core rows do not sum to the Immune total."), " An ",
      "immune cell whose identity was not resolved is still immune."),
    shiny::tags$p(shiny::tags$strong("Pancreas compartments."), " Core membership follows the ",
      "exported islet outline. The peri-islet expansion is non-core cells at 0 < distance ",
      "\u2264 20 \u00b5m, derived from the signed distance and not the drawn annotation. Other ",
      "pancreas is every remaining eligible cell. Fractions are immune counts over all cells in ",
      "that compartment.")
  )
}

#' One line, with the methods a click away. Replaces five standing paragraphs.
.ie_pln_claim <- function(modal_link = "s3_show_pln_methods") {
  shiny::tags$p(
    class = "ie-card-note", style = "font-size: 12px",
    "Compositional only \u00b7 counts, not a gradient \u00b7 one modality at a time",
    shiny::tags$span(" \u00b7 "),
    shiny::actionLink(modal_link, "how this is measured \u2192", style = "font-size: 12px")
  )
}

# The artifact's pLN sections that belong to donors outside this domain's
# cohort. They are measured and they exist; what they are not is selectable
# here, and the reason is a domain boundary. Saying so beats letting them
# vanish between a union and an intersect.
.ie_pln_outside_cohort_note <- function(layer, cohort) {
  donors <- .ie_distance_pln_outside_cohort(layer, cohort)
  if (!length(donors)) return(NULL)
  shiny::tags$p(
    class = "ie-card-note", role = "status",
    paste0("withheld(pln_donor_outside_active_domain_cohort): ", length(donors),
           " further pLN section(s) — donor ", paste(donors, collapse = ", "),
           " — are present in the artifact but not in this domain's cohort, ",
           "so they are not offered above. They are not missing data.")
  )
}

.ie_pln_table <- function(rows) {
  data.frame(
    Section = rows$roi,
    Tissue = rows$tissue,
    Modality = rows$modality,
    Status = ifelse(is.na(rows$donor_status), "not recorded", rows$donor_status),
    check.names = FALSE
  )
}

# rows are the artifact's own pLN rows for one donor; the reader has already
# refused any that carry a distance, a bin or a compartment.
#' As with `.ie_compartment_pln_panel`: the caller emits the boundary claim once.
.ie_distance_pln_panel <- function(layer, donor) {
  claim <- NULL
  if (!ie_is_measured(layer$value)) {
    return(shiny::tagList(claim, .ie_distance_note(layer$value, "guard_s3_pln")))
  }
  rows <- .ie_distance_pln_rows(layer, donor)
  if (is.null(rows)) {
    exported <- is.data.frame(layer$pln) && nrow(layer$pln)
    value <- ie_avail("unavailable", if (exported) "pln_section_not_exported_for_donor"
                      else "guard_pln_export_absent")
    return(shiny::tagList(claim, .ie_distance_note(
      value, "guard_pln_export_absent",
      artifact = "distance_profiles/data.parquet (roi != panc)",
      builder = "scripts/senior/build_distance_profiles.py")))
  }
  measured <- rows[rows$state == "measured", , drop = FALSE]
  absent <- rows[rows$state != "measured", , drop = FALSE]
  notes <- lapply(seq_len(nrow(absent)), function(i) {
    row <- absent[i, , drop = FALSE]
    shiny::tagList(
      shiny::tags$p(class = "ie-card-note",
                    paste0(row$roi, " (", row$tissue, ", ", row$modality, "): ")),
      .ie_distance_note(ie_avail(row$state, row$reason), "guard_s3_pln",
                        artifact = "distance_profiles/data.parquet (roi != panc)",
                        builder = "scripts/senior/build_distance_profiles.py")
    )
  })
  # Modality is a column, not a footnote. A donor can carry a PhenoCycler and a
  # Xenium section of the same node, typed by different lanes against different
  # panels; their counts are not interchangeable and must never be read as one
  # measurement, or added together. The denominators differ too, which is why
  # "Section cells" is repeated on every row rather than stated once.
  composition <- if (nrow(measured)) {
    shiny::tagList(
      shiny::tags$table(
        class = "ie-pln-composition",
        shiny::tags$tr(lapply(c("Modality", "Call", "Resolved by", "Cells",
                                "Section cells", "Fraction"), shiny::tags$th)),
        lapply(order(measured$modality, measured$resolved_by != "all",
                     measured$immune_subtype), function(i) {
          row <- measured[i, , drop = FALSE]
          shiny::tags$tr(
            shiny::tags$td(row$modality),
            shiny::tags$td(row$immune_subtype), shiny::tags$td(row$resolved_by),
            shiny::tags$td(format(row$n_immune, big.mark = ",", trim = TRUE)),
            shiny::tags$td(format(row$n_cells, big.mark = ",", trim = TRUE)),
            shiny::tags$td(sprintf("%.2f%%", 100 * row$immune_fraction))
          )
        })
      ),
    )
  } else NULL
  shiny::tagList(claim, composition, notes,
                 shiny::tags$table(
                   class = "ie-pln-sections",
                   shiny::tags$tr(lapply(names(.ie_pln_table(rows)), shiny::tags$th)),
                   lapply(seq_len(nrow(rows)), function(i) {
                     entry <- .ie_pln_table(rows)[i, , drop = FALSE]
                     shiny::tags$tr(lapply(as.character(unlist(entry)), shiny::tags$td))
                   })
                 ))
}

card_s3_compartments_ui <- function(ns) {
  shiny::div(class = "ie-card", style = "overflow: visible;",
             shiny::tags$h3(class = "ie-card-title", "Immune cells: pancreas compartments and pLN"),
             shiny::uiOutput(ns("s3_compartments_controls")),
             shiny::fluidRow(
               shiny::column(8, shiny::uiOutput(ns("s3_compartments_body")),
                             shiny::tags$h4("Compartment gradient — immune composition"),
                             shiny::uiOutput(ns("s3_compartments_gradient"))),
               shiny::column(4, shiny::tags$h4("Pancreatic lymph node (pLN)"),
                             shiny::uiOutput(ns("s3_compartments_pln")))
             ))
}

card_s3_compartments_server <- function(input, output, session, prepared = NULL,
                                       cohort = NULL, active_tab = NULL, tab_key = NULL, ...) {
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                             identical(active_tab(), tab_key))
  # The standing methods text lives in one modal rather than above every panel. Registered here
  # because the shared phenotype-rules loop is keyed by id and does not know about this card.
  shiny::observeEvent(input$s3_show_pln_methods, {
    shiny::showModal(.ie_pln_methods_modal())
  }, ignoreInit = TRUE)
  layer <- shiny::reactive({
    shiny::req(active())
    .ie_distance_read()
  })
  # The active domain's cohort, and only it. A pLN section belonging to a donor
  # this domain does not contain is named in the pLN column as withheld, not
  # offered here -- see .ie_distance_donors for why widening this list was a
  # no-op that only looked like it reached them.
  donors <- shiny::reactive(.ie_distance_donors(layer(), cohort))
  donor <- shiny::reactive({
    value <- input$s3_compartments_donor
    if (is.null(value) || !value %in% donors()) value <- head(donors(), 1L)
    value
  })
  selected <- shiny::reactive({
    pair <- input$s3_compartments_pair
    if (is.null(pair)) pair <- "Immune | all"
    .ie_distance_select(layer(), donor(), pair)
  })
  output$s3_compartments_controls <- shiny::renderUI({
    .ie_distance_controls(session$ns, "s3_compartments", layer(), donors(),
                          input$s3_compartments_donor, input$s3_compartments_pair)
  })
  # The compartment gradient artifact, read once per activation. It is a
  # separate export from distance_profiles and is allowed to be absent: when it
  # is, the panels below say so by name instead of going blank.
  gradient <- shiny::reactive({
    shiny::req(active())
    .ie_compartment_read()
  })
  output$s3_compartments_pln <- shiny::renderUI({
    shiny::req(active())
    # The compositional pLN panel is the compartment gradient's. The section
    # table from distance_profiles stays beneath it: it names which section was
    # measured on which modality, which the composition alone does not say.
    shiny::tagList(.ie_pln_claim(session$ns("s3_show_pln_methods")),
                   .ie_compartment_pln_panel(gradient(), donor()),
                   .ie_distance_pln_panel(layer(), donor()),
                   .ie_pln_outside_cohort_note(layer(), cohort))
  })
  output$s3_compartments_gradient <- shiny::renderUI({
    shiny::req(active())
    .ie_compartment_gradient_panel(gradient(), donor())
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
      shiny::tags$p(class = "ie-card-note", style = "font-size: 12px", "Resolved by: ",
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
                    caption = paste("Counts: immune / eligible cells. Pancreas only;",
                                    "the pLN column beside this plot is compositional",
                                    "and shares no axis with it.")) +
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
  invisible(list(active = active, selected = selected, plot = plot,
                 data = compartment_data, donor = donor, gradient = gradient))
}
