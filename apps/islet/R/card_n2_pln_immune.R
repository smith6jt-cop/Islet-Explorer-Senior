# ---------------------------------------------------------------------------
# N2 -- pLN immune composition per donor, PhenoCycler and Xenium SIDE BY SIDE.
#
# The two lanes are never pooled and never added. They are different sections of
# different blocks, typed by different methods against different panels: one
# names immune classes from PhenoCycler clusters, the other partitions the
# Xenium pLN release. A cell counted on one side is not the same cell on the
# other, and there is no cell identifier that crosses -- the PhenoCycler
# object_id is a QuPath UUID and a Xenium cell id is unique only within its own
# bundle. The only join either side makes is on donor and disease status.
#
# Sources, all read as they stand:
#   <ISLET_LAYERS_DIR>/grid_composition  scripts/senior/build_grid_composition.py --
#       the one layer that already reads the pLN Xenium release AND the
#       PhenoCycler pLN typing, per donor, with a state and a reason for every
#       cell of the grid it cannot fill. This card reads its pLN rows only.
#   <ISLET_LAYERS_DIR>/pln               the node layer, for how many cells the
#       node annotations hold on each slide -- a different denominator from the
#       one the typing lane assessed, and shown as such.
#   data/development/integration/concordance/pln_lineage_summary.csv
#       the measured cross-modality verdict, carried verbatim when present.
#
# WHAT THIS CARD REFUSES TO DO.
#
# * It does not flatter either side and does not make them agree. The protein
#   pLN typing is under repair; whatever the layer says today is what is shown,
#   including the share of the section it resolves as not immune. Nothing here
#   is corrected toward the other modality, and neither modality is the target.
# * It never pools the two modalities into one number, one bar or one
#   denominator. Every count is labelled with the lane that assessed it.
# * It never computes a distance, adjacency or neighbourhood between the
#   pancreas and the lymph node. There are no islets in a lymph node, so a
#   signed distance to one is undefined there rather than zero.
# * It never turns an absence into a zero, and never writes a donor id, a donor
#   count or a share into a literal: all of them are read from the artifacts.

.n2_builder <- "scripts/senior/build_grid_composition.py"
.n2_tissue <- "pln"
.n2_verdict_rel <- "data/development/integration/concordance/pln_lineage_summary.csv"
.n2_verdict_manifest_rel <- "data/development/integration/concordance/manifest.json"

.n2_layers_dir <- function() ie_layers_dir()

.n2_note <- function(reason, state = "unavailable", artifact = NULL, builder = .n2_builder) {
  ie_require_measured(ie_avail(state, reason = reason), guard = "n2_pln_immune",
                      artifact = if (is.null(artifact) || is.na(artifact)) .n2_layers_dir() else artifact,
                      builder = builder)
}

.n2_note_ui <- function(note) shiny::tags$p(class = "ie-card-note", role = "status", format(note))

# A MEASURED statement. ie_require_measured() correctly returns nothing for a
# measured value, so a measurement says itself rather than being dressed as an
# absence.
.n2_fact <- function(text) shiny::tags$p(class = "ie-card-note", role = "status", text)

# The grid layer's contract is enforced in ONE place, card_g1_grid.R's reader:
# both cards must refuse the same artifact for the same reasons, and a second
# copy of those checks would be a second contract that could drift. Guarded, so
# this card names the missing reader instead of failing if it is ever sourced
# alone.
.n2_read_grid <- function(layers_dir, active_domain_root = NULL) {
  dir <- file.path(layers_dir, "grid_composition")
  if (!exists("card_g1_grid_read", mode = "function")) {
    return(list(value = ie_avail("unavailable", "grid_composition_reader_not_loaded", source = dir)))
  }
  data <- card_g1_grid_read(dir, active_domain_root)
  if (!is.null(data$note)) {
    # The reader's own reason is carried through unchanged, so the two cards can
    # never disagree about why the layer was refused. Only the guard differs,
    # which is what lets a reader trace the message to the card that showed it.
    detail <- if (length(data$note$lines)) data$note$lines[[1L]] else "no reason recorded"
    parts <- strsplit(detail, " \u2014 ", fixed = TRUE)[[1L]]
    reason <- if (length(parts) > 1L) paste(parts[-1L], collapse = " \u2014 ") else detail
    return(list(value = ie_avail("unavailable", reason, source = dir)))
  }
  presence <- data$presence[data$presence$tissue == .n2_tissue, , drop = FALSE]
  composition <- data$composition[data$composition$tissue == .n2_tissue, , drop = FALSE]
  if (!nrow(presence)) {
    return(list(value = ie_avail("blank", "grid_composition_carries_no_lymph_node_rows",
                                 source = dir)))
  }
  presence$donor_id <- as.character(presence$donor_id)
  composition$donor_id <- as.character(composition$donor_id)
  list(value = ie_avail("measured", "pln_composition_present", n = nrow(presence), source = dir),
       presence = presence, composition = composition, manifest = data$manifest, dir = dir)
}

# The node layer's per-donor cell totals: how many cells the node annotations
# hold on that slide, which is NOT the number the typing lane assessed. Reuses
# card_n1_sections.R's reader so one contract governs that artifact.
.n2_read_nodes <- function(layers_dir) {
  dir <- file.path(layers_dir, "pln")
  if (!exists(".n1_read_nodes", mode = "function")) {
    return(list(value = ie_avail("unavailable", "pln_node_reader_not_loaded", source = dir)))
  }
  layer <- .n1_read_nodes(dir)
  if (!ie_is_measured(layer$value)) return(layer)
  rows <- layer$data[layer$data$state == "measured", , drop = FALSE]
  if (!nrow(rows)) {
    return(list(value = ie_avail("unavailable", "no_usable_node_geometry_in_the_node_layer",
                                 source = layer$value$source)))
  }
  totals <- stats::aggregate(list(n_cells_in_nodes = rows$n_cells),
                             list(donor_id = as.character(rows$donor_id)), sum, na.rm = TRUE)
  list(value = layer$value, totals = totals)
}

# The measured cross-modality verdict, carried rather than recomputed. Its
# thresholds were fixed before it was measured and its donor count is its own --
# it is not the donor count of the layer above, and the two are printed
# separately for that reason.
.n2_read_verdict <- function(path = NULL, manifest_path = NULL) {
  resolve <- function(rel, given) {
    if (!is.null(given)) return(given)
    if (exists(".ie_pln_repo_file", mode = "function")) return(.ie_pln_repo_file(rel))
    if (file.exists(rel)) return(rel)
    NA_character_
  }
  csv <- resolve(.n2_verdict_rel, path)
  absent <- function(reason, state = "unavailable") {
    list(value = ie_avail(state, reason, source = if (is.na(csv)) .n2_verdict_rel else csv))
  }
  if (is.na(csv) || !file.exists(csv)) return(absent("cross_modality_verdict_not_exported"))
  data <- tryCatch(utils::read.csv(csv, check.names = FALSE), error = function(e) NULL)
  if (is.null(data)) return(absent("cross_modality_verdict_unreadable", "failed"))
  required <- c("tissue", "lineage", "n_donors", "median_share_protein_pct",
                "median_share_rna_pct", "median_diff_pp", "verdict", "verdict_reason")
  if (!all(required %in% names(data))) return(absent("cross_modality_verdict_schema_invalid", "failed"))
  data <- data[data$tissue == .n2_tissue, , drop = FALSE]
  if (!nrow(data)) return(absent("cross_modality_verdict_carries_no_lymph_node_rows", "blank"))
  sources <- NULL
  manifest_file <- resolve(.n2_verdict_manifest_rel, manifest_path)
  if (!is.na(manifest_file) && file.exists(manifest_file) &&
      requireNamespace("jsonlite", quietly = TRUE)) {
    manifest <- tryCatch(jsonlite::fromJSON(manifest_file, simplifyVector = TRUE),
                         error = function(e) NULL)
    if (!is.null(manifest)) sources <- manifest[[.n2_tissue]]
  }
  list(value = ie_avail("measured", "cross_modality_verdict_present", n = nrow(data), source = csv),
       data = data, sources = sources)
}

card_n2_pln_immune_read <- function(layers_dir = .n2_layers_dir(), active_domain_root = NULL,
                                    verdict_path = NULL, verdict_manifest_path = NULL) {
  list(
    grid = .n2_read_grid(layers_dir, active_domain_root),
    nodes = .n2_read_nodes(layers_dir),
    verdict = .n2_read_verdict(verdict_path, verdict_manifest_path),
    layers_dir = layers_dir
  )
}

.n2_count <- function(x) {
  if (length(x) != 1L || is.na(x)) return("unavailable")
  format(x, big.mark = ",", trim = TRUE)
}

# Disease-stage order comes from the DOMAIN, never from a literal: the two datasets declare
# their own levels, and which of those levels a BUILD actually carries is a separate question
# that has already changed -- the run cohort had no Aab+ donor before 2026-09-15 and has three
# now. A literal would have been wrong on one side of that date whichever way it was written.
.n2_status_order <- function(values) {
  observed <- unique(as.character(values[!is.na(values) & nzchar(as.character(values))]))
  levels <- if (exists("domain_group_levels", mode = "function")) {
    tryCatch(as.character(domain_group_levels()), error = function(e) character())
  } else character()
  c(levels[levels %in% observed], sort(setdiff(observed, levels)))
}

.n2_donor_order <- function(presence) {
  status <- vapply(unique(presence$donor_id), function(donor) {
    values <- presence$disease_status[presence$donor_id == donor]
    values <- values[!is.na(values) & nzchar(as.character(values))]
    if (length(values)) as.character(values)[[1L]] else NA_character_
  }, character(1))
  order_levels <- .n2_status_order(status)
  rank <- match(status, order_levels, nomatch = length(order_levels) + 1L)
  names(status)[order(rank, names(status))]
}

# One lane's own buckets for one donor, as a table. Counts are summed within a
# bucket -- a producer may write a bucket more than once -- and never across
# modalities. Returns NULL when this lane measured nothing for this donor, so a
# table output is mounted only where it will have rows in it.
.n2_bucket_table <- function(layer, donor, modality) {
  if (!ie_is_measured(layer$grid$value)) return(NULL)
  composition <- layer$grid$composition
  rows <- composition[composition$donor_id == donor & composition$modality == modality &
                        composition$state == "measured" &
                        is.finite(composition$frac_of_assessed), , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  rows <- stats::aggregate(list(n_cells = rows$n_cells, frac = rows$frac_of_assessed),
                           list(lineage_common = rows$lineage_common,
                                bucket_kind = rows$bucket_kind,
                                is_immune = rows$is_immune), sum)
  rows <- rows[order(rows$bucket_kind != "lineage", -rows$frac), , drop = FALSE]
  data.frame(
    Bucket = rows$lineage_common, Kind = rows$bucket_kind,
    Immune = ifelse(is.na(rows$is_immune) | !nzchar(rows$is_immune), "\u2014", rows$is_immune),
    Cells = vapply(rows$n_cells, .n2_count, character(1)),
    `% of cells this lane assessed` = sprintf("%.2f%%", 100 * rows$frac),
    check.names = FALSE, stringsAsFactors = FALSE)
}

# The output id carrying one lane's table. Derived from the modality the layer
# declares, so a lane added upstream gets its own output without this file
# changing.
.n2_table_id <- function(modality) paste0("n2_table_", gsub("[^A-Za-z0-9_]", "_", modality))

# One modality's own panel: its own denominator, its own state, its own counts.
# Nothing in it is shared with the panel beside it.
.n2_modality_panel <- function(layer, donor, modality, ns = identity) {
  presence <- layer$grid$presence
  row <- presence[presence$donor_id == donor & presence$modality == modality, , drop = FALSE]
  if (!nrow(row)) {
    return(shiny::tagList(
      shiny::tags$h5(modality),
      .n2_note_ui(.n2_note("this_lane_has_no_row_for_this_donor", "blank"))))
  }
  header <- shiny::tagList(
    shiny::tags$h5(modality),
    .n2_fact(sprintf("%s — %s", row$method[[1L]],
                     ie_state_label(row$state[[1L]], row$n_cells_assessed[[1L]]))))
  if (!identical(row$state[[1L]], "measured")) {
    conflict <- row$citation_conflict[[1L]]
    return(shiny::tagList(
      header,
      .n2_note_ui(.n2_note(row$reason[[1L]], row$state[[1L]], row$source[[1L]])),
      if (!is.na(conflict) && nzchar(conflict))
        .n2_note_ui(.n2_note(conflict, "unavailable", row$source[[1L]])) else NULL))
  }
  composition <- layer$grid$composition
  table <- .n2_bucket_table(layer, donor, modality)
  if (is.null(table)) {
    return(shiny::tagList(header,
                          .n2_note_ui(.n2_note("no_measured_bucket_for_this_lane_and_donor"))))
  }
  holes <- composition[composition$donor_id == donor & composition$modality == modality &
                         composition$state != "measured", , drop = FALSE]
  shiny::tagList(
    header,
    shiny::tableOutput(ns(.n2_table_id(modality))),
    .n2_fact(sprintf(paste0("Denominator: %s cells this lane assessed on this donor's section. ",
                            "Cells with no nucleus, excluded before anything was counted: %s."),
                     .n2_count(row$n_cells_assessed[[1L]]),
                     .n2_count(row$n_no_nucleus_excluded[[1L]]))),
    # The buckets this lane did not measure, grouped by state so every reason is
    # carried verbatim and the artifact is named once rather than once per
    # bucket. Five near-identical guard notes bury the counts above them.
    if (nrow(holes)) lapply(sort(unique(holes$state)), function(state) {
      reasons <- unique(holes$reason[holes$state == state & nzchar(holes$reason)])
      if (!length(reasons)) reasons <- "no reason recorded"
      .n2_note_ui(.n2_note(sprintf("%d bucket(s) this lane did not measure -- %s",
                                   sum(holes$state == state), paste(reasons, collapse = "; ")),
                           state, row$source[[1L]]))
    }) else NULL)
}

# What the node annotations hold on this donor's PhenoCycler slide, beside what
# the typing lane assessed. They are different numbers and the difference is not
# an error: one is every cell inside a node annotation, the other is the cells
# that lane could type.
.n2_node_note <- function(layer, donor) {
  nodes <- layer$nodes
  if (!ie_is_measured(nodes$value)) {
    return(.n2_note_ui(.n2_note(nodes$value$reason, nodes$value$state, nodes$value$source,
                                "core/phenocycler/scripts/pln/build_pln_layer.py")))
  }
  row <- nodes$totals[nodes$totals$donor_id == donor, , drop = FALSE]
  if (!nrow(row)) {
    return(.n2_note_ui(.n2_note("no_lymph_node_annotation_for_this_donor_in_the_node_layer",
                                "blank", nodes$value$source,
                                "core/phenocycler/scripts/pln/build_pln_layer.py")))
  }
  .n2_fact(sprintf(paste0("The node annotations on this donor's PhenoCycler slide hold %s cells, ",
                          "placed by point-in-polygon on the QuPath object_id. That is the ",
                          "section's own total, not a denominator of either lane above."),
                   .n2_count(row$n_cells_in_nodes[[1L]])))
}

# THE MODALITY GAP, as a measured quantity with its own donor count.
#
# Paired within donor: a donor measured on one side only contributes to neither
# the median nor the count. The two lanes have different denominators, so the
# difference is in percentage points of each lane's own assessed cells -- it is
# a difference between two shares, not a difference between two cell counts.
.n2_gap_table <- function(layer) {
  if (!ie_is_measured(layer$grid$value)) return(NULL)
  composition <- layer$grid$composition
  rows <- composition[composition$state == "measured" &
                        is.finite(composition$frac_of_assessed), , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  modalities <- sort(unique(rows$modality))
  if (length(modalities) != 2L) return(NULL)
  out <- list()
  per_donor <- function(sub, modality) {
    side <- sub[sub$modality == modality, , drop = FALSE]
    if (!nrow(side)) return(NULL)
    stats::aggregate(list(frac = side$frac_of_assessed), list(donor_id = side$donor_id), sum)
  }
  for (lineage in sort(unique(rows$lineage_common))) {
    sub <- rows[rows$lineage_common == lineage, , drop = FALSE]
    left <- per_donor(sub, modalities[[1L]])
    right <- per_donor(sub, modalities[[2L]])
    if (is.null(left) || is.null(right)) next
    shared <- intersect(left$donor_id, right$donor_id)
    if (!length(shared)) next
    a <- left$frac[match(shared, left$donor_id)]
    b <- right$frac[match(shared, right$donor_id)]
    diff <- 100 * (b - a)
    entry <- data.frame(
      Bucket = lineage, `Donors measured in both` = length(shared),
      check.names = FALSE, stringsAsFactors = FALSE)
    entry[[paste0("Median % (", modalities[[1L]], ")")]] <- 100 * stats::median(a)
    entry[[paste0("Median % (", modalities[[2L]], ")")]] <- 100 * stats::median(b)
    entry[[paste0("Median gap, pp (", modalities[[2L]], " − ", modalities[[1L]], ")")]] <-
      stats::median(diff)
    entry[["Gap IQR, pp"]] <- stats::IQR(diff)
    out[[length(out) + 1L]] <- entry
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Each column is formatted as what it is: a count is a count, a share carries a
# per cent sign, and only the gap carries a sign -- a "+11.0" donor count would
# read as a gap.
.n2_gap_cell <- function(name, value) {
  if (!is.numeric(value)) return(as.character(value))
  if (startsWith(name, "Median %")) return(sprintf("%.1f%%", value))
  if (startsWith(name, "Median gap")) return(sprintf("%+.1f", value))
  if (startsWith(name, "Gap IQR")) return(sprintf("%.1f", value))
  format(value, big.mark = ",", trim = TRUE)
}

.n2_gap_panel <- function(layer) {
  table <- .n2_gap_table(layer)
  if (is.null(table)) {
    return(.n2_note_ui(.n2_note("no_bucket_is_measured_in_both_modalities_for_any_donor", "blank")))
  }
  shiny::tagList(
    shiny::tags$table(
      class = "ie-pln-composition",
      shiny::tags$tr(lapply(names(table), shiny::tags$th)),
      lapply(seq_len(nrow(table)), function(i) {
        shiny::tags$tr(lapply(seq_along(table), function(j) {
          shiny::tags$td(.n2_gap_cell(names(table)[[j]], table[[j]][[i]]))
        }))
      })),
    .n2_fact(paste0("Paired within donor, over the donors each row names -- a donor measured on ",
                    "one side only is in neither the median nor the count. Each share is of the ",
                    "cells its own lane assessed, so the gap is a difference between two shares ",
                    "and never a difference between two cell counts. The gap is shown as it ",
                    "falls: neither lane is corrected toward the other, and neither is the ",
                    "target of the other.")))
}

# The carried verdict. Its thresholds were fixed before it was measured, its
# donor count is its own, and the releases it was measured against are named --
# a verdict measured against one release is not evidence about another.
.n2_verdict_panel <- function(layer) {
  verdict <- layer$verdict
  if (!ie_is_measured(verdict$value)) {
    return(.n2_note_ui(.n2_note(verdict$value$reason, verdict$value$state, verdict$value$source,
                                "scripts/senior/modality_concordance.py")))
  }
  data <- verdict$data
  sources <- verdict$sources
  source_note <- if (is.null(sources) || is.null(sources[["protein_source"]])) {
    .n2_note_ui(.n2_note("cross_modality_verdict_does_not_record_its_sources", "blank",
                         verdict$value$source, "scripts/senior/modality_concordance.py"))
  } else {
    layer_roots <- NULL
    grid_manifest <- layer$grid$manifest
    if (!is.null(grid_manifest) && !is.null(grid_manifest[["grid"]])) {
      grid <- grid_manifest[["grid"]]
      if (is.data.frame(grid)) {
        pln <- grid[grid$tissue == .n2_tissue, , drop = FALSE]
        if (nrow(pln)) layer_roots <- stats::setNames(as.character(pln$root), pln$modality)
      }
    }
    shiny::tagList(
      .n2_fact(sprintf("Measured against %s (protein) and %s (RNA).",
                       sources[["protein_source"]], sources[["rna_source"]])),
      if (!is.null(layer_roots))
        .n2_fact(sprintf(paste0("The composition above is read from %s. Where those differ, the ",
                                "verdict and the composition were measured on different builds ",
                                "and are two statements, not one."),
                         paste(sprintf("%s: %s", names(layer_roots), layer_roots),
                               collapse = "; ")))
        else NULL)
  }
  shiny::tagList(
    shiny::tags$table(
      class = "ie-pln-composition",
      shiny::tags$tr(lapply(c("Lineage", "Donors", "Median % protein", "Median % RNA",
                              "Median gap, pp", "Verdict", "Why"), shiny::tags$th)),
      lapply(seq_len(nrow(data)), function(i) {
        shiny::tags$tr(
          shiny::tags$td(data$lineage[[i]]),
          shiny::tags$td(format(data$n_donors[[i]])),
          shiny::tags$td(sprintf("%.1f%%", data$median_share_protein_pct[[i]])),
          shiny::tags$td(sprintf("%.1f%%", data$median_share_rna_pct[[i]])),
          shiny::tags$td(sprintf("%+.1f", data$median_diff_pp[[i]])),
          shiny::tags$td(data$verdict[[i]]),
          shiny::tags$td(data$verdict_reason[[i]]))
      })),
    .n2_fact(paste0("Carried from the concordance export, not recomputed here. Its shares are ",
                    "of each lane's resolved immune cells, which is a different denominator ",
                    "from the table above, and its donor count is its own.")),
    source_note)
}

card_n2_pln_immune_ui <- function(ns) {
  shiny::div(
    class = "ie-card", id = ns("n2_pln_immune"), style = "overflow: visible;",
    shiny::tags$h3(class = "ie-card-title", "pLN immune composition — each modality on its own"),
    shiny::uiOutput(ns("n2_body"))
  )
}

card_n2_pln_immune_server <- function(
    input, output, session, prepared = NULL, cohort = NULL,
    layers_dir = ie_layers_dir(),
    active_domain_root = NULL, verdict_path = NULL, verdict_manifest_path = NULL,
    active_tab = NULL, tab_key = NULL, ...) {
  ns <- session$ns
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                              identical(active_tab(), tab_key))

  layer <- shiny::reactive({
    shiny::req(active())
    card_n2_pln_immune_read(layers_dir = layers_dir, active_domain_root = active_domain_root,
                            verdict_path = verdict_path,
                            verdict_manifest_path = verdict_manifest_path)
  })

  roster <- shiny::reactive({
    if (!exists(".ie_pln_roster", mode = "function")) return(NULL)
    .ie_pln_roster(cohort)
  })

  all_donors <- shiny::reactive({
    data <- layer()
    if (!ie_is_measured(data$grid$value)) return(character())
    .n2_donor_order(data$grid$presence)
  })

  donors <- shiny::reactive({
    everyone <- all_donors()
    allowed <- roster()
    if (is.null(allowed) || !length(everyone)) return(everyone)
    keep <- everyone[everyone %in% allowed]
    if (length(keep)) keep else everyone
  })

  donor <- shiny::reactive({
    choices <- donors()
    shiny::req(length(choices) > 0L)
    selected <- input$n2_donor
    if (is.null(selected) || !selected %in% choices) choices[[1L]] else selected
  })

  modalities <- shiny::reactive({
    data <- layer()
    if (!ie_is_measured(data$grid$value)) return(character())
    sort(unique(data$grid$presence$modality))
  })

  output$n2_body <- shiny::renderUI({
    data <- layer()
    if (!ie_is_measured(data$grid$value)) {
      return(shiny::tagList(
        .n2_note_ui(.n2_note(data$grid$value$reason, data$grid$value$state,
                             data$grid$value$source)),
        shiny::tags$h4("The carried cross-modality verdict"),
        .n2_verdict_panel(data)))
    }
    shiny::tagList(
      shiny::tags$p(
        paste0("Each lane below typed a different section of a different block against a ",
               "different panel, and neither is a correction of the other. They are shown side ",
               "by side and are never pooled, never added and never placed on one denominator. ",
               "No distance, adjacency or neighbourhood is computed between a lymph node and a ",
               "pancreas anywhere on this card: there are no islets in a lymph node, so a ",
               "distance to one is undefined there rather than zero.")),
      shiny::selectInput(ns("n2_donor"), "Donor", choices = donors(),
                         selected = shiny::isolate(donor())),
      shiny::uiOutput(ns("n2_withheld")),
      shiny::uiOutput(ns("n2_modalities")),
      shiny::tags$h4("What the two lanes differ by"),
      shiny::uiOutput(ns("n2_gap")),
      shiny::tags$h4("The carried cross-modality verdict"),
      shiny::uiOutput(ns("n2_verdict")))
  })

  output$n2_withheld <- shiny::renderUI({
    outside <- setdiff(all_donors(), donors())
    if (!length(outside)) return(shiny::tagList())
    .n2_note_ui(.n2_note(sprintf(
      paste0("%d donor(s) carry a lymph node row in this layer and are not in this domain's ",
             "cohort roster, so they are not selectable here: %s. They are measured and they ",
             "exist; they are named rather than dropped in silence."),
      length(outside), paste(sort(outside), collapse = ", ")), "withheld",
      layer()$grid$value$source))
  })

  output$n2_modalities <- shiny::renderUI({
    data <- layer()
    names <- modalities()
    shiny::req(length(names) > 0L)
    width <- max(3L, as.integer(floor(12 / length(names))))
    shiny::tagList(
      shiny::fluidRow(lapply(names, function(modality) {
        shiny::column(width, .n2_modality_panel(data, donor(), modality, ns))
      })),
      .n2_node_note(data, donor()))
  })

  # A table renderer per lane the layer declares, registered the way card_g1_grid
  # registers its per-cell plots: a lane added upstream draws without this file
  # changing. The output is mounted by the panel above only when the lane
  # measured something for the selected donor, so no empty table is ever drawn.
  shiny::observeEvent(modalities(), {
    for (name in modalities()) {
      local({
        this <- name
        output[[.n2_table_id(this)]] <- shiny::renderTable({
          table <- .n2_bucket_table(layer(), donor(), this)
          shiny::req(is.data.frame(table), nrow(table) > 0L)
          table
        }, rownames = FALSE)
      })
    }
  }, ignoreNULL = TRUE)

  output$n2_gap <- shiny::renderUI(.n2_gap_panel(layer()))
  output$n2_verdict <- shiny::renderUI(.n2_verdict_panel(layer()))

  invisible(list(active = active, layer = layer, donor = donor, donors = donors,
                 gap = shiny::reactive(.n2_gap_table(layer()))))
}
