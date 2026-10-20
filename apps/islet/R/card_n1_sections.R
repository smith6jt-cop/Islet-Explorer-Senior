# ---------------------------------------------------------------------------
# N1 -- which annotation is the node on each pLN slide.
#
# The node is the pLN's structural unit: it plays the part the islet plays in
# the pancreas. WHICH annotation is the node is a fact about the acquisition,
# recorded by hand in core/phenocycler/data/integration/pln_node_annotations.csv,
# and it is not something a builder or this card may infer. The card shows that
# record beside what the pLN layer MEASURED from the same annotations
# (<ISLET_LAYERS_DIR>/pln, built by core/multiplex5000/scripts/pln/build_pln_layer.py)
# and beside the sections the sections layer carries for the same pLN
# (<ISLET_LAYERS_DIR>/sections, built by scripts/senior/build_sections_layer.py --
# which until now had no reader anywhere in the app).
#
# WHAT THIS CARD REFUSES TO DO.
#
# * It never hides an exception. A node annotation with no area bounds no
#   region, an image whose only annotation is a whole-tissue outline gives a
#   per-slide number wearing a node's name, and a donor the node layer carries
#   that the ingest root does not cannot be joined to any measurement. All three
#   are named here, with their counts, taken from the artifacts rather than
#   written into a string.
# * It never writes a donor id, a cohort size or a node count into a literal.
#   Every one of those is read out of the record or the layer. A cohort fact
#   written as a literal is a domain leak the test suite cannot see.
# * It never turns an absence into a zero. A node with no usable geometry has no
#   cell count -- not zero cells -- and a section with no recorded pixel size
#   reads unavailable, because reading its coordinates at a guessed pixel size
#   would rescale every distance derived from them.
# * It never computes anything spatial between the pancreas and the pLN. They
#   are different blocks; this card only says which annotation bounds the node.
#
# The sections layer is joined only when it records the domain it was built
# from. The legacy CODEX build and the run cohort share five donors and no
# vocabulary, so a matching donor id is not a licence to read.

.n1_builder_layer <- "core/multiplex5000/scripts/pln/build_pln_layer.py"
.n1_builder_sections <- "scripts/senior/build_sections_layer.py"
.n1_record_rel <- "core/phenocycler/data/integration/pln_node_annotations.csv"

.n1_layers_dir <- function() ie_layers_dir()

.n1_note <- function(reason, state = "unavailable", artifact = NULL, builder = .n1_builder_layer) {
  ie_require_measured(ie_avail(state, reason = reason), guard = "n1_sections",
                      artifact = if (is.null(artifact) || is.na(artifact)) .n1_layers_dir() else artifact,
                      builder = builder)
}

.n1_note_ui <- function(note) shiny::tags$p(class = "ie-card-note", role = "status", format(note))

# A MEASURED statement, which is not a guard note. ie_require_measured() returns
# nothing for a measured value -- correctly, since there is nothing to guard --
# so a measured exception says itself in plain text rather than being dressed as
# an absence. A hole and a finding must never print the same way.
.n1_fact <- function(text) shiny::tags$p(class = "ie-card-note", role = "status", text)

# A data.frame rendered as the card's own table. Inline rather than a
# tableOutput: this card is one output, so every note and every table it has to
# show is in the one place a reader -- and the registry suite -- looks.
.n1_table_ui <- function(table, class) {
  shiny::tags$table(
    class = class,
    shiny::tags$tr(lapply(names(table), shiny::tags$th)),
    lapply(seq_len(nrow(table)), function(i) {
      shiny::tags$tr(lapply(as.character(unlist(table[i, , drop = FALSE])), shiny::tags$td))
    }))
}

#' Locate a repository-committed record from the running app's directory.
#'
#' The acquisition records are committed to the repository, not to a domain data
#' root: they describe which annotation was drawn on which image, which is true
#' of the acquisition whichever domain is loaded. `domain_data_root()` is still
#' the only route to scientific data, and nothing here reads any.
#'
#' Shared with card_n2_pln_immune.R, which reads the concordance export the same
#' way. Returns NA_character_ when the file is not found, so the caller can say
#' so by name instead of failing.
.ie_pln_repo_file <- function(rel) {
  root <- Sys.getenv("ISLET_REPO_ROOT", unset = "")
  candidates <- c(if (nzchar(root)) file.path(root, rel), rel,
                  file.path("..", rel), file.path("../..", rel),
                  file.path("../../..", rel), file.path("../../../..", rel))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) return(NA_character_)
  normalizePath(hit[[1L]], mustWork = FALSE)
}

#' The donor roster of the active domain's cohort, or NULL when there is none.
#'
#' Shared with card_n2_pln_immune.R. A donor the artifact carries that the
#' roster does not is named as withheld rather than dropped in silence: a donor
#' that vanishes with no note reads as a donor that was never acquired.
.ie_pln_roster <- function(cohort) {
  if (is.null(cohort)) return(NULL)
  value <- NULL
  if (is.list(cohort) && is.function(cohort$donor_ids)) {
    value <- tryCatch(as.character(cohort$donor_ids()), error = function(e) NULL)
  }
  if (is.null(value)) value <- tryCatch(as.character(cohort), error = function(e) character())
  value <- unique(value[grepl("^[0-9]+$", value)])
  if (length(value)) value else NULL
}

# A hand-recorded flag is read strictly. "True"/"False" are the only values the
# record writes, and anything else is not silently taken for FALSE.
.n1_flag <- function(x) {
  v <- tolower(trimws(as.character(x)))
  ifelse(v %in% c("true", "false"), v == "true", NA)
}

.n1_read_record <- function(path = .ie_pln_repo_file(.n1_record_rel)) {
  absent <- function(reason, state = "unavailable") {
    list(value = ie_avail(state, reason, source = if (is.na(path)) .n1_record_rel else path))
  }
  if (is.na(path) || !file.exists(path)) return(absent("pln_node_record_not_found"))
  data <- tryCatch(utils::read.csv(path, colClasses = c(donor_id = "character"),
                                   comment.char = "#", check.names = FALSE),
                   error = function(e) NULL)
  if (is.null(data)) return(absent("pln_node_record_unreadable", "failed"))
  required <- c("image_file", "donor_id", "roi", "qupath_entry", "node_classes",
                "n_nodes", "n_cells", "excluded", "recorded_by", "recorded_on",
                "excluded_reason")
  if (!all(required %in% names(data))) return(absent("pln_node_record_schema_invalid", "failed"))
  if (!nrow(data)) return(absent("pln_node_record_is_empty"))
  data$excluded_flag <- .n1_flag(data$excluded)
  if (anyNA(data$excluded_flag)) return(absent("pln_node_record_exclusion_not_readable", "failed"))
  # The record's own rule: a blank required cell is refused, never guessed. The
  # row stays visible and says what is missing.
  data$node_classes <- trimws(as.character(data$node_classes))
  data$recorded_classes <- lapply(data$node_classes, function(value) {
    if (!nzchar(value)) return(character())
    trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
  })
  data$n_classes_recorded <- vapply(data$recorded_classes, length, integer(1))
  list(value = ie_avail("measured", "pln_node_record_present", n = nrow(data), source = path),
       data = data)
}

.n1_read_nodes <- function(dir) {
  path <- file.path(dir, "nodes.parquet")
  absent <- function(reason, state = "unavailable") list(value = ie_avail(state, reason, source = path))
  if (!file.exists(path)) return(absent("pln_node_layer_not_built"))
  if (!requireNamespace("arrow", quietly = TRUE)) return(absent("arrow_reader_unavailable"))
  data <- tryCatch(as.data.frame(arrow::read_parquet(path)), error = function(e) NULL)
  if (is.null(data)) return(absent("pln_node_layer_unreadable", "failed"))
  required <- c("donor_id", "image", "qupath_entry", "node_index", "node_name",
                "annotation_geom", "delimiter_kind", "area_um2", "state", "reason",
                "n_cells", "node_scope", "node_scope_reason", "qupath_check",
                "delta_vs_qupath")
  if (!all(required %in% names(data))) return(absent("pln_node_layer_schema_invalid", "failed"))
  if (!nrow(data)) return(absent("pln_node_layer_is_empty"))
  data$donor_id <- as.character(data$donor_id)
  if (anyNA(data$state) || any(!data$state %in% IE_STATES)) {
    return(absent("pln_node_state_invalid", "failed"))
  }
  # A node the layer could not measure must not carry a count: a missing
  # measurement is unavailable, never zero.
  holes <- data[data$state != "measured", , drop = FALSE]
  if (nrow(holes) && (any(!is.na(holes$n_cells)) || any(!is.na(holes$area_um2)))) {
    return(absent("unmeasured_node_carries_a_count", "failed"))
  }
  if (nrow(holes) && (anyNA(holes$reason) || any(!nzchar(holes$reason)))) {
    return(absent("unmeasured_node_has_no_reason", "failed"))
  }
  list(value = ie_avail("measured", "pln_node_layer_present", n = nrow(data), source = path),
       data = data)
}

# The layer's own record of which donors the ingest root holds. A donor this
# layer carries that the ingest root does not cannot be joined to any
# measurement -- an object_id join drops it in silence, which is exactly the
# failure this panel exists to make loud.
.n1_read_coverage <- function(dir) {
  path <- file.path(dir, "MANIFEST.json")
  absent <- function(reason, state = "unavailable") list(value = ie_avail(state, reason, source = path))
  if (!file.exists(path)) return(absent("pln_layer_manifest_absent"))
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(absent("json_reader_unavailable"))
  manifest <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(manifest)) return(absent("pln_layer_manifest_unreadable", "failed"))
  coverage <- manifest[["ingested_cells_coverage"]]
  if (is.null(coverage)) return(absent("pln_layer_manifest_records_no_ingest_coverage"))
  list(value = ie_avail("measured", "pln_layer_ingest_coverage_present", source = path),
       coverage = coverage, manifest = manifest)
}

.n1_read_sections <- function(dir, active_domain_root = NULL) {
  path <- file.path(dir, "data.parquet")
  absent <- function(reason, state = "unavailable") list(value = ie_avail(state, reason, source = path))
  if (!file.exists(path)) return(absent("sections_layer_not_built"))
  if (!requireNamespace("arrow", quietly = TRUE)) return(absent("arrow_reader_unavailable"))
  table <- tryCatch(arrow::read_parquet(path, as_data_frame = FALSE), error = function(e) NULL)
  if (is.null(table)) return(absent("sections_layer_unreadable", "failed"))
  recorded <- table$schema$metadata[["source_domain_root"]]
  root <- if (is.null(active_domain_root)) {
    if (exists("domain_data_root", mode = "function")) domain_data_root() else NULL
  } else active_domain_root
  if (is.null(root)) return(absent("active_domain_root_unknown"))
  if (is.null(recorded) || !nzchar(recorded)) return(absent("sections_layer_domain_not_recorded"))
  canonical <- function(p) normalizePath(p, mustWork = FALSE)
  if (!identical(canonical(recorded), canonical(root))) {
    return(absent("sections_layer_belongs_to_a_different_domain"))
  }
  data <- tryCatch(as.data.frame(table), error = function(e) NULL)
  if (is.null(data)) return(absent("sections_layer_unreadable", "failed"))
  required <- c("donor_id", "roi", "modality", "image_file", "frame")
  if (!all(required %in% names(data))) return(absent("sections_layer_schema_invalid", "failed"))
  if (!nrow(data)) return(absent("sections_layer_is_empty"))
  data$donor_id <- as.character(data$donor_id)
  list(value = ie_avail("measured", "sections_layer_present", n = nrow(data), source = path),
       data = data)
}

card_n1_sections_read <- function(layers_dir = .n1_layers_dir(),
                                  record_path = .ie_pln_repo_file(.n1_record_rel),
                                  active_domain_root = NULL) {
  list(
    record = .n1_read_record(record_path),
    nodes = .n1_read_nodes(file.path(layers_dir, "pln")),
    coverage = .n1_read_coverage(file.path(layers_dir, "pln")),
    sections = .n1_read_sections(file.path(layers_dir, "sections"), active_domain_root),
    layers_dir = layers_dir
  )
}

.n1_count <- function(x) {
  if (length(x) != 1L || is.na(x)) return("unavailable")
  format(x, big.mark = ",", trim = TRUE)
}

# One row per recorded pLN image: what the maintainer recorded, and what the
# layer measured from the same annotations. The two are shown side by side
# rather than reconciled, so a disagreement is visible rather than resolved
# behind the reader's back.
.n1_record_table <- function(layer) {
  record <- layer$record$data
  nodes <- if (ie_is_measured(layer$nodes$value)) layer$nodes$data else NULL
  rows <- lapply(seq_len(nrow(record)), function(i) {
    row <- record[i, , drop = FALSE]
    mine <- if (is.null(nodes)) nodes else
      nodes[nodes$donor_id == row$donor_id & nodes$image == row$image_file, , drop = FALSE]
    usable <- if (is.null(mine)) NA_integer_ else sum(mine$state == "measured")
    measured_cells <- if (is.null(mine) || !nrow(mine) || !any(mine$state == "measured")) NA_real_
      else sum(mine$n_cells[mine$state == "measured"], na.rm = TRUE)
    scope <- if (is.null(mine) || !nrow(mine)) NA_character_ else
      paste(sort(unique(stats::na.omit(mine$node_scope))), collapse = ", ")
    data.frame(
      Donor = row$donor_id,
      Section = row$roi,
      `QuPath entry` = row$qupath_entry,
      `Annotation that is the node` = if (nzchar(row$node_classes)) row$node_classes
        else "unavailable(node_classes_not_recorded)",
      `Nodes recorded` = row$n_nodes,
      `Nodes with usable geometry` = if (is.na(usable)) "unavailable" else as.character(usable),
      `Cells recorded` = .n1_count(row$n_cells),
      `Cells assigned to a node` = .n1_count(measured_cells),
      Scope = if (is.na(scope) || !nzchar(scope)) "unavailable" else scope,
      check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Every exception, derived from the artifacts and named. Nothing here is a list
# of donors typed into this file: each block asks the record or the layer which
# rows it applies to, so a new one appears without this card changing and a
# resolved one disappears.
.n1_exception_notes <- function(layer) {
  notes <- list()
  record <- if (ie_is_measured(layer$record$value)) layer$record$data else NULL

  if (!is.null(record)) {
    # An image whose node is delimited by something other than a drawn pLN
    # region. The record carries the maintainer's confirmation that for such an
    # image that annotation IS the node; it is shown, not assumed away.
    other <- vapply(record$recorded_classes, function(classes) {
      length(classes) > 0L && !all(grepl("^pLN", classes))
    }, logical(1))
    for (i in which(other)) {
      notes[[length(notes) + 1L]] <- .n1_fact(sprintf(
        paste0("donor %s (%s): the node is delimited by %s, not by a drawn pLN region. The ",
               "record carries that as the node for this image, with %s cells; it is the ",
               "recorded exception, not an inference."),
        record$donor_id[[i]], record$image_file[[i]],
        paste(record$recorded_classes[[i]], collapse = ", "),
        .n1_count(record$n_cells[[i]])))
    }
    blank <- which(record$n_classes_recorded == 0L)
    for (i in blank) {
      notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
        "donor %s (%s): no node annotation is recorded for this image, and none is guessed",
        record$donor_id[[i]], record$image_file[[i]]),
        artifact = layer$record$value$source, builder = NULL))
    }
    disagree <- which(record$n_classes_recorded != record$n_nodes)
    for (i in disagree) {
      notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
        paste0("donor %s (%s): the record names %d annotation(s) and counts %s. The two ",
               "columns of one hand-recorded row disagree."),
        record$donor_id[[i]], record$image_file[[i]],
        record$n_classes_recorded[[i]], format(record$n_nodes[[i]])), "failed",
        artifact = layer$record$value$source, builder = NULL))
    }
  }

  if (ie_is_measured(layer$nodes$value)) {
    nodes <- layer$nodes$data
    holes <- nodes[nodes$state != "measured", , drop = FALSE]
    for (i in seq_len(nrow(holes))) {
      notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
        paste0("donor %s, node %s on %s: %s (%s). It bounds no region, so it holds no cells ",
               "-- not zero cells -- and nothing on this slide may be attributed to it."),
        holes$donor_id[[i]], holes$node_name[[i]], holes$image[[i]],
        holes$reason[[i]], holes$annotation_geom[[i]]),
        holes$state[[i]], artifact = layer$nodes$value$source))
    }
    slide <- nodes[nodes$state == "measured" & !is.na(nodes$node_scope) &
                     nodes$node_scope == "slide", , drop = FALSE]
    for (i in seq_len(nrow(slide))) {
      notes[[length(notes) + 1L]] <- .n1_fact(sprintf(
        paste0("donor %s, node %s on %s: %s Pooling it with per-node values from other ",
               "slides mixes the two."),
        slide$donor_id[[i]], slide$node_name[[i]], slide$image[[i]],
        slide$node_scope_reason[[i]]))
    }
    differs <- nodes[nodes$state == "measured" & !is.na(nodes$qupath_check) &
                       nodes$qupath_check == "differs", , drop = FALSE]
    for (i in seq_len(nrow(differs))) {
      notes[[length(notes) + 1L]] <- .n1_fact(sprintf(
        paste0("donor %s, node %s: point-in-polygon membership and QuPath's own child count ",
               "differ by %s cell(s) out of %s. Membership is the geometric one, carried on ",
               "the QuPath object_id; the QuPath count is a cross-check, never the source."),
        differs$donor_id[[i]], differs$node_name[[i]],
        format(differs$delta_vs_qupath[[i]]), .n1_count(differs$n_cells[[i]])))
    }
  } else {
    notes[[length(notes) + 1L]] <- .n1_note_ui(
      .n1_note(layer$nodes$value$reason, layer$nodes$value$state, layer$nodes$value$source))
  }

  if (ie_is_measured(layer$coverage$value)) {
    coverage <- layer$coverage$coverage
    missing <- as.character(coverage[["in_this_layer_but_not_ingested"]])
    rec <- if (ie_is_measured(layer$record$value)) layer$record$data else NULL
    for (donor in missing) {
      cells <- if (!is.null(rec) && donor %in% rec$donor_id)
        rec$n_cells[rec$donor_id == donor][[1L]] else NA_real_
      excluded <- if (!is.null(rec) && donor %in% rec$donor_id)
        rec$excluded_flag[rec$donor_id == donor][[1L]] else NA
      notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
        paste0("donor %s is in the node layer with %s cells and is not in the ingested cell ",
               "root the layer names (%s). The node record's own exclusion flag for it reads ",
               "%s. Until the ingest is re-run, an object_id join to this donor returns ",
               "nothing, and nothing measured on its cells can be shown -- the cells are not ",
               "missing, the join is."),
        donor, .n1_count(cells), coverage[["root"]],
        if (isTRUE(excluded)) "excluded=True" else if (isFALSE(excluded)) "excluded=False"
        else "unavailable"),
        artifact = layer$coverage$value$source,
        builder = "core/multiplex5000/scripts/pln/ingest_pln_cells.py"))
    }
    stranded <- as.character(coverage[["ingested_but_not_in_this_layer"]])
    for (donor in stranded) {
      notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
        paste0("donor %s has ingested cells and no node record in this layer, so no cell of ",
               "it can be placed in a node."), donor),
        artifact = layer$coverage$value$source))
    }
  } else {
    notes[[length(notes) + 1L]] <- .n1_note_ui(
      .n1_note(layer$coverage$value$reason, layer$coverage$value$state, layer$coverage$value$source))
  }
  notes
}

# Which sections exist for the pLN of each recorded donor. The sections layer
# spans the modalities; the pLN annotation record is a PhenoCycler acquisition
# record, so the two carry different donor sets on purpose and both differences
# are named.
.n1_sections_panel <- function(layer, donors = NULL) {
  if (!ie_is_measured(layer$sections$value)) {
    return(.n1_note_ui(.n1_note(layer$sections$value$reason, layer$sections$value$state,
                                layer$sections$value$source, .n1_builder_sections)))
  }
  data <- layer$sections$data
  rows <- data[data$roi != "panc", , drop = FALSE]
  if (!nrow(rows)) {
    return(.n1_note_ui(.n1_note("sections_layer_carries_no_lymph_node_section", "blank",
                                layer$sections$value$source, .n1_builder_sections)))
  }
  if (!is.null(donors)) rows <- rows[rows$donor_id %in% donors, , drop = FALSE]
  if (!nrow(rows)) {
    return(.n1_note_ui(.n1_note("no_lymph_node_section_for_the_donors_on_this_card", "blank",
                                layer$sections$value$source, .n1_builder_sections)))
  }
  pixel <- if ("pixel_size_um" %in% names(rows)) {
    ifelse(is.na(rows$pixel_size_um), "unavailable", sprintf("%.4f", rows$pixel_size_um))
  } else rep("unavailable", nrow(rows))
  serial <- if ("xenium_serial" %in% names(rows)) {
    ifelse(is.na(rows$xenium_serial) | !nzchar(rows$xenium_serial), "unavailable",
           as.character(rows$xenium_serial))
  } else rep("unavailable", nrow(rows))
  table <- data.frame(
    Donor = rows$donor_id, Section = rows$roi, Modality = rows$modality,
    Image = ifelse(is.na(rows$image_file) | !nzchar(rows$image_file), "unavailable",
                   basename(as.character(rows$image_file))),
    Serial = serial, `Pixel size (µm)` = pixel, Frame = rows$frame,
    check.names = FALSE, stringsAsFactors = FALSE)
  # Lane C's section_relationship, when the layer carries it. It is a dataset
  # attribute -- same slide, same block, different block -- and never a distance
  # gate: no distance crosses between the pancreas and a lymph node here or
  # anywhere else.
  relationship <- if ("section_relationship" %in% names(rows)) {
    table$`Section relationship` <- ifelse(
      is.na(rows$section_relationship) | !nzchar(as.character(rows$section_relationship)),
      "unavailable", as.character(rows$section_relationship))
    NULL
  } else {
    .n1_note_ui(.n1_note("sections_layer_records_no_section_relationship", "blank",
                         layer$sections$value$source, .n1_builder_sections))
  }
  shiny::tagList(.n1_table_ui(table, "ie-pln-sections"), relationship)
}

.n1_coverage_note <- function(layer, donors, record_donors) {
  if (!ie_is_measured(layer$sections$value)) return(NULL)
  rows <- layer$sections$data
  rows <- rows[rows$roi != "panc", , drop = FALSE]
  no_section <- setdiff(donors, unique(rows$donor_id))
  no_record <- setdiff(unique(rows$donor_id), record_donors)
  notes <- list()
  if (length(no_section)) {
    notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
      paste0("%d donor(s) with a recorded pLN node carry no lymph node section in the sections ",
             "layer: %s. The node is recorded; the section that would carry it is not."),
      length(no_section), paste(sort(no_section), collapse = ", ")), "blank",
      layer$sections$value$source, .n1_builder_sections))
  }
  if (length(no_record)) {
    notes[[length(notes) + 1L]] <- .n1_note_ui(.n1_note(sprintf(
      paste0("%d donor(s) carry a lymph node section in the sections layer and no recorded node ",
             "annotation: %s. Which annotation is the node is recorded per PhenoCycler image, so ",
             "a section from another modality has no row to read."),
      length(no_record), paste(sort(no_record), collapse = ", ")), "blank",
      layer$sections$value$source, .n1_builder_sections))
  }
  notes
}

card_n1_sections_ui <- function(ns) {
  shiny::div(
    class = "ie-card", id = ns("n1_sections"), style = "overflow: visible;",
    shiny::tags$h3(class = "ie-card-title", "Which annotation is the node"),
    shiny::uiOutput(ns("n1_body"))
  )
}

card_n1_sections_server <- function(
    input, output, session, prepared = NULL, cohort = NULL,
    layers_dir = ie_layers_dir(),
    record_path = NULL, active_domain_root = NULL, active_tab = NULL, tab_key = NULL, ...) {
  ns <- session$ns
  active <- shiny::reactive(is.null(active_tab) || is.null(tab_key) ||
                              identical(active_tab(), tab_key))

  layer <- shiny::reactive({
    shiny::req(active())
    card_n1_sections_read(
      layers_dir = layers_dir,
      record_path = if (is.null(record_path)) .ie_pln_repo_file(.n1_record_rel) else record_path,
      active_domain_root = active_domain_root)
  })

  roster <- shiny::reactive(.ie_pln_roster(cohort))

  record_donors <- shiny::reactive({
    data <- layer()
    if (!ie_is_measured(data$record$value)) return(character())
    sort(unique(data$record$data$donor_id))
  })

  donors <- shiny::reactive({
    all_donors <- record_donors()
    allowed <- roster()
    if (is.null(allowed) || !length(all_donors)) return(all_donors)
    keep <- all_donors[all_donors %in% allowed]
    if (length(keep)) keep else all_donors
  })

  # ONE output for the whole card. The card has no control of its own, so there
  # is nothing to keep out of the re-render, and a note that lands in a child
  # output is a note the card's own output does not carry -- which is how a
  # reader, and a test that reads one output, sees a card with nothing to say
  # instead of a card saying why.
  # The record, filtered to the donors this domain's cohort can select. NULL
  # rather than an empty frame when there is nothing to show, so the table
  # output below is mounted only where it will have rows in it.
  record_table <- shiny::reactive({
    data <- layer()
    if (!ie_is_measured(data$record$value)) return(NULL)
    rows <- .n1_record_table(data)
    selectable <- donors()
    if (length(selectable)) rows <- rows[rows$Donor %in% selectable, , drop = FALSE]
    if (nrow(rows)) rows else NULL
  })

  output$n1_record_table <- shiny::renderTable({
    rows <- record_table()
    shiny::req(is.data.frame(rows), nrow(rows) > 0L)
    rows
  }, rownames = FALSE)

  output$n1_body <- shiny::renderUI({
    data <- layer()
    record <- data$record
    selectable <- donors()
    outside <- setdiff(record_donors(), selectable)
    table <- if (!ie_is_measured(record$value)) {
      .n1_note_ui(.n1_note(record$value$reason, record$value$state, record$value$source, NULL))
    } else if (is.null(record_table())) {
      .n1_note_ui(.n1_note("no_recorded_pln_node_for_this_cohort", "blank",
                           record$value$source, NULL))
    } else {
      shiny::tableOutput(ns("n1_record_table"))
    }
    exceptions <- .n1_exception_notes(data)
    if (!length(exceptions)) {
      exceptions <- list(.n1_fact(paste0(
        "No exception is recorded or measured on these slides: every recorded node annotation ",
        "bounds a region, and the node layer and the ingested cell root carry the same donors.")))
    }
    shiny::tagList(
      # Collapsed by default. The measured state below -- the record table, the named exceptions,
      # the withheld donors -- stays visible; only the standing explanation of HOW membership is
      # derived is behind the toggle, which is the rule test_app_text_is_collapsed_by_default pins.
      ie_explain(
        label = "How node membership is decided",
        shiny::tags$p(
          paste0("The node is the lymph node's structural unit, and which annotation delimits it ",
                 "is a fact about the acquisition rather than something a builder infers. The ",
                 "record below is hand-recorded per image; the measured columns beside it come ",
                 "from the node layer, which places each cell by point-in-polygon against the ",
                 "annotation and carries the QuPath object_id verbatim. QuPath's Parent column is ",
                 "never the source of membership: it reads the same whole-tissue annotation on ",
                 "every row of this export."))),
      shiny::tags$h4("The node on each pLN slide"),
      table,
      if (length(outside)) .n1_note_ui(.n1_note(sprintf(
        paste0("%d donor(s) have a recorded pLN node and are not in this domain's cohort roster, ",
               "so they are not shown above: %s. They are named rather than dropped in silence."),
        length(outside), paste(sort(outside), collapse = ", ")), "withheld",
        record$value$source, NULL)) else NULL,
      shiny::tags$h4("The exceptions, named"),
      exceptions,
      shiny::tags$h4("Sections carrying these lymph nodes"),
      .n1_sections_panel(data, if (length(selectable)) selectable else NULL),
      .n1_coverage_note(data, selectable, record_donors())
    )
  })

  invisible(list(active = active, layer = layer, donors = donors,
                 record_donors = record_donors))
}
