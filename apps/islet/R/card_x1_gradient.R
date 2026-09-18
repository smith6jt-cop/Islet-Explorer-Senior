# Small layer artifacts only. Both cards use these readers; no cell-tree scan
# or tissue-coordinate transform belongs in a Shiny reactive.
#
# The artifact carries two arms and they are not the same kind of measurement.
# roi == "panc" is spatial: signed distance to the islet edge, in the
# signed_distance_um frame. Every other roi is a lymph node section cut from a
# DIFFERENT BLOCK than the pancreas on the same slide; it shares a pixel grid
# with the pancreas and nothing anatomically. Those rows are compositional only,
# joined on donor_id and donor_status, and the reader rejects the whole artifact
# if one of them arrives carrying a distance, a bin or a compartment.

IE_PANCREAS_ROI <- "panc"

# Every source path the artifact records, and whether it is scoped to a domain.
#
# source_domain_root alone is not enough. It says which domain the BUILDER
# believed it was serving, and a builder that resolved its cells somewhere else
# still records the right domain root beside the wrong cells -- a right-looking
# islet drawn from the other cohort's measurements, which is the exact leak
# tests/testthat/test-domain-data-paths.R exists to catch. So each recorded path
# is checked, not just the one that names the domain.
#
# Domain-scoped sources are the ones feeding the pancreas distance axis. They
# must live in the active domain's tree: the domain root, or beside it under the
# same parent. For run20 that parent holds app/, cells/ and islet_membership/
# together, and the legacy CODEX root (data/app_data) is nowhere beneath it, so
# cells borrowed from the other cohort are refused by the source they came from
# rather than passing on a correct-looking source_domain_root.
#
# What this does NOT do is enumerate every configured domain and check the path
# against each: ieCore hands out the ACTIVE domain only. So the rule is
# containment in the active tree, which is what the two real roots differ by.
IE_DISTANCE_DOMAIN_SOURCES <- c("source_cells", "source_membership")

# Not domain-scoped, and not silently exempt either. The Xenium ingest root and
# the immune vocabulary are shared across domains by construction -- neither is a
# cohort measurement and neither feeds the distance axis -- so they are checked
# for being recorded and absolute rather than for domain containment.
#
# Each is required exactly when something depends on it, because an empty string
# is how the builder honestly records "this arm produced nothing":
#   source_xen_cells          -- required once any pLN row exists; the rows came
#                                from there whatever state they are in.
#   source_immune_vocabulary  -- required once a pLN row is MEASURED, which is
#                                the only case where a number was derived from
#                                the vocabulary. Rows that are blank or withheld
#                                never consulted it, and demanding it then would
#                                refuse the pancreas arm over a pLN-arm gap.
IE_DISTANCE_SHARED_SOURCES <- c("source_xen_cells", "source_immune_vocabulary")

.ie_path_within <- function(path, tree) {
  path <- normalizePath(path, mustWork = FALSE)
  tree <- normalizePath(tree, mustWork = FALSE)
  identical(path, tree) || startsWith(path, paste0(tree, .Platform$file.sep))
}

# Returns NULL when every recorded source path is acceptable, otherwise the
# reason naming the offending key.
.ie_distance_check_sources <- function(metadata, domain_root, has_pln_rows,
                                       has_measured_pln) {
  required_when <- c(source_xen_cells = has_pln_rows,
                     source_immune_vocabulary = has_measured_pln)
  read <- function(key) {
    value <- metadata[[key]]
    if (is.null(value)) return(NA_character_)
    as.character(value)
  }
  tree <- dirname(normalizePath(domain_root, mustWork = FALSE))
  for (key in IE_DISTANCE_DOMAIN_SOURCES) {
    value <- read(key)
    if (is.na(value) || !nzchar(value)) return(paste0("distance_profile_", key, "_not_recorded"))
    if (!startsWith(value, .Platform$file.sep)) {
      return(paste0("distance_profile_", key, "_not_absolute"))
    }
    if (!.ie_path_within(value, tree)) {
      return(paste0("distance_profile_", key, "_outside_domain"))
    }
  }
  for (key in IE_DISTANCE_SHARED_SOURCES) {
    value <- read(key)
    needed <- isTRUE(required_when[[key]])
    if (is.na(value) || !nzchar(value)) {
      # Absent or empty is how the builder records "nothing here depended on
      # this". A contradiction only when something did.
      if (needed) return(paste0("distance_profile_", key, "_not_recorded"))
      next
    }
    if (!startsWith(value, .Platform$file.sep)) {
      return(paste0("distance_profile_", key, "_not_absolute"))
    }
  }
  NULL
}

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
    required <- c("donor_id", "donor_status", "roi", "tissue", "modality",
                  "distance_bin_um", "n_cells", "n_immune",
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
    n_rows <- nrow(data)
    if (anyNA(data$roi) || any(!nzchar(data$roi))) {
      return(absent("distance_profile_roi_not_recorded", "failed"))
    }
    # A different block cannot be put on the pancreas distance axis. Refusing
    # the artifact is the only safe response: a pLN row that carries a bin has
    # already had a meaningless number computed for it upstream.
    pln <- data[data$roi != IE_PANCREAS_ROI, , drop = FALSE]
    if (nrow(pln) && (any(pln$frame != "none") || any(!is.na(pln$distance_bin_um)) ||
                      any(!is.na(pln$compartment)))) {
      return(absent("distance_profile_pln_carries_distance", "failed"))
    }
    unresolved <- pln[pln$state != "measured", , drop = FALSE]
    if (nrow(unresolved) && (anyNA(unresolved$reason) || any(!nzchar(unresolved$reason)))) {
      return(absent("distance_profile_pln_reason_missing", "failed"))
    }
    # Checked once the pLN rows are known, because whether a shared source has
    # to be recorded depends on whether that arm produced anything.
    bad_source <- .ie_distance_check_sources(
      table$schema$metadata, domain_data_root(),
      has_pln_rows = nrow(pln) > 0L,
      has_measured_pln = nrow(pln) > 0L && any(pln$state == "measured"))
    if (!is.null(bad_source)) return(absent(bad_source, "failed"))
    data <- data[data$roi == IE_PANCREAS_ROI, , drop = FALSE]
    if (!nrow(data)) return(absent("distance_profile_pancreas_arm_absent"))
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
    # n describes the rows these cards draw on the distance axis -- the pancreas
    # arm -- not the artifact's row count. The two differ by the compositional
    # pLN rows, which share no axis with the plot and are counted separately;
    # reporting the total made the card claim measurements it never drew.
    list(value = ie_avail("measured", n = nrow(data), source = path), data = data,
         pln = pln, n_artifact_rows = n_rows, n_pln_rows = nrow(pln), ledger = ledger,
         bin_width = table$schema$metadata[["bin_width_um"]],
         pln_semantics = table$schema$metadata[["pln_semantics"]])
  }, error = function(e) absent(paste0("distance_profile_read_failed: ", conditionMessage(e)), "failed"))
}

# The donors these cards may offer: the artifact's pancreas arm, narrowed to the
# active domain's cohort.
#
# The pLN rows are deliberately NOT unioned in. The artifact's pLN arm reaches
# donors the active domain does not contain -- it takes its donor_status from the
# integration manifest, so a section with no pancreas in this cohort still gets a
# row. Offering those donors here would put a donor the domain's cohort has never
# heard of into a domain-scoped selector, and every other thing keyed on the
# selection (status ordering, donor colour, drill-down) has no record of it.
#
# An earlier version did union them and then intersected with the cohort anyway,
# which removed exactly the donors the union had just added: code that looked
# like it reached pLN-only donors and did not. They are not silently dropped now
# -- .ie_distance_pln_outside_cohort names them so the card can say they are
# withheld and why.
.ie_distance_donors <- function(layer, cohort) {
  if (!ie_is_measured(layer$value)) return(character())
  donors <- unique(layer$data$donor_id)
  if (!is.null(cohort) && is.function(cohort$donor_ids)) {
    donors <- intersect(donors, .ie_distance_norm(cohort$donor_ids()))
  }
  sort(donors)
}

# Donors carrying a pLN section in the artifact that the active domain's cohort
# does not contain. Named rather than dropped: the sections exist and were
# measured, and the reason they are not offered is a domain boundary, not an
# absence of data.
.ie_distance_pln_outside_cohort <- function(layer, cohort) {
  if (!ie_is_measured(layer$value) || !is.data.frame(layer$pln) || !nrow(layer$pln)) {
    return(character())
  }
  donors <- unique(layer$pln$donor_id)
  if (!is.null(cohort) && is.function(cohort$donor_ids)) {
    donors <- setdiff(donors, .ie_distance_norm(cohort$donor_ids()))
  } else {
    donors <- character()
  }
  sort(donors)
}

# The pLN rows for one donor. Never joined to the pancreas arm by anything but
# donor_id and donor_status, and never given a distance.
.ie_distance_pln_rows <- function(layer, donor) {
  if (!ie_is_measured(layer$value) || !is.data.frame(layer$pln)) return(NULL)
  if (is.null(donor) || !length(donor)) return(NULL)
  rows <- layer$pln[layer$pln$donor_id == donor, , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  rows
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

#: Past this the gradient is reporting the rest of the section, not the peri-islet environment.
#: The maintainer set 500 um on 2026-09-18; the slider only narrows it, never widens it, so no
#: reading of this card can reach a distance the cohort was not asked about.
X1_MAX_DIST_UM <- 500


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
      shiny::sliderInput(session$ns("x1_max_dist"), "Maximum distance from islet edge (\u00b5m)",
                         min = 50, max = X1_MAX_DIST_UM, value = X1_MAX_DIST_UM,
                         step = 50, width = "320px"),
      shiny::plotOutput(session$ns("x1_gradient_plot")),
      shiny::tags$p(class = "ie-card-note", style = "font-size: 12px",
                    "Resolved by: ",
                    if (length(value$provenance)) paste(value$provenance, collapse = ", ")
                    else "no immune calls observed", "."),
      ie_explain(
        shiny::tags$p("Existing immune calls / all QC-eligible, non-excluded cells in each bin. ",
                      "Signed distance is negative inside the islet and positive outside; zero is ",
                      "retained and has its own bin."),
        label = "How the gradient is measured"),
      shiny::downloadButton(session$ns("x1_gradient_download"), "Download plot")
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
    # Clip BEFORE the fraction is read, not with xlim afterwards: a scale limit drops rows
    # silently and ggplot would still draw a line across the gap.
    cap <- suppressWarnings(as.numeric(input$x1_max_dist))
    if (!length(cap) || !is.finite(cap)) cap <- X1_MAX_DIST_UM
    cap <- min(cap, X1_MAX_DIST_UM)
    data <- data[is.finite(data$distance_bin_um) & abs(data$distance_bin_um) <= cap, , drop = FALSE]
    ie_attach_frame(data[order(data$distance_bin_um), ], ie_frame_of(source))
  })
  plot <- shiny::reactive({
    data <- plot_data()
    reference <- ie_attach_frame(data.frame(distance_bin_um = 0), ie_frame("signed_distance_um"))
    ie_assert_same_frame(data, reference, "peri-islet immune gradient plot")
    ggplot2::ggplot(data, ggplot2::aes(distance_bin_um, immune_fraction)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linetype = "dashed") +
      ggplot2::geom_line() + ggplot2::geom_point() +
      # Linear in micrometres. pseudo_log compressed 0-500 um into the left third of the panel
      # and spent the rest on distances this card no longer shows, so the peri-islet band --
      # the thing being read -- was the part that was squashed.
      ggplot2::scale_x_continuous(
        breaks = scales::breaks_width(100),
        limits = c(min(data$distance_bin_um, na.rm = TRUE),
                   max(data$distance_bin_um, na.rm = TRUE)),
        expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::scale_y_continuous(labels = scales::label_percent()) +
      ggplot2::labs(x = "Signed distance to islet edge (µm)", y = "Immune fraction of eligible cells",
                    title = paste("Donor", data$donor_id[[1L]], "—", data$donor_status[[1L]]),
                    subtitle = .ie_distance_pair(selected()$data)[[1L]],
                    caption = paste("Bin labels are outer edges; width", layer()$bin_width,
                                    "µm. Zero has its own bin. Existing typing, no refitted thresholds.")) +
      .ie_distance_theme()
  })
  output$x1_gradient_plot <- shiny::renderPlot(plot())
  output$x1_gradient_download <- shiny::downloadHandler(
    filename = function() "peri-islet-immune-gradient.pdf",
    content = function(file) ggplot2::ggsave(file, plot = plot(), width = 8, height = 5)
  )
  invisible(list(active = active, selected = selected, plot = plot, data = plot_data))
}
