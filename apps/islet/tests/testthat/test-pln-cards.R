# The two pLN cards: card_n1_sections.R and card_n2_pln_immune.R.
#
# What is asserted here, and why each assertion exists:
#
# * Both halves of each pair exist and the Lymph Node registry reaches them. A
#   registry entry with no card is the defect that shipped twice on these tabs
#   -- `ledger` and, before this, `n1_sections` -- and it renders unavailable on
#   every donor while looking like data.
# * Every missing input produces a NAMED note rather than an error or a blank
#   panel, on the card's OWN output. A note that lands only in a child output is
#   a note the card does not carry.
# * The exceptions n1 exists to surface are surfaced from the artifacts: a node
#   annotation that bounds no region, a node that is really the whole slide, a
#   membership count that disagrees with QuPath's own, and a donor the node
#   layer holds that the ingested cell root does not.
# * n2 keeps the two modalities apart. Its gap is paired within donor and
#   carries its own donor count, and no number on it is the sum of the two
#   lanes.
#
# The fixtures are synthetic and their donor ids belong to no cohort, so nothing
# here can pass by reproducing a real result. No measured value of the real
# cohort is a target of any assertion.

.pln_card_ids <- c("n1_sections", "n2_pln_immune")

.pln_env <- function() {
  env <- new.env(parent = globalenv())
  if (requireNamespace("shiny", quietly = TRUE)) {
    for (name in c("NS", "tagList", "div", "tags", "moduleServer", "reactive")) {
      env[[name]] <- getExportedValue("shiny", name)
    }
  }
  # Every card the Lymph Node registry mounts, plus what those cards read:
  # card_n2_pln_immune.R reads the grid layer through card_g1_grid.R's reader
  # rather than keeping a second copy of that contract, card_s3_compartments.R
  # reads card_x1_gradient.R's, and the placeholder the registry falls back to
  # is defined once in mod_crossmodal_ui.R.
  for (file in c("card_x1_gradient.R", "card_s3_compartments.R", "card_g1_grid.R",
                 "card_n1_sections.R", "card_n2_pln_immune.R",
                 "mod_crossmodal_ui.R", "mod_lymphnode_ui.R", "mod_lymphnode_server.R")) {
    sys.source(.find_source(file.path("R", file)), envir = env)
  }
  env
}

.pln_html <- function(value) {
  if (is.list(value) && !is.null(value$html)) return(as.character(value$html))
  paste(as.character(value), collapse = " ")
}

.pln_missing <- function(prefix) file.path(tempfile(prefix), "never-written")

# ---------------------------------------------------------------------------
# Fixtures. Synthetic throughout; the donor ids are outside every cohort.

.pln_n1_fixture <- function(dir, root, with_domain = TRUE) {
  layers <- file.path(dir, "layers")
  dir.create(file.path(layers, "pln"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(layers, "sections"), recursive = TRUE, showWarnings = FALSE)
  record <- file.path(dir, "pln_node_annotations.csv")
  writeLines(c(
    "# Synthetic record. These donor ids are in no cohort.",
    paste("image_file,donor_id,roi,qupath_entry,node_classes,n_nodes,n_cells,excluded",
          "recorded_by,recorded_on,excluded_reason", sep = ","),
    "A.qptiff,9001,pln_1,1,pLN_top;pLN_bottom,2,300,False,tester,2026-09-18,",
    "B.qptiff,9002,pln_1,2,Tissue,1,112,False,tester,2026-09-18,",
    "C.qptiff,9003,pln_1,3,pLN_top;pLN_bottom,2,50,False,tester,2026-09-18,"), record)

  nodes <- data.frame(
    donor_id = c("9001", "9001", "9002", "9003", "9003"),
    image = c("A.qptiff", "A.qptiff", "B.qptiff", "C.qptiff", "C.qptiff"),
    qupath_entry = c(1L, 1L, 2L, 3L, 3L),
    node_index = c(0L, 1L, 0L, 0L, 1L),
    node_name = c("pLN_top", "pLN_bottom", "Tissue", "pLN_top", "pLN_bottom"),
    annotation_geom = c("Polygon", "Polygon", "Polygon", "LineString", "LineString"),
    delimiter_kind = c("rectangle", "rectangle", "outline", "none", "none"),
    area_um2 = c(1e6, 2e6, 3e6, NA_real_, NA_real_),
    state = c("measured", "measured", "measured", "unavailable", "unavailable"),
    reason = c("", "", "",
               "annotation_is_linestring_with_zero_area_it_bounds_no_region",
               "annotation_is_linestring_with_zero_area_it_bounds_no_region"),
    n_cells = c(200, 100, 112, NA_real_, NA_real_),
    node_scope = c("node", "node", "slide", NA_character_, NA_character_),
    node_scope_reason = c("one of 2 usable node(s) on this image", "one of 2 usable node(s) on this image",
                          "the only usable node on this image, holds all of its assigned cells.",
                          "not_applicable_node_bounds_no_region",
                          "not_applicable_node_bounds_no_region"),
    qupath_check = c("reproduced", "reproduced", "differs", "not_applicable", "not_applicable"),
    delta_vs_qupath = c(0, 0, -8, NA_real_, NA_real_),
    stringsAsFactors = FALSE)
  arrow::write_parquet(nodes, file.path(layers, "pln", "nodes.parquet"))

  jsonlite::write_json(list(ingested_cells_coverage = list(
    root = "/synthetic/cells_pln",
    donors_in_that_root = c("9001", "9002"),
    in_this_layer_but_not_ingested = c("9003"),
    ingested_but_not_in_this_layer = character(),
    agrees = FALSE)), file.path(layers, "pln", "MANIFEST.json"), auto_unbox = TRUE)

  sections <- data.frame(
    donor_id = c("9001", "9001", "9002"),
    roi = c("panc", "pln_1", "pln_1"),
    modality = c("xenium", "xenium", "xenium"),
    image_file = c("dir/x.ome.tif", "dir/x.ome.tif", "dir/y.ome.tif"),
    xenium_serial = c("0001", "0001", "0002"),
    pixel_size_um = c(NA_real_, NA_real_, NA_real_),
    frame = c("xenium_um:9001:panc", "xenium_um:9001:pln_1", "xenium_um:9002:pln_1"),
    stringsAsFactors = FALSE)
  table <- arrow::arrow_table(sections)
  if (with_domain) table$metadata <- list(source_domain_root = root)
  arrow::write_parquet(table, file.path(layers, "sections", "data.parquet"))
  list(layers = layers, record = record)
}

.pln_n2_fixture <- function(dir, root) {
  layers <- file.path(dir, "layers")
  grid <- file.path(layers, "grid_composition")
  dir.create(grid, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(
    layer = "grid_composition", source_domain_root = root,
    grid = data.frame(
      tissue = c("pln", "pln"), modality = c("protein", "RNA"),
      root = c("/synthetic/typing/pln/pheno", "/synthetic/xenium/release_pln"),
      stringsAsFactors = FALSE)),
    file.path(grid, "manifest.json"), auto_unbox = TRUE)

  # A pancreas row of the same layer. It must never reach this card: the two
  # tissues are different blocks and this card reads the lymph node only.
  presence <- data.frame(
    tissue = "pln",
    modality = c("protein", "RNA", "protein", "RNA"),
    donor_id = c("9001", "9001", "9002", "9002"),
    state = c("measured", "measured", "unavailable", "measured"),
    reason = c("", "",
               "the clustering track produced no cluster for this donor; its export carries no intensity measurement",
               ""),
    n_cells_assessed = c(1000L, 500L, NA_integer_, 400L),
    n_no_nucleus_excluded = c(0L, 7L, NA_integer_, 5L),
    n_files = c(1L, 1L, 0L, 1L),
    source = c("/synthetic/typing/pln/pheno", "/synthetic/xenium/release_pln",
               "/synthetic/typing/pln/pheno", "/synthetic/xenium/release_pln"),
    method = c("cluster-named immune classes", "immune partition with a T/NK split",
               "cluster-named immune classes", "immune partition with a T/NK split"),
    citation_conflict = "",
    disease_status = c("ND", "ND", "T1D", "T1D"),
    stringsAsFactors = FALSE)
  presence <- rbind(presence, data.frame(
    tissue = "pancreas", modality = "protein", donor_id = "9003", state = "measured",
    reason = "", n_cells_assessed = 999999L, n_no_nucleus_excluded = 3L, n_files = 1L,
    source = "/synthetic/typing/cells", method = "subpanel classification",
    citation_conflict = "", disease_status = "ND", stringsAsFactors = FALSE))
  utils::write.csv(presence, file.path(grid, "presence.csv"), row.names = FALSE)

  bucket <- function(modality, donor, lineage, kind, n, denom, immune, state = "measured",
                     reason = "") {
    data.frame(tissue = "pln", modality = modality, donor_id = donor,
               lineage_common = lineage, bucket_kind = kind, n_cells = n,
               denom_assessed = denom,
               frac_of_assessed = if (is.na(n)) NA_real_ else n / denom,
               state = state, reason = reason, pheno_resolution = "subclass",
               xenium_resolution = "yes", is_immune = immune, resolution_conflict = "",
               disease_status = if (donor == "9001") "ND" else "T1D",
               stringsAsFactors = FALSE)
  }
  composition <- rbind(
    bucket("protein", "9001", "T_NK", "lineage", 100, 1000, "yes"),
    bucket("protein", "9001", "Myeloid", "lineage", 200, 1000, "yes"),
    bucket("protein", "9001", "Non_immune", "negative", 700, 1000, ""),
    bucket("RNA", "9001", "T_NK", "lineage", 250, 500, "yes"),
    bucket("RNA", "9001", "Myeloid", "lineage", 50, 500, "yes"),
    bucket("RNA", "9001", "Non_immune", "negative", 200, 500, ""),
    bucket("protein", "9002", "T_NK", "lineage", NA_real_, NA_real_, "yes", "blank",
           "the clustering track produced no cluster for this donor"),
    bucket("RNA", "9002", "T_NK", "lineage", 200, 400, "yes"),
    bucket("RNA", "9002", "Myeloid", "lineage", 40, 400, "yes"),
    bucket("RNA", "9002", "Non_immune", "negative", 160, 400, ""))
  pancreas <- bucket("protein", "9003", "T_NK", "lineage", 111111, 999999, "yes")
  pancreas$tissue <- "pancreas"
  composition <- rbind(composition, pancreas)
  utils::write.csv(composition, file.path(grid, "composition.csv"), row.names = FALSE)

  verdict <- file.path(dir, "pln_lineage_summary.csv")
  utils::write.csv(data.frame(
    tissue = c("pancreas", "pln"), lineage = c("T_NK", "T_NK"), n_donors = c(9L, 11L),
    median_share_protein_pct = c(10.5, 31.1), median_share_rna_pct = c(12.5, 50.8),
    median_diff_pp = c(2.0, -22.7), verdict = c("CONCORDANT", "MAY NOT CARRY A CROSS-MODALITY CLAIM"),
    verdict_reason = c("within the threshold fixed before measurement",
                       "median offset exceeds the threshold fixed before measurement"),
    stringsAsFactors = FALSE), verdict, row.names = FALSE)
  manifest <- file.path(dir, "concordance_manifest.json")
  jsonlite::write_json(list(pln = list(protein_source = "/synthetic/typing/pln/pheno",
                                       rna_source = "/synthetic/xenium/release_pln_older")),
                       manifest, auto_unbox = TRUE)
  list(layers = layers, verdict = verdict, manifest = manifest)
}

# ---------------------------------------------------------------------------

test_that("both pLN cards ship a ui/server pair the Lymph Node registry reaches", {
  env <- .pln_env()
  ids <- vapply(env$IE_LYMPHNODE_CARDS, `[[`, "", "id")
  for (id in .pln_card_ids) {
    expect_true(id %in% ids, info = paste(id, "is registered"))
    expect_equal(sum(ids == id), 1L, info = paste(id, "is registered exactly once"))
  }
  # Derived over the whole registry: a card added here without a builder fails
  # in this suite rather than in front of a reader.
  for (id in ids) {
    expect_true(is.function(env[[paste0("card_", id, "_ui")]]),
                info = paste0("card_", id, "_ui"))
    expect_true(is.function(env[[paste0("card_", id, "_server")]]),
                info = paste0("card_", id, "_server"))
  }
  for (id in .pln_card_ids) {
    args <- formals(env[[paste0("card_", id, "_server")]])
    # The tab key is passed in, never guessed: a card that guessed its own key
    # shipped with a key matching nothing and could not draw.
    expect_identical(tail(names(args), 3L), c("active_tab", "tab_key", "..."), info = id)
    expect_null(args$active_tab, info = id)
    expect_null(args$tab_key, info = id)
  }
})

test_that("the Lymph Node tab mounts a slot for each pLN card", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("htmltools")
  env <- .pln_env()
  html <- htmltools::renderTags(env$lymphnode_ui("lymphnode"))$html
  for (id in .pln_card_ids) {
    expect_match(html, paste0('data-card-id="', id, '"'), fixed = TRUE)
    # Not the placeholder: the pair exists, so the registry must reach it.
    expect_false(grepl(paste0("unavailable(not_built): no card_", id, "_ui()"), html, fixed = TRUE))
  }
})

test_that("neither pLN card carries a pancreas-to-lymph-node distance", {
  # A lymph node has no islet in it, so a signed distance to one is undefined
  # there rather than zero. The check is on the artifact columns that would
  # carry such a number, not on the English word.
  for (file in c("card_n1_sections.R", "card_n2_pln_immune.R")) {
    source_text <- paste(readLines(.find_source(file.path("R", file)), warn = FALSE),
                         collapse = "\n")
    for (column in c("signed_dist_islet_um", "distance_bin_um", "cell_distance",
                     "neighbourhood_")) {
      expect_false(grepl(column, source_text, fixed = TRUE),
                   info = paste(file, "must not read", column))
    }
  }
})

test_that("card_n1_sections_read names each absent input instead of failing", {
  env <- .pln_env()
  layer <- env$card_n1_sections_read(layers_dir = .pln_missing("n1-layers-"),
                                     record_path = .pln_missing("n1-record-"),
                                     active_domain_root = tempdir())
  for (part in c("record", "nodes", "coverage", "sections")) {
    value <- layer[[part]]$value
    expect_s3_class(value, "ie_avail")
    expect_false(ie_is_measured(value), info = part)
    expect_true(nzchar(value$reason), info = part)
    expect_true(value$state %in% IE_STATES, info = part)
  }
})

test_that("card_n2_pln_immune_read names each absent input instead of failing", {
  env <- .pln_env()
  layer <- env$card_n2_pln_immune_read(layers_dir = .pln_missing("n2-layers-"),
                                       active_domain_root = tempdir(),
                                       verdict_path = .pln_missing("n2-verdict-"),
                                       verdict_manifest_path = .pln_missing("n2-manifest-"))
  for (part in c("grid", "nodes", "verdict")) {
    value <- layer[[part]]$value
    expect_s3_class(value, "ie_avail")
    expect_false(ie_is_measured(value), info = part)
    expect_true(nzchar(value$reason), info = part)
  }
})

test_that("each pLN card degrades to a named note on its own output when its layer is absent", {
  skip_if_not_installed("shiny")
  env <- .pln_env()
  withr::local_envvar(c(ISLET_LAYERS_DIR = .pln_missing("missing-pln-layers-")))
  cards <- list(
    n1_sections = list(output = "n1_body",
                       args = list(layers_dir = .pln_missing("n1-layers-"),
                                   record_path = .pln_missing("n1-record-"),
                                   active_domain_root = tempdir())),
    n2_pln_immune = list(output = "n2_body",
                         args = list(layers_dir = .pln_missing("n2-layers-"),
                                     active_domain_root = tempdir(),
                                     verdict_path = .pln_missing("n2-verdict-"),
                                     verdict_manifest_path = .pln_missing("n2-manifest-"))))
  for (id in names(cards)) {
    spec <- cards[[id]]
    run_card <- function(input, output, session) {
      card <- do.call(env[[paste0("card_", id, "_server")]],
                      c(list(input = input, output = output, session = session), spec$args))
    }
    shiny::testServer(run_card, {
      session$flushReact()
      expect_true(card$active(), info = id)
      html <- .pln_html(output[[spec$output]])
      expect_true(nzchar(html), info = id)
      # The guard is named so the message can be traced to the code that
      # emitted it, and the artifact is named so the reader knows what is
      # missing rather than only that something is.
      expect_match(html, paste0("[", id, "]"), fixed = TRUE, info = id)
      expect_match(html, "Artifact:", fixed = TRUE, info = id)
      # Every panel of the card behaves the same way: a named note, never an
      # error and never an empty frame that reads as "nothing was found".
      if (identical(id, "n2_pln_immune")) {
        expect_null(card$gap())
        expect_match(.pln_html(output$n2_gap), "[n2_pln_immune]", fixed = TRUE)
        expect_match(.pln_html(output$n2_verdict), "Artifact:", fixed = TRUE)
      }
    })
  }
})

test_that("n1 names every exception it measures rather than hiding it", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("arrow")
  env <- .pln_env()
  dir <- tempfile("n1-fixture-")
  dir.create(dir)
  root <- normalizePath(dir, mustWork = FALSE)
  fixture <- .pln_n1_fixture(dir, root)
  run_card <- function(input, output, session) {
    card <- env$card_n1_sections_server(input, output, session,
                                        layers_dir = fixture$layers,
                                        record_path = fixture$record,
                                        active_domain_root = root)
  }
  shiny::testServer(run_card, {
    session$flushReact()
    html <- .pln_html(output$n1_body)
    # The node annotation that bounds no region: unusable, and it holds no
    # cells -- not zero cells.
    expect_match(html, "annotation_is_linestring_with_zero_area_it_bounds_no_region", fixed = TRUE)
    # The whole-slide annotation wearing a node's name: the layer's own scope
    # reason is carried, and what it means for pooling is said out loud.
    expect_match(html, "the only usable node on this image", fixed = TRUE)
    expect_match(html, "Pooling it with per-node values from other slides mixes the two.",
                 fixed = TRUE)
    # Membership is geometric on the object_id; QuPath's own count is the
    # cross-check, and where the two differ the difference is stated.
    expect_match(html, "differ by -8 cell(s)", fixed = TRUE)
    # The donor the node layer holds that the ingested cell root does not.
    expect_match(html, "/synthetic/cells_pln", fixed = TRUE)
    expect_match(html, "ingest_pln_cells.py", fixed = TRUE)
    # A node with no usable geometry has no assigned-cell count at all.
    expect_match(html, "unavailable", fixed = TRUE)
    # Lane C's attribute is absent from this layer, so it is named, not faked.
    expect_match(html, "sections_layer_records_no_section_relationship", fixed = TRUE)
    # The section table is read, which is the first reader this layer has had.
    expect_match(html, "y.ome.tif", fixed = TRUE)
    expect_identical(card$record_donors(), c("9001", "9002", "9003"))

    # The record itself is a drawn table, not prose: the recorded node classes
    # beside what the layer measured from them.
    recorded <- .pln_html(output$n1_record_table)
    expect_match(recorded, "pLN_top;pLN_bottom", fixed = TRUE)
    expect_match(recorded, "Tissue", fixed = TRUE)
    # The donor whose annotations bound no region has no assigned-cell count at
    # all -- not a zero.
    expect_match(recorded, "unavailable", fixed = TRUE)
    expect_false(grepl(">0<", recorded, fixed = TRUE))
  })
})

test_that("n1 refuses the sections layer by name when it records no domain", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("arrow")
  env <- .pln_env()
  dir <- tempfile("n1-nodomain-")
  dir.create(dir)
  root <- normalizePath(dir, mustWork = FALSE)
  fixture <- .pln_n1_fixture(dir, root, with_domain = FALSE)
  layer <- env$card_n1_sections_read(layers_dir = fixture$layers,
                                     record_path = fixture$record,
                                     active_domain_root = root)
  expect_false(ie_is_measured(layer$sections$value))
  expect_identical(layer$sections$value$reason, "sections_layer_domain_not_recorded")
  # The five donors the two domains share mean a matching donor id is not a
  # licence to read, so a mismatch is refused the same way.
  mismatch <- env$card_n1_sections_read(layers_dir = fixture$layers,
                                        record_path = fixture$record,
                                        active_domain_root = tempfile("other-domain-"))
  expect_false(ie_is_measured(mismatch$sections$value))
})

test_that("n2 shows each lane on its own and states the gap with its own donor count", {
  skip_if_not_installed("shiny")
  env <- .pln_env()
  dir <- tempfile("n2-fixture-")
  dir.create(dir)
  root <- normalizePath(dir, mustWork = FALSE)
  fixture <- .pln_n2_fixture(dir, root)
  run_card <- function(input, output, session) {
    card <- env$card_n2_pln_immune_server(input, output, session,
                                          layers_dir = fixture$layers,
                                          active_domain_root = root,
                                          verdict_path = fixture$verdict,
                                          verdict_manifest_path = fixture$manifest)
  }
  shiny::testServer(run_card, {
    session$flushReact()
    panels <- .pln_html(output$n2_modalities)
    # Both lanes are drawn, each with its OWN denominator stated on it.
    expect_match(panels, "protein", fixed = TRUE)
    expect_match(panels, "RNA", fixed = TRUE)
    expect_match(panels, "1,000 cells this lane assessed", fixed = TRUE)
    expect_match(panels, "500 cells this lane assessed", fixed = TRUE)

    # Each lane's buckets are a table of its own, mounted per modality, so
    # neither lane's counts can land in the other's column.
    protein <- .pln_html(output$n2_table_protein)
    rna <- .pln_html(output$n2_table_RNA)
    expect_match(protein, "10.00%", fixed = TRUE)
    expect_match(rna, "50.00%", fixed = TRUE)
    both <- paste(panels, protein, rna)
    # NOT POOLED: 1,000 + 500 is a number this card must never produce, and
    # neither is 100 + 250 T_NK cells.
    expect_false(grepl("1,500", both, fixed = TRUE))
    expect_false(grepl("350", both, fixed = TRUE))
    # NOT ACROSS TISSUES: the pancreas rows of the same layer are another
    # tissue, another block and another section, and they reach neither the
    # donor list nor any panel here.
    expect_false("9003" %in% card$donors())
    expect_false(grepl("999,999", both, fixed = TRUE))
    expect_false(grepl("111,111", both, fixed = TRUE))

    gap <- card$gap()
    expect_true(is.data.frame(gap))
    expect_true(all(gap$`Donors measured in both` == 1L))
    row <- gap[gap$Bucket == "T_NK", , drop = FALSE]
    expect_equal(nrow(row), 1L)
    shares <- sort(unlist(row[grepl("^Median %", names(row))]))
    expect_equal(unname(shares), c(10, 50))
    expect_equal(abs(row[[grep("^Median gap", names(row))]]), 40)

    # The donor whose protein lane produced nothing: its reason is shown, and
    # nothing is imputed for it. It is also absent from the paired count above.
    session$setInputs(n2_donor = "9002")
    session$flushReact()
    other <- .pln_html(output$n2_modalities)
    expect_match(other, "no cluster for this donor", fixed = TRUE)

    verdict <- .pln_html(output$n2_verdict)
    expect_match(verdict, "MAY NOT CARRY A CROSS-MODALITY CLAIM", fixed = TRUE)
    # The carried verdict brings its own donor count, which is not the layer's.
    expect_match(verdict, ">11<", fixed = TRUE)
    # And it names what it was measured against, so a verdict measured on one
    # release is not read as evidence about another.
    expect_match(verdict, "/synthetic/xenium/release_pln_older", fixed = TRUE)
    expect_match(verdict, "/synthetic/xenium/release_pln", fixed = TRUE)
    # The pancreas row of the same export is not a lymph node statement.
    expect_false(grepl("CONCORDANT", verdict, fixed = TRUE))
  })
})
