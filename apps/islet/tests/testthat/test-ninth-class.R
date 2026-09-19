# The ninth entry in the run cohort's display vocabulary is an ABSENCE, not a cell type.
#
# `Unresolved` is what the refinement pass writes -- producer label
# `Unresolved (no lineage evidence)` -- on a cell it demoted out of Exocrine. Exocrine is the
# complement class, so there is no complement beneath it and the cell is left with no lineage
# evidence at all.
#
# Everything here guards the SAME failure from a different side: an absence being read as a
# ninth population. It has to be drawn, so it has to be in `display_order` and in every
# palette; being in `display_order` is exactly what makes every consumer that walks that list
# treat it as a cell type. So the vocabulary carries it, and two declarations carry the fact
# that it is not one -- `absence: true` in the DOMAIN, and a colour that is achromatic while
# every cell type's colour has a hue.

.run20 <- function() ie_load_domain(.find_repo_file("apps/islet/config/run20_default.yml"))

#: The eight classes as they stood before 2026-09-18, in order. A literal on purpose: this is
#: the thing that must NOT have changed, so deriving it from the file under test would prove
#: nothing.
.EIGHT_CLASSES <- c("Beta_cell", "Alpha_cell", "Delta_cell", "Exocrine",
                    "Blood_Vessel", "Lymphatic", "Neural", "Immune")

.is_achromatic <- function(hex) {
  rgb <- grDevices::col2rgb(hex)
  apply(rgb, 2L, function(ch) max(ch) - min(ch)) == 0L
}


test_that("the eight classes are untouched and the absence is appended last", {
  domain <- .run20()
  order <- domain_phenotype_order(domain)

  expect_length(order, 9L)
  # Not `setequal`: the order of the eight is part of the display contract.
  expect_identical(order[seq_along(.EIGHT_CLASSES)], .EIGHT_CLASSES)
  expect_identical(order[[9L]], "Unresolved")

  colors <- domain_phenotype_colors(domain)
  expect_setequal(names(colors), order)
  # The eight kept their colours; the absence did not take one of them.
  expect_identical(unname(colors[.EIGHT_CLASSES]),
                   c("#E63946", "#F4A261", "#E9C46A", "#A8DADC",
                     "#BC4749", "#84A98C", "#B5838D", "#3A5A8C"))
})


test_that("the absence is declared, never inferred from its spelling", {
  domain <- .run20()

  absences <- ie_absence_phenotypes(domain_phenotype_residuals(domain))
  expect_identical(absences, "Unresolved")

  # The declaration carries a reason. An absence without a recorded reason is
  # indistinguishable from an oversight -- the same rule ie_avail() enforces.
  entry <- Filter(function(e) identical(e[["name"]], "Unresolved"),
                  domain_phenotype_residuals(domain))[[1L]]
  expect_true(isTRUE(entry[["absence"]]))
  expect_true(nzchar(entry[["description"]]))
  expect_match(entry[["description"]], "Unresolved (no lineage evidence)", fixed = TRUE)

  # THE OTHER DIRECTION, which a string match would get wrong: the legacy CODEX vocabulary
  # declares `Structural` and `Unknown` as residuals, and neither is a recorded non-call. They
  # are real upstream display classes, so nothing here may promote them to absences.
  codex <- ie_load_domain(.find_repo_file("apps/islet/config/islet_default.yml"))
  expect_true(length(domain_phenotype_residuals(codex)) > 0L)
  expect_identical(ie_absence_phenotypes(domain_phenotype_residuals(codex)), character(0))
})


test_that("the absence colour is ieCore's, and reads as absence rather than as a cell type", {
  domain <- .run20()
  colors <- domain_phenotype_colors(domain)

  # One authority. YAML cannot call R, so the hex is repeated in the config; this is the only
  # thing stopping the two drifting.
  expect_identical(unname(colors[["Unresolved"]]), ie_absence_class_color())
  expect_identical(ie_absence_class_color(), unname(ie_state_palette()[["uncalled"]]))

  # No hue at all, while every one of the eight carries one. That is the visual statement, and
  # it is asserted rather than left to a reviewer's eye.
  expect_true(.is_achromatic(colors[["Unresolved"]]))
  expect_false(any(.is_achromatic(colors[.EIGHT_CLASSES])))

  # It cannot be mistaken for a class, or for one of the other availability states.
  expect_false(colors[["Unresolved"]] %in% colors[.EIGHT_CLASSES])
  others <- ie_state_palette()[setdiff(names(ie_state_palette()), "uncalled")]
  expect_false(ie_absence_class_color() %in% others)
})


test_that("every selectable palette covers the run cohort vocabulary, absence included", {
  domain <- .run20()
  previous <- ie_domain()
  on.exit(ie_set_domain(previous), add = TRUE)
  ie_set_domain(domain)

  env <- new.env(parent = globalenv())
  sys.source(.find_repo_file("apps/islet/R/drilldown_helpers.R"), envir = env)
  expected <- domain_phenotype_order(domain)
  expect_true(all(vapply(
    env$PHENOTYPE_PALETTES,
    function(palette) all(expected %in% names(palette)),
    logical(1)
  )))
  # The alternate palettes are keyed on the CODEX spellings and fall back to the DOMAIN for
  # anything they do not name, so the absence must be the SAME colour in all four. A palette
  # that recoloured it would be inventing a ninth cell type per palette.
  expect_true(all(vapply(
    env$PHENOTYPE_PALETTES,
    function(palette) identical(unname(palette[["Unresolved"]]), ie_absence_class_color()),
    logical(1)
  )))
})


test_that("00_globals.R asks the DOMAIN instead of matching the label", {
  domain <- .run20()
  previous <- ie_domain()
  on.exit(ie_set_domain(previous), add = TRUE)
  ie_set_domain(domain)

  # Evaluate the REAL expressions out of 00_globals.R rather than restating them: the file
  # cannot be sourced in the suite (it attaches shiny, plotly and sf), but the three
  # definitions this lane added are self-contained.
  wanted <- c("PHENOTYPE_ABSENCE_CLASSES", "is_absence_phenotype", "phenotype_cell_type_order")
  env <- new.env(parent = globalenv())
  found <- character(0)
  for (expr in parse(.find_repo_file("apps/islet/R/00_globals.R"))) {
    if (!is.call(expr) || !identical(as.character(expr[[1L]]), "<-")) next
    name <- as.character(expr[[2L]])
    if (length(name) != 1L || !name %in% wanted) next
    found <- c(found, name)
    eval(expr, envir = env)
  }
  expect_setequal(found, wanted)

  expect_identical(env$PHENOTYPE_ABSENCE_CLASSES, "Unresolved")
  expect_identical(env$is_absence_phenotype(c("Beta_cell", "Unresolved", "Immune")),
                   c(FALSE, TRUE, FALSE))
  # Total, not a lookup: an unknown class is a cell type, so a caller that forgets to guard
  # gets the old behaviour instead of a silent reclassification.
  expect_false(env$is_absence_phenotype("SomethingNew"))
  expect_identical(env$phenotype_cell_type_order(), .EIGHT_CLASSES)
})


test_that("a producer's honest non-call is not read as a measured identity", {
  # This is the leak on the availability side. The test that shipped matched the exact string
  # "Unavailable" and nothing else, so every label below read `measured`.
  got <- ie_identity_state(c(
    "Beta_cell",
    "Unresolved (no lineage evidence)",   # protein refinement
    "Uncalled_no_naming_evidence",        # the RNA analogue
    "Uncalled_no_enrichment",
    "Uncalled_below_cluster_qc",
    "Uncalled_ambiguous",                 # pLN release
    "Immune_uncalled",                    # the `*_uncalled` spelling, same file
    "Unavailable_no_nucleus",
    "Unavailable",
    "Non_immune",                         # a measured negative IS an answer
    "",
    NA_character_
  ))

  expect_identical(as.character(got), c(
    "measured",
    "uncalled", "uncalled", "uncalled", "uncalled", "uncalled", "uncalled",
    "unavailable", "unavailable",
    "measured",
    "blank", "blank"
  ))
  expect_identical(levels(got), IE_STATES)
})


test_that("the sixth availability state is complete wherever the five were", {
  expect_true("uncalled" %in% IE_STATES)
  # Appended, not inserted: these are factor levels and a card legend walks them in order.
  expect_identical(IE_STATES[seq_len(5L)],
                   c("measured", "blank", "unavailable", "withheld", "failed"))

  # Every state has a label and a colour, or a card renders an untraceable blank.
  for (state in IE_STATES) expect_true(nzchar(ie_state_label(state)))
  expect_setequal(names(ie_state_palette()), IE_STATES)
  expect_false(anyDuplicated(ie_state_palette()) > 0L)

  # It is an absence, so it obeys the reason rule the other absences obey.
  expect_error(ie_avail("uncalled"), "requires a reason")
  value <- ie_avail("uncalled", reason = "no lineage evidence")
  expect_false(ie_is_measured(value))
  expect_match(ie_state_text(value), "no lineage evidence", fixed = TRUE)
  # A hole and a measurement must never render the same string.
  expect_false(identical(ie_state_text(value),
                         ie_state_text(ie_avail("measured", n = 0L))))
})
