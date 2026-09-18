# A selector offers only what it can draw.
#
# The maintainer's rule, 2026-09-18: "Only donors that have data ready should be available to choose
# in the app ... The same goes for similar situations elsewhere."
#
# The post-Xenium IF overlay offered all 25 donors of its availability layer and rendered a named
# note for 23 of them. An unavailable state is the right thing to record and the wrong thing to put
# in a dropdown.

test_that("ie_selectable_donors offers measured donors and nothing else", {
  data <- data.frame(
    donor_id = c("6436", "6442", "6539", "115", "6591"),
    state    = c("measured", "measured", "unavailable", "failed", "blank"),
    stringsAsFactors = FALSE)

  expect_identical(ie_selectable_donors(data), c("6436", "6442"))

  # `blank` is NOT offered: it means the lane ran and found nothing, which is a result to record,
  # not something to plot. `failed` is not offered either -- 115's registration missed its gate at
  # 7.198 um against 3.67, and offering it would invite reading an overlay that is known wrong.
  expect_false(any(c("115", "6591") %in% ie_selectable_donors(data)))

  # Opting a state in is explicit, never a default.
  expect_identical(ie_selectable_donors(data, states = c("measured", "blank")),
                   c("6436", "6442", "6591"))
})

test_that("ie_selectable_donors is intersected with the active domain's roster", {
  data <- data.frame(donor_id = c("6436", "6442", "112"),
                     state = rep("measured", 3), stringsAsFactors = FALSE)
  # 112 is measured in the artifact and outside this domain. The two domains share five donors and
  # no vocabulary, so a shared artifact must not smuggle one across.
  cohort <- list(donor_ids = function() c("6436", "6442"))
  expect_identical(ie_selectable_donors(data, cohort = cohort), c("6436", "6442"))
})

test_that("what is left out is stated once, with its states counted", {
  data <- data.frame(
    donor_id = c("6436", "6442", "6539", "115", "6591"),
    state    = c("measured", "measured", "unavailable", "unavailable", "failed"),
    stringsAsFactors = FALSE)
  offered <- ie_selectable_donors(data)
  note <- ie_selectable_note(data, offered)

  expect_false(is.null(note))
  text <- gsub("<[^>]*>", "", as.character(note))
  expect_true(grepl("3 of 5 donors are not offered", text, fixed = TRUE))
  expect_true(grepl("2 unavailable", text, fixed = TRUE))
  expect_true(grepl("1 failed", text, fixed = TRUE))
  # Narrowing a selector must not read as data that never existed.
  expect_true(grepl("recorded, not missing", text, fixed = TRUE))

  # Nothing withheld, nothing said.
  all_measured <- data.frame(donor_id = c("a", "b"), state = c("measured", "measured"),
                             stringsAsFactors = FALSE)
  expect_null(ie_selectable_note(all_measured, ie_selectable_donors(all_measured)))
})

test_that("an empty or malformed table is refused rather than guessed at", {
  expect_identical(ie_selectable_donors(NULL), character())
  expect_identical(ie_selectable_donors(data.frame()), character())
  expect_error(ie_selectable_donors(data.frame(donor_id = "a", stringsAsFactors = FALSE)),
               "state", fixed = TRUE)
})
