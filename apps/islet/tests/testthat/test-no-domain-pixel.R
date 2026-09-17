# PIXEL_SIZE_UM is a DOMAIN presentation constant (a rounded display value -- see
# apps/R/ieCore/R/ie_frames.R's note on ie_um_per_px()) and the spec forbids it in new files:
# a frame's scale is a property of the image a producer recorded, never a constant borrowed
# from whichever cohort happens to be loaded. This test enforces that ban for the two files
# this worker wrote.
#
# THIS TEST IS ALLOWLIST-SCOPED ON PURPOSE. There are 10 live PIXEL_SIZE_UM uses in the repo
# today (measured on 2026-09-16; `segmentation_helpers.R` is an explicit pre-existing
# exemption per the brief and is not to be reused as precedent for a new file). A repo-wide
# version of this test fails on arrival for reasons that have nothing to do with the files
# this worker owns.

ALLOWLISTED_FILES <- c(
  "apps/R/ieCore/R/ie_layer_io.R",
  "apps/R/ieCore/R/ie_holds.R"
)

.pixel_size_um_pattern <- function() "\\bPIXEL_SIZE_UM\\b"

#' Lines in `path` (or literal R source given via `text=`) that reference PIXEL_SIZE_UM
#' outside of a `#`-comment.
.find_pixel_size_um <- function(path = NULL, text = NULL) {
  lines <- if (!is.null(text)) text else readLines(path, warn = FALSE)
  stripped <- sub("#.*$", "", lines)
  hit <- grepl(.pixel_size_um_pattern(), stripped, perl = TRUE)
  which(hit)
}

test_that("neither allowlisted file uses PIXEL_SIZE_UM", {
  for (rel in ALLOWLISTED_FILES) {
    path <- .find_repo_file(rel)
    hits <- .find_pixel_size_um(path)
    expect_length(hits, 0)
    if (length(hits) > 0) {
      cat(sprintf("%s: PIXEL_SIZE_UM at line(s) %s\n",
                  rel, paste(hits, collapse = ", ")))
    }
  }
})

test_that("the scanner actually catches a PIXEL_SIZE_UM use, demonstrating the other direction", {
  offending <- c(
    "px_to_um <- function(px) {",
    "  px * PIXEL_SIZE_UM",
    "}"
  )
  expect_true(length(.find_pixel_size_um(text = offending)) > 0)

  # A frame's own scale accessor (ie_um_per_px(), what this repo actually uses) must not
  # be flagged -- the ban is on the DOMAIN constant by name, not on the concept of a scale.
  clean <- c(
    "px_to_um <- function(px, section_row) {",
    "  px * ie_um_per_px(section_row)",
    "}"
  )
  expect_length(.find_pixel_size_um(text = clean), 0)

  # A comment merely mentioning the constant must not be flagged either.
  commented <- c("# never read PIXEL_SIZE_UM here")
  expect_length(.find_pixel_size_um(text = commented), 0)
})

test_that("re-measuring the repo-wide count shows why this test is scoped, not why it is weak", {
  ie_layer_io_path <- .find_repo_file("apps/R/ieCore/R/ie_layer_io.R")
  repo_root <- sub("/apps/R/ieCore/R/ie_layer_io\\.R$", "", ie_layer_io_path)
  r_files <- list.files(repo_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  r_files <- r_files[!grepl("/\\.venv/|/renv/library/", r_files)]
  total_hits <- sum(vapply(r_files, function(f) length(.find_pixel_size_um(f)), integer(1)))
  expect_gt(total_hits, 0)
})
