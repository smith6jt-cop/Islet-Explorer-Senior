# Guard against a bare `return(NULL)` creeping into the file this worker owns.
#
# THIS TEST IS ALLOWLIST-SCOPED ON PURPOSE. There are 161 bare `return(NULL)` uses across 15
# files in this repo today (measured on 2026-09-16 -- re-measure with the command in the
# second test below). A repo-wide version of this test fails on arrival and tells the reader
# nothing about the file this worker actually wrote. Scoping it to `ie_layer_io.R` is what
# makes "a bare return(NULL) in a new ieCore file is a defect" a claim this suite can actually
# enforce, rather than a claim it can only assert and immediately contradict.
#
# The allowlist held a second entry, the hold framework, until 2026-09-22. That file is
# CODE_TIER1 in scripts/retire/dispositions.py and the ledger owns its fate; a path literal
# here is a live reference by the ledger's own check, so the entry left before the file does.
#
# The rule itself: `ie_layer_io.R` returns an `ie_avail(...)` or an `ie_guard_note` (or NULL
# from ie_require_layer(), which is itself produced by ie_require_measured() -- a call, not a
# literal) for every "not found" case. A bare `return(NULL)` in it would be the exact defect
# the availability contract exists to prevent: a silent NULL standing in for an absence that
# should carry a state and a reason.

ALLOWLISTED_FILES <- c(
  "apps/R/ieCore/R/ie_layer_io.R"
)

#' Bare `return(NULL)` -- not `return(NULL_something)`, not `# return(NULL)` inside a string,
#' just the literal call, on its own, optionally inside other syntax on the same line.
.bare_return_null_pattern <- function() "\\breturn\\s*\\(\\s*NULL\\s*\\)"

#' Lines in `path` (or literal R source given via `text=`) that trip the bare-return(NULL)
#' pattern outside of a `#`-comment.
.find_bare_return_null <- function(path = NULL, text = NULL) {
  lines <- if (!is.null(text)) text else readLines(path, warn = FALSE)
  # Strip a trailing `# ...` comment before matching, so a comment that merely
  # DISCUSSES `return(NULL)` (as this very file does, above) is not a hit.
  stripped <- sub("#.*$", "", lines)
  hit <- grepl(.bare_return_null_pattern(), stripped, perl = TRUE)
  which(hit)
}

test_that("no allowlisted file contains a bare return(NULL)", {
  for (rel in ALLOWLISTED_FILES) {
    path <- .find_repo_file(rel)
    hits <- .find_bare_return_null(path)
    expect_length(hits, 0)
    if (length(hits) > 0) {
      cat(sprintf("%s: bare return(NULL) at line(s) %s\n",
                  rel, paste(hits, collapse = ", ")))
    }
  }
})

test_that("the scanner actually catches a bare return(NULL), demonstrating the other direction", {
  # Same detector as above, run against synthetic source that is never written to disk. This
  # is the "fails if you add the thing back" half of the acceptance check: it proves the
  # test would have caught the violation had it landed in an allowlisted file, without
  # requiring this suite to actually break the tree to prove it.
  offending <- c(
    "f <- function(x) {",
    "  if (is.null(x)) return(NULL)",
    "  x",
    "}"
  )
  expect_true(length(.find_bare_return_null(text = offending)) > 0)

  # A guard-note return (what ie_layer_io.R actually does) must NOT be flagged.
  clean <- c(
    "f <- function(x) {",
    "  if (is.null(x)) return(ie_require_measured(x, guard = 'g'))",
    "  x",
    "}"
  )
  expect_length(.find_bare_return_null(text = clean), 0)

  # A comment merely mentioning the pattern must not be flagged either.
  commented <- c("# never do this: return(NULL)")
  expect_length(.find_bare_return_null(text = commented), 0)
})

test_that("re-measuring the repo-wide count shows why this test is scoped, not why it is weak", {
  # Not an assertion on the count itself (it will keep moving as other workers land code) --
  # just proof that a repo-wide version of this rule would presently fail, which is the
  # justification the brief asks this test to carry visibly rather than only in prose.
  ie_layer_io_path <- .find_repo_file("apps/R/ieCore/R/ie_layer_io.R")
  repo_root <- sub("/apps/R/ieCore/R/ie_layer_io\\.R$", "", ie_layer_io_path)
  r_files <- list.files(repo_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  r_files <- r_files[!grepl("/\\.venv/|/renv/library/", r_files)]
  total_hits <- sum(vapply(r_files, function(f) length(.find_bare_return_null(f)), integer(1)))
  expect_gt(total_hits, 0)
})
