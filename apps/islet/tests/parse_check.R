#!/usr/bin/env Rscript
# Syntax-parse every R source file in the Shiny app.
#
# Catches syntax errors without needing the app's data files or heavy runtime
# deps (sf, reticulate, anndata) — those are only required to actually launch the
# app, and the data artifacts are gitignored (not present in CI or a fresh clone).
# parse() checks syntax only; it does not source/evaluate, so no library or data
# is loaded.
#
#   Rscript apps/islet/tests/parse_check.R    # from the repo root
#   Rscript tests/parse_check.R                  # from apps/islet

root <- if (dir.exists("apps/islet/R")) {
  "apps/islet"
} else if (dir.exists("R")) {
  "."
} else {
  stop("Cannot locate the Shiny app R/ directory from wd = ", getwd())
}

iecore <- if (dir.exists("apps/R/ieCore/R")) {
  "apps/R/ieCore/R"
} else {
  file.path("..", "R", "ieCore", "R")
}
if (!dir.exists(iecore)) stop("Cannot locate apps/R/ieCore/R from wd = ", getwd())

files <- c(
  list.files(iecore, pattern = "[.][Rr]$", full.names = TRUE),
  list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE),
  file.path(root, "app.R")
)
files <- files[file.exists(files)]
if (length(files) == 0) stop("No R files found to parse under ", root)

ok <- TRUE
for (f in files) {
  res <- tryCatch({ parse(f); TRUE },
                  error = function(e) {
                    message("PARSE ERROR in ", f, ": ", conditionMessage(e))
                    FALSE
                  })
  if (isTRUE(res)) message("ok: ", f)
  ok <- ok && isTRUE(res)
}

if (!ok) quit(status = 1L)
message(sprintf("All %d R source files parsed cleanly.", length(files)))
