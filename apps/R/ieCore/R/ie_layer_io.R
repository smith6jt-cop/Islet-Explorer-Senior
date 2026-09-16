# Reading layer artifacts -- the small, pre-aggregated parquet a builder writes under
# `$ISLET_LAYERS_DIR`, so a card never scans the big source trees at render time.
#
# Nothing under `$ISLET_LAYERS_DIR` exists yet (see docs/ISLET_LAYERS_SPEC.md). Every card
# that reads a layer therefore hits the "directory not there" case on this tree today, and
# the repo rule for a missing measurement applies here too: a card that cannot find its
# artifact renders a NAMED NOTE, never a blank panel. This file is the one place that turns
# "directory absent" into that note, so no card re-derives the absence logic itself.

#' The root builders write layer artifacts under
#'
#' `$ISLET_LAYERS_DIR`, defaulting to `/data/islet_layers`. There is deliberately no other
#' default and no search path: a layer root is configuration, not something to discover by
#' walking the filesystem.
#' @export
ie_layers_dir <- function() {
  d <- Sys.getenv("ISLET_LAYERS_DIR", unset = "/data/islet_layers")
  if (is.na(d) || !nzchar(d)) d <- "/data/islet_layers"
  d
}

#' Path to one named layer directory under [ie_layers_dir()]
#' @export
ie_layer_path <- function(layer) {
  if (length(layer) != 1L || is.na(layer) || !nzchar(layer)) {
    stop("ie_layer_path(): layer must be one non-empty string", call. = FALSE)
  }
  file.path(ie_layers_dir(), layer)
}

#' Availability state for a layer directory
#'
#' `measured` only when the directory exists AND holds at least one file. An existing-but-empty
#' directory is not a measurement -- it is what `mkdir -p` leaves behind before a builder has
#' written anything -- so it reports the same `unavailable` state as a directory that was never
#' created, with a reason that says which.
#' @export
ie_layer_state <- function(layer) {
  path <- ie_layer_path(layer)
  if (!dir.exists(path)) {
    return(ie_avail("unavailable",
                     reason = paste0("layer directory not found: ", path),
                     source = path))
  }
  files <- list.files(path, recursive = TRUE, all.files = FALSE, no.. = TRUE)
  if (length(files) == 0L) {
    return(ie_avail("unavailable",
                     reason = paste0("layer directory is empty: ", path),
                     source = path))
  }
  ie_avail("measured", n = length(files), source = path)
}

#' Guard a card on a layer directory existing and holding files
#'
#' Returns `NULL` when the layer is present, otherwise a named [ie_guard_note] built by
#' [ie_require_measured()] -- never `NULL` for a missing layer. `guard` and `builder` are
#' passed straight through so the note names both the check that failed and the command that
#' would fix it.
#' @export
ie_require_layer <- function(layer, guard, builder = NULL) {
  state <- ie_layer_state(layer)
  ie_require_measured(state, guard = guard, artifact = ie_layer_path(layer), builder = builder)
}
