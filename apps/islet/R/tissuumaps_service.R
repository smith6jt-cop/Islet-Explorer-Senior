# ---------------------------------------------------------------------------
# TissUUmaps tile server: make sure one is running, without starting a herd
#
# The Viewer tab's TissUUmaps backend iframes a server that has to be running
# somewhere. On this workstation the durable answer is the systemd --user unit
# `islet-tissuumaps-workstation` (deploy/islet-tissuumaps-workstation.service).
# This file is the FALLBACK for a development server, and for the case where the
# unit is not installed.
#
# WHY THE LOCKING IS NOT DECORATION. shiny-server runs one R worker PER SESSION.
# The obvious implementation -- "if the port is closed, spawn one" -- was measured
# with six concurrent workers against a cold port and spawned SIX servers, because
# each checked the port before any child had bound it. Taking a lock does not fix
# that on its own either: releasing the lock as soon as the child process is
# launched has the same race, since launching is not binding. The winner must hold
# the lock until the server actually ANSWERS.
#
# It must also never block app startup. Every path here is wrapped so that a
# failure leaves the Viewer tab showing its ordinary "no projects" message rather
# than taking the whole app down.
# ---------------------------------------------------------------------------

#' Is a TissUUmaps server answering at this base URL?
#'
#' `/ping` is the health endpoint. The `.tmap` routes are NOT usable as one: a
#' missing project returns a full ~31 KB HTML body, so a body-length or
#' content-type check would call a 404 healthy.
tissuumaps_alive <- function(base, timeout_s = 2) {
  url <- paste0(sub("/+$", "", base), "/ping")
  isTRUE(tryCatch({
    con <- url(url, open = "rb")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    length(readBin(con, "raw", 1L)) >= 0L
  }, error = function(e) FALSE, warning = function(w) FALSE))
}

#' Port out of a base URL, for the spawn path. NA when the base is a path such
#' as "/tissuumaps" -- in that case something in front of us owns the service and
#' we must not try to start one.
.tissuumaps_port <- function(base) {
  m <- regmatches(base, regexpr("(?<=:)[0-9]+(?=/|$)", base, perl = TRUE))
  if (!length(m)) NA_integer_ else as.integer(m)
}

#' Ensure a TissUUmaps server is up, and return TRUE if one is answering.
#'
#' @param base        base URL, default from resolve_tissuumaps_base()
#' @param project_dir slide root. MUST be the directory holding the .tmap files,
#'   because their layers reference `channels/<case>/...` and
#'   `regions/<case>.geojson` relatively against it.
#' @param wait_s      how long to wait for the spawned server to answer. The
#'   lock is held for this whole period on purpose -- see the header.
ensure_tissuumaps <- function(base = resolve_tissuumaps_base(),
                              project_dir = resolve_tissuumaps_project_dir(),
                              wait_s = 20) {
  tryCatch({
    # 1. Somebody is already serving -- the common case once the systemd unit is
    #    installed, and the only case we want in production.
    if (tissuumaps_alive(base)) return(TRUE)

    # 2. A path-shaped base means a reverse proxy owns the route. Starting a
    #    server here would bind a port nothing forwards to.
    port <- .tissuumaps_port(base)
    if (is.na(port)) return(FALSE)

    bin <- file.path(normalizePath("../..", mustWork = FALSE), ".venv", "bin", "tissuumaps_server")
    if (!file.exists(bin) || is.null(project_dir) || !dir.exists(project_dir)) return(FALSE)
    cfg <- file.path(normalizePath("../..", mustWork = FALSE), "deploy", "tissuumaps.cfg")

    # THE LOCK IS A DIRECTORY, IN A SHARED PLACE. Two things this got wrong first
    # time round, both caught by running the race rather than reasoning about it:
    #   * `file(path, open = "xb")` is NOT a valid R connection mode -- R warns
    #     "Invalid argument" and opens nothing, so every worker believed it had
    #     LOST the lock and none of them ever spawned.
    #   * `tempdir()` is per-PROCESS. A lock created there is invisible to every
    #     other worker, so it can never exclude anything.
    # `dir.create()` is atomic and returns FALSE without creating when the
    # directory exists, which is the R idiom for exactly this.
    lock <- file.path(Sys.getenv("TMPDIR", "/tmp"),
                      sprintf("islet-tissuumaps-%d.lock", port))
    # A worker killed mid-spawn would otherwise leave the lock forever. Anything
    # older than twice the wait cannot still be starting a server.
    if (dir.exists(lock)) {
      age <- difftime(Sys.time(), file.info(lock)$mtime, units = "secs")
      if (is.finite(age) && age > 2 * wait_s) unlink(lock, recursive = TRUE)
    }
    won <- isTRUE(suppressWarnings(dir.create(lock, showWarnings = FALSE)))

    if (won) {
      on.exit(unlink(lock, recursive = TRUE), add = TRUE)
      args <- c("-l", "127.0.0.1", "-p", port, "-r")
      # -r and the cfg both matter: without READ_ONLY the server exposes a
      # WRITABLE file browser over the whole slide root.
      if (file.exists(cfg)) args <- c(args, "-c", cfg)
      args <- c(args, project_dir)
      system2(bin, args = shQuote(args), stdout = FALSE, stderr = FALSE, wait = FALSE)
    }

    # 3. Winner and losers alike wait here. The winner still holds the lock, so
    #    no second worker can reach the spawn above until the server answers or
    #    the timeout expires.
    deadline <- Sys.time() + wait_s
    while (Sys.time() < deadline) {
      if (tissuumaps_alive(base)) return(TRUE)
      Sys.sleep(0.4)
    }
    FALSE
  }, error = function(e) FALSE)
}
