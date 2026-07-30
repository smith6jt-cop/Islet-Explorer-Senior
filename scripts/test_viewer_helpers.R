#!/usr/bin/env Rscript
# Unit tests for the TissUUmaps viewer helpers (base R only; no Shiny needed).
# Run from the repo root:  Rscript scripts/test_viewer_helpers.R

# Exercise the new TissUUmaps helpers without booting Shiny: parse
# viewer_helpers.R and eval only the function definitions we care about.
project_root <- NULL

exprs <- parse("app/shiny_app/R/viewer_helpers.R")
wanted <- c("resolve_viewer_backend", "resolve_tissuumaps_base",
            "resolve_tissuumaps_project_dir", "list_available_tmap_projects",
            "build_tissuumaps_iframe_url")
for (e in exprs) {
  if (is.call(e) && identical(as.character(e[[1]]), "<-") &&
      as.character(e[[2]]) %in% wanted) eval(e)
}

ok <- TRUE
check <- function(label, got, want) {
  pass <- identical(got, want)
  ok <<- ok && pass
  cat(sprintf("%-58s %s\n", label, if (pass) "PASS" else
    sprintf("FAIL (got %s, want %s)", deparse(got), deparse(want))))
}

# --- backend selection ------------------------------------------------------
Sys.unsetenv("ISLET_VIEWER_BACKEND")
check("backend default", resolve_viewer_backend(), "avivator")
Sys.setenv(ISLET_VIEWER_BACKEND = "tissuumaps")
check("backend=tissuumaps", resolve_viewer_backend(), "tissuumaps")
Sys.setenv(ISLET_VIEWER_BACKEND = "  TissUUmaps ")
check("backend case/space tolerant", resolve_viewer_backend(), "tissuumaps")
Sys.setenv(ISLET_VIEWER_BACKEND = "vitessce")
check("backend unknown -> avivator",
      suppressWarnings(resolve_viewer_backend()), "avivator")
Sys.setenv(ISLET_VIEWER_BACKEND = "")
check("backend empty -> avivator", resolve_viewer_backend(), "avivator")

# --- base URL ---------------------------------------------------------------
Sys.unsetenv("TISSUUMAPS_URL")
check("base default same-origin", resolve_tissuumaps_base(), "/tissuumaps")
Sys.setenv(TISSUUMAPS_URL = "http://127.0.0.1:5100/")
check("base strips trailing slash", resolve_tissuumaps_base(), "http://127.0.0.1:5100")

# --- iframe URL -------------------------------------------------------------
Sys.unsetenv("TISSUUMAPS_PROJECT_PATH")
check("iframe url", build_tissuumaps_iframe_url("6539"),
      "http://127.0.0.1:5100/6539.tmap")
check("iframe url strips .tmap", build_tissuumaps_iframe_url("6539.tmap"),
      "http://127.0.0.1:5100/6539.tmap")
check("iframe url NULL project", build_tissuumaps_iframe_url(NULL), NULL)
check("iframe url empty project", build_tissuumaps_iframe_url(""), NULL)
Sys.setenv(TISSUUMAPS_PROJECT_PATH = "senior/projects")
check("iframe url with ?path=", build_tissuumaps_iframe_url("6539"),
      "http://127.0.0.1:5100/6539.tmap?path=senior%2Fprojects")
Sys.unsetenv("TISSUUMAPS_PROJECT_PATH")

# --- project listing --------------------------------------------------------
tmp <- file.path(tempdir(), "tmaps")
dir.create(tmp, showWarnings = FALSE)
for (f in c("6539.tmap", "6450.tmap", "115.tmap", "notes.txt")) {
  file.create(file.path(tmp, f))
}
Sys.setenv(TISSUUMAPS_PROJECT_DIR = tmp)
check("lists .tmap stems, sorted, no stragglers",
      list_available_tmap_projects(), c("115", "6450", "6539"))
Sys.setenv(TISSUUMAPS_PROJECT_DIR = file.path(tmp, "does-not-exist"))
check("missing dir -> empty (no OME-TIFFs here)",
      list_available_tmap_projects(), character(0))
Sys.unsetenv("TISSUUMAPS_PROJECT_DIR")

cat(if (ok) "\nALL PASS\n" else "\nFAILURES\n")
quit(status = if (ok) 0 else 1)
