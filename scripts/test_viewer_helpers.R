#!/usr/bin/env Rscript
# Unit tests for the TissUUmaps viewer helpers (base R only; no Shiny needed).
# Run from the repo root:  Rscript scripts/test_viewer_helpers.R

# Load every helper function, including their dependencies, without running
# startup code or registering Shiny resource paths. Availability is base R too.
helpers <- new.env(parent = globalenv())
sys.source("apps/R/ieCore/R/ie_availability.R", envir = helpers)
for (e in parse("apps/islet/R/viewer_helpers.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
    eval(e, envir = helpers)
  }
}
# Keep this test independent of installed data readers and local cohort data,
# while exercising the real default argument and cohort-filtering path.
helpers$read_viewer_cohort <- function() {
  ids <- c("112", "115", "6450", "6539")
  list(donor_ids = ids, availability = helpers$ie_avail("measured", n = length(ids)))
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
check("backend default", helpers$resolve_viewer_backend(), "avivator")
Sys.setenv(ISLET_VIEWER_BACKEND = "tissuumaps")
check("backend=tissuumaps", helpers$resolve_viewer_backend(), "tissuumaps")
Sys.setenv(ISLET_VIEWER_BACKEND = "  TissUUmaps ")
check("backend case/space tolerant", helpers$resolve_viewer_backend(), "tissuumaps")
Sys.setenv(ISLET_VIEWER_BACKEND = "vitessce")
check("backend unknown -> avivator",
      suppressWarnings(helpers$resolve_viewer_backend()), "avivator")
Sys.setenv(ISLET_VIEWER_BACKEND = "")
check("backend empty -> avivator", helpers$resolve_viewer_backend(), "avivator")

# --- base URL ---------------------------------------------------------------
Sys.unsetenv("TISSUUMAPS_URL")
check("base default same-origin", helpers$resolve_tissuumaps_base(), "/tissuumaps")
Sys.setenv(TISSUUMAPS_URL = "http://127.0.0.1:5100/")
check("base strips trailing slash", helpers$resolve_tissuumaps_base(), "http://127.0.0.1:5100")

# --- iframe URL -------------------------------------------------------------
Sys.unsetenv("TISSUUMAPS_PROJECT_PATH")
check("iframe url", helpers$build_tissuumaps_iframe_url("6539"),
      "http://127.0.0.1:5100/6539.tmap")
check("iframe url strips .tmap", helpers$build_tissuumaps_iframe_url("6539.tmap"),
      "http://127.0.0.1:5100/6539.tmap")
check("iframe url NULL project", helpers$build_tissuumaps_iframe_url(NULL), NULL)
check("iframe url empty project", helpers$build_tissuumaps_iframe_url(""), NULL)
Sys.setenv(TISSUUMAPS_PROJECT_PATH = "senior/projects")
check("iframe url with ?path=", helpers$build_tissuumaps_iframe_url("6539"),
      "http://127.0.0.1:5100/6539.tmap?path=senior%2Fprojects")
Sys.unsetenv("TISSUUMAPS_PROJECT_PATH")

# --- project listing --------------------------------------------------------
tmp <- tempfile("viewer-helpers-")
app_dir <- file.path(tmp, "apps", "islet")
project_dir <- file.path(tmp, "projects")
dir.create(app_dir, recursive = TRUE)
dir.create(project_dir)
helpers$project_root <- tmp
old_wd <- getwd()
setwd(app_dir)
for (f in c("6539.tmap", "6450.tmap", "115.tmap", "0112.tmap", "7003.tmap", "notes.txt")) {
  file.create(file.path(project_dir, f))
}
Sys.setenv(TISSUUMAPS_PROJECT_DIR = project_dir)
check("lists sorted cohort .tmap stems, preserves padded ids",
      helpers$list_available_tmap_projects(), c("0112", "115", "6450", "6539"))
Sys.setenv(TISSUUMAPS_PROJECT_DIR = file.path(tmp, "does-not-exist"))
check("missing dir -> empty (no OME-TIFFs here)",
      helpers$list_available_tmap_projects(), character(0))
Sys.unsetenv("TISSUUMAPS_PROJECT_DIR")
setwd(old_wd)
unlink(tmp, recursive = TRUE)

cat(if (ok) "\nALL PASS\n" else "\nFAILURES\n")
quit(status = if (ok) 0 else 1)
