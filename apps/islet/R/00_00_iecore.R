# Load the monorepo-local shared package sources before Shiny auto-sources the
# remaining application files. The deployed app already relies on the monorepo
# layout for ../../data/app_data, so this does not introduce a new layout rule.

.iecore_candidates <- unique(c(
  Sys.getenv("IECORE_R_PATH", unset = ""),
  file.path("..", "R", "ieCore", "R"),
  file.path("apps", "R", "ieCore", "R"),
  file.path("..", "..", "apps", "R", "ieCore", "R")
))
.iecore_candidates <- .iecore_candidates[nzchar(.iecore_candidates)]
.iecore_hits <- .iecore_candidates[dir.exists(.iecore_candidates)]
if (!length(.iecore_hits)) {
  stop("Cannot locate apps/R/ieCore/R; set IECORE_R_PATH", call. = FALSE)
}
.iecore_r_dir <- normalizePath(.iecore_hits[[1L]], mustWork = TRUE)
for (.iecore_file in sort(list.files(.iecore_r_dir, pattern = "[.]R$", full.names = TRUE))) {
  sys.source(.iecore_file, envir = globalenv())
}

.domain_candidates <- unique(c(
  Sys.getenv("ISLET_DOMAIN_CONFIG", unset = ""),
  file.path("config", "islet_default.yml"),
  file.path("apps", "islet", "config", "islet_default.yml"),
  file.path("..", "..", "apps", "islet", "config", "islet_default.yml")
))
.domain_candidates <- .domain_candidates[nzchar(.domain_candidates)]
.domain_hits <- .domain_candidates[file.exists(.domain_candidates)]
if (!length(.domain_hits)) {
  stop("Cannot locate apps/islet/config/islet_default.yml; set ISLET_DOMAIN_CONFIG", call. = FALSE)
}
DOMAIN <- ie_load_domain(.domain_hits[[1L]])
ie_set_domain(DOMAIN)
message("[DOMAIN] Loaded ", domain_app_name(), " from ", domain_config_path())

rm(.iecore_candidates, .iecore_hits, .iecore_r_dir, .iecore_file,
   .domain_candidates, .domain_hits)
