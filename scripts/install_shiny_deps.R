pkgs <- c(
  "shiny","readxl","dplyr","stringr","tidyr","ggplot2","plotly","broom","base64enc","cluster","sf",
  # Phase 16 statistics: mixed-effects model + post-hoc contrasts. Required by
  # utils_stats.R (lmerTest::lmer) and mod_statistics_server.R (emmeans::emmeans).
  # If missing in /usr/local/lib/R/site-library the shiny-server worker fails
  # to initialize with HTTP 500 since user `shiny` can't see personal user libs.
  "lmerTest","emmeans"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
message("Installed/verified: ", paste(pkgs, collapse=", "))

# Install Bioconductor packages (slingshot and friends)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("slingshot", "SingleCellExperiment", "DelayedMatrixStats")
try({
  BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
  message("Installed/verified (Bioconductor): ", paste(bioc_pkgs, collapse=", "))
}, silent = TRUE)

# destiny (DiffusionMap) and tradeSeq from Bioconductor; lmtest from CRAN
try({
  BiocManager::install(c("destiny","tradeSeq"), ask = FALSE, update = FALSE)
  message("Installed/verified (Bioconductor): destiny, tradeSeq")
}, silent = TRUE)

if (!requireNamespace("lmtest", quietly = TRUE)) {
  install.packages("lmtest", repos = "https://cloud.r-project.org")
  message("Installed/verified: lmtest")
}

# Install vitessceR (from CRAN if available; fallback to GitHub)
install_vitesscer <- function() {
  if (requireNamespace("vitessceR", quietly = TRUE)) {
    message("vitessceR already installed")
    return(invisible(TRUE))
  }
  ok <- FALSE
  # Try CRAN first
  try({
    install.packages("vitessceR", repos = "https://cloud.r-project.org")
    ok <- requireNamespace("vitessceR", quietly = TRUE)
  }, silent = TRUE)
  if (!ok) {
    # Fallback to GitHub
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    try({
      remotes::install_github("vitessce/vitessceR")
      ok <- requireNamespace("vitessceR", quietly = TRUE)
    }, silent = TRUE)
  }
  if (ok) message("Installed/verified: vitessceR") else warning("vitessceR install did not complete; see logs")
  invisible(ok)
}

# Only install vitessceR if requested explicitly (not required for the embedded viewer)
if (identical(tolower(Sys.getenv("INSTALL_VITESSCER", "0")), "1")) {
  install_vitesscer()
}
