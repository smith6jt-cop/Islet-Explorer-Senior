# deploy/install_runtime_R.R — OPTIONAL. Only run if the deployment must include the
# GitHub-only packages that environment.yml omits for v1. The app runs fine without them:
#   * vitessceR — optional embedded-viewer path; degrades gracefully if absent.
#   * rdeck     — WebGL scatter; only engages when ISLET_USE_RDECK=TRUE (default FALSE).
#
# Run inside the activated conda env:
#   /pubapps/<project>/prod/conda/bin/Rscript deploy/install_runtime_R.R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("vitessce/vitessceR", upgrade = "never")
remotes::install_github("qfes/rdeck",         upgrade = "never")
