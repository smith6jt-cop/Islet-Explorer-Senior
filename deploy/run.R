# deploy/run.R — launcher for UF RC PubApps (invoked by the systemd ExecStart).
#
# RC's model runs the app directly (no shiny-server): a single R process listens on
# the RC-assigned port and RC's load balancer reverse-proxies https://<name>.rc.ufl.edu
# to it. The systemd unit sets WorkingDirectory to apps/islet/, so getwd() there makes
# the app's relative data paths (file.path("..","..","data","app_data",...)) resolve.
library(shiny)
options(shiny.sanitize.errors = TRUE)
runApp(
  appDir         = getwd(),
  launch.browser = FALSE,
  host           = "0.0.0.0",
  port           = as.integer(Sys.getenv("ISLET_APP_PORT", "8080"))
)
