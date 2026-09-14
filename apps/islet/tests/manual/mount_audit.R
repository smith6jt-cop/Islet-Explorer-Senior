#!/usr/bin/env Rscript
# Every `output$x <-` in a module must have a matching `ns("x")` SOMEWHERE it can
# actually be reached from -- its own UI file, OR another server-side renderUI
# (several cards here are built that way: spatial's `neighborhood_cards` mounts
# ~17 outputs from inside the server, and statistics builds Section 4 the same
# way). Searching only the UI file reports those as unmounted, which buries the
# real hits.
#
# Nothing else catches this class: parse() is happy, the testthat suite is happy,
# and the session smoke test is happy, because an unmounted output simply never
# renders. It is how the infiltration card's plot, table, note and download
# button vanished behind a section heading with nothing under it.
#
#   cd apps/islet && Rscript tests/manual/mount_audit.R
#
# ACCEPTED BASELINE (recorded 2026-09-14, before the modernist port): exactly TWO hits, both
# pre-existing and unrelated to any styling --
#   trajectory  traj_scatter_container    (mod_trajectory_server.R)
#   spatial     tissue_scatter_container  (mod_spatial_server.R)
# Both are the Phase-3 rdeck backend-switch renderUIs. The UI files still mount the legacy
# traj_scatter / tissue_scatter directly, so the WebGL path is unreachable. Whether rdeck should
# be reachable at all is a functionality question, not a styling one; until it is answered, TWO
# is the pass mark and any third hit is the port's fault.
# ai_assistant added: Senior keeps the chat card behind ISLET_ENABLE_AI where Junior
# deleted it, so its outputs must be audited too.
mods <- c("plot", "trajectory", "statistics", "spatial", "inflammation", "viewer",
          "ai_assistant")
all_src <- paste(unlist(lapply(list.files("R", "[.]R$", full.names = TRUE), readLines)),
                 collapse = "\n")
app_src <- if (file.exists("app.R")) paste(readLines("app.R"), collapse = "\n") else ""
bad <- 0
for (m in mods) {
  sv <- sprintf("R/mod_%s_server.R", m)
  if (!file.exists(sv)) next
  s <- readLines(sv)
  outs <- unique(sub(".*output\\$([A-Za-z0-9_.]+).*", "\\1",
                     grep("output\\$[A-Za-z0-9_.]+ *<-", s, value = TRUE)))
  mounted <- function(o) {
    grepl(sprintf('ns("%s")', o), all_src, fixed = TRUE) ||
    grepl(sprintf('"%s"', o), app_src, fixed = TRUE)
  }
  miss <- outs[!vapply(outs, mounted, TRUE)]
  if (length(miss)) {
    cat(sprintf("  %-14s NOT MOUNTED ANYWHERE: %s\n", m, paste(miss, collapse = ", ")))
    bad <- bad + length(miss)
  }
}
if (bad == 0) cat("  every server output is reachable\n") else quit(status = 1)
