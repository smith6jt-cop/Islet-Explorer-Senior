viewer_ui <- function(id) {
  ns <- NS(id)
  div(style = "width: 100%;",
    # Phenotype-rules pill moved to the rail's Help entry, which reaches it from
    # every tab rather than from the five that happened to mount a copy.
    # In-place iframe src updater. The server pushes the new Avivator URL here
    # instead of letting renderUI recreate the <iframe>. Recreating the iframe
    # (a WebGL canvas) on every reactive tick causes heavy flicker under
    # Firefox WebRender on macOS (Mozilla bug #1555544). We update the src of
    # the existing node and dedup against the last value via a data-relsrc
    # attribute so spurious messages don't trigger a reload.
    tags$script(HTML("
Shiny.addCustomMessageHandler('islet_update_viewer_src', function(msg){
  if (!msg || !msg.id) return;
  var f = document.getElementById(msg.id);
  if (!f) return;
  if (f.getAttribute('data-relsrc') === msg.src) return;
  f.setAttribute('data-relsrc', msg.src);
  f.src = msg.src;
});
")),
    uiOutput(ns("image_selector")),
    uiOutput(ns("viewer_frame"))
  )
}
