viewer_ui <- function(id) {
  ns <- NS(id)
  div(style = "width: 100%;",
    phenotype_rules_button("show_phenotype_rules_viewer"),
    uiOutput(ns("local_image_picker"))
  )
}
