# Lymph Node tab -- server.
#
# A THIN REGISTRY. See mod_lymphnode_ui.R for why. This file dispatches to each
# registered card's server and must not grow card logic of its own.

#' Lymph Node tab server.
#'
#' @param id module id
#' @param prepared reactive holding the prepared dataset
#' @param cohort the run cohort (see ie_cohort.R)
#' @param active_tab reactive returning the selected tab key
lymphnode_server <- function(id, prepared, cohort = NULL, active_tab = reactive("")) {
  moduleServer(id, function(input, output, session) {
    for (card in IE_LYMPHNODE_CARDS) {
      fn <- paste0("card_", card$id, "_server")
      if (exists(fn, mode = "function")) {
        # local() so `card` is captured per iteration rather than by reference:
        # without it every card server would see the last element of the list.
        local({
          this <- card
          get(paste0("card_", this$id, "_server"), mode = "function")(
            input, output, session,
            prepared = prepared, cohort = cohort, active_tab = active_tab
          )
        })
      }
    }
    invisible(NULL)
  })
}
