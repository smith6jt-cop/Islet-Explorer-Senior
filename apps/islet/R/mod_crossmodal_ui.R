# Cross-Modal tab -- UI.
#
# A THIN REGISTRY, on purpose. Increments 2, 3 and 5 each add a card here, and
# they are built by different agents in parallel. If the cards lived in this file
# those agents would all edit it and collide. Instead each card ships as its own
# `card_*.R` declaring `card_<id>_ui(ns)` and `card_<id>_server(input, output,
# session, ...)`, registers itself in IE_CROSSMODAL_CARDS below, and this file
# stops changing.
#
# Adding a card is one line in the registry plus one new file. Nothing else.

#' Cards mounted on the Cross-Modal tab, in display order.
#'
#' `id` names the pair of functions `card_<id>_ui` / `card_<id>_server`.
#' `title` is what the card shows when it cannot render -- so the note has a
#' subject even before the builder exists.
IE_CROSSMODAL_CARDS <- list(
  list(id = "ledger",   title = "Donor × capability"),
  list(id = "x1_gradient", title = "Peri-islet immune gradient"),
  list(id = "x2_overlay",  title = "Post-Xenium IF overlay"),
  list(id = "x3_identity", title = "IF hormone identity")
)

#' Cross-Modal tab UI.
#'
#' @param id module id
crossmodal_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      class = "ie-card-stack",
      lapply(IE_CROSSMODAL_CARDS, function(card) {
        fn <- paste0("card_", card$id, "_ui")
        if (exists(fn, mode = "function")) {
          get(fn, mode = "function")(ns)
        } else {
          # The card has no builder yet. Say which one, and say it in the place
          # the card will occupy, rather than rendering an empty panel that
          # reads as "no data".
          ie_card_placeholder_ui(ns(card$id), card$title, fn)
        }
      })
    )
  )
}

#' A named note standing in for a card whose implementation does not exist yet.
#'
#' Not a blank panel: a blank panel is indistinguishable from a card that ran and
#' found nothing, which is the confusion the five availability states exist to
#' prevent.
ie_card_placeholder_ui <- function(id, title, expected_function) {
  div(
    class = "ie-card ie-card-placeholder",
    tags$h3(class = "ie-card-title", title),
    tags$p(
      class = "ie-card-note",
      sprintf("unavailable(not_built): no %s() is defined yet.", expected_function)
    )
  )
}
