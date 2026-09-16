# Lymph Node tab -- UI.
#
# A THIN REGISTRY, on purpose. Increments 2, 3 and 5 each add a card here, and
# they are built by different agents in parallel. If the cards lived in this file
# those agents would all edit it and collide. Instead each card ships as its own
# `card_*.R` declaring `card_<id>_ui(ns)` and `card_<id>_server(input, output,
# session, ...)`, registers itself in IE_LYMPHNODE_CARDS below, and this file
# stops changing.
#
# Adding a card is one line in the registry plus one new file. Nothing else.

#' Cards mounted on the Lymph Node tab, in display order.
#'
#' `id` names the pair of functions `card_<id>_ui` / `card_<id>_server`.
#' `title` is what the card shows when it cannot render -- so the note has a
#' subject even before the builder exists.
IE_LYMPHNODE_CARDS <- list(
  list(id = "n1_sections",     title = "Donor 112 · two sections, two frames"),
  list(id = "s3_compartments", title = "Immune cells: islet · exocrine · pLN")
)

#' Lymph Node tab UI.
#'
#' @param id module id
lymphnode_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(
      class = "ie-card-stack",
      lapply(IE_LYMPHNODE_CARDS, function(card) {
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

# ie_card_placeholder_ui() is defined once, in mod_crossmodal_ui.R.
