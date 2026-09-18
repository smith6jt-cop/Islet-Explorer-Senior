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

#' The nav_panel_hidden value this registry mounts into. See mod_crossmodal_ui.R.
IE_LYMPHNODE_TAB <- "LymphNode"

#' Cards mounted on the Lymph Node tab, in display order.
#'
#' `id` names the pair of functions `card_<id>_ui` / `card_<id>_server`.
#' `title` is what the card shows when it cannot render -- so the note has a
#' subject even before the builder exists.
# `n1_sections` was once registered with no card_n1_sections.R and no server, so the first card on
# this tab rendered unavailable on every donor -- the same defect the Cross-Modal `ledger` entry
# had. It is registered again here because the card now exists: card_n1_sections.R declares both
# halves of the pair, as does card_n2_pln_immune.R. A registry entry is a promise that a card
# exists; unregistered is honest, "unavailable" over a card nobody wrote is not.
#
# Display order is the reading order: which annotation is the node, then what that node is made of
# in each modality separately, then the compartments the pancreas and the node are read in.
IE_LYMPHNODE_CARDS <- list(
  list(id = "n1_sections",     title = "Which annotation is the node"),
  list(id = "n2_pln_immune",   title = "pLN immune composition — each modality on its own"),
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
        ui <- if (exists(fn, mode = "function")) {
          get(fn, mode = "function")(ns)
        } else {
          # The card has no builder yet. Say which one, and say it in the place
          # the card will occupy, rather than rendering an empty panel that
          # reads as "no data".
          ie_card_placeholder_ui(ns(card$id), card$title, fn)
        }
        htmltools::tagAppendAttributes(ui, `data-card-id` = card$id)
      })
    )
  )
}

# ie_card_placeholder_ui() is defined once, in mod_crossmodal_ui.R.
