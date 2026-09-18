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

#' The nav_panel_hidden value this registry mounts into.
#'
#' Passed to every card server, so a card never hardcodes a tab key of its own. A card that
#' guessed its own key is how card_s1_size_composition shipped with `tab_key = "spatial"` --
#' which matches nothing in IE_TAB_KEYS, so `req(identical(active_tab(), tab_key))` was never
#' true and the card could not draw. Checked against IE_TAB_KEYS in test-card-registry.R.
IE_CROSSMODAL_TAB <- "CrossModal"

#' Cards mounted on the Cross-Modal tab, in display order.
#'
#' `id` names the pair of functions `card_<id>_ui` / `card_<id>_server`.
#' `title` is what the card shows when it cannot render -- so the note has a
#' subject even before the builder exists.
# `ledger` was registered here with no card_ledger.R and no card_ledger_server anywhere, so the
# first card on the tab rendered unavailable on every donor. A registry entry is a promise that a
# card exists; unregistered is honest, "unavailable" is not. Removed 2026-09-18 -- rebuild the
# entry when the card is written, not before.
#
# `s2_trajectory` moved to the Trajectory tab, where a trajectory belongs. It is a within-donor
# islet ordering, not a cross-modality comparison, and it was on this tab only because the
# registry that first mounted it was this one.
IE_CROSSMODAL_CARDS <- list(
  list(id = "x1_gradient", title = "Peri-islet immune gradient"),
  list(id = "x2_overlay",  title = "Post-Xenium IF overlay"),
  list(id = "x3_identity", title = "IF hormone identity"),
  # S1 and S2 existed as card_*.R files but were listed in NO registry, so they were dead
  # code: six cards on disk, four mounted. Found by worker I-islet, which filed NEEDS.md
  # rather than edit a registry it did not own.
  list(id = "s1_size_composition", title = "Islet size · composition"),
  # The four-cell grid -- pancreas and pLN, protein and RNA. The app reads one
  # of the four; this card reads the layer that presents all four at donor
  # level, and names the state of every cell it cannot show.
  list(id = "g1_grid",             title = "Both tissues, both modalities")
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
