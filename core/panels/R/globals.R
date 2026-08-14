#' Project-wide paths and options.
#'
#' Read at app start; downstream modules pull paths from `app_paths`.

app_paths <- list(
  # Constant — the 10x Xenium Prime 5K Human Pan-Tissue gene list.
  reference_5k     = "data/reference_5k",
  # 10x's other pre-designed Xenium panels (hLung, hBrain, hIO, etc.);
  # populated by scripts/fetch_10x_panels.R.
  reference_panels = "data/reference_panels",
  # Tissue-agnostic subpanel biology definitions (no detection_pct etc).
  subpanels_shared = "data/subpanels_shared",
  # Per-tissue inputs (subpanels, audit, optional custom panel).
  tissues_root     = "data/tissues",
  cache            = "cache"
)

app_milestones <- list(
  overview       = "M2",
  panel_browser  = "M2",
  load_xenium    = "M3",
  panel_compare  = "M4",
  cluster        = "M5",
  clustree       = "M6",
  subcluster     = "M7",
  marker         = "M8",
  export         = "M9"
)

#' Placeholder card used by every module before its milestone lands.
placeholder_card <- function(title, milestone) {
  bslib::card(
    bslib::card_header(title),
    bslib::card_body(
      shiny::p(shiny::strong(title), " — coming in ", milestone, "."),
      shiny::p("This is an M1 skeleton; the module's controls and outputs ",
               "are added in the milestone above per CLAUDE.md §10.")
    )
  )
}
