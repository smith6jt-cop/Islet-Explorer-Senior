# ---------- Spatial Tab module UI ----------
# Exports: spatial_selection_ui(id), spatial_ui(id)
# Layout: Controls sidebar (col-2) + Tissue Scatter (col-6) + Leiden Panel (col-4)
# + 3 neighborhood analysis card rows below (A: Infiltration, B: Enrichment, C: Proximity)

#' Selection block for the shell's frozen control panel.
#'
#' These controls used to sit in a col-2 card at the left of the tab, so Spatial
#' was the one tab whose selectors were somewhere else. Everything keeps its
#' module namespace and its input id; only where it renders changes.
#'
#' Dropped on the way: the two palette selectInputs. The rail footer owns both
#' now and `ie_sync_selects()` keeps every writer in step, so a second copy here
#' would be a third writer of a value that already has two.
spatial_selection_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "ie-aside-head", "Spatial tissue map"),
    uiOutput(ns("donor_selector")),
    selectInput(ns("color_by"), "Color by",
                choices = c("Phenotype" = "phenotype",
                            "Leiden Cluster" = "leiden"),
                selected = "phenotype", width = "100%"),
    uiOutput(ns("leiden_res_selector")),
    radioButtons(ns("region_filter"), "Show cells",
                 choices = c("All" = "all",
                             "Core + Peri" = "core_peri",
                             "Core Only" = "core"),
                 selected = "all", inline = FALSE),
    checkboxInput(ns("color_background"), "Color background cells", value = FALSE),
    uiOutput(ns("phenotype_filter")),
    # Donor status is NOT here. It read `choices = c("ND", "Aab+", "T1D")` --
    # the other cohort's vocabulary printed as a control for whichever dataset
    # is loaded -- and it was one of three unsynced copies of the same filter.
    # The cohort panel in the top bar owns it for every tab now.
    downloadButton(ns("download_spatial"), "Download CSV", class = "btn-block")
  )
}

spatial_ui <- function(id) {
  ns <- NS(id)

  # Helper: numbered section heading (matches statistics tab pattern)
  # `section_heading()` was defined locally here (and in three other
  # files) with the same hardcoded blue gradient inline. Replaced by
  # ie_section_heading() in R/ui_modernist_bits.R.
  section_heading <- ie_section_heading

  tagList(
    # Phenotype-rules pill moved to the rail's Help entry, which reaches it from
    # every tab rather than from the five that happened to mount a copy.
    # ==== ROW 1 (existing): Controls + Scatter + Leiden ====
    fluidRow(style = "display: flex; flex-wrap: wrap;",
      # The col-2 control card that stood here is now spatial_selection_ui(),
      # mounted in the shell's frozen panel by app.R. Same module namespace and
      # the same input ids -- a relocation, not a rewrite.

      # ---- Center: Tissue Scatter ----
      # 7 + 5 rather than 6 + 4: with the control column moved to the frozen
      # panel the row was two-thirds full and a twelfth of it sat empty.
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 8px;",
            h5("Tissue Scatter Plot", style = "margin: 0;"),
            actionButton(ns("scatter_reset_zoom"), "Reset Zoom",
                         style = "font-size: 12px; padding: 2px 10px",
                         class = "btn-sm btn-outline-secondary")
          ),
          tags$small("Drag to zoom, double-click to reset.", style = "display: block; margin-bottom: 5px"),
          plotOutput(ns("tissue_scatter"), height = "700px",
                     dblclick = ns("scatter_dblclick"),
                     brush = brushOpts(id = ns("scatter_brush"), resetOnNew = TRUE)),
          div(style = "text-align: right; margin-top: 8px",
            downloadButton(ns("dl_tissue_scatter"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      ),

      # ---- Right: Leiden Panel ----
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%",
          h5("Islet-Level UMAPs"),
          uiOutput(ns("leiden_not_available")),
          fluidRow(
            column(6,
              h5("Leiden Clustering", style = "font-size: 13px; margin-top: 0"),
              plotlyOutput(ns("leiden_umap"), height = "300px")
            ),
            column(6,
              h5("Donor Status", style = "font-size: 13px; margin-top: 0"),
              plotOutput(ns("spatial_umap_donor"), height = 300),
              div(style = "text-align: right; margin-top: 6px;",
                downloadButton(ns("dl_spatial_umap_donor"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            )
          ),
          hr(style = "margin: 8px 0;"),
          h5("Cluster Composition", style = "font-size: 15px"),
          plotlyOutput(ns("cluster_composition"), height = "350px")
        )
      )
    ),

    # ==== ROW 2: Card A — Immune Infiltration Overview ====
    uiOutput(ns("neighborhood_cards"))
  )
}
