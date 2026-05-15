# ---------- Spatial Tab module UI ----------
# Exports: spatial_ui(id)
# Layout: Controls sidebar (col-2) + Tissue Scatter (col-6) + Leiden Panel (col-4)
# + 3 neighborhood analysis card rows below (A: Infiltration, B: Enrichment, C: Proximity)

spatial_ui <- function(id) {
  ns <- NS(id)

  # Helper: numbered section heading (matches statistics tab pattern)
  section_heading <- function(step, title, subtitle) {
    div(style = "margin-bottom: 14px; margin-top: 22px; padding-bottom: 8px; border-bottom: 2px solid #d0e0f0;",
      div(style = "display: flex; align-items: baseline; gap: 10px;",
        span(step,
             style = paste0("display: inline-block; background: linear-gradient(135deg, #4477AA, #5599CC);",
                            " color: white; font-weight: 700; font-size: 14px; padding: 2px 10px;",
                            " border-radius: 12px; min-width: 28px; text-align: center;")),
        span(title, style = "font-weight: 700; font-size: 18px;")
      ),
      tags$small(subtitle, style = "color: #777; display: block; margin-top: 4px;")
    )
  }

  tagList(
    # ==== ROW 1 (existing): Controls + Scatter + Leiden ====
    fluidRow(style = "display: flex; flex-wrap: wrap;",
      # ---- Left sidebar: Controls ----
      column(2,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%;",
          h5("Spatial Tissue Map", style = "margin-top: 0;"),
          uiOutput(ns("donor_selector")),
          hr(style = "margin: 8px 0;"),
          selectInput(ns("color_by"), "Color by",
                      choices = c("Phenotype" = "phenotype",
                                  "Leiden Cluster" = "leiden"),
                      selected = "phenotype"),
          uiOutput(ns("leiden_res_selector")),
          hr(style = "margin: 8px 0;"),
          radioButtons(ns("region_filter"), "Show Cells",
                       choices = c("All" = "all",
                                   "Core + Peri" = "core_peri",
                                   "Core Only" = "core"),
                       selected = "all", inline = FALSE),
          checkboxInput(ns("color_background"), "Color background cells", value = FALSE),
          hr(style = "margin: 8px 0;"),
          uiOutput(ns("phenotype_filter")),
          hr(style = "margin: 8px 0;"),
          checkboxGroupInput(ns("groups"), "Donor Status",
                             choices = c("ND", "Aab+", "T1D"),
                             selected = c("ND", "Aab+", "T1D"), inline = FALSE),
          hr(style = "margin: 8px 0;"),
          h5("Palettes", style = "font-size: 14px; margin-bottom: 8px;"),
          selectInput("spatial_palette", "Phenotype",
                      choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                      selected = "High Contrast"),
          selectInput("spatial_donor_palette", "Donor Status",
                      choices = names(DONOR_PALETTES),
                      selected = "Bright"),
          hr(style = "margin: 8px 0;"),
          downloadButton(ns("download_spatial"), "Download CSV", style = "width: 100%; font-size: 13px;")
        )
      ),

      # ---- Center: Tissue Scatter ----
      column(6,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; height: 100%;",
          div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 8px;",
            h5("Tissue Scatter Plot", style = "margin: 0;"),
            actionButton(ns("scatter_reset_zoom"), "Reset Zoom",
                         style = "font-size: 12px; padding: 2px 10px;",
                         class = "btn-sm btn-outline-secondary")
          ),
          tags$small("Drag to zoom, double-click to reset.", style = "color: #888; display: block; margin-bottom: 5px;"),
          plotOutput(ns("tissue_scatter"), height = "700px",
                     dblclick = ns("scatter_dblclick"),
                     brush = brushOpts(id = ns("scatter_brush"), resetOnNew = TRUE))
        )
      ),

      # ---- Right: Leiden Panel ----
      column(4,
        div(class = "card", style = "padding: 15px; margin-bottom: 15px; overflow: visible; height: 100%;",
          h5("Islet-Level UMAPs"),
          uiOutput(ns("leiden_not_available")),
          fluidRow(
            column(6,
              h5("Leiden Clustering", style = "font-size: 13px; margin-top: 0;"),
              plotlyOutput(ns("leiden_umap"), height = "300px")
            ),
            column(6,
              h5("Donor Status", style = "font-size: 13px; margin-top: 0;"),
              plotOutput(ns("spatial_umap_donor"), height = 300)
            )
          ),
          hr(style = "margin: 8px 0;"),
          h5("Cluster Composition", style = "font-size: 15px;"),
          plotlyOutput(ns("cluster_composition"), height = "350px")
        )
      )
    ),

    # ==== ROW 2: Card A — Immune Infiltration Overview ====
    uiOutput(ns("neighborhood_cards"))
  )
}
