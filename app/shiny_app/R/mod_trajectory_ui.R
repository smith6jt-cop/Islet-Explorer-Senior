trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(style = "background-color: #e8f4fd; border: 1px solid #b8daff; border-radius: 5px; padding: 10px 15px; margin-bottom: 15px; color: #004085; font-size: 13px;",
      "Click any point on the pseudotime scatter plot to view the segmentation map for that islet."
    ),
    div(style = "display: flex; gap: 20px; align-items: center; margin-bottom: 8px; flex-wrap: wrap;",
      tags$div(
        title = GLOBAL_OUTLIER_TOOLTIP,
        checkboxInput("traj_remove_outliers",
                      "Exclude outliers (z > 3 per donor status)",
                      value = TRUE,
                      width = "320px")
      ),
      checkboxInput(ns("show_outlier_table"), "Show excluded-outlier table", value = FALSE)
    ),
    uiOutput(ns("traj_status")),
    fluidRow(
      # Left card: Scatterplot and heatmaps (narrowed from col-8)
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          fluidRow(
            column(3,
              uiOutput(ns("traj_feature_selector")),
              selectInput(ns("traj_show_trend"), "Trend lines:",
                         choices = c("None" = "none",
                                     "Overall" = "overall",
                                     "By Donor Status" = "by_donor"),
                         selected = "by_donor"),
              selectInput(ns("pt_mode"), "Pseudotime mode:",
                         choices = c("Core only" = "core",
                                     "Core + Peri" = "combined"),
                         selected = "core")
            ),
            column(3,
              selectInput(ns("traj_color_by"), "Color points by:",
                         choices = c("Donor Status" = "donor_status",
                                     "Donor ID" = "donor_id",
                                     "Leiden (res 0.3)" = "leiden_0.3",
                                     "Leiden (res 0.5)" = "leiden_0.5",
                                     "Leiden (res 0.8)" = "leiden_0.8",
                                     "Leiden (res 1.0)" = "leiden_1.0"),
                         selected = "donor_status"),
              selectInput(ns("traj_point_size"), "Point size by:",
                         choices = c("Cell Count" = "total_cells",
                                     "Islet Diameter" = "islet_diam_um",
                                     "Uniform" = "uniform"),
                         selected = "total_cells")
            ),
            column(3,
              sliderInput(ns("traj_alpha"), "Point transparency:",
                         min = 0.1, max = 1.0, value = 0.6, step = 0.05),
              sliderInput(ns("traj_point_size_slider"), "Point size:",
                         min = 0.5, max = 5.0, value = 3.0, step = 0.1),
              numericInput(ns("min_cells"), "Min cells/islet",
                           value = 1, min = 1, max = 200, step = 1)
            ),
            column(3,
              div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #f9f9f9; min-height: 100px;",
                h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-weight: bold; color: #333; font-size: 15px;"),
                uiOutput(ns("traj_legend"))
              ),
              selectInput("traj_donor_palette", "Donor Status Colors",
                          choices = names(DONOR_PALETTES),
                          selected = "Paul Tol (default)")
            )
          ),
          uiOutput(ns("traj_outlier_badge")),
          plotlyOutput(ns("traj_scatter"), height = 650),
          br(),
          plotlyOutput(ns("traj_heatmap"), height = 100),
          hr(),
          h5("Multi-Feature Heatmap", style = "color: #2c5aa0;"),
          fluidRow(
            column(8, uiOutput(ns("multi_heatmap_marker_selector"))),
            column(4, sliderInput(ns("multi_heatmap_nbins"), "Number of bins:", min = 10, max = 40, value = 20, step = 5))
          ),
          plotOutput(ns("traj_multi_heatmap"), height = "auto")
        )
      ),
      # Right column: UMAPs stacked + segmentation viewer below
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
          fluidRow(
            column(6,
              h5("UMAP: Donor Status", style = "font-size: 13px; margin-top: 0;"),
              plotOutput(ns("traj_umap_donor"), height = 500)
            ),
            column(6,
              h5("UMAP: Selected Feature", style = "font-size: 13px; margin-top: 0;"),
              plotOutput(ns("traj_umap_feature"), height = 500)
            )
          )
        ),
        uiOutput(ns("segmentation_viewer_panel"))
      )
    ),
    fluidRow(column(12, uiOutput(ns("traj_outlier_info"))))
  )
}
