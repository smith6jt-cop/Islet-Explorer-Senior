trajectory_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Pill and hint strip moved to the rail's Help entries -- see IE_TAB_TIPS.
    # `traj_remove_outliers` and `traj_outlier_z` are gone: they were
    # non-namespaced MIRRORS of the Plot and Spatial copies, kept in step by
    # observers in app.R. The cohort panel in the top bar owns the rule for
    # every tab. The excluded-outlier TABLE is a display choice and stays.
    div(style = "display: flex; gap: 20px; align-items: center; margin-bottom: 8px; flex-wrap: wrap",
      checkboxInput(ns("show_outlier_table"), "Show excluded-outlier table", value = FALSE)
    ),
    uiOutput(ns("traj_status")),
    fluidRow(
      # Left card: Scatterplot and heatmaps (narrowed from col-8)
      column(7,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px",
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
                         min = 0, max = 0.95, value = 0.4, step = 0.05),
              sliderInput(ns("traj_point_size_slider"), "Point size:",
                         min = 0.5, max = 5.0, value = 3.0, step = 0.1),
              # Min cells/islet moved to the cohort panel. It was one of four
              # unsynced copies, and THIS one bound to `total_cells` -- core
              # only -- where the others used core+peri, so the same number
              # meant a different subset on this tab than on any other.
            ),
            column(3,
              div(class = "ie-legend-box",
                h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-weight: bold; font-size: 15px"),
                uiOutput(ns("traj_legend"))
              )
              # `traj_donor_palette` moved to the rail footer, which is the only
              # writer of either palette now. It was the third copy of one
              # selectInput, kept in step by observers that each had to name the
              # other two.
            )
          ),
          uiOutput(ns("traj_outlier_badge")),
          plotlyOutput(ns("traj_scatter"), height = 650),
          br(),
          plotlyOutput(ns("traj_heatmap"), height = 100),
          hr(),
          h5("Multi-Feature Heatmap"),
          fluidRow(
            column(8, uiOutput(ns("multi_heatmap_marker_selector"))),
            column(4, sliderInput(ns("multi_heatmap_nbins"), "Number of bins:", min = 10, max = 40, value = 20, step = 5))
          ),
          plotOutput(ns("traj_multi_heatmap"), height = "auto"),
          div(style = "text-align: right; margin-top: 8px",
            downloadButton(ns("dl_multi_heatmap"), "Download PNG",
                           class = "btn-sm btn-outline-primary"))
        )
      ),
      # Right column: UMAPs stacked + segmentation viewer below
      column(5,
        div(class = "card", style = "padding: 15px; margin-bottom: 20px",
          # STACKED, not side by side. These are coord_fixed() plots -- a UMAP has
          # no meaningful aspect other than 1:1 -- and this card sits in 5 of 12
          # columns, so two of them abreast got ~190px of width each against
          # 500px of height. Most of that width went to the legend and the data
          # collapsed into a dot in the corner with the axis titles overlapping.
          # One per row gets the full column, and the height comes down to match
          # so each is about square.
          fluidRow(
            column(12,
              h5("UMAP: Donor Status", style = "font-size: 13px; margin-top: 0"),
              plotOutput(ns("traj_umap_donor"), height = 320),
              div(style = "text-align: right; margin-top: 6px",
                downloadButton(ns("dl_umap_donor"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            )
          ),
          fluidRow(
            column(12,
              h5("UMAP: Selected Feature", style = "font-size: 13px; margin-top: 12px"),
              plotOutput(ns("traj_umap_feature"), height = 320),
              div(style = "text-align: right; margin-top: 6px",
                downloadButton(ns("dl_umap_feature"), "Download PNG",
                               class = "btn-sm btn-outline-primary"))
            )
          )
        ),
        uiOutput(ns("segmentation_viewer_panel"))
      )
    ),
    fluidRow(column(12, uiOutput(ns("traj_outlier_info"))))
  )
}
