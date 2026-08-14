#' Panel Browser module — selectizeInput across the subpanels + optional
#' custom panel + the joined 5K audit reference.
#'
#' Layout: sidebar with a selectize for the primary panel, a second
#' selectize for an optional comparison panel, a download button, and a
#' filter for `exclude_recommended`. Main panel: a gene `DT`, a `plotly`
#' detection-percentage scatter (rendered when the active tissue has at
#' least two `reference_runs`), and a 2-set comparison summary.

panel_browser_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shinyWidgets_or_shiny_selectize(ns("panel"),  "Primary panel",
                                      choices = NULL, multiple = FALSE),
      shinyWidgets_or_shiny_selectize(ns("compare"), "Compare with (optional)",
                                      choices = NULL, multiple = FALSE),
      shiny::checkboxInput(ns("hide_excluded"),
                           "Hide genes with exclude_recommended = yes",
                           value = FALSE),
      shiny::checkboxInput(ns("show_gene_labels"),
                           "Label genes on scatter (line + name)",
                           value = FALSE),
      shiny::p(shiny::strong("Tip: "),
               "select two panels to see set-overlap stats below the table."),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Custom panel source",
          shiny::uiOutput(ns("custom_upload_intro")),
          shiny::fileInput(ns("custom_upload"),
                           "Replace custom panel CSV",
                           accept = c(".csv", "text/csv")),
          shiny::actionButton(ns("custom_reset"),
                              "Reset to default",
                              class = "btn-sm btn-outline-secondary"),
          shiny::uiOutput(ns("custom_status"))
        )
      )
    ),
    bslib::card(
      bslib::card_header(shiny::textOutput(ns("primary_caption"),
                                           inline = TRUE)),
      DT::DTOutput(ns("genes"))
    ),
    bslib::card(
      bslib::card_header(shiny::textOutput(ns("scatter_caption"),
                                           inline = TRUE)),
      shiny::uiOutput(ns("scatter_help")),
      plotly::plotlyOutput(ns("scatter"), height = "420px")
    ),
    bslib::card(
      bslib::card_header("Overlap with comparison panel"),
      shiny::uiOutput(ns("compare_summary"))
    )
  )
}

# Use shiny's built-in selectize; we wrap via this thin shim so the call
# site can be swapped to shinyWidgets later if M10 polish wants it.
shinyWidgets_or_shiny_selectize <- function(id, label, choices, multiple) {
  shiny::selectizeInput(id, label, choices = choices, multiple = multiple,
                        options = list(allowEmptyOption = TRUE,
                                       placeholder = "(none)"))
}

panel_browser_server <- function(id, panels, panels_default, app_state) {
  shiny::moduleServer(id, function(input, output, session) {

    custom_slot <- custom_panel_slot_key()

    # Handle a user upload of a replacement custom-panel CSV. Enrichment
    # always uses the on-disk defaults (panels_default()), not the
    # currently-active override, so re-uploads are idempotent.
    shiny::observeEvent(input$custom_upload, {
      f <- input$custom_upload
      if (is.null(f)) return()
      tryCatch({
        result <- upload_custom_panel(
          csv_path    = f$datapath,
          panels      = panels_default(),
          source_name = f$name
        )
        app_state$custom_panel_override <- result$df
        app_state$custom_panel_status   <- result
      }, error = function(e) {
        app_state$custom_panel_override <- NULL
        app_state$custom_panel_status   <- NULL
        shinyWidgets::sendSweetAlert(
          session = session,
          title   = "Custom panel upload failed",
          text    = conditionMessage(e),
          type    = "error",
          btn_labels = "OK"
        )
      })
    }, ignoreNULL = TRUE)

    shiny::observeEvent(input$custom_reset, {
      app_state$custom_panel_override <- NULL
      app_state$custom_panel_status   <- NULL
      shiny::updateTextInput(session, "custom_upload", value = "")
    }, ignoreInit = TRUE)

    # Tissue-aware upload intro
    output$custom_upload_intro <- shiny::renderUI({
      p <- panels()
      default_lbl <- p$tissue$manifest$custom_panel$display_name
      if (is.null(default_lbl)) {
        shiny::p(shiny::em(
          "Upload a CSV to attach a custom panel for this session. ",
          "A `gene` column is required. The active tissue (",
          shiny::strong(p$tissue$display_name),
          ") has no default custom panel; uploaded genes are enriched ",
          "from the tissue's 5K audit table."))
      } else {
        shiny::p(shiny::em(
          "Upload a CSV to replace the default custom panel (",
          shiny::strong(default_lbl), ") for this session. ",
          "A `gene` column is required; missing canonical columns are ",
          "filled in from the default panel and the tissue's 5K audit ",
          "for any matching genes."))
      }
    })

    output$custom_status <- shiny::renderUI({
      st <- app_state$custom_panel_status
      p  <- panels()
      default_lbl <- p$tissue$manifest$custom_panel$display_name
      if (is.null(st)) {
        cust <- p$custom
        if (is.null(cust) || nrow(cust) == 0L) {
          return(shiny::tags$p(shiny::tags$small(
            shiny::strong("Custom panel: "),
            "no default for this tissue.")))
        }
        return(shiny::tags$p(shiny::tags$small(
          shiny::strong("Custom panel: "), "default ",
          if (!is.null(default_lbl)) sprintf("(%s, %d genes)",
                                             default_lbl, nrow(cust))
          else sprintf("(%d genes)", nrow(cust)))))
      }
      shiny::tags$div(
        shiny::tags$p(shiny::tags$small(
          shiny::strong("Custom panel: "), "uploaded ",
          shiny::tags$code(st$source_name),
          sprintf(" (%d genes)", st$n_genes))),
        shiny::tags$ul(
          class = "small mb-0",
          shiny::tags$li(sprintf("%d enriched from default panel",
                                 st$n_enriched_from_custom)),
          shiny::tags$li(sprintf("%d enriched from 5K audit",
                                 st$n_enriched_from_5k)),
          shiny::tags$li(sprintf("%d unmatched (no reference data)",
                                 st$n_unmatched))
        )
      )
    })

    # Prefix used to namespace the 10x pre-designed panel keys so they
    # don't collide with tissue subpanel filenames if someone names a
    # subpanel "hLung".
    ref_panel_key <- function(panel_id) paste0("ref10x:", panel_id)

    # Flat vector of every selectable panel value (subpanels + custom +
    # tissue's joined 5K audit + 10x pre-designed panels).
    panel_choices <- shiny::reactive({
      p <- panels()
      out <- sort(p$meta$subpanel_keys)
      has_custom <- (!is.null(p$custom) && nrow(p$custom) > 0L) ||
                    !is.null(app_state$custom_panel_status)
      if (has_custom) out <- c(out, custom_slot)
      out <- c(out, "xenium5k_in_audit")
      mf <- p$reference_panels_manifest
      if (!is.null(mf) && nrow(mf)) {
        out <- c(out, ref_panel_key(mf$panel_id))
      }
      out
    })

    # Optgroup-structured choices for the selectize dropdown:
    #   "Tissue subpanels"        -> sorted subpanel keys
    #   "Custom panel"            -> the (renamed) custom slot
    #   "Tissue 5K audit"         -> xenium5k_in_audit
    #   "10x pre-designed panels" -> hLung, hBrain, ... labeled by display_name
    panel_choices_grouped <- shiny::reactive({
      p <- panels()
      st <- app_state$custom_panel_status
      groups <- list()
      subp_keys <- sort(p$meta$subpanel_keys)
      if (length(subp_keys)) {
        groups[["Tissue subpanels"]] <- stats::setNames(subp_keys, subp_keys)
      }
      has_custom <- (!is.null(p$custom) && nrow(p$custom) > 0L) ||
                    !is.null(st)
      if (has_custom) {
        groups[["Custom panel"]] <- stats::setNames(custom_slot,
                                                    custom_panel_label(st, p))
      }
      groups[["Tissue 5K audit"]] <- c("xenium5k_in_audit" = "xenium5k_in_audit")
      mf <- p$reference_panels_manifest
      if (!is.null(mf) && nrow(mf)) {
        vals   <- ref_panel_key(mf$panel_id)
        species_tag <- ifelse(is.na(mf$species) | !nzchar(mf$species),
                              "", paste0(mf$species, ", "))
        labels <- sprintf("%s (%s%d genes)",
                          mf$display_name, species_tag, mf$n_genes)
        groups[["Reference panels"]] <- stats::setNames(vals, labels)
      }
      groups
    })

    shiny::observe({
      grouped <- panel_choices_grouped()
      flat    <- panel_choices()
      cur_panel   <- shiny::isolate(input$panel)
      cur_compare <- shiny::isolate(input$compare)
      sel_panel <- if (!is.null(cur_panel) && cur_panel %in% flat) cur_panel
                   else flat[1]
      sel_compare <- if (!is.null(cur_compare) &&
                         (cur_compare %in% flat || cur_compare == "(none)"))
                       cur_compare else "(none)"
      shiny::updateSelectizeInput(session, "panel",
                                  choices = grouped,
                                  server  = TRUE,
                                  selected = sel_panel)
      shiny::updateSelectizeInput(session, "compare",
                                  choices = c(list(`(none)` = c("(none)" = "(none)")),
                                              grouped),
                                  server  = TRUE,
                                  selected = sel_compare)
    })

    # When a user upload lands, jump to the custom slot.
    shiny::observeEvent(app_state$custom_panel_status, {
      shiny::updateSelectizeInput(session, "panel", selected = custom_slot)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # External nav from Overview: a subpanel key arrives via app_state
    shiny::observeEvent(app_state$selected_subpanel, {
      sel <- app_state$selected_subpanel
      if (is.null(sel) || !nzchar(sel)) return()
      if (sel %in% panel_choices()) {
        shiny::updateSelectizeInput(session, "panel", selected = sel)
      }
    }, ignoreNULL = TRUE)

    get_panel_df <- function(key) {
      p <- panels()
      if (is.null(key) || !nzchar(key) || identical(key, "(none)")) return(NULL)
      if (identical(key, custom_slot))     return(p$custom)
      if (identical(key, "xenium5k_in_audit")) return(p$xenium5k)
      if (startsWith(key, "ref10x:")) {
        pid <- sub("^ref10x:", "", key)
        return(p$reference_panels[[pid]])
      }
      p$subpanels[[key]]
    }

    primary_df <- shiny::reactive({
      df <- get_panel_df(input$panel)
      if (is.null(df)) return(NULL)
      if (isTRUE(input$hide_excluded) &&
          "exclude_recommended" %in% names(df)) {
        df <- df[!(df$exclude_recommended %in% c("yes", "TRUE", "true")), ,
                 drop = FALSE]
      }
      df
    })

    compare_df <- shiny::reactive({
      get_panel_df(input$compare)
    })

    output$primary_caption <- shiny::renderText({
      df <- primary_df()
      if (is.null(df)) return("Select a panel.")
      sprintf("%s — %d genes",
              panel_display_label(input$panel,
                                  app_state$custom_panel_status,
                                  panels()),
              nrow(df))
    })

    output$genes <- DT::renderDT({
      df <- primary_df()
      shiny::req(df)
      DT::datatable(
        df,
        rownames = FALSE,
        filter   = "top",
        options  = list(pageLength = 25, scrollX = TRUE,
                        autoWidth = FALSE),
        class    = "stripe hover compact nowrap"
      )
    })

    # Tissue-aware scatter caption / help — names the actual reference runs.
    output$scatter_caption <- shiny::renderText({
      runs <- panels()$reference_runs
      if (length(runs) >= 2L) {
        sprintf("Detection across reference Xenium runs (%s vs %s)",
                runs[1], runs[2])
      } else if (length(runs) == 1L) {
        sprintf("Detection in reference Xenium run (%s)", runs[1])
      } else {
        "Detection across reference Xenium runs"
      }
    })

    output$scatter_help <- shiny::renderUI({
      runs <- panels()$reference_runs
      if (length(runs) >= 2L) {
        shiny::tagList(
          shiny::p("Each point is a gene in the primary panel. ",
                   "x = ",
                   shiny::code(sprintf("detection_pct_%s", runs[1])),
                   ", y = ",
                   shiny::code(sprintf("detection_pct_%s", runs[2])),
                   ", colour = ",
                   shiny::code("log2_detection_ratio"),
                   ". Diagonal is y = x."),
          shiny::p(shiny::em("Hover over a point to see the gene name ",
                             "and per-run detection %. Toggle ",
                             shiny::strong("Label genes on scatter"),
                             " in the sidebar to print gene names next ",
                             "to each point."))
        )
      } else if (length(runs) == 1L) {
        shiny::p("Histogram of ",
                 shiny::code(sprintf("detection_pct_%s", runs[1])),
                 " across genes in the primary panel.")
      } else {
        shiny::p(shiny::em(
          "The active tissue declares no reference_runs in its manifest; ",
          "detection columns are not available."))
      }
    })

    output$scatter <- plotly::renderPlotly({
      df <- primary_df()
      shiny::req(df)
      runs <- panels()$reference_runs

      empty_plot <- function(msg) {
        plotly::plotly_empty(type = "scatter", mode = "markers") |>
          plotly::layout(annotations = list(list(
            text = msg, showarrow = FALSE,
            x = 0.5, y = 0.5, xref = "paper", yref = "paper")))
      }

      if (length(runs) == 0L) {
        return(empty_plot("No reference_runs configured for this tissue."))
      }

      run_x_col <- sprintf("detection_pct_%s", runs[1])
      run_y_col <- if (length(runs) >= 2L)
        sprintf("detection_pct_%s", runs[2]) else NA_character_

      if (!run_x_col %in% names(df) ||
          (!is.na(run_y_col) && !run_y_col %in% names(df))) {
        return(empty_plot("Detection columns not present for this panel."))
      }

      # Single-run fallback: histogram of the one column we have.
      if (is.na(run_y_col)) {
        return(plotly::plot_ly(
          x = df[[run_x_col]], type = "histogram", nbinsx = 40
        ) |> plotly::layout(
          xaxis = list(title = run_x_col),
          yaxis = list(title = "genes")))
      }

      # Two-run scatter (log10-spaced ticks; original % preserved on hover).
      df$lx <- log10(1 + df[[run_x_col]])
      df$ly <- log10(1 + df[[run_y_col]])
      tick_pct  <- c(0, 1, 3, 10, 30, 100)
      tick_vals <- log10(1 + tick_pct)
      tick_text <- as.character(tick_pct)
      ax_max    <- log10(101)

      ratio_col <- "log2_detection_ratio"
      has_ratio <- ratio_col %in% names(df) &&
                   any(is.finite(df[[ratio_col]]))
      marker <- if (has_ratio) {
        r <- df[[ratio_col]]
        q <- stats::quantile(abs(r), 0.95, na.rm = TRUE)
        cmax_q <- min(5, max(2, q))
        list(size = 6, color = r,
             cmin = -cmax_q, cmax = cmax_q, cmid = 0, cauto = FALSE,
             colorscale = "RdBu", reversescale = TRUE,
             colorbar = list(title = sprintf("log2 %s/%s",
                                             runs[2], runs[1])))
      } else {
        list(size = 6, color = "#1f77b4")
      }

      hover_text <- sprintf("%s<br>%s: %.2f%%<br>%s: %.2f%%",
                            df$gene, runs[1], df[[run_x_col]],
                            runs[2], df[[run_y_col]])

      fig <- plotly::plot_ly(
        df,
        x = ~lx, y = ~ly,
        type = "scatter", mode = "markers",
        marker = marker,
        text   = hover_text,
        hoverinfo = "text"
      ) |>
        plotly::layout(
          xaxis = list(title    = run_x_col,
                       tickvals = tick_vals, ticktext = tick_text,
                       range    = c(0, ax_max), zeroline = FALSE),
          yaxis = list(title    = run_y_col,
                       tickvals = tick_vals, ticktext = tick_text,
                       range    = c(0, ax_max), zeroline = FALSE,
                       scaleanchor = "x", scaleratio = 1),
          shapes = list(list(type = "line",
                             x0 = 0, y0 = 0, x1 = ax_max, y1 = ax_max,
                             line = list(dash = "dot", color = "gray")))
        )

      if (isTRUE(input$show_gene_labels) && nrow(df) > 0L) {
        # Plotly slows to a crawl past a few hundred annotations; cap at
        # 300 most-detectable genes (max of the two runs) so toggling the
        # label on for the full 5K panel doesn't freeze the browser.
        label_cap <- 300L
        sub <- df
        if (nrow(sub) > label_cap) {
          ord <- order(pmax(sub[[run_x_col]], sub[[run_y_col]]),
                       decreasing = TRUE)
          sub <- sub[ord[seq_len(label_cap)], , drop = FALSE]
        }
        # Fan labels out across four horizontal tracks: assigning by
        # y-rank so vertically-adjacent labels always land in different
        # tracks instead of stacking on top of each other. A tiny
        # vertical jitter further separates labels that fall into the
        # same track at nearby y.
        n_sub <- nrow(sub)
        y_rank <- integer(n_sub)
        y_rank[order(sub$ly)] <- seq_len(n_sub)
        track  <- (y_rank - 1L) %% 4L
        xoff_tracks <- c(160, 70, -70, -160)
        yoff_tracks <- c(-8,  4, -4,  8)
        xoff <- xoff_tracks[track + 1L]
        yoff <- yoff_tracks[track + 1L]
        fig <- fig |> plotly::add_annotations(
          x          = sub$lx,
          y          = sub$ly,
          text       = sub$gene,
          showarrow  = TRUE,
          arrowhead  = 1,
          arrowsize  = 0.5,
          arrowwidth = 0.5,
          arrowcolor = "rgba(80,80,80,0.55)",
          ax         = xoff,
          ay         = yoff,
          font       = list(size = 9, color = "#333"),
          standoff   = 3,
          bgcolor    = "rgba(255,255,255,0.75)",
          borderpad  = 1
        )
      }

      fig
    })

    output$compare_summary <- shiny::renderUI({
      a <- primary_df()
      b <- compare_df()
      if (is.null(a) || is.null(b)) {
        return(shiny::p(shiny::em("Select a comparison panel in the sidebar.")))
      }
      ga <- a$gene; gb <- b$gene
      both <- intersect(ga, gb)
      a_only <- setdiff(ga, gb); b_only <- setdiff(gb, ga)
      st <- app_state$custom_panel_status
      p  <- panels()
      lbl_a <- panel_display_label(input$panel,   st, p)
      lbl_b <- panel_display_label(input$compare, st, p)
      shiny::tagList(
        shiny::tags$table(
          class = "table table-sm",
          shiny::tags$thead(shiny::tags$tr(
            shiny::tags$th("Set"), shiny::tags$th("Genes"))),
          shiny::tags$tbody(
            shiny::tags$tr(shiny::tags$td(lbl_a, " only"),
                           shiny::tags$td(length(a_only))),
            shiny::tags$tr(shiny::tags$td("Both"),
                           shiny::tags$td(length(both))),
            shiny::tags$tr(shiny::tags$td(lbl_b, " only"),
                           shiny::tags$td(length(b_only))),
            shiny::tags$tr(shiny::tags$td(shiny::strong("Union")),
                           shiny::tags$td(shiny::strong(length(union(ga, gb))))))
        ),
        shiny::p(shiny::strong("Shared genes: "),
                 if (length(both) == 0) "(none)"
                 else paste(head(sort(both), 60), collapse = ", "),
                 if (length(both) > 60) shiny::span(
                   shiny::em(sprintf(" ... and %d more", length(both) - 60))))
      )
    })
  })
}
