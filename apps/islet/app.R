# Islet Explorer - Modular App Entrypoint
# All utility functions and module code live in R/ (auto-sourced by Shiny)
# Package loading and constants live in global.R (auto-sourced by Shiny)

# ---------- UI ----------

ui <- bslib::page_fillable(
  padding = 0, gap = 0,

  # bs_theme(version = 5) UNMODIFIED, on purpose, and this is the one genuinely
  # behavioural change the design port makes. Two reasons it is the right one:
  #
  #   * The Modernist sheets were authored against Bootstrap 5. On Bootstrap 3
  #     a third of modernist-shiny.css addressed markup Shiny was not emitting
  #     -- most visibly `.nav-tabs .nav-link`, which is the BS5 form, so the tab
  #     bar was never styled by this app at all and showed Bootstrap's own blue.
  #   * bslib serves an UNMODIFIED BS5 theme from a precompiled bundle (~0.15 s,
  #     no Sass run). Any customisation -- primary=, base_font=, enable-rounded=
  #     -- flips it to a real compile, cold once per shiny-server worker. The
  #     whole design system here is CSS custom properties, so customising
  #     bs_theme would buy nothing and cost the app's startup time.
  #
  # Shiny emits identical classes under 3 and 5 for every input this app uses
  # (selectInput / radioButtons / checkboxGroupInput / actionButton / helpText),
  # and bslib's BS5 bundle ships a BS3 compatibility layer that still styles
  # .btn-default, .checkbox-inline, .well, .form-group, .help-block,
  # .control-label and .irs--shiny.
  theme = bslib::bs_theme(version = 5),

  useShinyjs(),
  tags$head(
    tags$title("Islet Explorer Senior"),
    # Favicons for browser tab, bookmarks, and mobile pinning. Browsers pick
    # the most appropriate sized PNG based on context; the ICO is the legacy
    # default location for tools that hard-fetch /favicon.ico. Apple touch
    # icon covers iOS Safari Add-to-Home-Screen.
    tags$link(rel = "icon", type = "image/png", sizes = "16x16",  href = "favicon-16.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "32x32",  href = "favicon-32.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "48x48",  href = "favicon-48.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "192x192", href = "favicon-192.png"),
    tags$link(rel = "icon", type = "image/png", sizes = "512x512", href = "favicon-512.png"),
    tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "apple-touch-icon.png"),
    tags$link(rel = "shortcut icon", type = "image/x-icon", href = "favicon.ico"),
    # The Modernist design system. Order matters: fonts, then styles.css (which
    # establishes the components), then modernist-shiny.css (which maps them onto
    # Shiny's own markup and rescues the native inputs styles.css hides).
    tags$link(rel = "preload", as = "font", type = "font/woff2",
              href = "fonts/Archivo-Variable.woff2", crossorigin = "anonymous"),
    tags$link(rel = "stylesheet", href = "fonts/archivo.css"),
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$link(rel = "stylesheet", href = "modernist-shiny.css"),
    tags$style(HTML("
    /* ---------------------------------------------------------------------
       What the Modernist stylesheets do NOT cover, and therefore stays here.
       Everything else that used to live in this block -- .card, buttons, form
       controls, wells, plotly figures, the sliders -- is now in
       www/styles.css and www/modernist-shiny.css, linked above.

       .ai-chat-*          LAYOUT ONLY. The chat card's skin moved to
                           modernist-shiny.css; what is left here is flex and sizing.
       .skip-link, focus   accessibility, added to Senior after the fork. The skip
                           link's own colours are in the sheet.

       REMOVED with the shell: `body.viewer-mode` / `body.trajectory-mode`, which
       were keyed to the Bootstrap-3 grid AND which nothing ever set -- grep found
       no writer of either class in any R or JS file; and `.equal-height-row`,
       whose row no longer exists.

       MOVED OUT 2026-09-14, to modernist-shiny.css's Senior section: .logo-header,
       every .sidebar rule, the tab bar, and the .ai-chat-* colours. They were the
       2019 blue skin and, being inline, they beat both linked sheets on source
       order -- so the app still rendered blue after the sheets were linked. The
       tab rules were also written as `.nav-tabs .nav-link`, the Bootstrap 5 form;
       Shiny renders Bootstrap 3, which is `.nav-tabs > li > a`, so that rule had
       never matched anything and Bootstrap's own blue showed through.

       Dropped deliberately: the global `body` and `.container-fluid` rules, so the
       sheet's --color-ground can take effect, and `accent-color` on native radios and
       checkboxes, which modernist-shiny.css:608-627 restyles properly.
       --------------------------------------------------------------------- */
.skip-link {
      position: absolute;
      left: 12px;
      top: -80px;
      z-index: 10000;
      padding: 10px 14px;
      font-weight: 700;
    }
.skip-link:focus { top: 12px; }
a:focus-visible, button:focus-visible, input:focus-visible,
    select:focus-visible, [tabindex]:focus-visible {
      outline: 3px solid #ffbf47 !important;
      outline-offset: 2px !important;
    }
.ai-chat-panel {
      padding: 15px;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
.ai-chat-logo {
      display: flex;
      justify-content: center;
      align-items: center;
      padding-bottom: 8px;
    }
.ai-chat-logo img {
      max-width: 100%;
      height: auto;
      max-height: 90px;
    }
.ai-chat-header {
      display: flex;
      justify-content: center;
      align-items: center;
      margin-top: 8px;
      margin-bottom: 4px;
    }
.ai-chat-history {
      flex: 1;
      overflow-y: auto;
      padding: 8px;
      display: flex;
      flex-direction: column;
      gap: 10px;
    }
.ai-chat-message {
      padding: 10px 12px;
      margin: 4px 0;
      max-width: 90%;
      word-wrap: break-word;
    }
.ai-chat-message.user {
      align-self: flex-end;
      margin-left: auto;
    }
.ai-chat-message.assistant {
      align-self: flex-start;
      margin-right: auto;
    }
.ai-chat-meta {
      font-weight: 600;
      font-size: 12px;
      display: block;
      margin-bottom: 4px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }
.js-plotly-plot .plotly .modebar {
      top: -28px !important;
    }
.ai-chat-logo img {
      max-width: 100%;
      height: auto;
      max-height: 120px;
    }
/* AI chat panel - match card heights */
    .ai-chat-panel-container {
      height: calc(100vh - 100px);
      display: flex;
      flex-direction: column;
      overflow-y: auto;
      margin-top: 20px;
    }
/* The four tissue x modality cells INSIDE the Cross-Modal grid card. They are
   drawn like a card and they are not one: `.ie-card-stack > .ie-card` is how
   tests/tab_smoke.R counts the cards a tab mounted, so four inner cells wearing
   .ie-card inside their own stack were counted as four extra cards and that tab
   failed its own harness with every card drawing correctly. Same box, its own
   name -- so the count is of registry cards only and the cells still read as
   cells. Kept here rather than in www/modernist-shiny.css only because this
   change did not own that sheet; it belongs beside .ie-card. */
    .ie-grid-cell {
      background: var(--color-bg);
      border: 2px solid var(--color-neutral-300);
      padding: 18px 20px 20px;
    }
"))),
  # The layout jQuery that lived here is GONE with the shell it served. It read
  # `$('#tabs li.active a').text()` -- the Bootstrap-3 tab bar, which no longer
  # exists -- and resized `.equal-height-panel` / `.main-content-panel`, two
  # classes the rail replaced. Under the flex shell nothing needs measuring: the
  # panels and the control panel are flex siblings and the conditionalPanel
  # takes the control panel out of the row by itself. Left as a note rather than
  # silently deleted, because "the sidebar stopped hiding on Viewer" would
  # otherwise look like a regression rather than a rule that moved to CSS.
  tags$a(class = "skip-link", href = "#main-content", "Skip to main content"),

  div(
    class = "ie-shell",

    # ---- Left rail: every destination, in one place ------------------------
    # Replaces the horizontal tabsetPanel. The VALUES are the old tab labels,
    # so `input$tabs` still reads "Plot" / "Trajectory" / ... and every
    # active-tab guard and conditionalPanel downstream is unchanged.
    railUI(
      "rail", IE_PIPELINE, IE_REFERENCE,
      brand = "Islet Explorer", logo = "logo.png",
      help = tagList(
        rail_help_link("help_rules", "Phenotype rules", "Marker definitions"),
        rail_help_link("help_tips", "Tips for this tab", "How to read it")
      ),
      footer = tagList(
        rail_palette_select("rail_phenotype_palette", "Phenotype palette",
                            names(PHENOTYPE_PALETTES), "High Contrast"),
        rail_palette_select("rail_donor_palette", "Donor palette",
                            names(DONOR_PALETTES), "Bright")
      )
    ),

    div(
      class = "ie-main",

      # ---- Frozen TOP bar: what is identical on every tab -------------------
      # The split between here and the right-hand panel is by SIMILARITY, not by
      # kind. FILTER is which islets are in scope and is the same question on
      # every tab, so it earns a fixed horizontal strip the eye stops seeing;
      # SELECTION is what the tab plots and stays in the right-hand panel, where
      # attention already is.
      #
      # Before this, donor status lived under three different ids, min-cells
      # under four with no syncing at all, and Trajectory and Inflammation had
      # no cohort filter whatsoever.
      #
      # Hidden on the Viewer, which is a read-only image surface and reads none
      # of it.
      conditionalPanel(
        condition = "input.tabs != 'Viewer'",
        div(class = "ie-topbar", cohortFilterUI("cohort"))
      ),

      div(
        class = "ie-body",

        # The panels and the control panel are SIBLINGS in a flex row, each
        # with its own scroller, not nested. That is what .equal-height-panel
        # gave both of them before: the Plot sidebar is ~1700px of controls
        # against a ~1000px viewport, so sharing a scroller makes the whole
        # document grow to the sidebar's height and the shell stops containing
        # anything.
        div(
          class = "ie-panels", id = "main-content", role = "main", tabindex = "-1",

          # navset_hidden mounts EVERY panel and toggles CSS visibility, just
          # as tabsetPanel did. That is why the drill-down duplicate-output
          # guards in mod_plot_server.R and mod_trajectory_server.R are still
          # required. Order matches IE_PIPELINE: navset_hidden selects its
          # FIRST panel by default and railUI checks its first row, so if the
          # two disagree the app opens on one tab with another lit in the rail.
          bslib::navset_hidden(
            id = "tabs",
            bslib::nav_panel_hidden("Plot",
              fluidRow(
                plot_main_ui("plot", extra_panel = if (AI_ENABLED) column(2, ai_assistant_ui("ai")) else NULL)
              )
            ),
            bslib::nav_panel_hidden("Trajectory", trajectory_ui("traj")),
            bslib::nav_panel_hidden("Statistics", statistics_ui("stats")),
            bslib::nav_panel_hidden("Spatial",    spatial_ui("spatial")),
            # Same module namespace -- a relocation of the three neighborhood
            # cards, not a second module.
            bslib::nav_panel_hidden("Neighborhoods", spatial_neighborhood_ui("spatial")),
            bslib::nav_panel_hidden("CrossModal", crossmodal_ui("crossmodal")),
            bslib::nav_panel_hidden("LymphNode",  lymphnode_ui("lymphnode")),
            # Hidden unless SHOW_INFLAMMATION_TAB=1 (see 00_globals.R).
            # navset_hidden drops NULL entries, so this mounts nothing when off
            # -- and IE_PIPELINE is gated on the same flag, so the rail can
            # never offer a destination the tabset does not mount.
            if (SHOW_INFLAMMATION_TAB) {
              bslib::nav_panel_hidden("Inflammation", inflammation_ui("inflam"))
            },
            bslib::nav_panel_hidden("Viewer", viewer_ui("viewer"))
          )
        ),

        # ---- Frozen control panel ------------------------------------------
        # The Data & Filters controls, in the same place on every tab that
        # reads them, instead of a column that appeared and vanished. It stays
        # OUTSIDE the navset because statistics_server consumes
        # plot_returns$raw_df_prezero, so these inputs must stay mounted and
        # visible on the Statistics tab.
        #
        # Under BS5 Shiny gives the conditionalPanel wrapper `display: contents`
        # (modernist-shiny.css), so `.ie-aside` is a direct flex item of
        # `.ie-body` and the whole subtree leaves the row when the condition is
        # false, rather than leaving a blank gutter.
        # Hidden on the Viewer, which is a read-only image surface and reads
        # none of it. Each tab mounts its OWN block, namespaced to its module,
        # so nothing about the modules' wiring changes -- only where the
        # controls render.
        conditionalPanel(
          condition = "input.tabs != 'Viewer'",
          tags$aside(
            class = "ie-aside", role = "complementary",
            `aria-label` = "Controls for the active tab",
            div(
              class = "ie-aside-block",
              conditionalPanel(
                condition = "input.tabs == 'Plot' || input.tabs == 'Statistics'",
                plot_sidebar_ui("plot")
              ),
              conditionalPanel(
                condition = "input.tabs == 'Spatial' || input.tabs == 'Neighborhoods'",
                spatial_selection_ui("spatial")
              )
            )
          )
        )
      )
    )
  )
)

# Fail loudly at startup if a rail key drifts from a nav_panel_hidden value. A
# drifted key does not error at runtime -- it silently greys out a panel, which
# is exactly the failure the active-tab guards exist to prevent.
stopifnot(setequal(
  IE_TAB_KEYS,
  c("Plot", "Trajectory", "Statistics", "Spatial", "Neighborhoods",
    "CrossModal", "LymphNode", "Viewer",
    if (SHOW_INFLAMMATION_TAB) "Inflammation")
))

# ---------- The rail's "Tips for this tab" modal ----------
#
# `ie_tab_tips_modal()` was CALLED at the bottom of this file and defined
# nowhere: it and IE_TAB_TIPS were dropped from R/00_globals.R by 4335647 -- a
# commit about TissUUmaps layer opacity -- while this call site stayed. The
# rail's second Help entry therefore threw `could not find function
# "ie_tab_tips_modal"` from EVERY tab, and tests/parse_check.R could not see it,
# because parsing a call proves nothing about the callee. Restored here, beside
# its only caller; test-card-registry.R now drives the rail's Help observers, so
# a second deletion fails a test instead of reaching a reader.
#
# NOTHING IN THIS MODAL DESCRIBES THE DATA. Every section is read back from
# something the app already declares, so none of it can drift into a claim about
# the cohort, the donors or the markers:
#
#   what the tab is         the rail's own `sub`, from IE_PIPELINE / IE_REFERENCE
#   what is mounted on it   that tab's card registry, where it has one
#   how its controls work   IE_TAB_TIPS below -- one line per control that exists
#                           in the module source, cited above each entry
#   what an absence means   IE_STATES and ie_state_label(), from ieCore
#
# A tab with no IE_TAB_TIPS entry still gets the other three sections, so the
# Help row never opens an empty modal and never invents a fourth.

IE_TAB_TIPS <- list(
  # checkboxInput(ns("show_points"), "Show individual points") -- R/mod_plot_ui.R;
  # the plotly_click observer and its own nudge message -- R/mod_plot_server.R;
  # selectInput(ns("plot_color_by"), "Color points by:", "Donor ID") -- R/mod_plot_ui.R.
  Plot = list(
    paste("Tick “Show individual points”, then click a point on the islet-size",
          "scatter to open that islet's segmentation map below. Without points shown",
          "there is nothing to click."),
    paste("With points shown, set “Color points by” to Donor ID; clicking a donor",
          "in the legend hides or shows that donor."),
    paste("The cards do not resize themselves. Use browser zoom (Ctrl +/−), and",
          "scroll inside a card to reach controls below its fold.")
  ),
  # observeEvent(event_data("plotly_click", source = ns("traj_scatter"))) -- R/mod_trajectory_server.R.
  Trajectory = list(
    paste("Click a point on the pseudotime scatter to open the segmentation map for",
          "that islet.")
  ),
  # plotOutput(dblclick = ns("scatter_dblclick"), brush = brushOpts(ns("scatter_brush")))
  # -- R/mod_spatial_ui.R; both handlers -- R/mod_spatial_server.R.
  Spatial = list(
    paste("Drag on the tissue scatter to zoom into the brushed region; double-click",
          "anywhere on it to go back to the whole section.")
  )
)

#' The titles of the cards a registry mounts on `key`, in display order.
#'
#' Read from the registries themselves rather than restated, so a card added to
#' or removed from a registry changes this modal with no edit here -- which is
#' the same rule tab_smoke.R and test-card-registry.R follow. A tab with no
#' registry gets character(0), not an empty bullet list.
.ie_tab_card_titles <- function(key) {
  registries <- list(
    list(tab = get0("IE_CROSSMODAL_TAB", ifnotfound = NULL),
         cards = get0("IE_CROSSMODAL_CARDS", ifnotfound = NULL)),
    list(tab = get0("IE_LYMPHNODE_TAB", ifnotfound = NULL),
         cards = get0("IE_LYMPHNODE_CARDS", ifnotfound = NULL))
  )
  for (reg in registries) {
    if (!is.null(reg$tab) && identical(reg$tab, key) && length(reg$cards)) {
      return(vapply(reg$cards, `[[`, "", "title"))
    }
  }
  character(0)
}

#' What every card on every tab means by a state it cannot draw.
#'
#' The states and their labels come from ieCore (IE_STATES, ie_state_label()) and are walked in
#' whatever order that vector declares, so this list cannot fall out of step with the notes the
#' cards actually render -- it picked up the sixth state, `uncalled`, on 2026-09-18 without an
#' edit here. The closing line is the repo's rule, not a reading of any artifact: an absence is
#' not a zero.
.ie_tab_states_ui <- function() {
  states <- get0("IE_STATES", ifnotfound = NULL)
  label <- get0("ie_state_label", mode = "function", ifnotfound = NULL)
  if (is.null(states) || is.null(label)) return(NULL)
  shiny::tagList(
    shiny::tags$p(shiny::tags$strong("When a card cannot show something"),
      " it names the state and the reason it carries, in the place the figure",
      " would have been:"),
    shiny::tags$ul(lapply(states, function(s)
      shiny::tags$li(shiny::tags$code(s), " — ", label(s)))),
    shiny::tags$p(style = "font-size: 12px",
      "A missing measurement is unavailable. It is not a zero, and a card will",
      " not draw one in its place.")
  )
}

#' The Help > "Tips for this tab" modal, for whichever tab is active.
#'
#' @param tab one of IE_TAB_KEYS, or NULL/"" before a tab has been selected.
ie_tab_tips_modal <- function(tab) {
  key <- if (is.null(tab) || !length(tab)) "" else as.character(tab)[1]
  if (is.na(key)) key <- ""
  entry <- Filter(function(d) identical(d$key, key), c(IE_PIPELINE, IE_REFERENCE))
  entry <- if (length(entry)) entry[[1]] else NULL
  label <- if (is.null(entry)) "This tab" else entry$label
  cards <- .ie_tab_card_titles(key)
  tips <- IE_TAB_TIPS[[key]]

  shiny::modalDialog(
    title = paste0(label, " — how to read it"),
    easyClose = TRUE, footer = shiny::modalButton("Close"),
    if (!is.null(entry))
      shiny::tags$p(shiny::tags$strong("This tab: "), entry$sub),
    if (length(cards)) shiny::tagList(
      shiny::tags$p(shiny::tags$strong("Cards on this tab, in order:")),
      shiny::tags$ul(lapply(cards, shiny::tags$li))
    ),
    if (length(tips)) shiny::tagList(
      shiny::tags$p(shiny::tags$strong("Working the controls:")),
      shiny::tags$ul(lapply(tips, shiny::tags$li))
    ),
    .ie_tab_states_ui(),
    shiny::tags$p(style = "font-size: 12px",
      "The marker rules behind every phenotype are in the Help entry above this one.")
  )
}


# ---------- Server ----------

server <- function(input, output, session) {

  # ---- Shared reactive state ----
  forced_image   <- reactiveVal(NULL)
  selected_islet <- reactiveVal(NULL)

  # ---- Data loading chain ----
  # Priority: H5AD (islet_explorer.h5ad) > Excel (master_results.xlsx)
  validate_file <- reactive({
    has_h5ad <- !is.null(h5ad_path) && file.exists(h5ad_path)
    has_excel <- file.exists(master_path)
    missing_message <- if (islet_is_production()) {
      "Application data are unavailable."
    } else {
      paste("No data found. Checked:", h5ad_path, "and", master_path)
    }
    shiny::validate(shiny::need(has_h5ad || has_excel, missing_message))
    TRUE
  })

  master <- reactive({
    req(validate_file())
    load_master_auto(h5ad_path = h5ad_path, excel_path = master_path)
  })

  prepared <- reactive({
    pd <- prep_data(master())
    try({
      audit_na(pd$targets_all, "targets_all")
      audit_na(pd$markers_all, "markers_all")
      audit_na(pd$comp, "comp")
    }, silent = TRUE)
    pd
  })

  # ---- Module servers ----

  # ---- The shared cohort filter ----
  # One definition of "which islets are in scope", mounted once and consumed by
  # every tab. `cohort$apply(df)` returns the subset AND attaches an audit of
  # which filters bound to which column, so a filter this frame cannot apply
  # surfaces as a banner rather than silently doing nothing.
  cohort <- cohortFilterServer("cohort", prepared)

  # Active tab reactive (used to prevent duplicate non-namespaced output IDs)
  active_tab <- reactive(input$tabs)

  # The rail drives the hidden tabset. Its values ARE the old tab labels, so
  # `input$tabs` keeps reading "Plot" / "Trajectory" / ... and every guard
  # downstream of active_tab() is unchanged by the shell.
  railServer("rail", "tabs")

  # Active phenotype palette (global, consumed by root renders + spatial module)
  # Writers: the rail footer and the drill-down panel's own selector. The
  # Spatial tab used to carry a third copy; it went with the control card that
  # moved into the frozen panel.
  palette_name <- reactiveVal("High Contrast")

  ie_sync_selects(input, session, palette_name,
                  c("rail_phenotype_palette", "drilldown_palette"))

  active_palette <- reactive({
    get_phenotype_palette(palette_name())
  })

  # Active donor color palette (global, consumed by Plot, Trajectory, Statistics, Spatial)
  donor_palette_name <- reactiveVal("Bright")
  DONOR_PALETTE_CHOICES <- names(DONOR_PALETTES)

  # Sync every donor-palette selector, the rail's included.
  # ONE writer left for each palette, plus the drill-down panel's own phenotype
  # selector. The Plot sidebar's, Trajectory's and Spatial's copies all went to
  # the rail footer. ie_sync_selects() stays because the drill-down copy is real
  # -- it sits inside a renderUI'd panel and users reach it there.
  ie_sync_selects(input, session, donor_palette_name, c("rail_donor_palette"))

  active_donor_colors <- reactive({
    get_donor_color_palette(donor_palette_name())
  })

  # GONE WITH THE COHORT PANEL: two reactiveVals and six observers that mirrored
  # `sidebar_remove_outliers` / `traj_remove_outliers` / `spatial_remove_outliers`
  # and the matching z-threshold sliders into each other. Three copies of one
  # value, each of which had to name the other two -- and Trajectory and
  # Inflammation still got no cohort filter at all. The panel owns the rule now;
  # modules read `cohort$outlier_enabled()` and `cohort$outlier_z()`.
  #
  # Unchanged by the move: Spatial Card A's IQR rule (1.5x) is still NOT routed
  # through this. Different scale, deliberately fixed for zero-inflated
  # phenotype enrichment.

  # Plot module: returns list(raw_df, summary_df, get_selection_description)
  plot_returns <- plot_server("plot", prepared, selected_islet, active_tab, active_donor_colors,
                              remove_outliers = cohort$outlier_enabled,
                              outlier_threshold = cohort$outlier_z,
                              cohort = cohort)

  # Trajectory module
  trajectory_server("traj", prepared, selected_islet, forced_image, active_tab, active_donor_colors,
                    remove_outliers = cohort$outlier_enabled,
                    outlier_threshold = cohort$outlier_z,
                    cohort = cohort)

  # Viewer module
  viewer_server("viewer", forced_image, reactive(input$tabs))

  # Both tabs are thin registries over card_*.R files; see mod_crossmodal_ui.R.
  crossmodal_server("crossmodal", prepared, cohort = cohort,
                    active_tab = reactive(input$tabs))
  lymphnode_server("lymphnode", prepared, cohort = cohort,
                   active_tab = reactive(input$tabs))

  # Statistics module: wired to Plot module's reactive outputs.
  # It takes `raw_df_prezero`, NOT `raw_df`, so the Plot sidebar's zero-exclusion
  # checkbox cannot silently move statistical results; Statistics applies its own
  # `stats_exclude_zero` and reports that choice in the CSV header. `summary_df`
  # is deliberately not passed — this module derives AUC itself rather than
  # integrating the Plot tab's binned summary.
  statistics_server("stats",
                    raw_df                  = plot_returns$raw_df_prezero,
                    get_selection_description = plot_returns$get_selection_description,
                    donor_colors             = active_donor_colors,
                    remove_outliers          = cohort$outlier_enabled,
                    cohort                   = cohort)

  # Spatial Neighborhood module
  spatial_server("spatial", prepared, active_palette, active_donor_colors,
                 remove_outliers = cohort$outlier_enabled,
                 outlier_threshold = cohort$outlier_z,
                 cohort = cohort)

  # Inflammation module. Gated on the same flag as its tabPanel (see
  # 00_globals.R) so the module is not wired up when the tab is hidden.
  if (SHOW_INFLAMMATION_TAB) {
    inflammation_server("inflam", prepared, active_donor_colors, cohort = cohort)
  }

  # AI Assistant module (self-contained; gated so public builds omit it — see AI_ENABLED)
  if (AI_ENABLED) ai_assistant_server("ai")

  # ---- Unified segmentation + drilldown plot (root-level, used by both Plot and Trajectory panels) ----
  islet_segmentation_view_plot <- reactive({
    req(selected_islet())
    info <- selected_islet()
    show_peri_bd <- if (!is.null(input$drilldown_show_peri_boundary)) input$drilldown_show_peri_boundary else TRUE
    show_structs <- if (!is.null(input$drilldown_show_structures)) input$drilldown_show_structures else TRUE
    show_cells <- if (!is.null(input$drilldown_show_cells)) input$drilldown_show_cells else TRUE

    if (show_cells) {
      cells <- load_islet_cells(info$case_id, info$islet_key)
      if (!is.null(cells)) {
        color_by <- input$drilldown_color_by %||% "phenotype"
        if (!islet_valid_scalar_choice(
          color_by, colnames(cells), "drilldown_color_invalid"
        )) color_by <- "phenotype"
        show_peri_cells <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
        return(tryCatch(
          render_islet_drilldown_plot(info, cells, color_by = color_by, show_peri = show_peri_cells,
                                      show_peri_boundary = show_peri_bd, show_structures = show_structs,
                                      palette = active_palette()),
          error = function(e) {
            message("event=drilldown_render_failed")
            ggplot2::ggplot() +
              ggplot2::annotate("text", x = 0.5, y = 0.5,
                                label = "Unable to render this view.",
                                size = 4, color = "red") +
              ggplot2::theme_void() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
          }
        ))
      }
    }
    # Fallback: boundaries only (no cells or cells unavailable)
    tryCatch(
      render_islet_segmentation_plot(info, show_peri_boundary = show_peri_bd,
                                      show_structures = show_structs),
      error = function(e) {
        message("event=segmentation_render_failed")
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = "Unable to render this view.",
                            size = 4, color = "red") +
          ggplot2::theme_void() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
      }
    )
  })
  output$islet_segmentation_view <- renderPlot({ islet_segmentation_view_plot() })
  output$dl_islet_segmentation_view <- downloadHandler(
    filename = function() {
      info <- selected_islet()
      stem <- if (!is.null(info)) paste0("islet_", info$case_id, "_", info$islet_key) else "islet_segmentation"
      islet_png_filename(stem)
    },
    content = function(file) {
      p <- islet_segmentation_view_plot(); if (is.null(p)) return()
      ggplot2::ggsave(file, plot = p,
                      width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                      dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
    }
  )

  islet_drilldown_summary_plot <- reactive({
    req(selected_islet())
    info <- selected_islet()
    cells <- load_islet_cells(info$case_id, info$islet_key)
    req(cells)
    show_peri <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
    if (!show_peri && "cell_region" %in% colnames(cells)) {
      cells <- cells[cells$cell_region == "core", , drop = FALSE]
    }
    tryCatch(render_drilldown_summary(cells, palette = active_palette()), error = function(e) NULL)
  })
  output$islet_drilldown_summary <- renderPlot({ islet_drilldown_summary_plot() })
  output$dl_islet_drilldown_summary <- downloadHandler(
    filename = function() {
      info <- selected_islet()
      stem <- if (!is.null(info)) paste0("islet_", info$case_id, "_", info$islet_key, "_composition") else "islet_composition"
      islet_png_filename(stem)
    },
    content = function(file) {
      p <- islet_drilldown_summary_plot(); if (is.null(p)) return()
      ggplot2::ggsave(file, plot = p,
                      width = ISLET_PNG_DIMS$width, height = ISLET_PNG_DIMS$height,
                      dpi = ISLET_PNG_DIMS$dpi, units = ISLET_PNG_DIMS$units)
    }
  )

  # ---- Marker-rule legend for the phenotypes present in the selected islet ----
  output$islet_drilldown_rules <- renderUI({
    req(selected_islet())
    info <- selected_islet()
    cells <- load_islet_cells(info$case_id, info$islet_key)
    req(cells, "phenotype" %in% colnames(cells))
    present <- sort(unique(as.character(cells$phenotype)))
    if (length(present) == 0) return(NULL)
    pal <- active_palette()
    rows <- lapply(present, function(nm) {
      colour <- if (!is.null(pal) && nm %in% names(pal)) pal[[nm]] else "#CCCCCC"
      rule_html <- format_phenotype_rule(nm, "html")
      tags$div(
        style = "display: flex; align-items: flex-start; gap: 8px; padding: 4px 0",
        tags$span(style = sprintf("display:inline-block; width:12px; height:12px; min-width:12px; margin-top:4px; background:%s; border:1px solid #888; border-radius:2px;", colour)),
        tags$div(
          tags$div(style = "font-weight:600; font-size:13px", nm),
          tags$div(style = "font-size:12px; line-height:1.35", rule_html)
        )
      )
    })
    tags$details(
      style = "margin-top: 10px; padding: 8px 10px; font-size: 13px",
      tags$summary(style = "cursor: pointer; font-weight: 600",
                   "Marker rules for phenotypes in this islet"),
      tags$div(style = "margin-top: 8px", rows)
    )
  })

  output$islet_drilldown_table <- renderTable({
    req(selected_islet())
    info <- selected_islet()
    cells <- load_islet_cells(info$case_id, info$islet_key)
    req(cells)
    show_peri <- if (!is.null(input$drilldown_show_peri)) input$drilldown_show_peri else TRUE
    if (!show_peri && "cell_region" %in% colnames(cells)) {
      cells <- cells[cells$cell_region == "core", , drop = FALSE]
    }
    if ("cell_region" %in% colnames(cells)) {
      counts <- as.data.frame(table(Region = cells$cell_region), stringsAsFactors = FALSE)
      colnames(counts) <- c("Region", "Cells")
      counts$Region <- ifelse(counts$Region == "core", "Core", "Peri-islet")
      rbind(counts, data.frame(Region = "Total", Cells = sum(counts$Cells)))
    } else {
      data.frame(Region = "All", Cells = nrow(cells))
    }
  }, striped = TRUE, spacing = "s")

  # ---- Phenotype Rules modal (shared across all tabs) ----
  # ONE observer now. Five tabs each mounted their own
  # `phenotype_rules_button("show_phenotype_rules_<tab>")` pill at the top of
  # the panel, and a sixth -- Inflammation's -- had no observer at all, so its
  # pill had always been inert. The rail's Help entry replaces all six and
  # reaches the modal from every tab, including the ones that never had a pill.
  # The INLINE links inside the Spatial module are separate and still live;
  # they are namespaced and registered in mod_spatial_server.R.
  observeEvent(input$help_rules, { showModal(phenotype_rules_modal_ui()) },
               ignoreInit = TRUE)

  # The rail's second Help entry. One observer for every tab, because the modal
  # reads the ACTIVE tab rather than being mounted per panel -- which is the
  # whole point of moving the hint strips off the panels.
  observeEvent(input$help_tips, {
    showModal(ie_tab_tips_modal(active_tab()))
  }, ignoreInit = TRUE)
}

shinyApp(ui, server)
