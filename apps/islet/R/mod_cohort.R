# ---------- Shared cohort filter module ----------
# Mounted ONCE, in the frozen right-hand panel, visible on every analysis tab.
# Thin reactive wrapper around the pure functions in R/cohort_filter.R.
#
# PORTED FROM JUNIOR (f827106). `GROUP_LEVELS` became `domain_group_levels()`
# throughout -- see the note at the top of R/cohort_filter.R for why a literal
# there would empty every panel on the other domain rather than error. The age
# and diameter bounds were ALREADY derived from the data rather than fixed, and
# that carries over unchanged: this cohort's ages are not the other's.
#
# NOT HERE, DELIBERATELY: the delivered mod_cohort.R also carried
# normalize_features(), which robust-z-scored every feature column and returned
# the transformed frame as cohort$islets(). Removed. (a) The Statistics tab
# reports the UNSTANDARDISED difference in the feature's own units with a Welch
# CI plus Glass's delta -- z-scoring would make every effect size a function of
# who else is in the cohort. (b) The Inflammation composition test needs cell
# counts. (c) Its zero-MAD branch returned merely-CENTRED values, so a frame
# could silently mix z-scores and centred values in one column -- and `prop_*`
# is mostly zero. The panel owns ROW SELECTION; it never transforms a value.

cohortFilterUI <- function(id) {
  ns <- NS(id)
  tagList(
    # RENDERED, not static: the domain declares the vocabulary but the loaded
    # dataset decides which of it exists. That is a fact about the BUILD, not about
    # either domain -- and it has already changed once. The run cohort carried no
    # Aab+ donor before 2026-09-15 and carries three now (6450, 6505, 6521; 3,966 of
    # 18,384 islets, measured 2026-09-18), so a control that hardcoded either state
    # would have been wrong on one side of that date. A static control also invites
    # the reader to filter on a category whose box can never contain an islet, which
    # with this filter is worse than cosmetic: unticking every PRESENT group empties
    # every panel while a phantom stays ticked.
    uiOutput(ns("groups_field")),
    uiOutput(ns("age_field")),
    uiOutput(ns("sex_field")),
    ie_field("Min cells / islet",
      numericInput(ns("min_cells"), NULL, value = 1, min = 1, max = 500,
                   step = 1, width = "100%")),
    uiOutput(ns("diam_field")),
    uiOutput(ns("aab_field")),
    div(title = GLOBAL_OUTLIER_TOOLTIP,
      ie_field("Outliers",
        ie_segmented(ns("outlier_enabled"),
                     c(Keep = "keep", Exclude = "exclude"), "keep"))),
    conditionalPanel(
      sprintf("input['%s'] == 'exclude'", ns("outlier_enabled")),
      ie_field("z-threshold",
        sliderInput(ns("outlier_z"), NULL, min = 0.5, max = 10, value = 3,
                    step = 0.5, ticks = FALSE, width = "100%"))),
    # NO ISLET/DONOR READOUT. It printed "N islets - M donors" in the corner of
    # the cohort bar. The counts are reported with the results, on the tabs that
    # compute them, where they mean something specific; a running total in the
    # chrome is a number with no question attached to it.
    #
    # The low-donor warning it also carried is NOT lost: the Statistics tab
    # gates its own tests on group coverage and says so there.
  )
}

#' @param prepared reactive returning prep_data() output (needs $comp).
cohortFilterServer <- function(id, prepared) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    donors <- reactive({
      pd <- prepared(); req(pd$comp)
      pd$comp[!duplicated(pd$comp$`Case ID`), , drop = FALSE]
    })

    # Ranges are DERIVED from the data, not hardcoded. The delivered module fixed
    # the age slider at 0-70; this cohort's ages are bimodal at 1.58-4.22 and
    # 14.56-25.06, so a fixed 0-70 slider makes the whole usable range a few
    # pixels wide.
    age_bounds <- reactive({
      a <- suppressWarnings(as.numeric(donors()$age)); a <- a[is.finite(a)]
      if (!length(a)) NULL else c(floor(min(a)), ceiling(max(a)))
    })
    diam_bounds <- reactive({
      pd <- prepared()
      d <- suppressWarnings(as.numeric(pd$comp$islet_diam_um))
      d <- d[is.finite(d)]
      if (!length(d)) NULL else c(floor(min(d)), ceiling(max(d)))
    })

    statuses_present <- reactive({
      pd <- prepared()
      cols <- cohort_resolve_cols(pd$comp)
      if (is.na(cols$status)) return(domain_group_levels())
      present <- unique(as.character(pd$comp[[cols$status]]))
      hit <- intersect(domain_group_levels(), present)
      if (length(hit)) hit else domain_group_levels()
    })

    output$groups_field <- renderUI({
      g <- statuses_present()
      ie_field("Disease groups",
        checkboxGroupInput(ns("groups"), NULL, choices = g, selected = g, inline = TRUE))
    })

    output$age_field <- renderUI({
      b <- age_bounds(); if (is.null(b)) return(NULL)
      ie_field("Donor age", sliderInput(ns("age"), NULL, min = b[1], max = b[2],
                                        value = b, step = 1, post = " yr",
                                        ticks = FALSE, width = "100%"))
    })

    output$diam_field <- renderUI({
      b <- diam_bounds(); if (is.null(b)) return(NULL)
      ie_field("Islet diameter", sliderInput(ns("diam"), NULL, min = b[1], max = b[2],
                                             value = b, step = 5, post = " µm",
                                             ticks = FALSE, width = "100%"))
    })

    # RE-RENDERED, never messaged. ie_segmented() emits label.ie-seg-opt, and
    # Shiny's RadioInputBinding.receiveMessage replaces options by removing
    # div.shiny-options-group / label.radio -- neither of which this markup has.
    # updateRadioButtons(choices=) therefore APPENDS a second Bootstrap option
    # set that shares the group name, unchecking the hand-rolled one, and the
    # control renders visually empty while input$<id> still returns a value.
    output$sex_field <- renderUI({
      g <- sort(unique(stats::na.omit(as.character(donors()$gender))))
      if (!length(g)) return(NULL)
      lab <- c(F = "Female", M = "Male")
      ch <- c(stats::setNames("__all__", "All"),
              stats::setNames(g, ifelse(g %in% names(lab), lab[g], g)))
      # "Sex" not "Gender" in every user-facing label; the data column stays `gender`.
      ie_field("Sex", ie_segmented(ns("sex"), ch, "__all__"))
    })

    output$aab_field <- renderUI({
      avail <- intersect(COHORT_AAB_COLS, colnames(donors()))
      if (!length(avail)) return(NULL)
      lab <- c(AAb_GADA = "GADA", AAb_IA2A = "IA2A", AAb_ZnT8A = "ZnT8A",
               AAb_IAA = "IAA", AAb_mIAA = "mIAA")
      ie_field("Autoantibodies (Aab+ only)",
        checkboxGroupInput(ns("aab_flags"), NULL,
          choices = stats::setNames(avail, unname(lab[avail])),
          selected = avail, inline = TRUE))
    })

    # DEBOUNCE 400 ms, not the delivered 150. Shiny's SliderInputBinding already
    # debounces 250 ms client-side, so anything below that merges nothing during
    # a drag and only adds latency. A full Plot flush (two ggplotly figures) is
    # ~530 ms server-side, so 400 ms coalesces a drag without letting the UI
    # queue. Statistics is unaffected either way -- stats_run() is an
    # eventReactive on the Run button.
    spec <- debounce(reactive({
      ab <- age_bounds(); db <- diam_bounds()
      list(
        groups    = input$groups %||% statuses_present(),
        age       = if (is.null(input$age)) ab else input$age,
        sex       = if (is.null(input$sex) || identical(input$sex, "__all__")) NULL else input$sex,
        aab_flags = input$aab_flags %||% COHORT_AAB_COLS,
        min_cells = { v <- suppressWarnings(as.numeric(input$min_cells))
                      if (length(v) == 1 && is.finite(v) && v >= 1) v else 1 },
        diam      = if (is.null(input$diam)) db else input$diam,
        outlier_enabled = identical(input$outlier_enabled, "exclude"),
        outlier_z = { v <- suppressWarnings(as.numeric(input$outlier_z))
                      if (length(v) == 1 && is.finite(v) && v >= 0.5) v else 3.0 }
      )
    }), 400)

    # The delivered readout used req(input$groups), so unticking every group made
    # it FREEZE on its last value -- the most confusing possible response to
    # "I have selected nothing".

    list(
      spec            = spec,
      apply           = function(df, cells_expr = "core_peri") cohort_apply(df, spec(), cells_expr),
      donor_ids       = function(tbl = NULL) cohort_donor_ids(tbl %||% donors(), spec()),
      groups          = reactive(spec()$groups),
      outlier_enabled = reactive(spec()$outlier_enabled),
      outlier_z       = reactive(spec()$outlier_z),
      provenance      = function(audit = NULL) cohort_provenance(spec(), audit)
    )
  })
}
