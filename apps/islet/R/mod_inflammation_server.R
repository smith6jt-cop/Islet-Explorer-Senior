# ---------- Inflammation tab module server ----------
# Exports: inflammation_server(id, prepared, donor_colors)
#
# Two cohort-level analyses adapted from Steenbuck et al. 2025 (Bodenmiller lab,
# bioRxiv 2025.03.05.641526), plus the nPOD donor table they need for context.
#
# Dependencies (auto-sourced):
#   load_cell_distance()                   -- spatial_helpers.R; the whole
#       `cell_distance` DuckDB view, cached per domain root. It was named here
#       and defined nowhere until 2026-09-18, so this tab could not be enabled.
#   PHENOTYPE_COLORS                       -- drilldown_helpers.R
#   theme_islet, plotly_islet_*, ISLET_*   -- 00_globals.R
#   write_islet_csv, islet_csv_header      -- 00_globals.R
#   format_p_value                         -- utils_stats.R
#
# NOTE: validate/need are called as shiny::validate(shiny::need(...)). jsonlite
# masks `validate` on this app's search path, so the bare form resolves to the
# wrong function and dies with "is.character(txt) is not TRUE". app.R:420 uses
# the namespaced form for the same reason.
#
# Both heavy computations are plain `reactive()`s so they run once per session
# per setting and are shared by the plot, the table and the CSV download — the
# rule learned the hard way in the Jul 2026 CSV pass: an export must consume the
# same reactive the plot does, never re-derive its own subset.

# PRE-SPECIFIED multiple-testing family for the infiltration score.
#
# These seven immune phenotypes are the ones the score exists to interrogate,
# and fixing the family in source is what makes the correction honest: it cannot
# be widened or narrowed after seeing a result. Size matters enormously at this
# cohort size -- with 5 donors per stage the best achievable Kruskal-Wallis p is
# 0.00193 (perfect separation of all 15 donors), so BH across 21 phenotypes has
# a floor of q = 0.041 while BH across these 7 has a floor of 0.014. The wider
# family left almost no headroom for a real effect to survive.
#
# Phenotypes outside this list can still be displayed, but they are reported
# UNCORRECTED and flagged as exploratory rather than being folded into the
# family after the fact.
#
# EVERY NAME IN IT IS THE LEGACY CODEX VOCABULARY, and so are the cohort numbers above -- 5
# donors per stage, 15 in total, 21 phenotypes. The PhenoCycler run cohort resolves immune as ONE
# label: its eight classes contain `Immune` and no subset of it, because no lane has established
# an immune subtyping for release (the second pass in core/multiplex5000/scripts/typing/
# immune_subtype.py never modifies cell_class, and no app artifact reads its output). Under that
# domain six of these seven names match nothing, the family collapses to `Immune` alone, and the
# floor figures quoted above do not apply. Do not restate this list as though it were a
# cross-domain immune vocabulary.
INFLAM_FDR_FAMILY <- c("CD8a Tcell", "CD4 Tcell", "B cell",
                       "T cell", "Macrophage", "APCs", "Immune")
INFLAM_DEFAULT_PHENOTYPES <- INFLAM_FDR_FAMILY

#' @param cohort the shared cohort filter from cohortFilterServer(), or NULL.
#'   NULL keeps this module's own inputs, which is what the testthat suite and
#'   any standalone use get; app.R always passes one.
inflammation_server <- function(id, prepared, donor_colors = reactive(DONOR_COLORS),
                                cohort = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # DOMAIN, not a literal: this is the other cohort's vocabulary, and it
    # orders every stage contrast in this module. The run cohort declares its
    # own levels in config/run20_default.yml.
    status_levels <- domain_group_levels()

    # =========================================================================
    # Shared data
    # =========================================================================

    # Tissue-wide single-cell table with the signed distance to the islet edge.
    cells <- reactive({
      df <- load_cell_distance()
      # NAME THE SOURCE THAT EXISTS. This said "data/cell_distance.csv not found
      # -- run scripts/compute_neighborhood_metrics.py", and neither the file nor
      # that output has been the source since the view moved to parquet: it is
      # <domain root>/parquet/cell_distance/, built by
      # scripts/senior/build_cell_distance_parquet.py and registered as the
      # `cell_distance` DuckDB view in R/00_globals.R.
      shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0,
                    paste("No cell_distance view for this domain -- build",
                          "parquet/cell_distance under the active domain\u2019s data root",
                          "with scripts/senior/build_cell_distance_parquet.py.")))
      df$donor_status <- factor(as.character(df$donor_status), levels = status_levels)
      df <- df[!is.na(df$donor_status), , drop = FALSE]
      # A CELL-level frame: it carries case_id and donor_status but no age, sex,
      # AAb, diameter or per-islet counts, so the panel's donor-scope filters
      # resolve against the islet table and land here as an id membership test.
      if (!is.null(cohort)) {
        keep <- cohort$donor_ids()
        cid <- intersect(c("case_id", "Case ID", "imageid"), colnames(df))
        if (length(keep) && length(cid)) {
          df <- df[as.character(df[[cid[1]]]) %in% keep, , drop = FALSE]
        }
      }
      df
    })

    # One row per donor from the nPOD metadata joined into comp by prep_data().
    # RENAMED from `cohort`, which now names the shared cohort FILTER passed in
    # as a parameter -- a local of the same name inside moduleServer() shadows
    # it, so the module would have silently kept its old unfiltered behaviour
    # while app.R believed it had wired the panel through.
    donor_table <- reactive({
      pd <- prepared()
      shiny::validate(shiny::need(!is.null(pd$comp) && nrow(pd$comp) > 0, "No donor data loaded."))
      cn <- colnames(pd$comp)
      want <- c("Case ID", "Donor Status", "age", "age_years_npod", "age_group",
                "t1d_duration_years", "gender", "BMI", "HbA1c_percent",
                "C_peptide_ng_mL", "C_peptide_censored", "n_aab",
                "AAb_GADA", "AAb_IA2A", "AAb_mIAA", "AAb_ZnT8A")
      have <- intersect(want, cn)
      d <- unique(pd$comp[, have, drop = FALSE])
      d <- d[!duplicated(d$`Case ID`), , drop = FALSE]
      # THE COHORT PANEL REACHES THIS TAB FOR THE FIRST TIME. Inflammation had
      # no donor filter at all: every card ran on all donors regardless of what
      # the rest of the app was scoped to, so a reader who had narrowed the
      # cohort on Plot saw an unnarrowed answer here with nothing saying so.
      # min-cells and diameter are islet-level and do not apply to a donor row.
      # min-cells and diameter are ISLET-level and a donor row has neither
      # column, so they are switched OFF rather than left to report themselves
      # "not applied" -- an audit banner about a filter that could never apply
      # to this frame is noise, and noise is how a real one gets ignored.
      if (!is.null(cohort)) {
        d <- cohort_apply(d, utils::modifyList(cohort$spec(),
                                               list(min_cells = 1, diam = NULL)))
      }
      d$`Donor Status` <- factor(as.character(d$`Donor Status`), levels = status_levels)
      d[order(d$`Donor Status`, suppressWarnings(as.numeric(d$age))), , drop = FALSE]
    })

    npod_available <- reactive("t1d_duration_years" %in% colnames(donor_table()))

    # =========================================================================
    # Card A: donor cohort table
    # =========================================================================

    output$cohort_note <- renderUI({
      d <- donor_table()
      if (!npod_available()) {
        return(tags$div(style = paste0("background: #fdecea; border: 1px solid #f5b7b1;",
                                       " border-radius: 4px; padding: 8px 12px; margin-bottom: 10px;",
                                       " font-size: 12px; color: #8a1c12;"),
          tags$strong("nPOD clinical metadata not loaded. "),
          "Only age and sex are available. Run ",
          tags$code("python scripts/build_donor_metadata.py"),
          " to generate ", tags$code("data/donor_metadata.csv"), "."))
      }
      dur <- suppressWarnings(as.numeric(d$t1d_duration_years))
      dur <- dur[is.finite(dur)]
      dur_txt <- if (length(dur)) {
        sprintf("All %d T1D donors are recent-onset (duration %s year%s).",
                length(dur),
                if (diff(range(dur)) == 0) sprintf("%.2g", dur[1])
                else sprintf("%.2g–%.2g", min(dur), max(dur)),
                if (length(dur) == 1) "" else "s")
      } else ""
      # The one place the age discrepancy is visible rather than buried.
      mism <- if (all(c("age", "age_years_npod") %in% colnames(d))) {
        a <- suppressWarnings(as.numeric(d$age)); b <- suppressWarnings(as.numeric(d$age_years_npod))
        d$`Case ID`[is.finite(a) & is.finite(b) & abs(a - b) > 0.01]
      } else character(0)
      tagList(
        tags$p(style = "font-size: 12px; margin: 0 0 8px 0; line-height: 1.5",
          tags$strong("Scope. "), dur_txt,
          " This cohort therefore corresponds to the Steenbuck et al. “Onset” group; ",
          "it contains no long-duration T1D donors. Disease duration has too little ",
          "variance to serve as a covariate — it is reported here as a descriptor."),
        if (length(mism)) {
          tags$p(style = "font-size: 11px; margin: 0 0 8px 0",
            tags$strong("Note: "),
            sprintf("donor %s has a different age in the H5AD (%s) than in the nPOD record (%s). ",
                    paste(mism, collapse = ", "),
                    paste(sprintf("%.2f", suppressWarnings(as.numeric(d$age))[d$`Case ID` %in% mism]), collapse = ", "),
                    paste(sprintf("%.2f", suppressWarnings(as.numeric(d$age_years_npod))[d$`Case ID` %in% mism]), collapse = ", ")),
            "The H5AD value is used for all modelling so existing results do not shift.")
        } else NULL
      )
    })

    cohort_display <- reactive({
      d <- donor_table()
      fmt <- function(x, digits = 1) ifelse(is.na(x), "—", sprintf(paste0("%.", digits, "f"), as.numeric(x)))
      aab_cols <- intersect(c("AAb_GADA", "AAb_IA2A", "AAb_mIAA", "AAb_ZnT8A"), colnames(d))
      profile <- if (length(aab_cols)) {
        vapply(seq_len(nrow(d)), function(i) {
          pos <- sub("^AAb_", "", aab_cols[vapply(aab_cols, function(cc) isTRUE(as.logical(d[[cc]][i])), logical(1))])
          if (!length(pos)) "—" else paste(pos, collapse = ", ")
        }, character(1))
      } else rep("—", nrow(d))

      out <- data.frame(
        `Case ID` = as.character(d$`Case ID`),
        Status = as.character(d$`Donor Status`),
        `Age (y)` = fmt(d$age, 2),
        Sex = if ("gender" %in% colnames(d)) as.character(d$gender) else "—",
        check.names = FALSE, stringsAsFactors = FALSE
      )
      if (npod_available()) {
        cp <- suppressWarnings(as.numeric(d$C_peptide_ng_mL))
        cens <- if ("C_peptide_censored" %in% colnames(d)) as.logical(d$C_peptide_censored) else rep(FALSE, nrow(d))
        out$`Age group` <- as.character(d$age_group)
        out$`T1D duration (y)` <- ifelse(is.na(d$t1d_duration_years), "n/a",
                                         fmt(d$t1d_duration_years, 2))
        out$BMI <- fmt(d$BMI, 1)
        out$`HbA1c (%)` <- fmt(d$HbA1c_percent, 1)
        # The censored value keeps its "<" so the table never implies a precision
        # the assay did not deliver.
        out$`C-peptide (ng/mL)` <- ifelse(cens %in% TRUE,
                                          paste0("<", fmt(cp, 2)), fmt(cp, 2))
        out$`AAb+` <- profile
        out$`n AAb` <- as.character(d$n_aab)
      }
      out
    })

    output$cohort_table <- renderTable({ cohort_display() },
      striped = TRUE, bordered = TRUE, spacing = "s", align = "l")

    output$dl_cohort <- downloadHandler(
      filename = function() paste0("inflammation_donor_cohort_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        write_islet_csv(file, donor_table(), islet_csv_header(
          "Donor cohort — nPOD clinical metadata",
          "Source: scripts/build_donor_metadata.py -> data/donor_metadata.csv",
          "age = H5AD value (canonical for modelling); age_years_npod = nPOD record",
          "t1d_duration_years is NA for non-diabetic donors (not applicable, not zero)"))
      }
    )

    # =========================================================================
    # Card B: edgeR differential abundance
    # =========================================================================

    # Per-donor x cell-type count matrix. Counts, not proportions: that is the
    # entire point of this card.
    abundance <- reactive({
      df <- cells()
      region <- input$edger_region %||% "core_peri"
      keep <- switch(region,
                     core_peri = df$cell_region %in% c("core", "peri"),
                     core      = df$cell_region == "core",
                     rep(TRUE, nrow(df)))
      sub <- df[keep, , drop = FALSE]
      shiny::validate(shiny::need(nrow(sub) > 0, "No cells in the selected region."))
      tab <- table(sub$phenotype, sub$case_id)
      mat <- matrix(as.integer(tab), nrow = nrow(tab),
                    dimnames = list(rownames(tab), colnames(tab)))
      meta <- unique(sub[, c("case_id", "donor_status")])
      meta <- meta[match(colnames(mat), as.character(meta$case_id)), , drop = FALSE]
      list(counts = mat,
           group = factor(as.character(meta$donor_status), levels = status_levels),
           region = region, n_cells = nrow(sub))
    })

    edger_fit <- reactive({
      shiny::validate(shiny::need(requireNamespace("edgeR", quietly = TRUE),
                    "The edgeR package is not installed."))
      ab <- abundance()
      grp <- ab$group
      shiny::validate(shiny::need(nlevels(droplevels(grp)) >= 2,
                    "Need at least 2 donor-status groups."))

      y <- edgeR::DGEList(counts = ab$counts, group = grp)
      keep <- edgeR::filterByExpr(y, group = grp)
      y <- y[keep, , keep.lib.sizes = FALSE]
      y <- edgeR::calcNormFactors(y)
      design <- stats::model.matrix(~ grp)
      y <- edgeR::estimateDisp(y, design)
      fit <- edgeR::glmQLFit(y, design)

      # `~ grp` with levels ND / Aab+ / T1D gives coef 2 = Aab+ vs ND and
      # coef 3 = T1D vs ND; T1D vs Aab+ needs an explicit contrast.
      cn <- colnames(design)
      spec <- switch(input$edger_contrast %||% "t1d_nd",
        aab_nd  = list(coef = 2, label = "Aab+ vs ND"),
        t1d_nd  = list(coef = 3, label = "T1D vs ND"),
        t1d_aab = list(contrast = c(0, -1, 1), label = "T1D vs Aab+"))
      qlf <- if (!is.null(spec$contrast)) {
        shiny::validate(shiny::need(length(spec$contrast) == ncol(design),
                      "Contrast does not match the design (a donor group may be missing)."))
        edgeR::glmQLFTest(fit, contrast = spec$contrast)
      } else {
        shiny::validate(shiny::need(spec$coef <= ncol(design),
                      "Contrast unavailable — a donor group is missing from this subset."))
        edgeR::glmQLFTest(fit, coef = spec$coef)
      }
      res <- edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table
      res$cell_type <- rownames(res)
      rownames(res) <- NULL
      list(res = res[order(res$FDR, res$PValue), , drop = FALSE],
           label = spec$label,
           dispersion = y$common.dispersion,
           n_tested = nrow(y), n_dropped = sum(!keep),
           dropped = rownames(ab$counts)[!keep],
           n_cells = ab$n_cells, group = grp)
    })

    output$edger_note <- renderUI({
      e <- edger_fit()
      tags$p(style = "font-size: 12px; margin: 4px 0 10px 0; line-height: 1.5",
        tags$strong(e$label), " — ",
        sprintf("%s cells from %d donors; %d cell type%s tested",
                formatC(e$n_cells, format = "d", big.mark = ","),
                length(e$group), e$n_tested, if (e$n_tested == 1) "" else "s"),
        if (e$n_dropped > 0) sprintf(", %d dropped by filterByExpr (%s)",
                                     e$n_dropped, paste(e$dropped, collapse = ", ")) else "",
        sprintf(". Common dispersion %.3f (BCV %.2f).", e$dispersion, sqrt(e$dispersion)),
        tags$br(),
        tags$span(
          "logFC is on a log2 scale after library-size normalisation, so it is a change in ",
          "abundance rather than in share of the tissue. FDR is Benjamini-Hochberg across ",
          "the cell types tested in this contrast."))
    })

    output$edger_table <- renderTable({
      e <- edger_fit()
      r <- e$res
      data.frame(
        `Cell type` = r$cell_type,
        `log2FC` = sprintf("%.2f", r$logFC),
        `logCPM` = sprintf("%.1f", r$logCPM),
        `p` = format_p_value(r$PValue),
        `FDR` = format_p_value(r$FDR),
        Sig = ifelse(r$FDR < 0.05, "*", ""),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lccccc")

    output$edger_volcano <- renderPlotly({
      e <- edger_fit()
      r <- e$res
      r$neglog10 <- -log10(pmax(r$FDR, .Machine$double.xmin))
      r$sig <- ifelse(r$FDR < 0.05,
                      ifelse(r$logFC > 0, "Up (FDR<0.05)", "Down (FDR<0.05)"), "NS")
      r$hover <- sprintf("%s<br>log2FC: %.2f<br>p: %s<br>FDR: %s",
                         r$cell_type, r$logFC,
                         format_p_value(r$PValue), format_p_value(r$FDR))
      g <- ggplot(r, aes(x = logFC, y = neglog10, colour = sig, text = hover)) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#999") +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "#999") +
        geom_point(size = 3, alpha = 0.85) +
        scale_colour_manual(values = c("Up (FDR<0.05)" = "#d62728",
                                       "Down (FDR<0.05)" = "#1f77b4",
                                       "NS" = "#bbbbbb"), name = NULL) +
        labs(x = paste0("log2 fold change — ", e$label),
             y = "-log10(FDR)") +
        theme_islet() +
        theme(legend.position = "bottom")
      # No in-plot point labels: they collide with each other and with the y-axis
      # ("3Beta cell") at every width tried, and the results table sits directly
      # beside this plot listing every cell type ordered by FDR. Identity comes
      # from the hover instead.

      fonts <- plotly_islet_fonts()
      ggplotly(g, tooltip = "text") %>%
        layout(font = fonts$global,
               # b covers the x-axis title AND the horizontal legend below it;
               # measured, not guessed (legend sat 16 px off-canvas at b = 90).
               margin = list(t = 20, b = 115, l = 70, r = 20),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               # y is a fraction of the PLOT AREA, which shrinks as the bottom
               # margin grows -- so enlarging the margin alone never fixed this.
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.14,
                             font = fonts$legend_text)) %>%
        plotly_islet_config("inflammation_edger_volcano")
    })

    output$dl_edger <- downloadHandler(
      filename = function() paste0("inflammation_edger_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        e <- tryCatch(edger_fit(), error = function(err) NULL)
        if (is.null(e)) return(write_islet_csv(file, NULL, "edgeR results unavailable"))
        write_islet_csv(file, e$res, islet_csv_header(
          paste0("edgeR differential abundance — ", e$label),
          "Negative-binomial GLM, quasi-likelihood F-test, library-size normalised",
          paste("Cells included:", e$n_cells, "| region:", input$edger_region %||% "core_peri"),
          paste("Cell types tested:", e$n_tested,
                "| dropped by filterByExpr:", e$n_dropped),
          sprintf("Common dispersion: %.4f (BCV %.3f)", e$dispersion, sqrt(e$dispersion)),
          "FDR = Benjamini-Hochberg across the cell types in this contrast"))
      }
    )

    # =========================================================================
    # Card C: infiltration score
    # =========================================================================

    output$infil_pheno_ui <- renderUI({
      df <- tryCatch(cells(), error = function(e) NULL)
      avail <- if (is.null(df)) character(0) else sort(unique(as.character(df$phenotype)))
      sel <- intersect(INFLAM_DEFAULT_PHENOTYPES, avail)
      if (!length(sel)) sel <- head(avail, 5)
      checkboxGroupInput(ns("infil_phenos"), "Phenotypes",
                         choices = avail, selected = sel, inline = TRUE)
    })

    # Steenbuck et al.: score 1 inside the islet, 1/d outside, normalised by the
    # total cell count of the image.
    #
    # THE DENOMINATOR IS THE ADAPTATION THAT MATTERS. Their images are 75 ROIs
    # per donor chosen as "islets with proximal exocrine tissue", so their
    # denominator is islet-local by construction. Ours are whole sections, and
    # how much of a section lies near an islet varies 4.1x to 23.2x between our
    # donors -- purely technical variance that lands straight in the denominator.
    # It is not confounded with disease stage (Kruskal-Wallis p = 0.53 on that
    # dilution factor), so it does not bias the direction; it inflates variance
    # and costs power. Measured cost of using the whole section: T cell KW p
    # 0.1845 -> 0.0081, Immune 0.1119 -> 0.0207, B cell 0.0652 -> 0.0155,
    # CD8a Tcell 0.0271 -> 0.0108. Three phenotypes cross from null to
    # significant purely on the denominator.
    #
    # Default is therefore the islet-proximal population (the ROI analogue).
    # Whole-section remains selectable for comparison, clearly labelled.
    #
    # The NUMERATOR needs no adaptation: 87.8% of the CD8a score mass already
    # comes from within 20 um of an islet even though only 16% of CD8a cells are
    # there, so the 1/d weight is doing its job. Capping it changes nothing.
    infiltration <- reactive({
      df <- cells()
      mode <- input$infil_denom %||% "isletperi"
      d <- suppressWarnings(as.numeric(df$signed_dist))
      # 1/d explodes as d -> 0; floor at 1 µm so a cell one micron outside the
      # edge cannot outweigh a cell inside it.
      score <- ifelse(!is.finite(d), NA_real_,
                      ifelse(d <= 0, 1, 1 / pmax(d, 1)))
      keep <- is.finite(score)
      in_denom <- switch(mode,
        isletperi = keep & df$cell_region %in% c("core", "peri"),
        near100   = keep & is.finite(d) & d <= 100,
        keep)
      # NUMERATOR SCOPE. Steenbuck et al. sum over the SAME population they
      # divide by (all cells in the ROI), so their score is a clean per-cell
      # mean. Restricting only the denominator, as we do by default, makes it a
      # ratio of two different populations: 29% of the numerator mass sits
      # outside the islet+peri denominator. That is a units mismatch, not a bias
      # -- but it is not the paper's statistic either, so it is stated rather
      # than assumed.
      #
      # Default is the unmatched form on purpose: matching to islet+peri confines
      # the 1/d gradient to the 20 um peri band, which guts the very weighting
      # the score exists for (a T cell 25 um out would contribute nothing).
      # Matching costs Immune (KW p 0.021 -> 0.075) and T cell (0.008 -> 0.014)
      # while slightly helping CD8a (0.011 -> 0.009).
      in_numer <- if (isTRUE(input$infil_match_numer)) in_denom else keep
      shiny::validate(shiny::need(sum(in_numer) > 0, "No cells in the numerator scope."))
      agg <- stats::aggregate(list(score = score[in_numer]),
                              by = list(case_id = df$case_id[in_numer],
                                        donor_status = df$donor_status[in_numer],
                                        phenotype = as.character(df$phenotype)[in_numer]),
                              FUN = sum)
      shiny::validate(shiny::need(sum(in_denom) > 0, "No cells in the chosen denominator."))
      tot <- stats::aggregate(list(n_cells = rep(1L, sum(in_denom))),
                              by = list(case_id = df$case_id[in_denom]), FUN = sum)
      out <- merge(agg, tot, by = "case_id")
      out$infiltration_score <- out$score / out$n_cells
      out$donor_status <- factor(as.character(out$donor_status), levels = status_levels)
      attr(out, "denom_mode") <- mode
      attr(out, "matched") <- isTRUE(input$infil_match_numer)
      out
    })

    # What this design can detect BEFORE any biology. At 5 donors per stage the
    # rank tests have a hard floor, and BH multiplies it by the family size, so a
    # null result may carry no information at all. Computed from the actual group
    # sizes rather than hard-coded, so it stays true if the cohort changes.
    power_ceiling <- reactive({
      ns_by_grp <- table(droplevels(factor(
        unique(infiltration()[, c("case_id", "donor_status")])$donor_status)))
      ns_by_grp <- as.integer(ns_by_grp[ns_by_grp > 0])
      if (length(ns_by_grp) < 2) return(NULL)
      # Perfect separation: groups occupy consecutive rank blocks.
      grp <- factor(rep(seq_along(ns_by_grp), ns_by_grp))
      best_kw <- tryCatch(stats::kruskal.test(seq_len(sum(ns_by_grp)), grp)$p.value,
                          error = function(e) NA_real_)
      n1 <- ns_by_grp[1]; n2 <- ns_by_grp[min(2, length(ns_by_grp))]
      best_pw <- tryCatch(stats::wilcox.test(seq_len(n1), n1 + seq_len(n2))$p.value,
                          error = function(e) NA_real_)
      m_kw <- length(INFLAM_FDR_FAMILY)
      list(best_kw = best_kw, best_pw = best_pw,
           m_kw = m_kw, best_kw_q = best_kw * m_kw,
           best_pw_q = best_pw * m_kw * 3,
           kw_detectable = is.finite(best_kw) && best_kw * m_kw < 0.05,
           pw_detectable = is.finite(best_pw) && best_pw * m_kw * 3 < 0.05,
           group_sizes = ns_by_grp)
    })

    infil_selected <- reactive({
      sel <- input$infil_phenos
      shiny::validate(shiny::need(length(sel) > 0, "Select at least one phenotype."))
      s <- infiltration()
      s[s$phenotype %in% sel, , drop = FALSE]
    })

    # Donor-level tests: Kruskal-Wallis across the three stages on donor scores,
    # then pairwise Wilcoxon.
    #
    # THE FAMILY IS PRE-SPECIFIED (INFLAM_FDR_FAMILY, 7 immune phenotypes) and
    # fixed in source, for two reasons:
    #
    #  * It cannot be widened or narrowed after seeing a result. Correcting
    #    across whatever the user has ticked would let unticking phenotypes
    #    shrink the correction and manufacture significance.
    #  * Family size is the binding constraint at this cohort size. BH across
    #    21 phenotypes has a floor of q = 0.041 (only a PERFECT separation of all
    #    15 donors could clear 0.05); across these 7 the floor is 0.014, which
    #    leaves real headroom. See power_ceiling().
    #
    # Phenotypes outside the family are still testable and displayable, but they
    # are reported UNCORRECTED and flagged, never folded in retrospectively.
    #
    # PAIRWISE CONTRASTS ARE NOT FDR-CORRECTED, deliberately. With 5 donors per
    # stage the best achievable Wilcoxon p is 0.0079, so BH over 7 phenotypes x
    # 3 contrasts has a floor of q = 0.167 -- no pairwise contrast could ever
    # reach q < 0.05 whatever the biology. Printing such a column would be a
    # guaranteed false negative dressed as a result. They are reported raw and
    # labelled exploratory instead.
    infil_tests_all <- reactive({
      s <- infiltration()
      phenos <- sort(unique(s$phenotype))
      rows <- lapply(phenos, function(ph) {
        sub <- s[s$phenotype == ph, , drop = FALSE]
        g <- split(sub$infiltration_score, droplevels(sub$donor_status))
        med <- vapply(status_levels, function(l)
          if (!is.null(g[[l]]) && length(g[[l]])) stats::median(g[[l]]) else NA_real_, numeric(1))
        # Omnibus is computable only with >= 2 groups present.
        kw <- if (length(g) >= 2) {
          tryCatch(stats::kruskal.test(sub$infiltration_score,
                                       droplevels(sub$donor_status))$p.value,
                   error = function(e) NA_real_)
        } else NA_real_
        pw <- function(a, b) {
          if (!is.finite(kw)) return(NA_real_)          # no omnibus -> no post-hoc
          if (is.null(g[[a]]) || is.null(g[[b]]) || length(g[[a]]) < 2 || length(g[[b]]) < 2)
            return(NA_real_)
          tryCatch(suppressWarnings(stats::wilcox.test(g[[a]], g[[b]], exact = FALSE)$p.value),
                   error = function(e) NA_real_)
        }
        data.frame(phenotype = ph,
                   in_family = ph %in% INFLAM_FDR_FAMILY,
                   ND = med[["ND"]], `Aab+` = med[["Aab+"]], T1D = med[["T1D"]],
                   p_kw = kw,
                   p_nd_aab = pw("ND", "Aab+"),
                   p_nd_t1d = pw("ND", "T1D"),
                   p_aab_t1d = pw("Aab+", "T1D"),
                   check.names = FALSE, stringsAsFactors = FALSE)
      })
      out <- do.call(rbind, rows)
      # BH within the pre-specified family only; everything else stays NA.
      out$q_kw <- NA_real_
      fam <- out$in_family
      if (any(fam)) out$q_kw[fam] <- stats::p.adjust(out$p_kw[fam], method = "BH")
      out[order(!out$in_family, out$p_kw, out$phenotype), , drop = FALSE]
    })

    # Display slice: the ticked phenotypes only. q-values come from the fixed
    # family above and do not change with the selection.
    infil_tests <- reactive({
      sel <- input$infil_phenos
      shiny::validate(shiny::need(length(sel) > 0, "Select at least one phenotype."))
      t <- infil_tests_all()
      t[t$phenotype %in% sel, , drop = FALSE]
    })

    output$infil_table <- renderTable({
      t <- infil_tests()
      fmt <- function(x) ifelse(is.na(x), "—", formatC(x, format = "g", digits = 3))
      data.frame(
        Phenotype = paste0(t$phenotype, ifelse(t$in_family, "", " \u2021")),
        ND = fmt(t$ND), `Aab+` = fmt(t$`Aab+`), T1D = fmt(t$T1D),
        `KW p` = format_p_value(t$p_kw),
        `KW q` = ifelse(is.na(t$q_kw), "n/a \u2021", format_p_value(t$q_kw)),
        # Uncorrected on purpose -- see infil_tests_all(): no pairwise q could
        # ever clear 0.05 at this cohort size, so a q column here would be a
        # guaranteed null presented as a finding.
        `ND vs Aab+ (raw)` = format_p_value(t$p_nd_aab),
        `ND vs T1D (raw)` = format_p_value(t$p_nd_t1d),
        `Aab+ vs T1D (raw)` = format_p_value(t$p_aab_t1d),
        check.names = FALSE, stringsAsFactors = FALSE
      )
    # 9 columns: Phenotype, ND, Aab+, T1D, KW p, KW q, and the 3 raw pairwise.
    # renderTable(align=) needs exactly one character per column and fails at
    # RENDER time, not parse time -- count them when the columns change.
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lcccccccc")

    output$infil_note <- renderUI({
      t_all <- infil_tests_all()
      pc <- power_ceiling()
      fam <- t_all[t_all$in_family, , drop = FALSE]
      n_sig <- sum(!is.na(fam$q_kw) & fam$q_kw < 0.05)
      shown_out <- sum(!infil_tests()$in_family)
      denom_lbl <- switch(input$infil_denom %||% "isletperi",
        isletperi = "islet + peri-islet cells",
        near100   = "cells within 100 µm of an islet",
        "all cells on the section")
      tagList(
        tags$p(style = "font-size: 11px; margin: 8px 0 0 0; line-height: 1.5",
          tags$strong("Testing. "),
          "Kruskal-Wallis across the three stages on donor-level scores. FDR is Benjamini-Hochberg ",
          "within a ", tags$strong("pre-specified family of the ", length(INFLAM_FDR_FAMILY),
                                   " immune phenotypes"),
          " fixed in source \u2014 not across whatever is ticked, and not across all 21 phenotypes. ",
          "Denominator: ", tags$strong(denom_lbl),
          if (isTRUE(input$infil_match_numer))
            " (numerator matched to it \u2014 the paper's per-cell mean). "
          else " (numerator over all cells \u2014 see Methods). ",
          "Pairwise columns are ", tags$strong("uncorrected"),
          " and exploratory; see the power note below for why a corrected pairwise column would be ",
          "meaningless here.",
          if (shown_out > 0) tags$span(
            tags$br(), "\u2021 ", shown_out, " displayed phenotype",
            if (shown_out == 1) " is" else "s are",
            " outside the pre-specified family, so no q is reported for ",
            if (shown_out == 1) "it" else "them", ".") else NULL),
        if (!is.null(pc) && is.finite(pc$best_kw)) {
          tags$div(style = paste0("margin-top: 8px; padding: 8px 10px; border-radius: 4px;",
                                  " font-size: 11px; line-height: 1.5;",
                                  if (n_sig == 0) " background: #fff8e1; border-left: 3px solid #f0ad4e; color: #6a5200;"
                                  else " background: #eef6ee; border-left: 3px solid #2e7d32; color: #1b4620;"),
            tags$strong("What this test can detect. "),
            sprintf("With %s donors per stage the best attainable Kruskal-Wallis p is %s \u2014 that is a ",
                    paste(pc$group_sizes, collapse = "/"), format_p_value(pc$best_kw)),
            "PERFECT separation of every donor. After FDR across the ", pc$m_kw,
            "-phenotype family the best attainable q is ", tags$strong(format_p_value(pc$best_kw_q)),
            if (pc$kw_detectable) ", so a real effect can clear \u03b1 = 0.05, but only a very strong one. "
            else ", so NO result could clear \u03b1 = 0.05 whatever the biology. ",
            tags$br(),
            sprintf("The best attainable pairwise Wilcoxon p is %s; corrected across %d contrasts that floors at q = %s, ",
                    format_p_value(pc$best_pw), pc$m_kw * 3, format_p_value(pc$best_pw_q)),
            tags$strong("which no biology could ever beat"),
            " \u2014 hence uncorrected pairwise values above.",
            if (n_sig == 0) tagList(tags$br(), tags$strong(
              "Nothing in the family survives FDR here. Read that as underpowered, not as no effect."))
            else tagList(tags$br(), sprintf("%d of %d family phenotypes survive FDR.", n_sig, nrow(fam))))
        } else NULL)
    })

    infil_plot_obj <- reactive({
      s <- infil_selected()
      s$phenotype <- factor(s$phenotype, levels = sort(unique(s$phenotype)))
      s$hover <- sprintf("%s<br>%s — donor %s<br>score: %.3g",
                         s$phenotype, s$donor_status, s$case_id, s$infiltration_score)
      pal <- donor_colors()
      g <- ggplot(s, aes(x = phenotype, y = infiltration_score,
                         colour = donor_status, text = hover)) +
        geom_point(position = position_dodge(width = 0.65), size = 2.8, alpha = 0.85) +
        stat_summary(aes(group = donor_status), fun = stats::median, geom = "crossbar",
                     position = position_dodge(width = 0.65), width = 0.5,
                     linewidth = 0.4, show.legend = FALSE) +
        scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
        labs(x = NULL, y = "Infiltration score (per donor)") +
        theme_islet() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 30, hjust = 1))
      if (isTRUE(input$infil_log)) {
        # Scores span orders of magnitude across phenotypes; pseudo_log keeps any
        # exact zeros visible instead of dropping them at -Inf.
        g <- g + scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-6, base = 10))
      }
      g
    })

    output$infil_plot <- renderPlotly({
      fonts <- plotly_islet_fonts()
      ggplotly(infil_plot_obj(), tooltip = "text") %>%
        layout(font = fonts$global,
               margin = list(t = 20, b = 145, l = 80, r = 20),
               xaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               yaxis = list(title = list(font = fonts$axis_title),
                            tickfont = fonts$axis_tick, automargin = TRUE),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.28,
                             font = fonts$legend_text)) %>%
        plotly_islet_config("inflammation_infiltration")
    })

    output$dl_infil <- downloadHandler(
      filename = function() paste0("inflammation_infiltration_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        s <- tryCatch(infil_selected(), error = function(e) NULL)
        t <- tryCatch(infil_tests_all(), error = function(e) NULL)
        if (is.null(s)) return(write_islet_csv(file, NULL, "No phenotypes selected"))
        hdr <- islet_csv_header(
          "Infiltration score (Steenbuck et al. 2025 definition)",
          "score = 1 inside the islet, 1/d outside (d = um from the islet edge, floored at 1)",
          paste("numerator scope:", if (isTRUE(input$infil_match_numer))
                "matched to the denominator (Steenbuck et al. per-cell mean)"
                else "all cells on the section (preserves the full 1/d gradient)"),
          paste("denominator:", switch(input$infil_denom %||% "isletperi",
                isletperi = "islet + peri-islet cells (ROI analogue, default)",
                near100   = "cells within 100 um of an islet",
                "all cells on the section")),
          "source: data/cell_distance.csv",
          paste("FDR family (pre-specified, fixed in source):",
                paste(INFLAM_FDR_FAMILY, collapse = "; ")),
          "in_family = FALSE rows are reported UNCORRECTED (q_kw is NA) and are exploratory",
          "pairwise Wilcoxon columns are UNCORRECTED by design: at 5 donors per stage the best",
          "  attainable pairwise p is 0.0079, so any BH correction floors above alpha = 0.05")
        # Two blocks in one file. Written through a single open connection:
        # write.csv/write.table with append = TRUE on a PATH re-opens the file and
        # emits R's "appending column names to file" warning (the reason
        # write_islet_csv() exists). Same pattern here for the second block.
        con <- base::file(file, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        writeLines(paste("#", hdr), con)
        writeLines("", con)
        writeLines("# Per-donor infiltration scores", con)
        utils::write.csv(s[, c("case_id", "donor_status", "phenotype",
                               "score", "n_cells", "infiltration_score")],
                         con, row.names = FALSE)
        if (!is.null(t)) {
          writeLines("", con)
          writeLines("# Donor-level tests (BH across all phenotypes present)", con)
          utils::write.csv(t, con, row.names = FALSE)
        }
      }
    )

    # =========================================================================
    # Methods
    # =========================================================================

    output$inflam_methods <- renderUI({
      tags$div(style = "padding: 10px; font-size: 13px; line-height: 1.6",
        tags$p(
          tags$strong("Source. "),
          "Both analyses follow Steenbuck N., Damond N., et al., ",
          tags$em("Imaging mass cytometry reveals early β-cell dysfunction and changes in ",
                  "immune signatures during type 1 diabetes progression in human pancreata"),
          ", bioRxiv 2025.03.05.641526 (Bodenmiller and Atkinson labs). Their cohort is ",
          "88 nPOD donors across five stages with 75 imaging ROIs each; ours is 15 donors ",
          "across three stages with whole-section CODEX images. Where that difference ",
          "forces an adaptation it is stated below."),
        tags$p(
          tags$strong("Differential abundance. "),
          "Per-donor cell-type counts are modelled with edgeR: ",
          tags$code("filterByExpr"), " → ", tags$code("calcNormFactors"), " → ",
          tags$code("estimateDisp"), " → ", tags$code("glmQLFit"), " → ",
          tags$code("glmQLFTest"), ", with Benjamini-Hochberg across cell types. ",
          "Quasi-likelihood F-tests are used because they hold their nominal error rate at ",
          "small replicate numbers, which matters at n = 5 donors per stage. ",
          tags$strong("Counts rather than proportions is the substantive choice: "),
          "proportions must sum to 1, so the loss of β cells in T1D inflates every other ",
          "cell type's proportion whether or not its absolute abundance moved. The ",
          "library-size offset removes that artefact."),
        tags$p(
          tags$strong("Infiltration score. "),
          "Each cell scores 1 if it lies inside an islet and 1/d otherwise, with d the ",
          "distance in µm from the islet edge; per-donor sums are divided by that donor's ",
          "total cell count. The score rewards accumulation in the peri-islet space and ",
          "penalises accumulation in distal tissue, which a plain density cannot ",
          "distinguish. ", tags$strong("Adaptations: "),
          "the paper normalises per imaging ROI and we normalise per donor section (this ",
          "cohort has no ROIs); and 1/d is floored at d = 1 µm so a cell one micron ",
          "outside the edge cannot outweigh a cell inside it."),
        tags$p(
          tags$strong("Multiple comparisons. "),
          "edgeR FDR is BH across the cell types in the displayed contrast. Infiltration ",
          "Kruskal-Wallis p-values are BH within a ", tags$strong("pre-specified family of 7 ",
          "immune phenotypes"), " (", paste(INFLAM_FDR_FAMILY, collapse = ", "),
          "), fixed in source. Phenotypes outside it can be displayed but are reported ",
          "uncorrected and marked \u2021. Neither card is corrected against the other, nor ",
          "across contrasts \u2014 fix a contrast before reading an FDR column."),
        tags$p(
          tags$strong("Why the infiltration pairwise columns are uncorrected. "),
          "At 5 donors per stage the rank tests have a hard floor: the best attainable ",
          "Kruskal-Wallis p is 0.0019 and the best attainable pairwise Wilcoxon p is 0.0079, ",
          "both requiring PERFECT separation of every donor. Benjamini-Hochberg multiplies ",
          "that floor by the family size, so a corrected pairwise column would bottom out at ",
          "q = 0.17 \u2014 ", tags$strong("no biology whatsoever could produce a significant ",
          "value in it"), ". Printing one would be a guaranteed false negative dressed as a ",
          "result, so the pairwise values are reported raw and labelled exploratory. The ",
          "omnibus family is small enough (7) to have genuine headroom: its floor is q = 0.014. ",
          "The card states these limits next to the table so a null reads as ",
          tags$em("underpowered"), ", not as ", tags$em("no effect"), "."),
        tags$p(
          tags$strong("Numerator scope, and where we depart from the paper. "),
          "Steenbuck et al. sum the weights over the same cells they divide by, so their score is a ",
          "per-cell mean over one population. Changing only the denominator makes ours a ratio of ",
          "two populations \u2014 29% of the numerator mass lies outside the islet+peri denominator. ",
          "We default to that unmatched form deliberately: matching to islet+peri would confine the ",
          tags$code("1/d"), " gradient to the 20 µm peri band and gut the weighting the score exists ",
          "for (a T cell 25 µm out would contribute nothing). The ",
          tags$em("Match numerator to denominator"), " box restores the paper's exact statistic; it ",
          "costs Immune (p 0.021 \u2192 0.075) and T cell (0.008 \u2192 0.014) and slightly helps ",
          "CD8a (0.011 \u2192 0.009). Neither form is wrong; they answer slightly different questions, ",
          "so the choice is exposed rather than hidden."),
        tags$p(
          tags$strong("Nearest-islet attribution (checked, not a confounder). "),
          "On a whole slide, ", tags$code("signed_dist"), " is the distance to the ",
          tags$em("nearest"), " islet, so in principle a peri-ductal or perivascular immune ",
          "aggregate could be credited as islet infiltration in a way a single-islet ROI makes ",
          "impossible. Tested: islet crowding and islet count do not differ by disease stage ",
          "(Kruskal-Wallis p = 0.83 and 0.53), and the score correlates weakly with islet geometry ",
          "(|rho| \u2264 0.35 against islet count, \u2264 0.26 against crowding). The correlations ",
          "that exist are negative \u2014 more islets means a larger denominator, not accidental ",
          "credit. The ambiguity is real but does not drive these results."),
        tags$p(
          tags$strong("Why the denominator is islet-proximal. "),
          "Steenbuck et al. image 75 islet-centred ROIs per donor, so their per-image ",
          "denominator is islet-local by construction. Our whole sections are not: the ratio ",
          "of all cells to islet-proximal cells ranges 4.1x to 23.2x across our 15 donors. ",
          "That variation is technical, and it is not confounded with disease stage ",
          "(Kruskal-Wallis p = 0.53), so it does not bias direction \u2014 it inflates variance ",
          "and costs power. Measured: normalising by the whole section instead of the ",
          "islet-proximal population moved T cell from p = 0.18 to 0.008, Immune 0.11 to 0.02, ",
          "B cell 0.065 to 0.016 and CD8a T cell 0.027 to 0.011. Three phenotypes crossed from ",
          "null to significant on the denominator alone, which is why the ROI analogue is the ",
          "default and the whole-section option is labelled as a comparison, not a baseline."),
        tags$p(
          tags$strong("What this tab does not do. "),
          "It does not adjust for age, sex or BMI — those adjustments live in the ",
          "Statistics tab's companion models on a per-islet feature. It does not use ",
          "disease duration, which has no usable variance in this cohort. And it does not ",
          "adjust for HbA1c or C-peptide: both are consequences of β-cell loss rather ",
          "than confounders of it, so conditioning on them would bias the disease effect ",
          "toward the null."))
    })
  })
}
