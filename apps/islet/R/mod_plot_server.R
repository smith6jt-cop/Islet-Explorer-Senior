# ---------- Plot module server ----------
# Exports: plot_server(id, prepared, selected_islet, active_tab)
#
# Dependencies (from utility files already sourced):
#   safe_left_join, add_islet_key, compute_diameter_um  -- utils_safe_join.R
#   bin_islet_sizes, summary_stats                      -- utils_stats.R / data_loading.R
#   islet_spatial_lookup, render_islet_segmentation_plot -- segmentation_helpers.R
#   SF_AVAILABLE, PIXEL_SIZE_UM                         -- global.R
#
# Parameters:
#   id              -- module namespace id
#   prepared        -- reactive returning list(targets_all, markers_all, comp)
#   selected_islet  -- reactiveVal for cross-module click-to-segmentation
#   active_tab      -- reactive returning current tab name (e.g. "Plot")

#' @param cohort the shared cohort filter from cohortFilterServer(), or NULL.
#'   NULL keeps this module's own sidebar inputs, which is what the testthat
#'   suite and any standalone use of the module get; app.R always passes one.
plot_server <- function(id, prepared, selected_islet, active_tab = reactive("Plot"),
                        donor_colors_reactive = reactive(DONOR_COLORS_BRIGHT),
                        remove_outliers = reactive(FALSE),
                        outlier_threshold = reactive(3.0),
                        cohort = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Construct selector vocabularies from server-owned data, then use the
    # same values for UI construction and browser-input validation.
    microenvironment_selection_choices <- function(pd) {
      classes <- pd$targets_all %>%
        dplyr::pull(class) %>% unique() %>% na.omit() %>% sort()
      list(
        "Core" = setNames(paste0(classes, "|islet_core"), classes),
        "Peri-Islet" = setNames(paste0(classes, "|islet_band"), classes),
        "Core + Peri" = setNames(paste0(classes, "|islet_union"), classes)
      )
    }

    composition_selection_choices <- function(pd) {
      base_choices <- c(
        "Ins_frac" = "Ins_any",
        "Glu_frac" = "Glu_any",
        "Stt_frac" = "Stt_any"
      )
      comp_cols <- colnames(pd$comp)
      # The phenotype measures are the prop_ columns whose stem is one of the domain's declared
      # phenotypes. Anything else with a prop_ prefix is a different kind of call -- the galectin
      # and macrophage classes come from a QuPath object classifier, not from the cell typing --
      # and listing the two together would present them as one vocabulary. They get their own
      # group instead, named for where they came from.
      .stem_of <- function(x, prefix) sub(paste0("^", prefix), "", x)
      phenotypes <- domain_phenotype_order()
      all_prop <- grep("^prop_", comp_cols, value = TRUE)
      all_peri <- grep("^peri_prop_", comp_cols, value = TRUE)
      prop_cols <- all_prop[.stem_of(all_prop, "prop_") %in% phenotypes]
      peri_prop_cols <- all_peri[.stem_of(all_peri, "peri_prop_") %in% phenotypes]
      other_prop_cols <- setdiff(all_prop, prop_cols)
      other_peri_cols <- setdiff(all_peri, peri_prop_cols)
      combined_pairs <- pd$combined_pairs
      combined_values <- if (!is.null(combined_pairs) && nrow(combined_pairs) > 0) {
        setNames(paste0("combined_prop_", combined_pairs$stem), combined_pairs$label)
      } else {
        character(0)
      }
      if (!length(prop_cols) && !length(peri_prop_cols)) return(base_choices)

      choices <- list("Core — Hormone Positivity (QuPath threshold)" = base_choices)
      if (length(prop_cols)) {
        choices[[sprintf("Core — Phenotype Proportions (%s)", domain_classifier_label())]] <-
          setNames(prop_cols, gsub("^prop_", "", prop_cols))
      }
      if (length(peri_prop_cols)) {
        choices[[sprintf("Peri-Islet — Phenotype Proportions (%s)", domain_classifier_label())]] <-
          setNames(peri_prop_cols, gsub("^peri_prop_", "", peri_prop_cols))
      }
      if (length(combined_values)) {
        choices[[sprintf("Core+Peri — Phenotype Proportions (%s)", domain_classifier_label())]] <- combined_values
      }
      if (length(other_prop_cols)) {
        choices[["Core — Galectin & macrophage (QuPath object classifier)"]] <-
          setNames(other_prop_cols, gsub("_", " ", .stem_of(other_prop_cols, "prop_")))
      }
      if (length(other_peri_cols)) {
        choices[["Peri-Islet — Galectin & macrophage (QuPath object classifier)"]] <-
          setNames(other_peri_cols, gsub("_", " ", .stem_of(other_peri_cols, "peri_prop_")))
      }
      choices
    }

    # The observer that narrowed this tab's donor-status checkboxes to the
    # statuses actually present is GONE with the checkbox group itself -- the
    # cohort panel owns donor status now and renders its own choices the same
    # way, from the loaded data rather than from a literal. Left as a note
    # because `updateCheckboxGroupInput` on an id that no longer exists is a
    # silent no-op, so dead sync code here would never announce itself.

    selection_values <- function(choices) {
      unname(unlist(choices, use.names = FALSE))
    }

    selection_vocabularies <- reactive({
      pd <- prepared()
      list(
        composition = composition_selection_choices(pd),
        microenvironment = microenvironment_selection_choices(pd)
      )
    })

    validated_mode <- reactive({
      value <- input$mode %||% "Cell Populations"
      if (!islet_valid_scalar_choice(
        value, c("Cell Populations", "Microenvironment"), "plot_mode_invalid"
      )) return("Cell Populations")
      value
    })

    validated_selection <- reactive({
      vocabularies <- selection_vocabularies()
      composition_values <- selection_values(vocabularies$composition)
      microenvironment_values <- selection_values(vocabularies$microenvironment)
      mode <- validated_mode()
      value <- input$which
      allowed <- if (identical(mode, "Microenvironment")) {
        microenvironment_values
      } else {
        composition_values
      }
      stale_allowed <- if (identical(mode, "Microenvironment")) {
        composition_values
      } else {
        microenvironment_values
      }

      valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
        nzchar(value) && value %in% allowed
      if (valid) return(value)

      # Dynamic selectInput replacement briefly retains the preceding mode's
      # valid value. Do not log that normal transition as an attack.
      stale <- is.character(value) && length(value) == 1L && !is.na(value) &&
        value %in% stale_allowed
      if (!is.null(value) && !stale) {
        islet_log_validation_failure("plot_selection_invalid")
      }
      if (length(allowed)) allowed[[1]] else NULL
    })

    validated_metric <- reactive({
      if (identical(validated_mode(), "Microenvironment")) return("Percentage")
      value <- input$metric %||% "Percentage"
      if (!islet_valid_scalar_choice(
        value, c("Percentage", "Counts", "Density"), "plot_metric_invalid"
      )) return("Percentage")
      value
    })

    validated_groups <- reactive({
      # The cohort panel owns donor status for every tab when it is mounted.
      if (!is.null(cohort)) return(cohort$groups())
      allowed <- domain_group_levels()
      value <- input$groups
      if (is.null(value) || !length(value)) return(allowed)
      if (!islet_valid_choice_vector(
        value, allowed, "plot_groups_invalid"
      )) return(allowed)
      value
    })

    validated_min_cells <- reactive({
      if (!is.null(cohort)) return(cohort$spec()$min_cells %||% 1)
      value <- input$min_cells %||% 1
      if (!islet_valid_bounded_number(
        value, 1, 200, "plot_min_cells_invalid"
      )) return(1)
      as.numeric(value)
    })

    # ---- Donor palette helper (consistent colors across scatter & distribution) ----
    get_donor_palette <- function(ids) {
      if (is.null(ids)) return(character(0))
      ids_chr <- sort(unique(as.character(ids)))
      n <- length(ids_chr)
      if (n == 0) return(character(0))
      if (n <= 12) {
        cols <- RColorBrewer::brewer.pal(min(12, max(3, n)), "Paired")
      } else {
        cols <- rainbow(n, s = 0.8, v = 0.8)
      }
      cols <- cols[seq_len(n)]
      names(cols) <- ids_chr
      cols
    }

    # ---- Pseudo-log break points (includes 0 + powers of 10) ----
    # Filters out breaks that are too close in pseudo-log transformed space
    # to prevent label overlap (e.g., 0 and 1 collide on large-range axes).
    pseudo_log_breaks <- function(base = 10) {
      function(limits) {
        max_val <- max(limits)
        if (max_val <= 0) return(0)
        max_pow <- ceiling(log(max_val, base))
        min_pow <- if (max_val < 1) floor(log(max_val, base)) else 0
        candidates <- sort(unique(c(0, base^seq(min_pow, max_pow))))
        candidates <- candidates[candidates <= max_val * 1.1]

        # Greedy filter: drop breaks whose transformed position is too close
        trans_pos <- asinh(candidates / 2) / log(base)
        total_range <- diff(range(trans_pos))
        if (total_range == 0) return(candidates)
        min_gap <- max(0.3, total_range * 0.15)

        keep <- logical(length(candidates))
        keep[1] <- TRUE
        last_kept <- trans_pos[1]
        for (i in seq_along(candidates)[-1]) {
          if (trans_pos[i] - last_kept >= min_gap) {
            keep[i] <- TRUE
            last_kept <- trans_pos[i]
          }
        }
        candidates[keep]
      }
    }

    # ---- Reactive storage for outlier table (replaces global <<- assignment) ----
    plot_outliers <- reactiveVal(NULL)

    # ===========================================================================
    # Dynamic UI outputs
    # ===========================================================================

    output$dynamic_selector <- renderUI({
      vocabularies <- selection_vocabularies()
      selected <- validated_selection()
      if (identical(validated_mode(), "Microenvironment")) {
        selectInput(
          ns("which"), "Structure & Region",
          choices = vocabularies$microenvironment,
          selected = selected
        )
      } else {
        choices <- vocabularies$composition
        pheno_label <- extract_phenotype_from_column(selected)
        hint <- if (!is.null(pheno_label) && !is.null(get_phenotype_rule(pheno_label))) {
          phenotype_hint_ui(pheno_label, show_modal_link = ns("show_phenotype_rules_inline"))
        }
        tagList(
          selectInput(ns("which"), "Population Measure", choices = choices,
                      selected = selected),
          hint,
          helpText(
            "Core = cells within the islet boundary. Peri-Islet = cells within a 20 µm zone around the islet. ",
            "Core+Peri = phenotype cells pooled across both zones, with the matching pooled denominator.",
            if (nzchar(domain_denominator_note())) tags$span(" ", domain_denominator_note())
          )
        )
      }
    })

    output$metric_selector <- renderUI({
      if (identical(validated_mode(), "Microenvironment")) {
        # Microenvironment exposes a single metric (% of region area
        # occupied by the structure). No radio shown — Nerve / Capillary /
        # Lymphatic are area-valued, so the area fraction is the only
        # measurement that makes biological sense.
        return(NULL)
      } else {
        # Composition — Density uses QuPath-measured region areas (core_region_um2,
        # peri_region_um2) for all three region variants. Peri-islet area is the
        # QuPath `_exp20um` annotation area, merged in via prep_data().
        w <- validated_selection() %||% ""
        if (startsWith(w, "combined_prop_")) {
          tagList(
            radioButtons(ns("metric"), "Metric",
                         choices = c("% of all (core+peri) cells"      = "Percentage",
                                     "Cells per islet (core+peri)"     = "Counts",
                                     "Cells per mm² of core+peri area" = "Density"),
                         selected = "Percentage", inline = TRUE),
            helpText(style = "font-size: 12px; margin-top: -4px",
                     "Per islet = phenotype cells in core + 20 µm peri-islet zone. ",
                     "Per mm² of core+peri area = phenotype cells ÷ (core_region_um2 + peri_region_um2), ",
                     "using QuPath-measured annotation areas.")
          )
        } else if (startsWith(w, "peri_prop_")) {
          tagList(
            radioButtons(ns("metric"), "Metric",
                         choices = c("% of peri cells"          = "Percentage",
                                     "Cells per islet (peri)"   = "Counts",
                                     "Cells per mm² of peri area" = "Density"),
                         selected = "Percentage", inline = TRUE),
            helpText(style = "font-size: 12px; margin-top: -4px",
                     "Per islet = phenotype cells tallied in the 20 µm peri-islet zone. ",
                     "Per mm² of peri area = peri cells ÷ peri_region_um2 (QuPath-measured ",
                     "`_exp20um` annotation area).")
          )
        } else {
          tagList(
            radioButtons(ns("metric"), "Metric",
                         choices = c("% of core cells"            = "Percentage",
                                     "Cells per islet (core)"     = "Counts",
                                     "Cells per mm² of core area" = "Density"),
                         selected = "Percentage", inline = TRUE),
            helpText(style = "font-size: 12px; margin-top: -4px",
                     "Per islet = phenotype cells tallied in the islet core. ",
                     "Per mm² of core area = core cells ÷ core_region_um2 ",
                     "(QuPath-measured islet annotation area).")
          )
        }
      }
    })

    output$aab_warning <- renderUI({
      df <- raw_df_base()
      if (!is.null(df) && nrow(df) > 0) {
        aabp_donors <- df %>% dplyr::filter(donor_status == "Aab+")
        if (nrow(aabp_donors) == 0 && "Aab+" %in% validated_groups()) {
          div(style = "font-weight: bold; margin-bottom: 10px",
              "\u26a0 No matching donors.")
        } else {
          NULL
        }
      } else {
        NULL
      }
    })

    # ---- Age & Gender filter UIs (visible only when demographics available) ----
    # `age_filter_ui` and `gender_filter_ui` are GONE, not orphaned: the cohort
    # panel renders both, from the same data-derived bounds this code used.
    # tests/manual/mount_audit.R is what caught them still being defined here
    # after their UI slots moved -- an unmounted renderUI never errors, it just
    # never renders, which is exactly the failure that audit exists to find.

    # ===========================================================================
    # Reactive data chain:  raw_df_base -> raw_df -> plot_df -> summary_df
    # ===========================================================================

    # ---- raw_df_base: per-islet dataset aligned with current selections ----
    raw_df_base <- reactive({
      pd <- prepared()
      mode <- validated_mode()
      groups <- validated_groups()
      w <- validated_selection()

      if (identical(mode, "Microenvironment")) {
        req(!is.null(w), nzchar(w))
        # Parse "className|islet_region" encoding
        parts <- strsplit(w, "\\|")[[1]]
        target_class <- parts[1]
        region_tag <- if (length(parts) >= 2) parts[2] else "islet_core"

        df <- pd$targets_all %>%
          dplyr::filter(`Donor Status` %in% groups,
                        class == target_class,
                        tolower(type) == region_tag) %>%
          dplyr::filter(is.finite(islet_diam_um))
        # Microenvironment value = % of region area occupied by the
        # structure (annotation area within the bounds of core / peri /
        # core+peri). Counts are not exposed — Nerve / Capillary / Lymphatic
        # are area-valued (a single nerve may segment into several pieces),
        # so object counts depend on segmentation morphology, not biology.
        df <- df %>%
          dplyr::mutate(value = {
            reg <- suppressWarnings(as.numeric(region_um2))
            a   <- suppressWarnings(as.numeric(area_um2))
            ifelse(is.finite(reg) & reg > 0 & is.finite(a), 100 * a / reg, NA_real_)
          })
      } else {
        # Composition
        metric <- validated_metric()
        df <- pd$comp %>%
          dplyr::filter(`Donor Status` %in% groups) %>%
          dplyr::filter(is.finite(islet_diam_um))
        # QuPath-measured annotation areas (µm²) joined by prep_data().
        # Density uses these directly rather than the π·(d/2)² estimate so
        # core / peri / core+peri are all measured on the same scale.
        core_area_um2 <- if ("core_region_um2" %in% colnames(df)) {
          suppressWarnings(as.numeric(df$core_region_um2))
        } else NA_real_
        peri_area_um2 <- if ("peri_region_um2" %in% colnames(df)) {
          suppressWarnings(as.numeric(df$peri_region_um2))
        } else NA_real_
        if (startsWith(w, "combined_prop_")) {
          stem <- sub("^combined_prop_", "", w)
          # Look up the original prop_* column via the pre-built pair table —
          # column naming diverges (prop_* keeps spaces/+, peri_prop_* sanitizes).
          pair <- if (!is.null(pd$combined_pairs)) {
            pd$combined_pairs[pd$combined_pairs$stem == stem, , drop = FALSE]
          } else NULL
          core_col <- if (!is.null(pair) && nrow(pair) >= 1) pair$core_col[1] else paste0("prop_", stem)
          peri_col <- if (!is.null(pair) && nrow(pair) >= 1) pair$peri_col[1] else paste0("peri_prop_", stem)
          core_prop <- if (core_col %in% colnames(df)) suppressWarnings(as.numeric(df[[core_col]])) else NA_real_
          peri_prop <- if (peri_col %in% colnames(df)) suppressWarnings(as.numeric(df[[peri_col]])) else NA_real_
          core_total <- suppressWarnings(as.numeric(df$total_cells_core))
          peri_total <- if ("total_cells_peri" %in% colnames(df)) suppressWarnings(as.numeric(df$total_cells_peri)) else NA_real_
          # Two absences meet here and they are not the same thing. A side with NO CELLS --
          # a total of 0, or a total AND a proportion both null -- has nothing there: an islet
          # with no peri-islet expansion still contributes its core, and vice versa, so that
          # side contributes 0 and the other stands alone; that is what this block has always
          # done. A side WITH cells and a null proportion, or a side with a MEASURED proportion
          # and a null total (blank or unparseable after as.numeric), is a missing measurement
          # (build_galectin_layer.py writes NaN for a blank count or a classifier disagreement,
          # and every galectin prop_/peri_prop_ pair reaches this branch through
          # combined_pairs), and the islet's Core+Peri value is NA -- unavailable, never 0 over
          # a kept denominator and never the other side alone. Until 2026-09-21 the first of
          # those was folded into "no cells" and this tab drew the islet at exactly 0; until
          # 2026-09-22 the second was, and the tab drew the islet at its other side's value.
          # The columns are refused where they are read, once, in ie_combine_core_peri() in
          # ieCore (a negative or infinite total is not a count and reads as null; a proportion
          # outside [0, 1] is not a proportion and the islet is NA with a reason), so nothing
          # below tests for those again; an NA n carries through Counts, Density and
          # Percentage alike, and raw_df_base()'s filter(!is.na(value)) then leaves the islet
          # off the plot rather than at 0.
          combined <- ie_combine_core_peri(core_prop, peri_prop, core_total, peri_total)
          combined_n <- combined$n
          combined_total <- combined$total
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = round(combined_n))
          } else if (metric == "Density") {
            ca <- ifelse(is.finite(core_area_um2), core_area_um2, 0)
            pa <- ifelse(is.finite(peri_area_um2), peri_area_um2, 0)
            combined_area <- ca + pa
            df <- df %>% dplyr::mutate(value = ifelse(combined_area > 0,
              combined_n / combined_area * 1e6, NA_real_))
          } else {
            df <- df %>% dplyr::mutate(value = ifelse(combined_total > 0,
              100.0 * combined_n / combined_total, NA_real_))
          }
        } else if (startsWith(w, "peri_prop_")) {
          raw_val <- suppressWarnings(as.numeric(df[[w]]))
          # peri_prop_ are relative to peri-islet cells, use total_cells_peri
          peri_total <- if ("total_cells_peri" %in% colnames(df)) {
            suppressWarnings(as.numeric(df$total_cells_peri))
          } else {
            suppressWarnings(as.numeric(df$cells_total))
          }
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = round(raw_val * peri_total))
          } else if (metric == "Density") {
            df <- df %>% dplyr::mutate(value = ifelse(is.finite(peri_area_um2) & peri_area_um2 > 0,
              raw_val * peri_total / peri_area_um2 * 1e6, NA_real_))
          } else {
            df <- df %>% dplyr::mutate(value = raw_val * 100.0)
          }
        } else if (startsWith(w, "prop_")) {
          raw_val <- suppressWarnings(as.numeric(df[[w]]))
          core_total <- suppressWarnings(as.numeric(df$total_cells_core))
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = round(raw_val * core_total))
          } else if (metric == "Density") {
            df <- df %>% dplyr::mutate(value = ifelse(is.finite(core_area_um2) & core_area_um2 > 0,
              raw_val * core_total / core_area_um2 * 1e6, NA_real_))
          } else {
            df <- df %>% dplyr::mutate(value = raw_val * 100.0)
          }
        } else if (w %in% c("immune_frac_peri", "immune_frac_core", "immune_ratio",
                             "cd8_to_macro_ratio", "tcell_density_peri")) {
          # Immune metrics — only immune_frac_peri/core are true fractions (×100 for %)
          raw_val <- suppressWarnings(as.numeric(df[[w]]))
          if (w %in% c("immune_frac_peri", "immune_frac_core") && metric == "Percentage") {
            df <- df %>% dplyr::mutate(value = raw_val * 100.0)
          } else {
            df <- df %>% dplyr::mutate(value = raw_val)
          }
        } else {
          # Hormone fractions (Ins_any, Glu_any, Stt_any) — QuPath threshold counts.
          # Must use cells_total (same QuPath/groovy pipeline), NOT total_cells_core
          # (DeepCell single-cell pipeline) — mismatched denominators can yield >100%.
          num <- suppressWarnings(as.numeric(df[[w]]))
          den <- suppressWarnings(as.numeric(df$cells_total))
          if (metric == "Counts") {
            df <- df %>% dplyr::mutate(value = num)
          } else if (metric == "Density") {
            df <- df %>% dplyr::mutate(value = ifelse(is.finite(core_area_um2) & core_area_um2 > 0,
              num / core_area_um2 * 1e6, NA_real_))
          } else {
            df <- df %>%
              dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0,
                                           100.0 * num / den, NA_real_))
          }
        }
      }

      out <- df %>%
        dplyr::mutate(donor_status = `Donor Status`) %>%
        dplyr::filter(!is.na(value))

      # ONE call where five hand-rolled filters stood: min-cells, diameter, age,
      # sex and the autoantibody rule. cohort_apply() carries the same semantics
      # -- including "an unchecked AAb flag EXCLUDES donors carrying it, and no
      # flags at all excludes every Aab+ donor", which was ported from this very
      # block -- and it attaches an audit saying which filter bound to which
      # column, so a filter this frame cannot serve surfaces instead of
      # silently doing nothing. Donor status is applied upstream, where the two
      # mode branches each build `df`.
      if (!is.null(cohort)) {
        out <- cohort$apply(out, cells_expr = "core_peri")
      } else {
        min_cells <- validated_min_cells()
        if (min_cells > 1 &&
            "total_cells_core" %in% colnames(out) &&
            "total_cells_peri" %in% colnames(out)) {
          tc_core <- suppressWarnings(as.numeric(out$total_cells_core))
          tc_peri <- suppressWarnings(as.numeric(out$total_cells_peri))
          tc_core <- ifelse(is.finite(tc_core), tc_core, 0)
          tc_peri <- ifelse(is.finite(tc_peri), tc_peri, 0)
          out <- out[is.finite(tc_core + tc_peri) & (tc_core + tc_peri) >= min_cells, , drop = FALSE]
        }
        if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
          diam_min <- as.numeric(input$diam_range[1])
          diam_max <- as.numeric(input$diam_range[2])
          if (is.finite(diam_min) && is.finite(diam_max)) {
            out <- out %>%
              dplyr::filter(is.finite(islet_diam_um) &
                            islet_diam_um >= diam_min &
                            islet_diam_um <= diam_max)
          }
        }
        if (!is.null(input$age_range) && length(input$age_range) == 2 && "age" %in% colnames(out)) {
          out <- out[is.finite(as.numeric(out$age)) &
                     as.numeric(out$age) >= input$age_range[1] &
                     as.numeric(out$age) <= input$age_range[2], , drop = FALSE]
        }
        if (!is.null(input$gender_filter) && "gender" %in% colnames(out)) {
          out <- out[as.character(out$gender) %in% input$gender_filter, , drop = FALSE]
        }
        flags <- input$aab_flags
        all_aab_cols <- COHORT_AAB_COLS
        others <- out %>% dplyr::filter(donor_status != "Aab+")
        aabp   <- out %>% dplyr::filter(donor_status == "Aab+")
        if (nrow(aabp) > 0) {
          if (is.null(flags) || length(flags) == 0) {
            aabp <- aabp[0, , drop = FALSE]
          } else {
            unchecked <- setdiff(all_aab_cols, flags)
            unchecked_avail <- intersect(unchecked, colnames(aabp))
            if (length(unchecked_avail) > 0) {
              exclude_mat <- as.data.frame(aabp[, unchecked_avail, drop = FALSE])
              for (cc in colnames(exclude_mat)) exclude_mat[[cc]] <- as.logical(exclude_mat[[cc]])
              has_unchecked <- rowSums(exclude_mat, na.rm = TRUE) > 0
              aabp <- aabp[!has_unchecked, , drop = FALSE]
            }
          }
        }
        out <- dplyr::bind_rows(others, aabp)
      }
      attr(out, "selection_used") <- w
      attr(out, "mode_used") <- mode
      out
    })

    # ---- raw_df: applies exclude_zero_top filter ----
    raw_df <- reactive({
      out <- raw_df_base()
      sel_used <- attr(out, "selection_used")
      mode_used <- attr(out, "mode_used")
      if (isTRUE(input$exclude_zero_top)) {
        out <- out %>% dplyr::filter(value != 0)
      }
      attr(out, "selection_used") <- sel_used
      attr(out, "mode_used") <- mode_used
      out
    })

    # ---- plot_df: bins data and marks outliers ----
    plot_df <- reactive({
      df <- raw_df()
      if (is.null(df) || !nrow(df)) {
        return(df %||% data.frame())
      }
      mode <- attr(df, "mode_used") %||% validated_mode()
      selection_label <- attr(df, "selection_used") %||% validated_selection()
      if (is.null(selection_label) || length(selection_label) == 0 ||
          is.na(selection_label[1]) || !nzchar(selection_label[1])) {
        selection_label <- if (identical(mode, "Cell Populations")) "Ins_any" else ""
      }
      selection_label <- selection_label[1]

      bw <- suppressWarnings(as.numeric(input$binwidth))
      if (length(bw) != 1 || !is.finite(bw) || bw <= 0) bw <- 50

      # Flag outliers using the shared helper. z-score per donor status,
      # threshold from the live slider. The exclusion flag controls whether
      # outliers are dropped from summary_stats() downstream; ghost points
      # always render.
      df <- flag_outliers(df, value_col = "value", group_col = "donor_status",
                          method = "zscore", threshold = outlier_threshold(),
                          enabled = isTRUE(remove_outliers()))
      summary <- summarize_outlier_filter(df)
      outliers <- df[df$is_outlier, , drop = FALSE]

      # Store outliers for the sidebar table
      plot_outliers(
        if (nrow(outliers) > 0) {
          data.frame(
            Mode          = mode,
            Selection     = selection_label,
            Case_ID       = outliers$`Case ID`,
            Islet         = outliers$islet_key,
            Donor_Status  = outliers$donor_status,
            Diameter_um   = round(outliers$islet_diam_um, 1),
            Value         = round(outliers$value, 3),
            Z_Score       = round(outliers$outlier_score, 2),
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      )

      attr(df, "outlier_summary") <- summary

      # Attach diameter bins
      df <- bin_islet_sizes(df, "islet_diam_um", bw)
      df
    })

    # ---- summary_df: mean/SE (or SD/IQR) per bin per group ----
    # Excludes flagged outliers when the global filter is ON; otherwise
    # all rows contribute to the binned summary.
    summary_df <- reactive({
      df <- plot_df()
      if (is.null(df) || !nrow(df)) return(data.frame())
      stat_choice <- input$stat %||% "mean_se"
      if (!stat_choice %in% c("mean_se", "mean_sd", "median_iqr")) {
        stat_choice <- "mean_se"
      }
      if ("is_outlier" %in% colnames(df) && isTRUE(remove_outliers())) {
        df <- df[!df$is_outlier | is.na(df$is_outlier), , drop = FALSE]
      }
      summary_stats(df,
                    group_cols = c("donor_status", "diam_bin", "diam_mid"),
                    value_col  = "value",
                    stat       = stat_choice)
    })

    # ===========================================================================
    # Outputs: scatter plot (plt), distribution (dist), dist_ui, outlier info
    # ===========================================================================

    # ---- Main scatter plot ----
    output$plt <- renderPlotly({
      sm <- summary_df()
      if (is.null(sm) || !nrow(sm)) {
        return(plotly_empty() %>% layout(title = "Loading data..."))
      }

      mode <- validated_mode()
      selection <- validated_selection() %||% ""
      metric <- validated_metric()

      grp_levels <- domain_group_levels()
      sm$donor_status <- factor(sm$donor_status, levels = grp_levels)

      color_by <- input$plot_color_by %||% "donor_status"

      # Summary lines always use donor_status colors
      color_map <- donor_colors_reactive()

      # Y-axis label
      ylab <- if (identical(mode, "Microenvironment")) {
        wsel <- selection
        parts <- strsplit(wsel, "\\|")[[1]]
        class_name <- if (length(parts) >= 1 && nzchar(parts[1])) parts[1] else "Structure"
        region_tag <- if (length(parts) >= 2) parts[2] else "islet_core"
        region_lbl <- switch(region_tag,
                             islet_core  = "core",
                             islet_band  = "peri-islet",
                             islet_union = "core+peri",
                             "region")
        paste0(class_name, " area (% of ", region_lbl, " area)")
      } else {
        wsel <- selection
        region_tag <- if (startsWith(wsel, "combined_prop_")) "core+peri"
                       else if (startsWith(wsel, "peri_prop_")) "peri"
                       else "core"
        switch(metric,
          "Counts"     = paste0("Cells per islet (", region_tag, ")"),
          "Density"    = paste0("Cells per mm\u00b2 (", region_tag, " area)"),
          "Percentage" = if (region_tag == "core+peri") "% of all (core+peri) cells"
                         else if (region_tag == "peri")    "% of peri cells"
                         else                              "% of core cells"
        )
      }

      # Title (with outlier badge appended)
      base_title <- if (identical(mode, "Microenvironment")) {
        # Extract class name from encoded "className|region" value
        wsel <- selection
        class_name <- strsplit(wsel, "\\|")[[1]][1]
        paste0(class_name, " vs Islet Size")
      } else {
        wsel <- selection
        if (startsWith(wsel, "combined_prop_")) {
          nm <- paste0(gsub("^combined_prop_", "", wsel), " core+peri proportion")
        } else if (startsWith(wsel, "peri_prop_")) {
          nm <- paste0(gsub("^peri_prop_", "", wsel), " peri-islet proportion")
        } else if (startsWith(wsel, "prop_")) {
          nm <- paste0(gsub("^prop_", "", wsel), " proportion")
        } else {
          nm <- switch(wsel,
                       Ins_any = "Insulin+ fraction",
                       Glu_any = "Glucagon+ fraction",
                       Stt_any = "Somatostatin+ fraction",
                       wsel)
        }
        paste0(nm, " vs Islet Size")
      }
      outlier_summary <- attr(plot_df(), "outlier_summary")
      title_text <- paste0(base_title, outlier_title_suffix(outlier_summary))

      # Pseudo-log transform: apply manually to avoid ggplotly axis corruption
      log_scale_on <- !is.null(input$log_scale) && input$log_scale
      plog <- function(x) asinh(x / 2) / log(10)
      orig_y_vals <- numeric(0)
      if (log_scale_on) {
        orig_y_vals <- c(sm$y, sm$ymin, sm$ymax)
        sm$y    <- plog(sm$y)
        sm$ymin <- plog(sm$ymin)
        sm$ymax <- plog(sm$ymax)
      }

      # Initialize base plot
      p <- ggplot(sm, aes(x = diam_mid, y = y, color = donor_status, group = donor_status))

      # Individual points layer FIRST (drawn underneath summary lines)
      show_ind_points <- isTRUE(input$show_points)
      donor_colors <- NULL
      donor_id_breaks <- NULL

      if (show_ind_points) {
        raw <- plot_df()
        raw$donor_status <- factor(raw$donor_status, levels = grp_levels)

        # Click key for segmentation viewer
        if ("Case ID" %in% names(raw) && "islet_key" %in% names(raw)) {
          raw$islet_click_key <- paste(raw$`Case ID`, raw$islet_key, sep = "|")
        } else {
          raw$islet_click_key <- NA_character_
        }

        # Transform y-values for pseudo-log if needed
        if (log_scale_on) {
          orig_y_vals <- c(orig_y_vals, raw$value)
          raw$value <- plog(raw$value)
        }

        # Separate normal vs outliers
        raw_normal   <- raw[!raw$is_outlier | is.na(raw$is_outlier), ]
        raw_outliers <- raw[raw$is_outlier & !is.na(raw$is_outlier), ]

        # Jitter width
        if (!is.null(input$log_scale_x) && input$log_scale_x) {
          jw <- 0.15
        } else {
          jw <- max(1, as.numeric(input$binwidth) * 0.35)
        }

        # Vertical jitter
        yr <- suppressWarnings(range(raw$value, na.rm = TRUE))
        ydiff <- if (all(is.finite(yr))) diff(yr) else 0
        is_counts <- identical(metric, "Counts")
        jh <- if (is_counts) {
          mx <- max(0.02 * ydiff, 0.5)
          if (!is.finite(mx) || mx <= 0) 0.5 else mx
        } else {
          mx <- 0.01 * ydiff
          if (!is.finite(mx) || mx < 0) 0 else mx
        }

        # Normal points layer. Jitter is seeded so visual positions are
        # stable across reactive re-renders (toggling zero-exclusion or the
        # outlier filter doesn't shuffle the cloud).
        if (color_by == "donor_id") {
          if ("Case ID" %in% colnames(raw_normal)) {
            raw_normal$donor_id <- factor(raw_normal$`Case ID`)
            donor_colors <- get_donor_palette(levels(raw_normal$donor_id))
            donor_id_breaks <- levels(raw_normal$donor_id)
            p <- p +
              geom_point(data = raw_normal,
                         aes(x = diam_mid, y = value, color = donor_id, key = islet_click_key),
                         position = position_jitter(width = jw, height = jh, seed = 17),
                         size = input$pt_size, alpha = 1 - input$pt_alpha,
                         inherit.aes = FALSE)
          } else {
            p <- p +
              geom_point(data = raw_normal,
                         aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                         position = position_jitter(width = jw, height = jh, seed = 17),
                         size = input$pt_size, alpha = 1 - input$pt_alpha,
                         inherit.aes = FALSE)
          }
        } else {
          p <- p +
            geom_point(data = raw_normal,
                       aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh, seed = 17),
                       size = input$pt_size, alpha = 1 - input$pt_alpha,
                       inherit.aes = FALSE)
        }

        # Ghost outliers: faint grey × markers so excluded points are still
        # visible. Drawn last (on top) so reviewers can see WHERE in the
        # distribution the exclusions land without those points pulling the
        # summary trend lines. Distinct seed so ghost markers don't sit
        # exactly on top of normal points of the same coordinate.
        if (nrow(raw_outliers) > 0) {
          p <- p +
            geom_point(data = raw_outliers,
                       aes(x = diam_mid, y = value, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh, seed = 31),
                       shape = 4,
                       size = input$pt_size * 1.1,
                       stroke = 0.9,
                       colour = "grey45",
                       alpha = 0.55,
                       inherit.aes = FALSE,
                       show.legend = FALSE)
        }
      }

      # Summary lines + error bars ON TOP of individual points
      p <- p +
        geom_line(alpha = ifelse(!is.null(input$add_smooth) && input$add_smooth == "LOESS", 0, 1),
                  linewidth = 1.2) +
        geom_point(shape = 21, aes(fill = donor_status), color = "black", stroke = 0.4, size = 2.5) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0)

      # Color scale (must include donor_id colors when applicable)
      if (!is.null(donor_colors)) {
        combined_colors <- c(color_map, donor_colors)
        p <- p + scale_color_manual(values = combined_colors,
                                     breaks = donor_id_breaks,
                                     name = "Donor ID")
      } else {
        p <- p + scale_color_manual(values = color_map, breaks = grp_levels, drop = FALSE)
      }
      # Fill scale for summary points (shape=21 outline)
      p <- p + scale_fill_manual(values = color_map, guide = "none")

      p <- p +
        labs(x = "Islet diameter (\u00b5m)", y = ylab, color = NULL, title = title_text) +
        theme_islet() +
        theme(legend.position = c(1, 1),
              legend.justification = c(1, 1),
              legend.direction = "vertical")

      # X-axis scale
      xmax <- suppressWarnings(max(sm$diam_mid, na.rm = TRUE))
      if (!is.finite(xmax)) xmax <- 300

      if (!is.null(input$log_scale_x) && input$log_scale_x) {
        max_pow <- ceiling(log2(xmax))
        major_breaks <- 2^(0:max_pow)
        p <- p + scale_x_continuous(trans = "log2", breaks = major_breaks)
      } else {
        major_breaks <- seq(0, ceiling(xmax / 50) * 50, by = 50)
        minor_breaks <- seq(0, ceiling(xmax / 10) * 10, by = 10)
        p <- p + scale_x_continuous(breaks = major_breaks, minor_breaks = minor_breaks)
      }

      # Y-axis log scale (manual pseudo-log — ggplotly can't handle trans objects)
      if (log_scale_on) {
        orig_range <- range(orig_y_vals[is.finite(orig_y_vals)])
        brks <- pseudo_log_breaks(10)(orig_range)
        p <- p + scale_y_continuous(
          breaks = plog(brks),
          labels = scales::label_comma()(brks)
        )
      }

      # Minor grid
      p <- p + theme(panel.grid.minor = element_line(size = 0.2, colour = "#eeeeee"))

      # LOESS smoothing overlay
      if (!is.null(input$add_smooth) && input$add_smooth == "LOESS") {
        p <- p + geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 0.9)
      }

      # Convert to plotly with click source (namespaced for module compatibility)
      gg <- ggplotly(p, source = ns("plot_scatter"))
      # Scrub the "colour: black" leak (from the shape-21 summary point's
      # color = "black" parameter) and any duplicate donor_status lines.
      gg <- strip_plotly_color_leak(gg)
      # Add thin black outlines to all marker traces
      for (i in seq_along(gg$x$data)) {
        tr <- gg$x$data[[i]]
        if (!is.null(tr$mode) && grepl("markers", tr$mode)) {
          gg$x$data[[i]]$marker$line <- list(color = "black", width = 0.5)
        }
      }
      # Clean plotly legend: strip "(ID, N)" -> "ID" from ggplotly trace naming
      for (i in seq_along(gg$x$data)) {
        nm <- gg$x$data[[i]]$name
        if (!is.null(nm)) {
          cleaned <- gsub("^\\((.+),\\s*\\d+\\)$", "\\1", nm)
          gg$x$data[[i]]$name <- cleaned
          gg$x$data[[i]]$legendgroup <- cleaned
        }
      }
      fonts <- plotly_islet_fonts()
      gg <- gg %>%
        layout(
          font = fonts$global,
          title = list(font = fonts$title),
          xaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(orientation = "h", x = 0, y = -0.25,
                        xanchor = "left", yanchor = "top",
                        font = fonts$legend_text,
                        title = list(font = fonts$legend_title)),
          margin = list(b = 100)) %>%
        plotly_islet_config("plot_scatter") %>%
        event_register("plotly_click") %>%
        islet_click_bridge(ns("plot_native_click"))

      gg
    })

    # ---- Distribution controls ----
    output$dist_ui <- renderUI({
      tagList(
        fluidRow(
          column(4,
            sliderInput(ns("dist_pt_size"), "Point size",
                        min = 0.3, max = 4.0, value = 3.0, step = 0.1),
            sliderInput(ns("dist_pt_alpha"), "Point transparency",
                        min = 0, max = 0.95, value = 0.4, step = 0.05)
          ),
          column(4,
            selectInput(ns("dist_color_by"), "Color points by:",
                        choices = c("Donor Status" = "donor_status",
                                    "Donor ID" = "donor_id"),
                        selected = "donor_status"),
            radioButtons(ns("dist_type"), "Plot type",
                         choices = c("Violin", "Box"),
                         selected = "Violin", inline = TRUE)
          ),
          column(4,
            checkboxInput(ns("exclude_zero_dist"), "Exclude zero values", value = FALSE),
            checkboxInput(ns("dist_show_points"), "Show individual points", value = TRUE),
            checkboxInput(ns("dist_log_scale"), "Log scale y-axis", value = FALSE)
          )
        )
      )
    })

    # ---- Distribution violin/box plot ----
    # Note: the right panel reads `raw_df_base()` (not `raw_df()`) so its
    # "Exclude zero values" checkbox is fully independent of the left panel's
    # equivalent \u2014 toggling one no longer affects the other. Outliers are
    # flagged on the dist-specific subset so the helper sees the same
    # distribution it is summarising, and ghost \u00d7 markers render at positions
    # consistent with the current zero-exclusion state.
    output$dist <- renderPlotly({
      rdf <- raw_df_base()
      if (is.null(rdf) || nrow(rdf) == 0) return(NULL)

      req(validated_mode(), validated_selection())

      mode <- validated_mode()
      selection <- validated_selection() %||% ""
      metric <- validated_metric()

      # Force reactivity on distribution-specific controls
      dist_color_by_val <- input$dist_color_by
      dist_show_points_val <- input$dist_show_points

      # Independent zero-exclusion for the Distribution panel.
      if (isTRUE(input$exclude_zero_dist)) {
        rdf <- rdf %>% dplyr::filter(value != 0)
      }

      color_by <- input$dist_color_by %||% "donor_status"

      # Factor ordering
      rdf$donor_status <- factor(rdf$donor_status, levels = domain_group_levels())

      # Create donor_id column before it's needed
      if (color_by == "donor_id" && "Case ID" %in% colnames(rdf)) {
        rdf$donor_id <- factor(rdf$`Case ID`)
      }

      # Flag outliers on the dist-specific subset so badges and ghost markers
      # reflect the zero-exclusion choice. The shared helper appends
      # `is_outlier` + `outlier_score`; ghost rows are split off and drawn as
      # a separate layer.
      rdf <- flag_outliers(rdf, value_col = "value", group_col = "donor_status",
                           method = "zscore", threshold = outlier_threshold(),
                           enabled = isTRUE(remove_outliers()))
      dist_summary <- summarize_outlier_filter(rdf)
      rdf_ghost <- rdf[rdf$is_outlier, , drop = FALSE]
      rdf <- rdf[!rdf$is_outlier | is.na(rdf$is_outlier), , drop = FALSE]
      if (nrow(rdf_ghost) > 0 && !is.null(rdf_ghost$donor_status)) {
        rdf_ghost$donor_status <- factor(
          rdf_ghost$donor_status,
          levels = domain_group_levels()
        )
      }

      # Y-axis label
      ylab <- if (identical(mode, "Microenvironment")) {
        wsel <- selection
        parts <- strsplit(wsel, "\\|")[[1]]
        class_name <- if (length(parts) >= 1 && nzchar(parts[1])) parts[1] else "Structure"
        region_tag <- if (length(parts) >= 2) parts[2] else "islet_core"
        region_lbl <- switch(region_tag,
                             islet_core  = "core",
                             islet_band  = "peri-islet",
                             islet_union = "core+peri",
                             "region")
        paste0(class_name, " area (% of ", region_lbl, " area)")
      } else {
        wsel <- selection
        region_tag <- if (startsWith(wsel, "combined_prop_")) "core+peri"
                       else if (startsWith(wsel, "peri_prop_")) "peri"
                       else "core"
        switch(metric,
          "Counts"     = paste0("Cells per islet (", region_tag, ")"),
          "Density"    = paste0("Cells per mm\u00b2 (", region_tag, " area)"),
          "Percentage" = if (region_tag == "core+peri") "% of all (core+peri) cells"
                         else if (region_tag == "peri")    "% of peri cells"
                         else                              "% of core cells"
        )
      }

      # Base violin or box
      g <- ggplot(rdf, aes(x = donor_status, y = value, fill = donor_status))

      if (!is.null(input$dist_type) && input$dist_type == "Box") {
        g <- g + geom_boxplot(
          outlier.shape  = NA,
          outlier.size   = 0,
          outlier.colour = NA,
          outlier.fill   = NA,
          outlier.alpha  = 0,
          outlier.stroke = 0,
          alpha = 0.8
        )
      } else {
        g <- g + geom_violin(trim = FALSE, alpha = 0.8)
      }

      g <- g +
        scale_fill_manual(values = donor_colors_reactive(),
                          guide = "none")

      # Individual points (kept = non-outliers). Jitter is seeded so the
      # visual positions are stable across reactive re-renders \u2014 toggling
      # other controls won't shuffle the cloud.
      if (isTRUE(input$dist_show_points)) {
        if (color_by == "donor_id" && "donor_id" %in% colnames(rdf)) {
          donor_colors <- get_donor_palette(levels(rdf$donor_id))
          g <- g +
            geom_jitter(aes(x = donor_status, y = value, color = donor_id),
                        position = position_jitter(width = 0.18, seed = 17),
                        alpha = 1 - ifelse(is.null(input$dist_pt_alpha), 0.5, input$dist_pt_alpha),
                        size  = ifelse(is.null(input$dist_pt_size), 1.5, input$dist_pt_size),
                        stroke = 0, inherit.aes = FALSE) +
            scale_color_manual(values = donor_colors,
                               breaks = levels(rdf$donor_id),
                               name = "Donor ID",
                               guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
        } else {
          g <- g +
            geom_jitter(aes(x = donor_status, y = value, color = donor_status),
                        position = position_jitter(width = 0.18, seed = 17),
                        alpha = 1 - ifelse(is.null(input$dist_pt_alpha), 0.75, input$dist_pt_alpha),
                        size  = ifelse(is.null(input$dist_pt_size), 0.7, input$dist_pt_size),
                        stroke = 0, inherit.aes = FALSE) +
            scale_color_manual(values = donor_colors_reactive(),
                               guide = "none")
        }
      }

      # Ghost \u00d7 layer: excluded outliers always rendered on the Distribution
      # panel so users can see WHERE in the violin/box they live. Drawn last
      # so they sit on top of the violin/box body.
      if (nrow(rdf_ghost) > 0) {
        ghost_size <- ifelse(is.null(input$dist_pt_size), 0.7, input$dist_pt_size) * 1.4
        g <- g + geom_jitter(
          data = rdf_ghost,
          aes(x = donor_status, y = value),
          inherit.aes = FALSE,
          position = position_jitter(width = 0.18, seed = 31),
          shape = 4,
          size = ghost_size,
          stroke = 0.9,
          colour = "grey45",
          alpha = 0.6,
          show.legend = FALSE
        )
      }

      title_text <- paste0("Distribution across donor groups",
                           outlier_title_suffix(dist_summary))
      g <- g +
        labs(x = "Donor Status", y = ylab, title = title_text) +
        theme_islet() +
        theme(legend.position = if (color_by == "donor_id") "right" else "none")

      if (!is.null(input$dist_log_scale) && input$dist_log_scale) {
        g <- g + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                                    breaks = pseudo_log_breaks(10),
                                    labels = scales::label_comma())
      }

      p <- ggplotly(g, tooltip = c("x", "y", "colour"))
      p <- strip_plotly_color_leak(p)

      # Suppress plotly's built-in box outlier markers
      if (!is.null(input$dist_type) && input$dist_type == "Box") {
        for (i in seq_along(p$x$data)) {
          if (!is.null(p$x$data[[i]]$type) && p$x$data[[i]]$type == "box") {
            p$x$data[[i]]$marker$opacity      <- 0
            p$x$data[[i]]$marker$size          <- 0
            p$x$data[[i]]$marker$color         <- "rgba(0,0,0,0)"
            p$x$data[[i]]$marker$outliercolor  <- "rgba(0,0,0,0)"
            p$x$data[[i]]$boxpoints            <- FALSE
          }
        }
      }

      # Add thin black outlines to jitter point traces
      for (i in seq_along(p$x$data)) {
        tr <- p$x$data[[i]]
        if (!is.null(tr$mode) && grepl("markers", tr$mode)) {
          p$x$data[[i]]$marker$line <- list(color = "black", width = 0.5)
        }
      }

      # Clean plotly legend: strip "(ID, N)" -> "ID"
      for (i in seq_along(p$x$data)) {
        nm <- p$x$data[[i]]$name
        if (!is.null(nm)) {
          cleaned <- gsub("^\\((.+),\\s*\\d+\\)$", "\\1", nm)
          p$x$data[[i]]$name <- cleaned
          p$x$data[[i]]$legendgroup <- cleaned
        }
      }
      fonts <- plotly_islet_fonts()
      p <- p %>%
        layout(
          showlegend = isTRUE(color_by == "donor_id"),
          font = fonts$global,
          title = list(font = fonts$title),
          xaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          yaxis = list(title = list(font = fonts$axis_title),
                       tickfont = fonts$axis_tick),
          legend = list(font = fonts$legend_text,
                        title = list(font = fonts$legend_title))
        ) %>%
        plotly_islet_config("plot_distribution")
      p
    })

    # ---- Outlier info table in sidebar ----
    output$plot_outlier_info <- renderUI({
      show_table <- isTRUE(input$show_outlier_table)
      outlier_data <- plot_outliers()
      filter_on <- isTRUE(remove_outliers())

      if (!show_table) return(NULL)
      if (is.null(outlier_data) || nrow(outlier_data) == 0) {
        return(tags$div(
          style = "margin-top: 12px; padding: 8px 10px; font-size: 12px",
          "No outliers detected for the current selection."
        ))
      }

      thr <- suppressWarnings(as.numeric(outlier_threshold()))
      if (!is.finite(thr)) thr <- 3
      header <- if (filter_on) {
        sprintf("\u26a0\ufe0f %d Outlier%s Excluded (|z| > %g within donor status)",
                nrow(outlier_data),
                ifelse(nrow(outlier_data) > 1, "s", ""), thr)
      } else {
        sprintf("\u24d8 %d Outlier%s Flagged (filter OFF \u2014 still in summary)",
                nrow(outlier_data),
                ifelse(nrow(outlier_data) > 1, "s", ""))
      }

      tags$div(
        style = "margin-top: 15px; padding: 10px",
        tags$h6(style = "margin-top: 0", header),
        tags$div(style = "max-height: 200px; overflow-y: auto",
          renderTable(outlier_data))
      )
    })

    # ===========================================================================
    # Click handler: scatter plot -> embedded segmentation panel
    # ===========================================================================

    # Parse a plotly "case_id|islet_key" key and open the drill-down via the
    # shared resolver. `key_null_msg` (when supplied) is shown if no key is
    # present — used by the plotly_click path to nudge users to enable points.
    open_islet_from_key <- function(click_key, key_null_msg = NULL) {
      if (is.null(click_key) || length(click_key) == 0 ||
          is.na(click_key[1]) || !nzchar(as.character(click_key[1]))) {
        if (!is.null(key_null_msg)) showNotification(key_null_msg, type = "message", duration = 3)
        return(invisible(FALSE))
      }
      key_parts <- strsplit(as.character(click_key[1]), "\\|")[[1]]
      if (length(key_parts) != 2L ||
          !islet_safe_path_segment(key_parts[1]) ||
          !islet_safe_path_segment(key_parts[2])) {
        islet_log_validation_failure("islet_key_invalid")
        showNotification("Could not parse islet identifier", type = "warning", duration = 3)
        return(invisible(FALSE))
      }
      select_islet_from_key(key_parts[1], key_parts[2], prepared, selected_islet)
    }

    # Legacy plotly_click path (works on Linux/Windows).
    observeEvent(event_data("plotly_click", source = ns("plot_scatter")), {
      click_data <- event_data("plotly_click", source = ns("plot_scatter"))
      if (is.null(click_data)) return()
      open_islet_from_key(click_data$key,
                          key_null_msg = "Click on individual points to view segmentation (enable 'Show individual points')")
    })

    # macOS-proof native-click bridge: islet_click_bridge() pushes the hovered
    # point's key here on a native click, bypassing plotly's click/drag
    # discrimination (which can suppress plotly_click on macOS).
    observeEvent(input$plot_native_click, {
      nc <- input$plot_native_click
      if (is.null(nc) || is.null(nc$key)) return()
      open_islet_from_key(nc$key)
    })

    # ---- Embedded segmentation viewer panel (below Distribution card) ----
    output$segmentation_viewer_panel <- renderUI({
      # Only render when Plot tab is active to avoid duplicate non-namespaced output IDs
      if (active_tab() != "Plot") return(NULL)
      info <- selected_islet()
      if (is.null(info)) return(NULL)

      # Check if single-cell drill-down data exists for this islet
      has_cells <- drilldown_available() &&
        !is.null(load_islet_cells(info$case_id, info$islet_key))

      # Marker choices for color-by dropdown
      marker_choices <- c("phenotype" = "phenotype")
      if (has_cells) {
        cells <- load_islet_cells(info$case_id, info$islet_key)
        marker_cols <- setdiff(colnames(cells),
                               c("X_centroid", "Y_centroid", "phenotype", "cell_region",
                                 "Cell.Area", "Cell Area", "Nucleus.Area", "Nucleus Area",
                                 EXCLUDED_MARKERS))
        if (length(marker_cols) > 0) {
          marker_choices <- c(marker_choices, setNames(marker_cols, marker_cols))
        }
      }

      # Build title: "Islet_N (case_id, status, age, gender)"
      title_detail <- info$case_id
      if (!is.null(info$donor_status)) title_detail <- paste0(title_detail, ", ", info$donor_status)
      if (!is.null(info$donor_age) && !is.na(info$donor_age)) title_detail <- paste0(title_detail, ", age ", round(info$donor_age))
      if (!is.null(info$donor_gender) && !is.na(info$donor_gender) && nzchar(info$donor_gender)) title_detail <- paste0(title_detail, ", ", info$donor_gender)

      div(class = "card", style = "padding: 20px; margin-bottom: 20px",
        # Header: title + close button
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px",
          h3(paste0(info$islet_key, " (", title_detail, ")"),
             style = "margin: 0; font-size: 22px"),
          actionButton(ns("clear_segmentation"), "Close", class = "btn btn-sm btn-outline-secondary")
        ),
        # Main content: controls sidebar (left) + plot (center) + composition (right)
        fluidRow(
          # Left sidebar: all controls
          column(2,
            div(style = "padding: 10px",
              if (has_cells) {
                tagList(
                  h5("Layers", style = "margin-top: 0; font-size: 15px"),
                  checkboxInput("drilldown_show_cells", "Single Cells", value = TRUE),
                  checkboxInput("drilldown_show_peri_boundary", "Peri Boundary", value = TRUE),
                  checkboxInput("drilldown_show_structures", "Structures", value = TRUE),
                  hr(style = "margin: 8px 0;"),
                  conditionalPanel(
                    condition = "input.drilldown_show_cells",
                    selectInput("drilldown_color_by", "Color by",
                                choices = marker_choices, selected = "phenotype"),
                    selectInput("drilldown_palette", "Phenotype Palette",
                                choices = c("Original", "High Contrast", "Colorblind Safe", "Maximum Distinction"),
                                selected = "High Contrast"),
                    checkboxInput("drilldown_show_peri", "Show peri-islet cells", value = TRUE)
                  )
                )
              } else {
                tagList(
                  h5("Legend", style = "margin-top: 0; margin-bottom: 10px; font-size: 15px"),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px",
                    div(style = "width: 24px; height: 3px; background-color: #0066CC; margin-right: 8px;"),
                    span("Islet", style = "font-size: 14px")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px",
                    div(style = "width: 24px; height: 3px; background-color: #00CCCC; margin-right: 8px;"),
                    span("+20\u00b5m", style = "font-size: 14px")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px",
                    div(style = "width: 24px; height: 3px; background-color: #CC00CC; margin-right: 8px;"),
                    span("Nerve", style = "font-size: 14px")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px",
                    div(style = "width: 24px; height: 3px; background-color: #CC0000; margin-right: 8px;"),
                    span("Capillary", style = "font-size: 14px")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px",
                    div(style = "width: 24px; height: 3px; background-color: #00AA00; margin-right: 8px;"),
                    span("Lymphatic", style = "font-size: 14px")
                  )
                )
              },
              hr(style = "margin: 8px 0;"),
              span(paste("Centroid:", round(info$centroid_x, 0), "/", round(info$centroid_y, 0), "\u00b5m"),
                   style = "font-size: 13px")
            )
          ),
          # Center: plot
          column(6,
            plotOutput("islet_segmentation_view", height = "650px"),
            div(style = "text-align: right; margin-top: 6px",
              downloadButton("dl_islet_segmentation_view", "Download PNG",
                             class = "btn-sm btn-outline-primary"))
          ),
          # Right: composition sidebar
          column(4,
            if (has_cells) {
              conditionalPanel(
                condition = "input.drilldown_show_cells",
                div(style = "padding: 15px",
                  h5("Cell Composition", style = "margin-top: 0; font-size: 17px"),
                  plotOutput("islet_drilldown_summary", height = "400px"),
                  div(style = "text-align: right; margin-top: 6px",
                    downloadButton("dl_islet_drilldown_summary", "Download PNG",
                                   class = "btn-sm btn-outline-primary")),
                  hr(),
                  tableOutput("islet_drilldown_table"),
                  uiOutput("islet_drilldown_rules")
                )
              )
            }
          )
        )
      )
    })

    # ---- Clear segmentation viewer ----
    observeEvent(input$clear_segmentation, {
      selected_islet(NULL)
    })

    # ===========================================================================
    # Download handlers
    # ===========================================================================

    output$dl_summary <- downloadHandler(
      filename = function() {
        paste0("summary_", gsub("[^0-9A-Za-z]+", "_", Sys.time()), ".csv")
      },
      content = function(file) {
        df <- summary_df()
        if (is.null(df) || nrow(df) == 0) df <- data.frame()

        descriptions <- get_selection_description()

        # Write descriptions as comments, then the data
        writeLines(paste("#", descriptions), file)
        writeLines("", file, sep = "\n")
        write.table(df, file, append = TRUE, sep = ",", row.names = FALSE)
      }
    )

    # ===========================================================================
    # Selection description helper (also returned for Statistics download)
    # ===========================================================================

    get_selection_description <- function() {
      desc <- c()
      desc <- c(desc, paste("Generated:", Sys.time()))
      desc <- c(desc, paste("Data Layer:", validated_mode()))
      desc <- c(desc, paste("Selection:", validated_selection() %||% "Unknown"))
      desc <- c(desc, paste("Metric:", validated_metric()))

      # COHORT SCOPE COMES FROM THE PANEL, in one block, generated from the spec
      # AND the audit. Hand-writing it here is how a header comes to report a
      # filter that was silently skipped: the old lines below printed the AAb
      # and diameter selections whether or not the frame had the columns to
      # apply them. cohort_provenance() marks a requested-but-unapplied filter
      # "*** NOT APPLIED ***" instead.
      if (!is.null(cohort)) {
        desc <- c(desc, cohort$provenance(attr(raw_df_base(), "cohort_audit")))
      } else {
        desc <- c(desc, paste("Donor Groups:", paste(validated_groups(), collapse = ", ")))
        if (length(input$aab_flags) > 0) {
          excluded <- setdiff(COHORT_AAB_COLS, input$aab_flags)
          desc <- c(desc, if (length(excluded) > 0)
            paste("AAb Filter: Excluding donors with", paste(excluded, collapse = ", "))
            else "AAb Filter: All Aab+ donors included")
        }
        if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
          desc <- c(desc, paste("Diameter Range:",
                                input$diam_range[1], "-", input$diam_range[2], "\u00b5m"))
        }
      }

      desc <- c(desc, paste("Statistic:", input$stat %||% "mean_se"))
      desc <- c(desc, paste("Bin Width:", input$binwidth %||% "50", "\u00b5m"))

      return(desc)
    }

    # ---- Inline "see all rules" link under the Population Measure hint ----
    observeEvent(input$show_phenotype_rules_inline,
                 { showModal(phenotype_rules_modal_ui()) }, ignoreInit = TRUE)

    # ===========================================================================
    # Return values for other modules (Statistics, etc.)
    # ===========================================================================

    # `raw_df_prezero` is `raw_df_base` — the per-islet frame with every sidebar
    # filter applied but BEFORE `input$exclude_zero_top`. The Statistics tab
    # consumes this one so the Plot sidebar's zero checkbox cannot silently
    # alter statistical results; Statistics applies its own `stats_exclude_zero`
    # and reports it truthfully in the CSV header.
    list(
      raw_df                  = raw_df,
      raw_df_prezero          = raw_df_base,
      summary_df              = summary_df,
      get_selection_description = get_selection_description
    )
  })
}
