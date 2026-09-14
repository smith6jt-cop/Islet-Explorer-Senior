# ---------- Data loading and wrangling ----------
# Extracted from app.R
# Dependencies: safe_left_join, add_islet_key, compute_diameter_um (from R/utils_safe_join.R)
#               master_path, h5ad_path (from 00_globals.R)
#
# Data loading priority: H5AD (islet_explorer.h5ad) > Excel (master_results.xlsx)
# H5AD contains groovy-derived data stored in .uns by scripts/build_h5ad_for_app.py

safe_read_sheet <- function(path, sheet) {
  # Increase guess_max to reduce type misguesses on sparse columns
  tryCatch(readxl::read_excel(path, sheet = sheet, guess_max = 100000), error = function(e) NULL)
}

load_master <- function(path = master_path) {
  list(
    markers = safe_read_sheet(path, "Islet_Markers"),
    targets = safe_read_sheet(path, "Islet_Targets"),
    comp    = safe_read_sheet(path, "Islet_Composition"),
    lgals3  = safe_read_sheet(path, "LGALS3")
  )
}

# OPTIMIZATION: Consolidate donor metadata extraction
get_donor_metadata <- function(master) {
  # Priority order: composition > targets > markers
  for (sheet in c("comp", "targets", "markers")) {
    df <- master[[sheet]]
    if (!is.null(df) && nrow(df) > 0) {
      required_cols <- c("Case ID", "Donor Status")

      # Determine which AAb columns exist
      aab_candidates <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
      aab_cols <- intersect(aab_candidates, names(df))

      if (all(required_cols %in% names(df))) {
        return(df %>%
          dplyr::select(dplyr::all_of(c(required_cols, aab_cols))) %>%
          dplyr::distinct())
      }
    }
  }
  NULL
}

prep_data <- function(master) {
  # Determine which autoantibody columns are available
  aab_cols_targets <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$targets))
  aab_cols_markers <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$markers))
  aab_cols_comp    <- intersect(c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"), colnames(master$comp))
  aab_cols_all     <- unique(c(aab_cols_targets, aab_cols_markers, aab_cols_comp))
  # Use new consolidated function
  donors_meta <- get_donor_metadata(master)
  # Islet size proxy per islet corefor diameter
  targets <- master$targets %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  core_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_core") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, core_region_um2 = region_um2) %>%
    dplyr::distinct()
  # Peri-islet (20 µm expansion) annotation area, measured by QuPath alongside core
  peri_area <- targets %>%
    dplyr::filter(tolower(type) == "islet_band") %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, peri_region_um2 = region_um2) %>%
    dplyr::distinct()
  # Diameter is computed ONLY from core area; if no core exists, diameter is NA
  size_area <- core_area %>%
    dplyr::mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)
  # Pre-built bundle of per-islet region areas + diameter for downstream joins
  region_areas <- size_area %>%
    { safe_left_join(., core_area, by = c("Case ID", "Donor Status", "islet_key"), context = "region_areas:core_area") } %>%
    { safe_left_join(., peri_area, by = c("Case ID", "Donor Status", "islet_key"), context = "region_areas:peri_area") }

  # Targets: keep all region types and later filter by user selection
  targets_all <- targets %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_targets), islet_key, type, class, area_um2, region_um2, area_density, count) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        tolower(type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(type)
      )
      # Density is already in um2, no conversion needed
    ) %>%
    { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "targets_all:size_area") }

  # Synthesize missing union rows for counts by summing core + band (do NOT fabricate union area)
  if (nrow(targets_all) > 0) {
    keys <- c("Case ID", "Donor Status", "islet_key", "class")
    # Identify which key combos already have a union row
    have_union <- targets_all %>%
      dplyr::filter(type == "islet_union") %>%
      dplyr::select(dplyr::all_of(keys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_rows <- targets_all %>% dplyr::filter(type == "islet_core") %>%
      dplyr::select(dplyr::all_of(c(keys, "count", "area_um2", "region_um2"))) %>%
      dplyr::rename(count_core = count, area_um2_core = area_um2, region_um2_core = region_um2)
    band_rows <- targets_all %>% dplyr::filter(type == "islet_band") %>%
      dplyr::select(dplyr::all_of(c(keys, "count", "area_um2", "region_um2"))) %>%
      dplyr::rename(count_band = count, area_um2_band = area_um2, region_um2_band = region_um2)
    union_missing <- core_rows %>%
      dplyr::inner_join(band_rows, by = keys) %>%
  { safe_left_join(., have_union, by = keys, context = "targets_union_missing:have_union") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = suppressWarnings(as.numeric(area_um2_core)) + suppressWarnings(as.numeric(area_um2_band)),
                          region_um2 = suppressWarnings(as.numeric(region_um2_core)) + suppressWarnings(as.numeric(region_um2_band)),
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band))) %>%
        dplyr::mutate(area_density = dplyr::if_else(is.finite(region_um2) & region_um2 > 0,
                                                     area_um2 / region_um2, NA_real_))
      # Attach diameter via size_area
  synth <- synth %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "targets_synth:size_area") }
      # Bind synthetic union rows to targets_all
      targets_all <- dplyr::bind_rows(targets_all, synth)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones
  if (!is.null(donors_meta)) {
    targets_all <- targets_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      { safe_left_join(., donors_meta, by = c("Case ID","Donor Status"), context = "targets_all:donors_meta") }
  }

  # Markers with fraction positive / mean intensity (include LGALS3 sheet).
  # Both sheets are optional: a dataset can carry per-islet composition without per-marker
  # positivity, which is the case whenever no lane has fitted a level for every marker. The
  # LGALS3 branch below was already guarded; this one was not, so a workbook missing
  # Islet_Markers took the whole app down instead of leaving the marker views empty.
  if (is.null(master$markers) || nrow(master$markers) == 0) {
    message("[DATA] Islet_Markers absent; marker views will be empty for this dataset")
    # Take the key columns' types from a sheet that does exist, so the downstream joins
    # against donor metadata line up -- an empty frame with the wrong type fails the join
    # just as hard as a NULL fails the filter.
    donor_key <- NULL
    for (sheet in c("comp", "targets")) {
      df <- master[[sheet]]
      if (!is.null(df) && all(c("Case ID", "Donor Status") %in% names(df))) {
        donor_key <- df[0, c("Case ID", "Donor Status"), drop = FALSE]
        break
      }
    }
    if (is.null(donor_key)) {
      donor_key <- tibble::tibble(`Case ID` = character(), `Donor Status` = character())
    }
    markers <- dplyr::bind_cols(
      tibble::as_tibble(donor_key),
      tibble::tibble(islet_key = character(), region_type = character(),
                     marker = character(), n_cells = numeric(),
                     pos_count = numeric(), pos_frac = numeric())
    )
  } else {
    markers <- master$markers %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  }
  markers_all <- markers %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::any_of(aab_cols_markers), islet_key, region_type, marker, n_cells, pos_count, pos_frac) %>%
    dplyr::mutate(
      region_type = dplyr::case_when(
        tolower(region_type) %in% c("islet_core", "core") ~ "islet_core",
        tolower(region_type) %in% c("islet_band", "band", "peri-islet", "peri_islet") ~ "islet_band",
        tolower(region_type) %in% c("islet_union", "union", "islet+20um", "islet_20um") ~ "islet_union",
        TRUE ~ tolower(region_type)
      )
    )
  # Add LGALS3 rows if available
  if (!is.null(master$lgals3) && nrow(master$lgals3) > 0) {
    g3 <- master$lgals3 %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key)) %>%
      mutate(marker = as.character(marker)) %>%
      dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(intersect(colnames(master$lgals3), c("AAb_GADA","AAb_IA2A","AAb_ZnT8A","AAb_IAA","AAb_mIAA"))), islet_key, region_type, marker, n_cells, pos_count, pos_frac)
    if (!is.null(g3) && nrow(g3) > 0) {
      markers_all <- bind_rows(markers_all, g3)
    }
  }
  # Ensure AAb flags are present for all rows, including synthetic ones (markers)
  if (!is.null(donors_meta)) {
    markers_all <- markers_all %>% dplyr::select(-dplyr::any_of(aab_cols_all)) %>%
      { safe_left_join(., donors_meta, by = c("Case ID","Donor Status"), context = "markers_all:donors_meta") }
  }
  markers_all <- markers_all %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_all:size_area") }

  # Synthesize missing union rows for markers by summing core + band counts (do NOT fabricate area)
  if (nrow(markers_all) > 0) {
    mkeys <- c("Case ID", "Donor Status", "islet_key", "marker")
    have_union_m <- markers_all %>%
      dplyr::filter(region_type == "islet_union") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_union = TRUE)

    core_m <- markers_all %>% dplyr::filter(region_type == "islet_core") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_core = n_cells, pos_core = pos_count)
    band_m <- markers_all %>% dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>%
      dplyr::rename(n_band = n_cells, pos_band = pos_count)
    union_m <- markers_all %>% dplyr::filter(region_type == "islet_union") %>% dplyr::select(dplyr::all_of(c(mkeys, "n_cells", "pos_count"))) %>% dplyr::rename(n_union = n_cells, pos_union = pos_count)
    union_missing_m <- core_m %>%
      dplyr::inner_join(band_m, by = mkeys) %>%
  { safe_left_join(., have_union_m, by = mkeys, context = "markers_union_missing:have_union_m") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing_m) > 0) {
      synth_m <- union_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_union", marker,
                          n_cells = suppressWarnings(as.numeric(n_core)) + suppressWarnings(as.numeric(n_band)),
                          pos_count = suppressWarnings(as.numeric(pos_core)) + suppressWarnings(as.numeric(pos_band))) %>%
        dplyr::mutate(pos_frac = ifelse(is.finite(n_cells) & n_cells > 0,
                                        suppressWarnings(as.numeric(pos_count)) / suppressWarnings(as.numeric(n_cells)),
                                        NA_real_)) %>%
  { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_synth_m:size_area") }
      markers_all <- dplyr::bind_rows(markers_all, synth_m)
    }

    # Backfill missing band rows when core and union exist but band is missing: band = union - core (counts only)
    have_band_m <- markers_all %>%
      dplyr::filter(region_type == "islet_band") %>%
      dplyr::select(dplyr::all_of(mkeys)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(.has_band = TRUE)
    band_missing_m <- core_m %>%
      dplyr::inner_join(union_m, by = mkeys) %>%
  { safe_left_join(., have_band_m, by = mkeys, context = "markers_band_missing:have_band_m") } %>%
      dplyr::filter(is.na(.has_band))
    if (nrow(band_missing_m) > 0) {
      synth_band <- band_missing_m %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, region_type = "islet_band", marker,
                          n_cells = suppressWarnings(as.numeric(n_union)) - suppressWarnings(as.numeric(n_core)),
                          pos_count = suppressWarnings(as.numeric(pos_union)) - suppressWarnings(as.numeric(pos_core))) %>%
        dplyr::mutate(
          n_cells = ifelse(is.finite(n_cells), n_cells, 0),
          pos_count = ifelse(is.finite(pos_count), pos_count, 0),
          n_cells = pmax(0, n_cells),
          pos_count = pmax(0, pos_count),
          pos_frac = ifelse(n_cells > 0, pos_count / n_cells, NA_real_)
        ) %>%
  { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "markers_synth_band:size_area") }
      markers_all <- dplyr::bind_rows(markers_all, synth_band)
    }
  }

  # Fill pos_frac when n_cells and pos_count are present but pos_frac is NA; also compute pos_pct (0..100)
  if (nrow(markers_all) > 0) {
    markers_all <- markers_all %>%
      dplyr::mutate(
        .n = suppressWarnings(as.numeric(n_cells)),
        .p = suppressWarnings(as.numeric(pos_count)),
        pos_frac = dplyr::coalesce(pos_frac, ifelse(is.finite(.n) & .n > 0 & is.finite(.p), .p/.n, NA_real_)),
        pos_pct = ifelse(is.finite(pos_frac), 100.0 * pos_frac, NA_real_)
      ) %>%
      dplyr::select(-.n, -.p)
  }

  # Composition by islet
  comp <- master$comp %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  comp <- comp %>% dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_comp), islet_key, cells_total, Ins_single, Glu_single, Stt_single,
                          Multi_Pos, Triple_Neg, Ins_any, Glu_any, Stt_any)
  comp <- comp %>% { safe_left_join(., region_areas, by = c("Case ID", "Donor Status", "islet_key"), context = "comp:region_areas") }

  # Merge phenotype proportions into comp (H5AD only; NULL from Excel path)
  if (!is.null(master$phenotypes) && nrow(master$phenotypes) > 0) {
    comp <- safe_left_join(comp, master$phenotypes,
                           by = c("Case ID", "islet_key"),
                           context = "comp:phenotypes")
    message("[prep-data] Merged ", ncol(master$phenotypes) - 2, " phenotype columns into comp")
  }

  # Merge neighborhood metrics into comp (H5AD only; NULL from Excel path)
  if (!is.null(master$neighborhood) && nrow(master$neighborhood) > 0) {
    comp <- safe_left_join(comp, master$neighborhood,
                           by = c("Case ID", "islet_key"),
                           context = "comp:neighborhood")
    message("[prep-data] Merged ", ncol(master$neighborhood) - 2, " neighborhood columns into comp")
  }

  # Merge donor demographics into all dataframes (H5AD only; NULL from Excel path)
  if (!is.null(master$donor_demographics) && nrow(master$donor_demographics) > 0) {
    targets_all <- safe_left_join(targets_all, master$donor_demographics, by = "Case ID", context = "targets_all:demographics")
    markers_all <- safe_left_join(markers_all, master$donor_demographics, by = "Case ID", context = "markers_all:demographics")
    comp <- safe_left_join(comp, master$donor_demographics, by = "Case ID", context = "comp:demographics")
    message("[prep-data] Merged donor demographics (age, gender)")
  }

  # Combined (core+peri) phenotype pair table. peri_prop_* columns sanitize
  # the phenotype name (spaces -> _, + -> plus) but prop_* columns keep the
  # original label, so a direct stem intersection misses ~half the pairs. Build
  # a lookup that resolves each peri-style safe-name back to its prop_* core
  # column by trying un-sanitization candidates (mirrors the Spatial tab's
  # phenotype_options() resolver). Consumed by mod_plot_server.R Cell
  # Populations selector and raw_df_base() Combined branch.
  # A DOMAIN may supply per-islet cell-type proportions as their own file rather than inside an
  # H5AD. That is the case for any dataset built from cell typing alone: the composition is the
  # result, and there is no trajectory H5AD to carry it. Declared as `data.islet_proportions`;
  # silently absent for the CODEX domain, which gets its prop_* columns from the H5AD .obs.
  #
  # This must run BEFORE combined_pairs below, which derives the core+peri lookup from
  # colnames(comp) and would otherwise see none of these columns.
  # Only when nothing upstream has already supplied them. The DuckDB path merges prop_* from
  # islets.parquet at :259 above; merging again here would suffix every column .x/.y and leave
  # the plot selectors pointing at names that no longer exist.
  already_have_props <- length(grep("^prop_", names(comp))) > 0
  props_path <- if (already_have_props) NULL else
    tryCatch(domain_data_file("islet_proportions"), error = function(e) NULL)
  if (!is.null(props_path) && file.exists(props_path)) {
    props <- tryCatch(readr::read_csv(props_path, show_col_types = FALSE),
                      error = function(e) NULL)
    if (!is.null(props) && nrow(props) > 0) {
      if ("case_id" %in% names(props)) props <- dplyr::rename(props, `Case ID` = case_id)
      drop <- intersect(c("combined_islet_id"), names(props))
      props <- props %>% dplyr::select(-dplyr::any_of(drop))
      # Match the key types comp carries, or the join silently drops every row.
      props$`Case ID` <- as(props$`Case ID`, class(comp$`Case ID`)[1])
      before <- nrow(comp)
      comp <- safe_left_join(comp, props, by = c("Case ID", "islet_key"),
                             context = "comp:islet_proportions")
      matched <- sum(!is.na(comp[[grep("^prop_", names(comp), value = TRUE)[1]]]))
      message("[DATA] islet proportions merged: ", matched, " of ", before, " islets, ",
              length(grep("^prop_|^peri_prop_", names(comp))), " proportion columns")
    }
  }

  # Demographics reach `comp` only through the H5AD path (the two obs reads above), so a dataset
  # built from cell typing alone -- which has no trajectory H5AD -- arrives with no age and no
  # sex, and the Plot tab's filters render empty. Attach them from the domain's own donor table
  # when it has one. The stored column names stay `age`/`gender`; the UI says "Sex".
  if (!all(c("age", "gender") %in% names(comp))) {
    donor_path <- file.path(domain_data_root(), "donor_metadata.csv")
    if (file.exists(donor_path)) {
      dm <- tryCatch(readr::read_csv(donor_path, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(dm) && all(c("case_id", "age", "sex") %in% names(dm))) {
        dm <- dm %>%
          dplyr::transmute(`Case ID` = case_id, age = as.numeric(age),
                           gender = as.character(sex)) %>%
          dplyr::distinct()
        dm$`Case ID` <- as(dm$`Case ID`, class(comp$`Case ID`)[1])
        comp <- safe_left_join(comp, dm, by = "Case ID", context = "comp:donor_metadata")
        message("[DATA] demographics attached for ",
                length(unique(dm$`Case ID`)), " donors")
      }
    }
  }

  combined_pairs <- (function() {
    cn <- colnames(comp)
    peri_cols <- grep("^peri_prop_", cn, value = TRUE)
    if (length(peri_cols) == 0) return(NULL)
    rows <- lapply(peri_cols, function(pc) {
      safe <- sub("^peri_prop_", "", pc)
      candidates <- unique(c(
        gsub("plus", "+", gsub("_", " ", safe), fixed = TRUE),
        gsub("_", " ", safe),
        gsub("plus", "+", safe, fixed = TRUE),
        safe
      ))
      core_col <- NA_character_; label <- safe
      for (cand in candidates) {
        guess <- paste0("prop_", cand)
        if (guess %in% cn) { core_col <- guess; label <- cand; break }
      }
      data.frame(stem = safe, core_col = core_col, peri_col = pc,
                 label = label, stringsAsFactors = FALSE)
    })
    out <- dplyr::bind_rows(rows)
    out <- out[!is.na(out$core_col), , drop = FALSE]
    out <- out[order(out$label), , drop = FALSE]
    rownames(out) <- NULL
    out
  })()
  if (!is.null(combined_pairs) && nrow(combined_pairs) > 0) {
    message("[prep-data] Built combined core+peri lookup for ", nrow(combined_pairs), " phenotypes")
  }

  message("[prep-final] size_area=", paste(class(size_area), collapse='/'),
    " targets_all=", paste(class(targets_all), collapse='/'),
    " markers_all=", paste(class(markers_all), collapse='/'),
    " comp=", paste(class(comp), collapse='/'))
  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all,
       comp = comp, combined_pairs = combined_pairs)
}

# Simple NA audit
audit_na <- function(df, label) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  na_cnt <- vapply(df, function(x) sum(is.na(x)), integer(1))
  total <- nrow(df)
  pct <- ifelse(total > 0, round(100 * na_cnt / total, 2), 0)
  msg <- paste0("[NA audit] ", label, ": ", paste(names(na_cnt), paste0(na_cnt, " (", pct, "%)"), sep = "=", collapse = "; "))
  message(msg)
}

bin_islet_sizes <- function(df, diam_col, width) {
  x <- suppressWarnings(as.numeric(df[[diam_col]]))
  max_x <- max(x, na.rm = TRUE)
  max_x <- ifelse(is.finite(max_x), max_x, 0)
  # Define numeric bins: [lo, hi)
  bin_lo <- floor(x / width) * width
  bin_hi <- bin_lo + width
  diam_mid <- bin_lo + width/2
  # Build a human-readable label
  diam_bin <- paste0("[", bin_lo, ", ", bin_hi, ")")
  df$diam_bin <- factor(diam_bin, levels = unique(diam_bin[order(bin_lo)]), ordered = TRUE)
  df$diam_mid <- diam_mid
  df
}

# ---------- H5AD-based data loading ----------
# Reconstructs the same data frames as load_master() + prep_data() from the
# enriched H5AD built by scripts/build_h5ad_for_app.py

#' Load data from enriched H5AD file (islet_explorer.h5ad)
#' Returns the same list structure as load_master(): list(markers, targets, comp, lgals3)
#' @param path Path to the enriched H5AD file
#' @return list with targets, markers, comp, lgals3 DataFrames, or NULL on failure
load_master_h5ad <- function(path) {
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("[H5AD] anndata package not available, falling back to Excel")
    return(NULL)
  }

  tryCatch({
    message("[H5AD] Loading enriched data")
    ad <- anndata::read_h5ad(path)

    # Cache ad$uns as R list once (0.2s) instead of per-key bridge crossings (~14s for 62 keys)
    uns <- ad$uns

    # Reconstruct targets/markers/comp from cached uns list
    targets <- reconstruct_groovy_df_from_list(uns, "targets")
    markers <- reconstruct_groovy_df_from_list(uns, "markers")
    comp    <- reconstruct_groovy_df_from_list(uns, "composition")

    # Split LGALS3 rows from markers (they share the same groovy storage)
    lgals3 <- NULL
    if (!is.null(markers) && "marker" %in% names(markers)) {
      lgals3_mask <- markers$marker == "LGALS3"
      if (any(lgals3_mask)) {
        lgals3 <- markers[lgals3_mask, , drop = FALSE]
        markers <- markers[!lgals3_mask, , drop = FALSE]
      }
    }

    # Extract phenotype proportions from .obs
    phenotype_df <- tryCatch({
      obs <- as.data.frame(ad$obs)
      prop_cols <- grep("^prop_", colnames(obs), value = TRUE)
      if (length(prop_cols) > 0 && "imageid" %in% colnames(obs) && "base_islet_id" %in% colnames(obs)) {
        phen <- obs[, c("imageid", "base_islet_id", prop_cols), drop = FALSE]
        phen$`Case ID` <- as.integer(as.character(phen$imageid))
        phen$islet_key <- gsub("^Islet_Islet_", "Islet_", as.character(phen$base_islet_id))
        phen[, c("Case ID", "islet_key", prop_cols), drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=h5ad_phenotype_extract_failed"); NULL })

    # Extract donor demographics (one row per donor)
    donor_demographics <- tryCatch({
      obs <- as.data.frame(ad$obs)
      if (all(c("imageid", "age", "gender") %in% colnames(obs))) {
        demo <- data.frame(
          `Case ID` = as.integer(as.character(obs$imageid)),
          age = as.numeric(obs$age),
          gender = as.character(obs$gender),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        demo[!duplicated(demo$`Case ID`), , drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=h5ad_demographics_extract_failed"); NULL })

    # Extract neighborhood metrics from .obs (added by compute_neighborhood_metrics.py → build_h5ad_for_app.py)
    neighborhood_df <- tryCatch({
      obs <- as.data.frame(ad$obs)
      nbr_cols <- grep("^peri_prop_|^peri_count_|^immune_|^cd8_|^tcell_|^enrich_z_|^min_dist_|^total_cells_peri|^total_cells_core|^immune_count_|^dpt_pseudotime|^leiden_|^sc_pseudotime|^sc_n_|^crossed_",
                       colnames(obs), value = TRUE)
      if (length(nbr_cols) > 0 && "imageid" %in% colnames(obs) && "base_islet_id" %in% colnames(obs)) {
        nbr <- obs[, c("imageid", "base_islet_id", nbr_cols), drop = FALSE]
        nbr$`Case ID` <- as.integer(as.character(nbr$imageid))
        nbr$islet_key <- gsub("^Islet_Islet_", "Islet_", as.character(nbr$base_islet_id))
        nbr[, c("Case ID", "islet_key", nbr_cols), drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=h5ad_neighborhood_extract_failed"); NULL })

    message("[H5AD] Loaded: targets=", if (!is.null(targets)) nrow(targets) else 0,
            " markers=", if (!is.null(markers)) nrow(markers) else 0,
            " comp=", if (!is.null(comp)) nrow(comp) else 0,
            " lgals3=", if (!is.null(lgals3)) nrow(lgals3) else 0,
            " phenotypes=", if (!is.null(phenotype_df)) ncol(phenotype_df) - 2 else 0,
            " demographics=", if (!is.null(donor_demographics)) nrow(donor_demographics) else 0,
            " neighborhood=", if (!is.null(neighborhood_df)) ncol(neighborhood_df) - 2 else 0)

    list(markers = markers, targets = targets, comp = comp, lgals3 = lgals3,
         phenotypes = phenotype_df, donor_demographics = donor_demographics,
         neighborhood = neighborhood_df)
  }, error = function(e) {
    message("event=h5ad_load_failed")
    NULL
  })
}

#' Reconstruct a groovy DataFrame from a pre-fetched .uns R list
#' @param uns_list R list from ad$uns (single bridge crossing)
#' @param sheet One of "targets", "markers", "composition"
#' @return data.frame or NULL
reconstruct_groovy_df_from_list <- function(uns_list, sheet) {
  prefix <- paste0("groovy_", sheet, "_")
  cols_key <- paste0("groovy_", sheet, "_columns")
  nrows_key <- paste0("groovy_", sheet, "_n_rows")

  if (is.null(uns_list[[cols_key]])) return(NULL)

  col_names <- uns_list[[cols_key]]
  n_rows <- as.integer(uns_list[[nrows_key]])

  if (n_rows == 0) return(NULL)

  df <- data.frame(row.names = seq_len(n_rows))
  for (col in col_names) {
    key <- paste0(prefix, col)
    vals <- uns_list[[key]]
    if (!is.null(vals)) {
      df[[col]] <- vals
    }
  }

  # Standardize Case ID to match Excel format (integer-like)
  if ("Case ID" %in% names(df)) {
    df[["Case ID"]] <- suppressWarnings(as.integer(df[["Case ID"]]))
  }

  df
}

#' Legacy: Reconstruct a groovy DataFrame via ad$uns (slow, ~14s for 62 keys)
#' Kept for backward compatibility but unused in normal path
reconstruct_groovy_df <- function(ad, sheet) {
  reconstruct_groovy_df_from_list(ad$uns, sheet)
}

#' Load master data from DuckDB-backed Parquet views (Phase 1)
#'
#' Rebuilds the same `list(markers, targets, comp, lgals3, phenotypes,
#' donor_demographics, neighborhood)` structure the legacy H5AD loader
#' returns, but pulls rows from the `islets`, `uns_targets`, `uns_markers`,
#' and `uns_composition` DuckDB views registered in `00_globals.R`.
#'
#' This function materialises the groovy sheets into in-memory data frames
#' (they are only ~5k-200k rows). The Parquet + DuckDB advantage is that the
#' *islet* table and per-donor *tissue* table can be queried lazily and only
#' the rows a module actually needs travel across the reticulate / pyarrow
#' boundary. See `mod_spatial_server.R` and `spatial_helpers.R` for the lazy
#' path.
#'
#' @param con A live `duckdb` connection (defaults to the session-shared `con`).
#' @return list matching `load_master_h5ad()`.
load_master_duckdb <- function(con = get0("con", envir = globalenv())) {
  if (is.null(con)) {
    message("[DUCKDB] connection not available")
    return(NULL)
  }
  if (!requireNamespace("DBI", quietly = TRUE)) {
    message("[DUCKDB] DBI not available")
    return(NULL)
  }

  tryCatch({
    message("[DUCKDB] Loading islets + groovy sheets from Parquet views")
    obs <- DBI::dbGetQuery(con, "SELECT * FROM islets")

    read_view <- function(view) {
      if (!DBI::dbExistsTable(con, view)) return(NULL)
      tryCatch(DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", view)),
               error = function(e) { message("event=duckdb_view_load_failed"); NULL })
    }

    targets <- read_view("uns_targets")
    markers <- read_view("uns_markers")
    comp    <- read_view("uns_composition")

    # Split LGALS3 rows out of the markers sheet (same convention as the H5AD path).
    lgals3 <- NULL
    if (!is.null(markers) && "marker" %in% names(markers)) {
      m <- markers$marker == "LGALS3"
      if (any(m, na.rm = TRUE)) {
        lgals3 <- markers[m, , drop = FALSE]
        markers <- markers[!m, , drop = FALSE]
      }
    }

    # Phenotype proportions from .obs (prop_* columns)
    phenotype_df <- tryCatch({
      prop_cols <- grep("^prop_", colnames(obs), value = TRUE)
      if (length(prop_cols) > 0 && "imageid" %in% colnames(obs) && "base_islet_id" %in% colnames(obs)) {
        phen <- obs[, c("imageid", "base_islet_id", prop_cols), drop = FALSE]
        phen$`Case ID` <- suppressWarnings(as.integer(as.character(phen$imageid)))
        phen$islet_key <- gsub("^Islet_Islet_", "Islet_", as.character(phen$base_islet_id))
        phen[, c("Case ID", "islet_key", prop_cols), drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=duckdb_phenotype_extract_failed"); NULL })

    # Donor demographics (one row per donor)
    donor_demographics <- tryCatch({
      if (all(c("imageid", "age", "gender") %in% colnames(obs))) {
        demo <- data.frame(
          `Case ID` = suppressWarnings(as.integer(as.character(obs$imageid))),
          age = suppressWarnings(as.numeric(obs$age)),
          gender = as.character(obs$gender),
          check.names = FALSE, stringsAsFactors = FALSE
        )
        demo[!duplicated(demo$`Case ID`), , drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=duckdb_demographics_extract_failed"); NULL })

    # Neighborhood metrics + Leiden + pseudotime (same regex as load_master_h5ad)
    neighborhood_df <- tryCatch({
      nbr_cols <- grep(
        "^peri_prop_|^peri_count_|^immune_|^cd8_|^tcell_|^enrich_z_|^min_dist_|^total_cells_peri|^total_cells_core|^immune_count_|^dpt_pseudotime$|^leiden_|^sc_pseudotime|^sc_n_|^crossed_",
        colnames(obs), value = TRUE
      )
      if (length(nbr_cols) > 0 && "imageid" %in% colnames(obs) && "base_islet_id" %in% colnames(obs)) {
        nbr <- obs[, c("imageid", "base_islet_id", nbr_cols), drop = FALSE]
        nbr$`Case ID` <- suppressWarnings(as.integer(as.character(nbr$imageid)))
        nbr$islet_key <- gsub("^Islet_Islet_", "Islet_", as.character(nbr$base_islet_id))
        nbr[, c("Case ID", "islet_key", nbr_cols), drop = FALSE]
      } else NULL
    }, error = function(e) { message("event=duckdb_neighborhood_extract_failed"); NULL })

    message("[DUCKDB] Loaded: targets=", if (!is.null(targets)) nrow(targets) else 0,
            " markers=", if (!is.null(markers)) nrow(markers) else 0,
            " comp=", if (!is.null(comp)) nrow(comp) else 0,
            " lgals3=", if (!is.null(lgals3)) nrow(lgals3) else 0,
            " phenotypes=", if (!is.null(phenotype_df)) ncol(phenotype_df) - 2 else 0,
            " demographics=", if (!is.null(donor_demographics)) nrow(donor_demographics) else 0,
            " neighborhood=", if (!is.null(neighborhood_df)) ncol(neighborhood_df) - 2 else 0)

    list(markers = markers, targets = targets, comp = comp, lgals3 = lgals3,
         phenotypes = phenotype_df, donor_demographics = donor_demographics,
         neighborhood = neighborhood_df)
  }, error = function(e) {
    message("event=duckdb_master_load_failed")
    NULL
  })
}

#' Try loading from DuckDB/Parquet first, then H5AD, then Excel
#' @param h5ad_path Path to enriched H5AD (or NULL to skip)
#' @param excel_path Path to master_results.xlsx
#' @return list(markers, targets, comp, lgals3) — same structure as load_master()
load_master_auto <- function(h5ad_path = NULL, excel_path = master_path) {
  # Phase 1: Prefer DuckDB-backed Parquet views when available.
  if (isTRUE(get0("USE_DUCKDB", envir = globalenv(), ifnotfound = FALSE))) {
    result <- load_master_duckdb()
    if (!is.null(result)) {
      message("[DATA] Using DuckDB Parquet source")
      return(result)
    }
    message("[DATA] DuckDB loader returned NULL; falling back to H5AD/Excel")
  }

  # Try H5AD next
  if (!is.null(h5ad_path) && file.exists(h5ad_path)) {
    result <- load_master_h5ad(h5ad_path)
    if (!is.null(result)) {
      message("[DATA] Using H5AD source")
      return(result)
    }
  }

  # Fall back to Excel
  if (file.exists(excel_path)) {
    message("[DATA] Using Excel source")
    return(load_master(excel_path))
  }

  stop("No data source found. Checked:\n  H5AD: ", h5ad_path, "\n  Excel: ", excel_path)
}

#' Compute Lamian-based pseudotime from AnnData
#' Uses Lamian's infer_tree_structure which accounts for multi-sample design
#' @param ad AnnData object (from anndata package)
#' @param features Character vector of feature names to use (NULL = use all)
#' @return Numeric vector of pseudotime values (same length as cells)
compute_lamian_pseudotime <- function(ad, features = NULL) {
  tryCatch({
    require(Lamian, quietly = TRUE)

    cat("Computing Lamian pseudotime with infer_tree_structure...\n")

    # Extract expression matrix (cells x genes)
    expr_mat <- as.matrix(ad$X)
    rownames(expr_mat) <- ad$obs_names
    colnames(expr_mat) <- ad$var_names

    # Define curated feature set for trajectory inference
    # Focus on hormones, immune markers, and spatial features that drive disease progression
    curated_features <- c(
      # Hormone markers (islet cell types)
      "INS", "GCG",
      # Immune infiltration markers
      "CD8a", "CD4", "HLADR", "CD163", "CD68",
      # Disease-associated markers
      "LGALS3", "BCatenin",
      # Spatial features (microenvironment)
      "Dist to Closest Lymphatic", "Dist to Closest Capillary", "Dist to Closest Nerve"
    )

    # Subset to requested features if provided, otherwise use curated set
    if (!is.null(features) && length(features) > 0) {
      available_features <- intersect(features, ad$var_names)
      if (length(available_features) == 0) {
        warning("None of the requested features found in AnnData, using curated features")
        available_features <- intersect(curated_features, ad$var_names)
      } else {
        cat(sprintf("Using %d of %d requested features\n", length(available_features), length(features)))
      }
    } else {
      # Use curated feature set by default
      available_features <- intersect(curated_features, ad$var_names)
      cat(sprintf("Using curated feature set: %d features\n", length(available_features)))
      cat(sprintf("  Features: %s\n", paste(available_features, collapse=", ")))
    }

    if (length(available_features) > 0) {
      feature_idx <- match(available_features, ad$var_names)
      expr_mat <- expr_mat[, feature_idx, drop = FALSE]
      colnames(expr_mat) <- ad$var_names[feature_idx]
    } else {
      stop("No valid features found for trajectory inference")
    }

    cat("Running PCA for dimensionality reduction...\n")
    # Compute PCA (Lamian expects cells x PCs)
    # Data is already z-scored, so no need for redundant centering/scaling
    n_pcs <- min(length(available_features) - 1, nrow(expr_mat) - 1)
    pca_result <- prcomp(expr_mat, rank. = n_pcs, center = FALSE, scale. = FALSE)
    pca_coords <- pca_result$x

    # Report variance explained
    var_explained <- summary(pca_result)$importance[2, ]
    cumvar <- cumsum(var_explained)
    n_pcs_80 <- if (any(cumvar >= 0.80)) which(cumvar >= 0.80)[1] else n_pcs
    cat(sprintf("  Using %d PCs (explains %.1f%% variance, %d PCs for 80%%)\n",
                n_pcs, cumvar[n_pcs] * 100, n_pcs_80))

    cat("Building cell annotation with sample IDs...\n")
    # Build cell annotation - CRITICAL: column 2 must be sample/donor ID
    # This allows Lamian to account for multi-sample structure
    cellanno <- data.frame(
      cell = ad$obs_names,
      sample = as.character(ad$obs$imageid),  # Donor/sample ID
      stringsAsFactors = FALSE
    )

    # Add cell type if available (optional)
    if ("donor_status" %in% colnames(ad$obs)) {
      cellanno$celltype <- as.character(ad$obs$donor_status)
    }

    cat("Running Lamian trajectory inference...\n")
    cat("  This accounts for donor/sample structure to avoid artificial grouping\n")
    cat("  Origin marker: INS (clusters with highest mean INS expression)\n")

    # Use Lamian's infer_tree_structure (designed for multi-sample data)
    # Transpose expression for Lamian (genes x cells expected)
    expr_mat_t <- t(expr_mat)

    # Determine reasonable max cluster number based on data size
    n_obs <- nrow(ad$obs)
    n_samples <- length(unique(cellanno$sample))
    # Use more clusters for better resolution: ~sqrt(n) or n/50, capped at 50
    max_clusters <- min(50, max(10, ceiling(sqrt(n_obs)), ceiling(n_obs / 50)))

    cat(sprintf("  Max clusters: %d (based on %d observations, %d samples)\n",
                max_clusters, n_obs, n_samples))

    res <- Lamian::infer_tree_structure(
      pca = pca_coords,
      cellanno = cellanno,
      expression = expr_mat_t,
      origin.marker = "INS",           # Use insulin as trajectory root - Lamian finds cluster with highest mean INS
      number.cluster = NA,              # Auto-determine cluster number
      max.clunum = max_clusters,        # Maximum clusters to consider (adaptive)
      kmeans.seed = 12345              # Reproducible clustering
    )

    # Extract pseudotime from Lamian result
    pseudotime <- res$pseudotime

    # Report what Lamian found
    if (!is.null(res$clusterRes)) {
      n_clusters <- length(unique(res$clusterRes))
      cat(sprintf("  Lamian identified %d clusters\n", n_clusters))
    }

    # Ensure it's a vector aligned with cells
    if (is.matrix(pseudotime)) {
      pseudotime <- as.numeric(pseudotime[, 1])
    }

    # Handle names if present
    if (!is.null(names(pseudotime))) {
      # Reorder to match ad$obs_names
      pseudotime <- pseudotime[ad$obs_names]
    }

    # Check pseudotime distribution before normalization
    cat(sprintf("  Raw pseudotime range: %.3f - %.3f (mean: %.3f, sd: %.3f)\n",
                min(pseudotime, na.rm = TRUE), max(pseudotime, na.rm = TRUE),
                mean(pseudotime, na.rm = TRUE), sd(pseudotime, na.rm = TRUE)))

    # Normalize to 0-1 range for consistency with PAGA
    pseudotime <- (pseudotime - min(pseudotime, na.rm = TRUE)) /
                  (max(pseudotime, na.rm = TRUE) - min(pseudotime, na.rm = TRUE))

    # REVERSE pseudotime direction (biological progression ND->Aab+->T1D)
    # Lamian origin.marker='INS' finds high-INS clusters, but trajectory may run backward
    # Reversal ensures: 0=healthy (high INS), 1=diseased (low INS)
    pseudotime <- 1 - pseudotime
    cat("  Pseudotime reversed for biological progression (0=healthy, 1=diseased)\n")

    cat(sprintf("Lamian pseudotime computed (range: %.3f - %.3f)\n",
                min(pseudotime, na.rm = TRUE), max(pseudotime, na.rm = TRUE)))
    cat(sprintf("  %d cells across %d samples\n",
                length(pseudotime), length(unique(cellanno$sample))))

    return(pseudotime)

  }, error = function(e) {
    warning(sprintf("Lamian pseudotime calculation failed: %s", e$message))
    cat("Error details:", e$message, "\n")
    return(rep(NA_real_, nrow(ad$obs)))
  })
}
