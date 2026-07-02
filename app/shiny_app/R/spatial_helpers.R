# ===== Spatial Tab Helper Functions =====
# Provides per-donor tissue-wide cell loading for the tissue scatter plot.
#
# Phase 1 (scaling): when `USE_DUCKDB` is TRUE and the `tissue` view is
# registered in DuckDB, donor loads become a partition-pruned SQL query
# against `data/parquet/tissue/case_id=<id>/*.parquet`. The fall-back path
# continues to read CSVs from `data/donors/` so the legacy single-dataset
# deployments keep working.

# ── Cache environment for per-donor tissue data ──────────────────────────
# Still used as an in-memory LRU so repeated brush interactions on the same
# donor do not round-trip to DuckDB. DuckDB's own buffer cache also helps,
# but an R-side cache skips the Arrow -> R conversion on repeat hits.
.donor_tissue_cache <- new.env(parent = emptyenv())

# Path to per-donor tissue CSVs (legacy fallback)
DONORS_DIR <- file.path("..", "..", "data", "app_data", "donors")

#' Check whether the `tissue` DuckDB view exists in the shared connection.
duckdb_tissue_available <- function() {
  if (!isTRUE(get0("USE_DUCKDB", envir = globalenv(), ifnotfound = FALSE))) return(FALSE)
  con <- get0("con", envir = globalenv(), ifnotfound = NULL)
  if (is.null(con)) return(FALSE)
  tryCatch(DBI::dbExistsTable(con, "tissue"), error = function(e) FALSE)
}

#' Check if per-donor tissue data is available (DuckDB OR CSV)
donor_tissue_available <- function() {
  if (duckdb_tissue_available()) return(TRUE)
  dir.exists(DONORS_DIR) && length(list.files(DONORS_DIR, pattern = "\\.csv$")) > 0
}

#' Get list of available donor image IDs
get_available_donors <- function() {
  if (duckdb_tissue_available()) {
    con <- get0("con", envir = globalenv(), ifnotfound = NULL)
    donors <- tryCatch(
      DBI::dbGetQuery(con, "SELECT DISTINCT case_id FROM tissue ORDER BY case_id")$case_id,
      error = function(e) NULL
    )
    if (!is.null(donors) && length(donors) > 0) return(as.character(donors))
  }
  files <- list.files(DONORS_DIR, pattern = "\\.csv$")
  gsub("\\.csv$", "", files)
}

#' Load per-donor tissue cells with caching
#'
#' Prefers DuckDB partition pruning when the `tissue` view is available;
#' otherwise reads the legacy CSV. The returned data.frame has the same
#' columns (`X_centroid`, `Y_centroid`, `phenotype`, `cell_region`,
#' `islet_name`) regardless of source.
#'
#' @param imageid Donor image ID (e.g., "6505")
#' @return data.frame or NULL if the donor is not available
load_donor_tissue <- function(imageid) {
  cache_key <- as.character(imageid)

  if (exists(cache_key, envir = .donor_tissue_cache)) {
    return(get(cache_key, envir = .donor_tissue_cache))
  }

  df <- NULL

  if (duckdb_tissue_available()) {
    con <- get0("con", envir = globalenv(), ifnotfound = NULL)
    df <- tryCatch({
      # Partition-pruned scan: DuckDB only touches the hive sub-directory for
      # this donor. The SELECT list excludes `case_id` (the partition key)
      # because downstream code relies on the historical 5-column layout.
      d <- DBI::dbGetQuery(con,
        "SELECT X_centroid, Y_centroid, phenotype, cell_region, islet_name
           FROM tissue
          WHERE case_id = ?",
        params = list(cache_key)
      )
      if (is.null(d) || nrow(d) == 0) NULL else d
    }, error = function(e) {
      message("[SPATIAL] DuckDB tissue query failed for ", imageid, ": ", e$message)
      NULL
    })
  }

  # Legacy CSV fallback (and dev boxes without Parquet files)
  if (is.null(df)) {
    fpath <- file.path(DONORS_DIR, paste0(imageid, ".csv"))
    if (!file.exists(fpath)) return(NULL)
    df <- tryCatch({
      d <- read.csv(fpath, stringsAsFactors = FALSE)
      if (nrow(d) == 0) NULL else d
    }, error = function(e) {
      message("[SPATIAL] Error loading donor tissue for ", imageid, ": ", e$message)
      NULL
    })
  }

  if (!is.null(df)) {
    assign(cache_key, df, envir = .donor_tissue_cache)
    message("[SPATIAL] Cached donor tissue for ", imageid, " (", nrow(df), " cells)")
  }

  df
}

## ===== Card C right plot: per-cell distance-to-islet, DuckDB-backed =====
# Reads ONE partition-pruned view registered in 00_globals.R:
#   cell_distance -> data/parquet/cell_distance/case_id=*
#     columns: donor_status, cell_region, phenotype, dist_islet (signed µm;
#              negative = inside islet, positive = outside), + `{m}_pos` booleans
#              for the 10 RESTORE-gated lineage markers, + case_id (hive).
# Marker positivity is the RESTORE `_pos` call (per-image calibrated threshold on
# REDSEA-corrected MFI), NOT a raw-intensity heuristic. Built by
# scripts/senior/build_cell_distance_parquet.py. See app/shiny_app/docs/spatial.md.

.cd_con <- function() get0("con", envir = globalenv(), ifnotfound = NULL)
cell_distance_view_ready <- function() {
  con <- .cd_con()
  if (is.null(con)) return(FALSE)
  tryCatch(DBI::dbExistsTable(con, "cell_distance"), error = function(e) FALSE)
}

.cd_columns <- function() {
  if (!cell_distance_view_ready()) return(character(0))
  tryCatch(colnames(DBI::dbGetQuery(.cd_con(), "SELECT * FROM cell_distance LIMIT 0")),
           error = function(e) character(0))
}

#' Distinct phenotype labels available in the cell_distance view (sorted).
cell_distance_phenotypes <- function() {
  if (!cell_distance_view_ready()) return(character(0))
  tryCatch(
    sort(DBI::dbGetQuery(.cd_con(),
      "SELECT DISTINCT phenotype FROM cell_distance WHERE phenotype IS NOT NULL ORDER BY phenotype"
    )$phenotype),
    error = function(e) character(0)
  )
}

#' The RESTORE-gated markers (one per `{m}_pos` column). Order preserved from the
#' parquet so the dropdown reads INS, GCG, SST, then immune/structural.
cell_distance_markers <- function() {
  sub("_pos$", "", grep("_pos$", .cd_columns(), value = TRUE))
}

# Distinct colours for the composition stack (Tableau-10; keyed by marker).
CELL_DIST_MARKER_COLORS <- c(
  INS = "#4E79A7", GCG = "#F28E2B", SST = "#E15759", CD3e = "#76B7B2",
  CD20 = "#59A14F", CD163 = "#EDC948", CD31 = "#B07AA1", SMA = "#FF9DA7",
  Vimentin = "#9C755F", Pan_Cytokeratin = "#BAB0AC"
)

# Shared SQL helpers.
.cd_esc <- function(x) gsub("'", "''", x, fixed = TRUE)
# NB: `dist_islet` is signed and often negative, so the bin offset MUST be
# parenthesised — `dist_islet - -70` becomes `dist_islet --70`, and `--` opens a
# SQL line comment that silently truncates the query. Always `(dist_islet - (%f))`.

#' Per-cell signed-distance density for the Card C right plot (Density view).
#'
#' @param mode "phenotype" or "marker".
#' @param group a phenotype label, a RESTORE marker name, or "__all__" (all cells).
#' @param statuses character subset of c("ND","Aab+","T1D").
#' @param n_bins number of distance bins.
#' @return NULL when the view is missing or no status is selected; otherwise
#'   list(hist = data.frame(donor_status, x, dens, n) | NULL, lo, hi, bw,
#'   n_by_status (named)). `hist` is NULL when the selection is empty.
query_cell_distance_hist <- function(mode, group, statuses, n_bins = 160L) {
  con <- .cd_con()
  if (is.null(con) || !cell_distance_view_ready()) return(NULL)
  statuses <- statuses[statuses %in% c("ND", "Aab+", "T1D")]
  if (length(statuses) == 0) return(NULL)
  if (is.null(group) || !nzchar(group)) group <- "__all__"
  q <- function(sql) tryCatch(DBI::dbGetQuery(con, sql),
                              error = function(e) { message("[CELL_DIST] ", conditionMessage(e)); NULL })
  status_sql <- paste(sprintf("'%s'", .cd_esc(statuses)), collapse = ",")

  is_marker <- identical(mode, "marker") && !identical(group, "__all__")
  # Guard the mode-switch transient: input$dist_mode can flip to "marker" while
  # input$dist_group still holds a phenotype label (the dropdown repopulates a
  # beat later). Treat an unknown group as all-cells until a real marker is set.
  if (is_marker && !(group %in% cell_distance_markers())) { is_marker <- FALSE; group <- "__all__" }
  cond <- if (is_marker) sprintf('"%s_pos"', group)
          else if (identical(group, "__all__")) NULL
          else sprintf("phenotype = '%s'", .cd_esc(group))
  wsel <- paste0("donor_status IN (", status_sql, ")", if (!is.null(cond)) paste0(" AND ", cond) else "")

  lim <- q(sprintf("SELECT quantile_cont(dist_islet,0.01) lo, quantile_cont(dist_islet,0.99) hi FROM cell_distance WHERE %s", wsel))
  if (is.null(lim) || !is.finite(lim$lo) || !is.finite(lim$hi) || lim$hi <= lim$lo) return(NULL)
  bw <- (lim$hi - lim$lo) / n_bins
  agg <- q(sprintf("WITH pop AS (SELECT donor_status, dist_islet FROM cell_distance WHERE %s AND dist_islet BETWEEN %f AND %f) SELECT donor_status, floor((dist_islet - (%f))/%f) AS binx, count(*) AS n FROM pop GROUP BY 1,2",
                   wsel, lim$lo, lim$hi, lim$lo, bw))
  if (is.null(agg) || nrow(agg) == 0)
    return(list(hist = NULL, lo = lim$lo, hi = lim$hi, bw = bw, n_by_status = integer(0)))

  agg$x <- lim$lo + (agg$binx + 0.5) * bw
  agg$dens <- agg$n / bw                     # cells per µm of distance
  list(hist = agg, lo = lim$lo, hi = lim$hi, bw = bw,
       n_by_status = tapply(agg$n, agg$donor_status, sum))
}

#' Marker-composition vs distance (Composition view): for each (status, distance
#' bin), the count of cells positive for each selected RESTORE marker. Returned
#' long, for grouped-stacked bars.
#'
#' @param markers character vector of RESTORE marker names.
#' @param statuses character subset of c("ND","Aab+","T1D").
#' @param n_bins number of (coarse) distance bins.
#' @return NULL if the view/inputs are empty; else list(comp = data.frame(
#'   donor_status, binx, marker, n, x) | NULL, lo, hi, bw, markers).
query_cell_distance_composition <- function(markers, statuses, n_bins = 16L) {
  con <- .cd_con()
  if (is.null(con) || !cell_distance_view_ready()) return(NULL)
  statuses <- statuses[statuses %in% c("ND", "Aab+", "T1D")]
  markers <- markers[markers %in% cell_distance_markers()]
  if (length(statuses) == 0 || length(markers) == 0) return(NULL)
  q <- function(sql) tryCatch(DBI::dbGetQuery(con, sql),
                              error = function(e) { message("[CELL_DIST] ", conditionMessage(e)); NULL })
  status_sql <- paste(sprintf("'%s'", .cd_esc(statuses)), collapse = ",")
  pos_or <- paste(sprintf('"%s_pos"', markers), collapse = " OR ")
  # X-range = 2–98 pct of distance among cells positive for ANY selected marker,
  # NOT all cells: 23M cells sit far out in tissue, so an all-cells range would
  # shove the whole inside-islet transition into the first bar. Clipping to the
  # selected markers' positive cells focuses the bins on where they actually live.
  lim <- q(sprintf("SELECT quantile_cont(dist_islet,0.02) lo, quantile_cont(dist_islet,0.98) hi FROM cell_distance WHERE donor_status IN (%s) AND (%s)", status_sql, pos_or))
  if (is.null(lim) || !is.finite(lim$lo) || !is.finite(lim$hi) || lim$hi <= lim$lo) return(NULL)
  # Bound to a near-islet window: structural markers (CD31, Vimentin, …) are
  # ubiquitous in bulk tissue and would otherwise stretch the axis to ~1.4 mm,
  # compressing the whole inside→outside boundary transition into one bar.
  lo <- max(lim$lo, -100); hi <- min(lim$hi, 300)
  if (hi <= lo) return(NULL)
  bw <- (hi - lo) / n_bins
  sums <- paste(sprintf('SUM(CASE WHEN "%s_pos" THEN 1 ELSE 0 END) AS "%s"', markers, markers), collapse = ", ")
  pos_cols <- paste(sprintf('"%s_pos"', markers), collapse = ", ")
  wide <- q(sprintf("WITH pop AS (SELECT donor_status, dist_islet, %s FROM cell_distance WHERE donor_status IN (%s) AND dist_islet BETWEEN %f AND %f) SELECT donor_status, floor((dist_islet - (%f))/%f) AS binx, %s FROM pop GROUP BY 1,2",
                    pos_cols, status_sql, lo, hi, lo, bw, sums))
  if (is.null(wide) || nrow(wide) == 0) return(list(comp = NULL, lo = lo, hi = hi, bw = bw, markers = markers))
  long <- do.call(rbind, lapply(markers, function(m) data.frame(
    donor_status = as.character(wide$donor_status), binx = wide$binx, marker = m,
    n = as.numeric(wide[[m]]), stringsAsFactors = FALSE)))
  long$x <- lo + (long$binx + 0.5) * bw
  list(comp = long, lo = lo, hi = hi, bw = bw, markers = markers)
}

#' Query cells inside an axis-aligned bounding box (for click-to-select).
#'
#' Uses DuckDB's partition pruning when available so a single brush on a
#' multi-million-cell donor stays sub-second. Falls back to an R-side filter
#' on the cached data.frame when DuckDB is not active.
#'
#' @param imageid Donor image ID.
#' @param xmin,xmax,ymin,ymax Bounding box in um (matches `X_centroid`/`Y_centroid`).
#' @return data.frame of matching cells or NULL.
query_tissue_bbox <- function(imageid, xmin, xmax, ymin, ymax) {
  if (duckdb_tissue_available()) {
    con <- get0("con", envir = globalenv(), ifnotfound = NULL)
    tryCatch({
      DBI::dbGetQuery(con,
        "SELECT X_centroid, Y_centroid, phenotype, cell_region, islet_name
           FROM tissue
          WHERE case_id = ?
            AND X_centroid BETWEEN ? AND ?
            AND Y_centroid BETWEEN ? AND ?",
        params = list(as.character(imageid), xmin, xmax, ymin, ymax)
      )
    }, error = function(e) {
      message("[SPATIAL] bbox query failed: ", e$message)
      NULL
    })
  } else {
    df <- load_donor_tissue(imageid)
    if (is.null(df)) return(NULL)
    keep <- df$X_centroid >= xmin & df$X_centroid <= xmax &
            df$Y_centroid >= ymin & df$Y_centroid <= ymax
    df[keep, , drop = FALSE]
  }
}
