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
DONORS_DIR <- file.path("..", "..", "data", "donors")

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
