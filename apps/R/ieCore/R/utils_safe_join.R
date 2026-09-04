# Safe wrapper around dplyr::left_join to handle non-data.frame inputs
safe_left_join <- function(x, y, by, context = "join") {
  # Ensure inputs are data.frames
  if (!inherits(x, "data.frame")) {
    x <- tryCatch(as.data.frame(x), error = function(e) {
      stop(sprintf("[safe_left_join] %s: Cannot coerce LHS to data.frame (class: %s)",
                   context, paste(class(x), collapse = "/")))
    })
  }
  if (!inherits(y, "data.frame")) {
    y <- tryCatch(as.data.frame(y), error = function(e) {
      stop(sprintf("[safe_left_join] %s: Cannot coerce RHS to data.frame (class: %s)",
                   context, paste(class(y), collapse = "/")))
    })
  }

  # Perform join with error context
  tryCatch(
    dplyr::left_join(x, y, by = by),
    error = function(e) {
      stop(sprintf("[safe_left_join] %s: %s", context, e$message))
    }
  )
}

add_unit_key <- function(df, domain = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)

  key_col <- domain_key_col(domain)
  key_regex <- domain_id_regex(domain)
  key_prefix <- domain_id_prefix(domain)

  if (!key_col %in% names(df)) df[[key_col]] <- NA_character_

  if ("region" %in% names(df)) {
    key_from_region <- stringr::str_extract(df$region, key_regex)
    df[[key_col]] <- dplyr::coalesce(df[[key_col]], key_from_region)
  }

  if ("name" %in% names(df)) {
    key_from_name <- stringr::str_extract(df$name, key_regex)
    only_digits  <- stringr::str_extract(df$name, "\\d+")
    fallback_name <- ifelse(!is.na(only_digits), paste0(key_prefix, "_", only_digits), NA_character_)
    df[[key_col]] <- dplyr::coalesce(df[[key_col]], key_from_name, fallback_name)
  }

  id_candidates <- unique(c(paste0(tolower(domain_unit_singular(domain)), "_id"), "islet_id"))
  id_col <- id_candidates[id_candidates %in% names(df)][1L]
  if (!is.na(id_col)) {
    id_str <- suppressWarnings(as.character(df[[id_col]]))
    id_digits <- stringr::str_extract(id_str, "\\d+")
    fallback_id <- ifelse(!is.na(id_digits), paste0(key_prefix, "_", id_digits), NA_character_)
    df[[key_col]] <- dplyr::coalesce(df[[key_col]], fallback_id)
  }

  df
}

# Compatibility API used by the islet app. Its behavior is now governed by the
# active DOMAIN rather than a filename constant.
add_islet_key <- function(df) add_unit_key(df)

compute_diameter_um <- function(area_um2) {
  area_um2 <- suppressWarnings(as.numeric(area_um2))
  ifelse(is.finite(area_um2) & area_um2 > 0, 2 * sqrt(area_um2 / pi), NA_real_)
}
