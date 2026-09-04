# Security-boundary helpers are dependency-free so Shiny can source them before
# package loading and application startup.

.ISLET_VALIDATION_ERROR_CODES <- c(
  "viewer_image_invalid",
  "spatial_donor_invalid",
  "plot_mode_invalid",
  "plot_groups_invalid",
  "plot_selection_invalid",
  "plot_metric_invalid",
  "plot_min_cells_invalid",
  "trajectory_mode_invalid",
  "trajectory_feature_invalid",
  "trajectory_ratio_invalid",
  "trajectory_region_invalid",
  "islet_key_invalid",
  "drilldown_color_invalid"
)

islet_flag_enabled <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    tolower(trimws(value)) %in% c("1", "true", "yes", "on")
}

islet_runtime_mode <- function() {
  tolower(trimws(Sys.getenv("ISLET_RUNTIME_MODE", "development")))
}

islet_is_production <- function() identical(islet_runtime_mode(), "production")

# Emit one fixed-schema line for a rejected browser-controlled value. The
# rejected value is deliberately not accepted as an argument, preventing log
# injection and accidental disclosure.
islet_log_validation_failure <- function(error_code) {
  if (!is.character(error_code) || length(error_code) != 1L ||
      is.na(error_code) || !(error_code %in% .ISLET_VALIDATION_ERROR_CODES)) {
    error_code <- "invalid_input"
  }
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  message(sprintf(
    paste0('{"ts":"%s","level":"WARN","event":"input_validation_failed",',
           '"component":"shiny","result":"rejected","error_code":"%s"}'),
    timestamp, error_code
  ))
  invisible(FALSE)
}

islet_valid_scalar_choice <- function(value, allowed, error_code,
                                      log_missing = FALSE) {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(value) && value %in% as.character(allowed)
  if (!valid && (log_missing || !is.null(value))) {
    islet_log_validation_failure(error_code)
  }
  valid
}

islet_valid_choice_vector <- function(value, allowed, error_code,
                                      allow_empty = FALSE) {
  valid <- is.character(value) && !anyNA(value) &&
    (allow_empty || length(value) > 0L) && all(value %in% as.character(allowed))
  if (!valid && !is.null(value)) islet_log_validation_failure(error_code)
  valid
}

islet_valid_bounded_number <- function(value, lower, upper, error_code) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= lower && value <= upper
  if (!valid && !is.null(value)) islet_log_validation_failure(error_code)
  valid
}

islet_valid_asset_basename <- function(value, allowed) {
  valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
    nzchar(value) && !grepl("[/\\\\]", value) &&
    grepl("[.]ome[.]tiff?$", value, ignore.case = TRUE) &&
    value %in% as.character(allowed)
  if (!valid && !is.null(value)) {
    islet_log_validation_failure("viewer_image_invalid")
  }
  valid
}

# Defense in depth for values that become a single filesystem path segment
# after first being checked against an authoritative application lookup.
islet_safe_path_segment <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", value) &&
    !(value %in% c(".", ".."))
}
