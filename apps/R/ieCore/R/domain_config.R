# Domain configuration engine extracted from LymphVizkit@820bb38. Configuration
# controls application vocabulary and display geometry only; it is never a
# source of cell identity, marker thresholds, or upstream scientific calls.

`%||%` <- function(x, y) if (is.null(x)) y else x

.iecore_state <- new.env(parent = emptyenv())

.domain_scalar <- function(value, label, type = "character") {
  valid <- length(value) == 1L && !is.null(value) && !is.na(value)
  if (type == "character") valid <- valid && is.character(value) && nzchar(value)
  if (type == "numeric") valid <- valid && is.numeric(value) && is.finite(value)
  if (!valid) stop("Invalid DOMAIN field: ", label, call. = FALSE)
  invisible(TRUE)
}

.domain_path <- function(domain, path, fallback = NULL) {
  value <- domain
  for (part in path) {
    if (!is.list(value) || is.null(value[[part]])) return(fallback)
    value <- value[[part]]
  }
  value
}

ie_domain_validate <- function(domain, source = "<memory>") {
  if (!is.list(domain)) stop("DOMAIN must be a mapping: ", source, call. = FALSE)

  .domain_scalar(.domain_path(domain, c("app", "name")), "app.name")
  .domain_scalar(.domain_path(domain, c("app", "title")), "app.title")
  .domain_scalar(.domain_path(domain, c("unit", "singular")), "unit.singular")
  .domain_scalar(.domain_path(domain, c("unit", "plural")), "unit.plural")
  .domain_scalar(.domain_path(domain, c("unit", "key_column")), "unit.key_column")
  .domain_scalar(.domain_path(domain, c("unit", "id_prefix")), "unit.id_prefix")
  .domain_scalar(.domain_path(domain, c("unit", "id_regex")), "unit.id_regex")

  pixel_size <- suppressWarnings(as.numeric(.domain_path(domain, "pixel_size_um")))
  .domain_scalar(pixel_size, "pixel_size_um", "numeric")
  if (pixel_size <= 0) stop("DOMAIN pixel_size_um must be positive", call. = FALSE)

  levels <- unlist(.domain_path(domain, c("grouping", "levels")), use.names = FALSE)
  if (!length(levels) || anyNA(levels) || any(!nzchar(levels)) || anyDuplicated(levels)) {
    stop("DOMAIN grouping.levels must be non-empty and unique", call. = FALSE)
  }
  palettes <- .domain_path(domain, c("grouping", "palettes"), list())
  if (!is.list(palettes) || !length(palettes)) {
    stop("DOMAIN grouping.palettes must be a non-empty mapping", call. = FALSE)
  }
  for (name in names(palettes)) {
    palette <- unlist(palettes[[name]], use.names = TRUE)
    if (is.null(names(palette)) || !all(levels %in% names(palette))) {
      stop("DOMAIN palette '", name, "' does not cover grouping.levels", call. = FALSE)
    }
  }

  order <- unlist(.domain_path(domain, c("phenotypes", "display_order")), use.names = FALSE)
  colors <- unlist(.domain_path(domain, c("phenotypes", "colors")), use.names = TRUE)
  if (!length(order) || anyNA(order) || any(!nzchar(order)) || anyDuplicated(order)) {
    stop("DOMAIN phenotypes.display_order must be non-empty and unique", call. = FALSE)
  }
  if (is.null(names(colors)) || !all(order %in% names(colors))) {
    stop("DOMAIN phenotypes.colors does not cover display_order", call. = FALSE)
  }

  regions <- .domain_path(domain, c("regions", "items"), list())
  canonicals <- vapply(regions, function(item) as.character(item$canonical %||% ""), character(1))
  if (!length(canonicals) || any(!nzchar(canonicals)) || anyDuplicated(canonicals)) {
    stop("DOMAIN regions.items canonicals must be non-empty and unique", call. = FALSE)
  }
  if (!any(vapply(regions, function(item) isTRUE(item$is_core), logical(1)))) {
    stop("DOMAIN regions.items must declare one core region", call. = FALSE)
  }

  attr(domain, "config_source") <- source
  domain
}

ie_load_domain <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop("DOMAIN configuration not found: ", paste(path, collapse = ""), call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The yaml package is required to load DOMAIN configuration", call. = FALSE)
  }
  normalized <- normalizePath(path, mustWork = TRUE)
  domain <- yaml::read_yaml(normalized)
  ie_domain_validate(domain, normalized)
}

ie_set_domain <- function(domain) {
  source <- attr(domain, "config_source") %||% "<memory>"
  .iecore_state$domain <- ie_domain_validate(domain, source)
  options(ieCore.domain = .iecore_state$domain)
  invisible(.iecore_state$domain)
}

ie_domain <- function(domain = NULL) {
  if (!is.null(domain)) return(domain)
  if (exists("domain", envir = .iecore_state, inherits = FALSE)) {
    return(.iecore_state$domain)
  }
  configured <- getOption("ieCore.domain")
  if (!is.null(configured)) return(configured)
  stop("No DOMAIN is active; call ie_set_domain(ie_load_domain(path))", call. = FALSE)
}

domain_config_path <- function(domain = NULL) attr(ie_domain(domain), "config_source")
domain_app_name <- function(domain = NULL) .domain_path(ie_domain(domain), c("app", "name"), "app")
domain_app_title <- function(domain = NULL) .domain_path(ie_domain(domain), c("app", "title"), domain_app_name(domain))

domain_unit_singular <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "singular"), "unit")
domain_unit_plural <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "plural"), "units")
domain_unit_Title <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "Title"), "Unit")
domain_key_col <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "key_column"), "unit_key")
domain_id_prefix <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "id_prefix"), "Unit")
domain_id_regex <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "id_regex"), "Unit_\\d+")
domain_peri_suffix <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "peri_suffix"), "_exp20um")
domain_pad_width <- function(domain = NULL) as.integer(.domain_path(ie_domain(domain), c("unit", "case_id_pad_width"), 4L))
domain_diam_col <- function(domain = NULL) .domain_path(ie_domain(domain), c("unit", "diam_column"), "unit_diam_um")
domain_pixel_size <- function(domain = NULL) as.numeric(.domain_path(ie_domain(domain), "pixel_size_um"))

domain_group_column <- function(domain = NULL) .domain_path(ie_domain(domain), c("grouping", "physical_column"), "Donor Status")
domain_group_label <- function(domain = NULL) .domain_path(ie_domain(domain), c("grouping", "display_label"), "Group")
domain_group_levels <- function(domain = NULL) unlist(.domain_path(ie_domain(domain), c("grouping", "levels")), use.names = FALSE)
domain_group_ordered <- function(domain = NULL) isTRUE(.domain_path(ie_domain(domain), c("grouping", "ordered"), TRUE))
domain_group_default_palette <- function(domain = NULL) .domain_path(ie_domain(domain), c("grouping", "default_palette"), "Bright")
domain_group_palettes <- function(domain = NULL) {
  palettes <- .domain_path(ie_domain(domain), c("grouping", "palettes"), list())
  lapply(palettes, function(palette) unlist(palette, use.names = TRUE))
}
domain_group_colors <- function(palette = domain_group_default_palette(), domain = NULL) {
  palettes <- domain_group_palettes(domain)
  chosen <- palettes[[palette]]
  if (is.null(chosen) && length(palettes)) chosen <- palettes[[1L]]
  chosen %||% stats::setNames(rep("#888888", length(domain_group_levels(domain))), domain_group_levels(domain))
}
domain_group_factor <- function(x, domain = NULL) {
  factor(as.character(x), levels = domain_group_levels(domain), ordered = domain_group_ordered(domain))
}

domain_region_items <- function(domain = NULL) .domain_path(ie_domain(domain), c("regions", "items"), list())
domain_region_canonicals <- function(domain = NULL) vapply(domain_region_items(domain), function(item) item$canonical, character(1))
domain_region_core <- function(domain = NULL) {
  core <- Filter(function(item) isTRUE(item$is_core), domain_region_items(domain))
  core[[1L]]$canonical
}
domain_region_display <- function(domain = NULL) {
  items <- domain_region_items(domain)
  stats::setNames(vapply(items, function(item) item$display, character(1)),
                  vapply(items, function(item) item$canonical, character(1)))
}
domain_region_geojson_classes <- function(domain = NULL) {
  classes <- vapply(domain_region_items(domain), function(item) {
    value <- item$geojson_class
    if (is.null(value) || is.na(value)) "" else as.character(value)
  }, character(1))
  classes[nzchar(classes)]
}
domain_expansion_class <- function(domain = NULL) .domain_path(ie_domain(domain), c("regions", "expansion", "geojson_class"), "ExpandedUnit")
domain_expansion_display <- function(domain = NULL) .domain_path(ie_domain(domain), c("regions", "expansion", "display"), "Peri-unit")
domain_expansion_um <- function(domain = NULL) as.numeric(.domain_path(ie_domain(domain), c("regions", "expansion", "distance_um"), 20))
domain_structures <- function(domain = NULL) .domain_path(ie_domain(domain), c("regions", "structures"), list())
domain_structure_classes <- function(domain = NULL) vapply(domain_structures(domain), function(item) item$class, character(1))

domain_excluded_markers <- function(domain = NULL) unlist(.domain_path(ie_domain(domain), c("markers", "excluded"), list()), use.names = FALSE)
domain_defining <- function(domain = NULL) .domain_path(ie_domain(domain), c("markers", "defining"), list())
domain_defining_markers <- function(domain = NULL) vapply(domain_defining(domain), function(item) item$marker, character(1))
domain_defining_single_cols <- function(domain = NULL) vapply(domain_defining(domain), function(item) item$single_col, character(1))
domain_defining_any_cols <- function(domain = NULL) vapply(domain_defining(domain), function(item) item$any_col, character(1))
domain_pseudotime_root <- function(domain = NULL) .domain_path(ie_domain(domain), c("markers", "pseudotime_root"), "")
domain_marker_panel <- function(domain = NULL) unlist(.domain_path(ie_domain(domain), c("markers", "panel"), list()), use.names = FALSE)

domain_phenotype_order <- function(domain = NULL) unlist(.domain_path(ie_domain(domain), c("phenotypes", "display_order"), list()), use.names = FALSE)
domain_phenotype_colors <- function(domain = NULL) unlist(.domain_path(ie_domain(domain), c("phenotypes", "colors"), list()), use.names = TRUE)
domain_phenotype_aliases <- function(domain = NULL) unlist(.domain_path(ie_domain(domain), c("phenotypes", "aliases"), list()), use.names = TRUE)
domain_phenotype_residuals <- function(domain = NULL) .domain_path(ie_domain(domain), c("phenotypes", "residuals"), list())
domain_feature <- function(name, default = FALSE, domain = NULL) isTRUE(.domain_path(ie_domain(domain), c("features", name), default))
