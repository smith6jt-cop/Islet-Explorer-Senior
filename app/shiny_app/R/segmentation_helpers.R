# ===== Islet Segmentation Viewer Helper Functions =====
# Extracted from app.R -- GeoJSON loading, spatial queries, and segmentation data helpers.
# References PIXEL_SIZE_UM and SF_AVAILABLE from global.R.

# Cache for loaded GeoJSON data (persists across sessions for performance)
geojson_cache <- new.env()

# islet_spatial_lookup is loaded at bottom of file via load_islet_spatial_lookup()

# Load GeoJSON for a specific case_id with caching
load_case_geojson <- function(case_id) {
  if (!SF_AVAILABLE) return(NULL)

  cache_key <- as.character(case_id)
  if (exists(cache_key, envir = geojson_cache)) {
    return(get(cache_key, envir = geojson_cache))
  }

  # Build candidate IDs: original + zero-padded variant for short numeric IDs (e.g. "112" -> "0112")
  candidate_ids <- case_id
  case_num <- suppressWarnings(as.integer(case_id))
  if (!is.na(case_num)) {
    padded <- sprintf("%04d", case_num)
    if (padded != case_id) candidate_ids <- c(candidate_ids, padded)
  }

  # Try each candidate: uncompressed first, then compressed
  json_path <- NULL
  for (cid in candidate_ids) {
    p <- file.path("..", "..", "data", "app_data", "json", paste0(cid, ".geojson"))
    if (file.exists(p)) { json_path <- p; break }
    p <- file.path("..", "..", "data", "app_data", "gson", paste0(cid, ".geojson.gz"))
    if (file.exists(p)) { json_path <- p; break }
  }

  if (is.null(json_path)) {
    message("[SEGMENTATION] GeoJSON not found for case ", case_id,
            " (tried: ", paste(candidate_ids, collapse = ", "), ")")
    return(NULL)
  }

  message("[SEGMENTATION] Loading GeoJSON for case ", case_id, " from ", json_path)

  geojson <- tryCatch({
    gj <- sf::st_read(json_path, quiet = TRUE)
    # Remove geographic CRS - these are actually image pixel/micron coordinates, not geographic
    # The GeoJSON export from QuPath uses pixel coordinates but sf assigns WGS84 by default
    sf::st_crs(gj) <- NA
    gj
  }, error = function(e) {
    message("[SEGMENTATION] Error reading GeoJSON: ", conditionMessage(e))
    NULL
  })

  if (!is.null(geojson)) {
    assign(cache_key, geojson, envir = geojson_cache)
    message("[SEGMENTATION] Cached GeoJSON for case ", case_id, " (", nrow(geojson), " features)")
  }

  return(geojson)
}

# Extract classification name from GeoJSON properties
# The classification field may be a JSON string like '{ "name": "Nerve", "color": [...] }'
get_classification_name <- function(geojson) {
  if (is.null(geojson) || nrow(geojson) == 0) return(character(0))

  if (!"classification" %in% names(geojson)) {
    return(rep(NA_character_, nrow(geojson)))
  }

  cls <- geojson$classification

  # Handle character vector (JSON strings)
  if (is.character(cls)) {
    sapply(cls, function(x) {
      if (is.na(x) || !nzchar(x)) return(NA_character_)
      parsed <- tryCatch(jsonlite::fromJSON(x), error = function(e) NULL)
      if (!is.null(parsed) && "name" %in% names(parsed)) {
        as.character(parsed$name)
      } else {
        NA_character_
      }
    }, USE.NAMES = FALSE)
  } else if (is.list(cls)) {
    # Handle list of classification objects
    sapply(cls, function(x) {
      if (is.list(x) && "name" %in% names(x)) x$name else NA_character_
    }, USE.NAMES = FALSE)
  } else if (is.data.frame(cls) && "name" %in% names(cls)) {
    cls$name
  } else {
    rep(NA_character_, nrow(geojson))
  }
}

# Get polygons within a bounding box around a centroid
# Input coordinates are in micrometers (from islet_spatial_lookup.csv)
# GeoJSON polygons are in pixel coordinates, so we convert um -> pixels
get_islet_region_polygons <- function(geojson, centroid_x_um, centroid_y_um, buffer_um = 200) {
  if (!SF_AVAILABLE || is.null(geojson) || nrow(geojson) == 0) {
    return(list(Islet = NULL, IsletExpanded = NULL, Nerve = NULL, Capillary = NULL, Lymphatic = NULL))
  }

  # Convert micrometers to pixels for GeoJSON query
  centroid_x_px <- centroid_x_um / PIXEL_SIZE_UM
  centroid_y_px <- centroid_y_um / PIXEL_SIZE_UM
  buffer_px <- buffer_um / PIXEL_SIZE_UM

  # Define bounding box limits in pixel coordinates
  xmin <- centroid_x_px - buffer_px
  xmax <- centroid_x_px + buffer_px
  ymin <- centroid_y_px - buffer_px
  ymax <- centroid_y_px + buffer_px

  # Create a bounding box polygon for filtering
  # Use st_intersects with a bbox polygon instead of st_crop (which has CRS issues)
  bbox_polygon <- tryCatch({
    # Create a simple polygon from bbox coordinates
    coords <- matrix(c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ), ncol = 2, byrow = TRUE)
    bbox_sf <- sf::st_sfc(sf::st_polygon(list(coords)))
    # Match CRS if available, otherwise treat as planar
    if (!is.na(sf::st_crs(geojson))) {
      sf::st_set_crs(bbox_sf, sf::st_crs(geojson))
    }
    bbox_sf
  }, error = function(e) {
    message("[SEGMENTATION] Bbox creation error: ", conditionMessage(e))
    NULL
  })

  if (is.null(bbox_polygon)) {
    return(list(Islet = NULL, IsletExpanded = NULL, Nerve = NULL, Capillary = NULL, Lymphatic = NULL))
  }

  # Find features that intersect with the bounding box
  intersects <- tryCatch({
    suppressWarnings(sf::st_intersects(geojson, bbox_polygon, sparse = FALSE)[, 1])
  }, error = function(e) {
    message("[SEGMENTATION] Intersection error: ", conditionMessage(e))
    rep(FALSE, nrow(geojson))
  })

  cropped <- geojson[intersects, ]

  if (nrow(cropped) == 0) {
    return(list(Islet = NULL, IsletExpanded = NULL, Nerve = NULL, Capillary = NULL, Lymphatic = NULL))
  }

  # Get classification names
  cls_names <- get_classification_name(cropped)

  # Extract by classification
  result <- list()
  for (cls in c("Islet", "IsletExpanded", "Nerve", "Capillary", "Lymphatic")) {
    matches <- grepl(paste0("^", cls, "$"), cls_names, ignore.case = TRUE)
    if (any(matches)) {
      result[[cls]] <- cropped[matches, ]
    } else {
      result[[cls]] <- NULL
    }
  }

  return(result)
}

# Load segmentation annotations from data directory
load_segmentation_data <- function() {
  seg_path <- file.path("..", "..", "data", "app_data", "annotations.tsv")
  if (!file.exists(seg_path)) {
    cat("[SEGMENTATION] Annotations file not found at", seg_path, "\n")
    return(NULL)
  }

  tryCatch({
    df <- read.delim(seg_path, stringsAsFactors = FALSE)
    cat("[SEGMENTATION] Loaded", nrow(df), "segmentation records\n")

    # Validate required columns
    required_cols <- c("Image", "Name", "Class", "Centroid X µm", "Centroid Y µm")
    if (!all(required_cols %in% colnames(df))) {
      warning("[SEGMENTATION] Missing required columns: ",
              paste(setdiff(required_cols, colnames(df)), collapse = ", "))
      return(NULL)
    }

    # Filter for islet-related classes
    df <- df[df$Class %in% c("Islet", "ExpandedIslet", "Nerve", "Lymphatic", "Capillary"), ]
    cat("[SEGMENTATION] Found", nrow(df), "annotations across classes:",
        paste(unique(df$Class), collapse = ", "), "\n")

    df
  }, error = function(e) {
    cat("[SEGMENTATION] Error loading annotations:", e$message, "\n")
    NULL
  })
}

# Load spatial lookup for trajectory zoom-to-islet
load_islet_spatial_lookup <- function() {
  lookup_path <- file.path("..", "..", "data", "app_data", "islet_spatial_lookup.csv")
  if (!file.exists(lookup_path)) {
    cat("[SPATIAL] islet_spatial_lookup.csv not found at", lookup_path, "\n")
    return(NULL)
  }

  tryCatch({
    df <- read.csv(lookup_path, stringsAsFactors = FALSE)
    cat("[SPATIAL] Loaded", nrow(df), "islet spatial records\n")

    # Validate required columns
    required_cols <- c("case_id", "islet_key", "centroid_x_um", "centroid_y_um")
    if (!all(required_cols %in% colnames(df))) {
      warning("[SPATIAL] Missing required columns: ",
              paste(setdiff(required_cols, colnames(df)), collapse = ", "))
      return(NULL)
    }

    df
  }, error = function(e) {
    cat("[SPATIAL] Error loading spatial lookup:", e$message, "\n")
    NULL
  })
}

discover_islet_assets <- function(image_dir = file.path("www", "local_images"),
                                  seg_dir = file.path("www", "LabelExports")) {
  # Discover base images
  if (!dir.exists(image_dir)) return(data.frame())
  image_files <- list.files(image_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE,
                            full.names = TRUE)
  if (!length(image_files)) return(data.frame())

  # Helper to extract the 4-digit sample id from image filename: ..._YYYY.ome.tiff
  extract_image_id <- function(fn) {
    b <- basename(fn)
    id <- sub(".*_([0-9]{4})\\.ome\\.tiff?$", "\\1", b, perl = TRUE, ignore.case = TRUE)
    if (!nzchar(id) || is.na(suppressWarnings(as.integer(id)))) {
      # Fallback to any 3-5 digit run in the name
      id <- stringr::str_extract(b, "[0-9]{3,5}")
    }
    if (is.na(id)) id <- NA_character_
    id
  }

  # Discover segmentation files once
  seg_df <- data.frame()
  if (dir.exists(seg_dir)) {
    seg_files <- list.files(seg_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE,
                            full.names = TRUE)
    if (length(seg_files)) {
      seg_bases <- basename(seg_files)
      # Sample id is first 4 digits at start of filename; label is after the first underscore
      seg_sid <- sub("^([0-9]{4}).*$", "\\1", seg_bases, perl = TRUE)
      seg_label_raw <- sub("^[0-9]{4}_", "", seg_bases, perl = TRUE)
      seg_label_raw <- sub("\\.ome\\.tif{1,2}f?$", "", seg_label_raw, perl = TRUE, ignore.case = TRUE)
      seg_label_disp <- gsub("_", " ", seg_label_raw)
      seg_label_disp <- gsub("([a-z])([A-Z])", "\\1 \\2", seg_label_disp, perl = TRUE)

      seg_df <- data.frame(
        seg_abs = seg_files,
        seg_rel = file.path("LabelExports", seg_bases),
        seg_base = seg_bases,
        seg_sample_id = seg_sid,
        seg_islet = stringr::str_extract(seg_bases, "Islet_[0-9]+"),
        seg_label = seg_label_raw,
        seg_display = seg_label_disp,
        stringsAsFactors = FALSE
      )
    }
  }

  # Build asset rows; include list-columns of segmentation rel paths and labels
  res <- lapply(image_files, function(img_abs) {
    img_base <- basename(img_abs)
    img_sample_id <- extract_image_id(img_base)
    # Derive case id (zero-padded when numeric)
    img_case_raw <- stringr::str_extract(img_base, "[0-9]{3,5}")
    img_case_num <- suppressWarnings(as.integer(img_case_raw))
    img_case <- if (!is.na(img_case_num)) sprintf("%04d", img_case_num) else (img_sample_id %||% img_case_raw)
    img_islet <- stringr::str_extract(img_base, "Islet_[0-9]+")

    seg_rel_list <- character(0)
    seg_label_list <- character(0)
    if (nrow(seg_df)) {
      cand <- seg_df
      # First filter by the stricter 4-digit sample id if available
      if (!is.na(img_sample_id) && nzchar(img_sample_id)) {
        cand <- cand[cand$seg_sample_id == img_sample_id, , drop = FALSE]
      }
      # If islet token present, prefer those matching it
      if (nrow(cand) && !is.na(img_islet) && nzchar(img_islet)) {
        iso <- cand[!is.na(cand$seg_islet) & cand$seg_islet == img_islet, , drop = FALSE]
        if (nrow(iso)) cand <- iso
      }
      if (nrow(cand)) {
        seg_rel_list <- cand$seg_rel
        seg_label_list <- cand$seg_display
      }
    }

    # Back-compat single seg columns: choose first if any
    seg_abs_first <- if (length(seg_rel_list)) file.path(seg_dir, basename(seg_rel_list[[1]])) else NA_character_
    seg_rel_first <- if (length(seg_rel_list)) seg_rel_list[[1]] else NA_character_

    # Return asset row with list columns
    data.frame(
      case_id = if (nzchar(img_case)) img_case else NA_character_,
      image_abs = img_abs,
      image_rel = file.path("local_images", img_base),
      image_name = img_base,
      islet_token = if (nzchar(img_islet)) img_islet else NA_character_,
      seg_abs = seg_abs_first,
      seg_rel = seg_rel_first,
      stringsAsFactors = FALSE
    ) -> row

    # Attach list-cols
    row$seg_rel_list <- I(list(seg_rel_list))
    row$seg_label_list <- I(list(seg_label_list))
    row
  })

  assets <- do.call(rbind, res)
  unique(assets)
}

choose_islet_asset <- function(assets, case_id, islet_key = NULL) {
  if (is.null(assets) || !nrow(assets)) return(NULL)
  if (is.null(case_id) || !nzchar(case_id)) return(assets[1, , drop = FALSE])

  case_vec <- unique(c(case_id, suppressWarnings(sprintf("%04d", as.integer(case_id)))))
  case_vec <- case_vec[nzchar(case_vec)]
  subset_assets <- assets[assets$case_id %in% case_vec, , drop = FALSE]

  if (!nrow(subset_assets)) {
    subset_assets <- assets[grepl(case_id, assets$image_name, fixed = TRUE), , drop = FALSE]
  }

  if (!is.null(islet_key) && nzchar(islet_key) && nrow(subset_assets)) {
    match_islet <- subset_assets[!is.na(subset_assets$islet_token) & subset_assets$islet_token == islet_key, , drop = FALSE]
    if (nrow(match_islet)) subset_assets <- match_islet
  }

  if (!nrow(subset_assets) && !is.null(islet_key) && nzchar(islet_key)) {
    subset_assets <- assets[grepl(islet_key, assets$image_name, fixed = TRUE) |
                              grepl(islet_key, assets$seg_rel, fixed = TRUE), , drop = FALSE]
  }

  if (!nrow(subset_assets)) return(NULL)
  subset_assets[1, , drop = FALSE]
}

# Improved helper function to get annotation info for an islet with class information
get_islet_annotations <- function(case_id, islet_key) {
  # Try spatial lookup first (faster and more reliable)
  if (!is.null(islet_spatial_lookup)) {
    case_id_char <- as.character(case_id)
    matches <- islet_spatial_lookup[
      as.character(islet_spatial_lookup$case_id) == case_id_char &
      islet_spatial_lookup$islet_key == islet_key,
    ]

    if (nrow(matches) > 0) {
      annotations <- list(
        centroid_x = as.numeric(matches$centroid_x_um[1]),
        centroid_y = as.numeric(matches$centroid_y_um[1]),
        area = as.numeric(matches$area_um2[1]),
        source = "spatial_lookup"
      )
      return(annotations)
    }
  }

  # Fallback to segmentation_data (lazy load on first access)
  if (!.seg_lazy$loaded) {
    .seg_lazy$data <- load_segmentation_data()
    .seg_lazy$loaded <- TRUE
  }
  segmentation_data <- .seg_lazy$data
  if (is.null(segmentation_data)) return(NULL)

  # Extract numeric islet ID from islet_key
  islet_numeric <- as.numeric(gsub(".*_(\\d+).*", "\\1", islet_key))
  if (is.na(islet_numeric)) return(NULL)

  # Convert case_id to character for matching
  case_id_char <- as.character(case_id)

  # Look for matching annotations across all classes
  matches <- segmentation_data[
    as.character(segmentation_data$Image) == case_id_char &
    grepl(paste0("Islet_", islet_numeric), segmentation_data$Name, ignore.case = TRUE),
  ]

  if (nrow(matches) == 0) return(NULL)

  # Safe numeric conversion
  safe_numeric <- function(x) {
    result <- suppressWarnings(as.numeric(as.character(x)))
    if (is.na(result)) 0 else result
  }

  # Return all matching annotations by class
  annotations <- list(
    centroid_x = safe_numeric(matches$`Centroid X µm`[1]),
    centroid_y = safe_numeric(matches$`Centroid Y µm`[1]),
    area = safe_numeric(matches$`Area µm^2`[1]),
    perimeter = safe_numeric(matches$`Perimeter µm`[1]),
    classes = unique(matches$Class),
    source = "annotations"
  )

  # Add class-specific counts if available
  for (cls in c("Islet", "ExpandedIslet", "Nerve", "Lymphatic", "Capillary")) {
    cls_matches <- matches[matches$Class == cls, ]
    if (nrow(cls_matches) > 0) {
      annotations[[paste0("has_", tolower(cls))]] <- TRUE
      annotations[[paste0(tolower(cls), "_count")]] <- nrow(cls_matches)
    }
  }

  annotations
}

# Deferred: annotations.tsv is 72 MB and only needed as a fallback
# when islet_spatial_lookup doesn't have the islet (never happens in practice)
.seg_lazy <- new.env(parent = emptyenv())
.seg_lazy$data <- NULL
.seg_lazy$loaded <- FALSE

# Load spatial lookup eagerly (small file, always needed)
islet_spatial_lookup <- load_islet_spatial_lookup()

# Resolve a clicked islet (case_id + islet_key) to the selected_islet payload
# and set it. Shared tail of the Plot and Trajectory click handlers (both the
# legacy plotly_click observers and the native-click bridge) so the centroid
# lookup, zero-padded case-id fallback, and donor-metadata join live in one
# place. `prepared` is the core data reactive; `selected_islet` is the shared
# reactiveVal. Returns TRUE on success; shows a notification and returns FALSE
# otherwise (suppress with notify = FALSE).
select_islet_from_key <- function(case_id, islet_key, prepared, selected_islet, notify = TRUE) {
  note <- function(msg, type = "warning", duration = 3) {
    if (notify) showNotification(msg, type = type, duration = duration)
  }
  if (!SF_AVAILABLE) {
    note("Segmentation requires the 'sf' package", "warning", 5); return(FALSE)
  }
  if (is.null(islet_spatial_lookup)) {
    note("Islet spatial lookup data not available", "warning", 5); return(FALSE)
  }
  if (is.null(case_id) || is.null(islet_key) ||
      is.na(case_id) || is.na(islet_key) ||
      !nzchar(as.character(case_id)) || !nzchar(as.character(islet_key))) {
    note("Could not identify clicked islet"); return(FALSE)
  }
  case_id <- as.character(case_id)
  islet_key <- as.character(islet_key)

  spatial_match <- islet_spatial_lookup[
    islet_spatial_lookup$case_id == case_id &
      islet_spatial_lookup$islet_key == islet_key, ]

  # Try zero-padded variant (data uses "112"; GeoJSON/lookup may use "0112")
  if (nrow(spatial_match) == 0) {
    case_id_num <- suppressWarnings(as.integer(case_id))
    if (!is.na(case_id_num)) {
      case_id_padded <- sprintf("%04d", case_id_num)
      spatial_match <- islet_spatial_lookup[
        islet_spatial_lookup$case_id == case_id_padded &
          islet_spatial_lookup$islet_key == islet_key, ]
      if (nrow(spatial_match) > 0) case_id <- case_id_padded
    }
  }

  if (nrow(spatial_match) == 0) {
    note(paste("Spatial data not found for", islet_key, "in case", case_id)); return(FALSE)
  }

  centroid_x <- spatial_match$centroid_x_um[1]
  centroid_y <- spatial_match$centroid_y_um[1]

  # Donor metadata from prepared()$comp
  pd <- prepared()
  donor_status_val <- NULL; donor_age <- NULL; donor_gender <- NULL
  if (!is.null(pd$comp) && "Case ID" %in% colnames(pd$comp)) {
    donor_row <- pd$comp[pd$comp$`Case ID` == suppressWarnings(as.integer(case_id)), , drop = FALSE]
    if (nrow(donor_row) > 0) {
      if ("Donor Status" %in% colnames(donor_row)) donor_status_val <- as.character(donor_row$`Donor Status`[1])
      if ("age" %in% colnames(donor_row)) donor_age <- donor_row$age[1]
      if ("gender" %in% colnames(donor_row)) donor_gender <- as.character(donor_row$gender[1])
    }
  }

  selected_islet(list(
    case_id      = case_id,
    islet_key    = islet_key,
    centroid_x   = centroid_x,
    centroid_y   = centroid_y,
    donor_status = donor_status_val,
    donor_age    = donor_age,
    donor_gender = donor_gender
  ))
  TRUE
}

# Reusable GeoJSON base plot builder
# Returns a ggplot with polygon layers + coord_sf (no crosshairs or title).
# Used by both render_islet_segmentation_plot() and render_islet_drilldown_plot().
build_segmentation_base_plot <- function(info, buffer_um = 150,
                                         show_peri_boundary = TRUE,
                                         show_structures = TRUE) {
  if (is.null(info)) return(NULL)

  geojson <- load_case_geojson(info$case_id)
  if (is.null(geojson)) return(NULL)

  polygons <- get_islet_region_polygons(geojson, info$centroid_x, info$centroid_y, buffer_um)

  centroid_x_px <- info$centroid_x / PIXEL_SIZE_UM
  centroid_y_px <- info$centroid_y / PIXEL_SIZE_UM
  buffer_px <- buffer_um / PIXEL_SIZE_UM

  xlim_range <- c(centroid_x_px - buffer_px, centroid_x_px + buffer_px)
  ylim_range <- c(centroid_y_px - buffer_px, centroid_y_px + buffer_px)

  colors <- c(
    "Islet" = "#0066CC", "IsletExpanded" = "#00CCCC",
    "Nerve" = "#CC00CC", "Capillary" = "#CC0000", "Lymphatic" = "#00AA00"
  )

  # Determine which layer classes to skip based on toggle params
  skip_classes <- character(0)
  if (!show_peri_boundary) skip_classes <- c(skip_classes, "IsletExpanded")
  if (!show_structures) skip_classes <- c(skip_classes, "Nerve", "Capillary", "Lymphatic")

  p <- ggplot2::ggplot() +
    theme_islet() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = ISLET_FONT_SIZES$title, face = "bold")
    ) +
    ggplot2::coord_sf(xlim = xlim_range, ylim = rev(ylim_range), expand = FALSE)

  # Filter Islet and IsletExpanded polygons to only the selected islet
  # (the bounding box query picks up neighboring islets)
  if (!is.null(polygons$Islet) && nrow(polygons$Islet) > 0 && !is.null(info$islet_key)) {
    islet_names <- if ("name" %in% names(polygons$Islet)) polygons$Islet$name
                   else if ("id" %in% names(polygons$Islet)) polygons$Islet$id
                   else NULL
    if (!is.null(islet_names)) {
      match_idx <- which(islet_names == info$islet_key)
      if (length(match_idx) == 0) match_idx <- grep(info$islet_key, islet_names, fixed = TRUE)
      if (length(match_idx) > 0) {
        polygons$Islet <- polygons$Islet[match_idx[1], , drop = FALSE]
      }
    }
  }
  if (!is.null(polygons$IsletExpanded) && nrow(polygons$IsletExpanded) > 0 && !is.null(info$islet_key)) {
    exp_names <- if ("name" %in% names(polygons$IsletExpanded)) polygons$IsletExpanded$name
                 else if ("id" %in% names(polygons$IsletExpanded)) polygons$IsletExpanded$id
                 else NULL
    if (!is.null(exp_names)) {
      exp_key <- paste0(info$islet_key, "_exp20um")
      match_idx <- which(exp_names == exp_key)
      if (length(match_idx) == 0) match_idx <- grep(info$islet_key, exp_names, fixed = TRUE)
      if (length(match_idx) > 0) {
        polygons$IsletExpanded <- polygons$IsletExpanded[match_idx[1], , drop = FALSE]
      }
    }
  }

  # Display names for legend
  display_names <- c(
    "Islet" = "Islet", "IsletExpanded" = "Peri Boundary",
    "Nerve" = "Nerve", "Capillary" = "Capillary", "Lymphatic" = "Lymphatic"
  )
  display_colors <- c(
    "Islet" = "#0066CC", "Peri Boundary" = "#00CCCC",
    "Nerve" = "#CC00CC", "Capillary" = "#CC0000", "Lymphatic" = "#00AA00"
  )

  # Combine all visible polygons into one sf with a structure_type column
  layer_order <- c("IsletExpanded", "Islet", "Lymphatic", "Capillary", "Nerve")
  poly_list <- list()
  for (cls in layer_order) {
    if (cls %in% skip_classes) next
    if (!is.null(polygons[[cls]]) && nrow(polygons[[cls]]) > 0) {
      pg <- sf::st_sf(
        structure_type = rep(display_names[cls], nrow(polygons[[cls]])),
        geometry = sf::st_geometry(polygons[[cls]])
      )
      poly_list[[cls]] <- pg
    }
  }

  if (length(poly_list) == 0) return(NULL)

  combined_sf <- do.call(rbind, poly_list)
  # Order factor so legend matches layer order
  visible_names <- display_names[layer_order[!layer_order %in% skip_classes]]
  visible_names <- visible_names[visible_names %in% unique(combined_sf$structure_type)]
  combined_sf$structure_type <- factor(combined_sf$structure_type, levels = visible_names)

  p <- p +
    ggplot2::geom_sf(data = combined_sf,
                     ggplot2::aes(colour = structure_type),
                     fill = NA, linewidth = 0.8, key_glyph = "path") +
    ggplot2::scale_colour_manual(values = display_colors, name = "Structures",
                                  drop = FALSE)

  p
}

# Reusable segmentation plot renderer
# Called by both Trajectory and Plot modules for their respective segmentation outputs
render_islet_segmentation_plot <- function(info, show_peri_boundary = TRUE, show_structures = TRUE) {
  if (is.null(info)) return(NULL)
  cat("[SEGMENTATION RENDER] Rendering for case=", info$case_id, ", islet=", info$islet_key, "\n")

  p <- build_segmentation_base_plot(info, show_peri_boundary = show_peri_boundary,
                                     show_structures = show_structures)
  if (is.null(p)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste("No segmentation data found for", info$islet_key,
                                        "in case", info$case_id),
                          size = 5, color = "gray50") +
        ggplot2::theme_void() +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
    )
  }

  centroid_x_px <- info$centroid_x / PIXEL_SIZE_UM
  centroid_y_px <- info$centroid_y / PIXEL_SIZE_UM

  # Add crosshairs and title (not in base plot)
  p <- p +
    ggplot2::geom_vline(xintercept = centroid_x_px, color = "gray50", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = centroid_y_px, color = "gray50", linetype = "dashed", linewidth = 0.3) +
    ggplot2::labs(
      title = paste(info$islet_key, "- Case", info$case_id),
      x = "X position (pixels)", y = "Y position (pixels)"
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "cm")
    )

  p
}
