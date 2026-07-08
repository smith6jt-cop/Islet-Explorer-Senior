library(shiny)
library(shinyjs)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(plotly)
library(broom)
library(jsonlite)

# Pixel size constant for coordinate conversion (micrometers per pixel)
# GeoJSON polygons use pixel coordinates, while islet_spatial_lookup.csv uses micrometers
PIXEL_SIZE_UM <- 0.3774  # micrometers per pixel

# Load sf package for GeoJSON/spatial operations (islet segmentation viewer)
SF_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (requireNamespace("sf", quietly = TRUE)) {
    library(sf)
    SF_AVAILABLE <- TRUE
    message("[SF] Package sf loaded for segmentation viewer")
  } else {
    message("[SF] Package 'sf' not found. Islet segmentation viewer will be disabled.")
    message("[SF] Install with: install.packages('sf')")
  }
})

# Check for httr2 package (required for AI assistant)
httr2_available <- requireNamespace("httr2", quietly = TRUE)
if (!httr2_available) {
  message("Package 'httr2' not found. The AI assistant will be disabled until you install it (install.packages('httr2')).")
}

# Load anndata package for trajectory analysis (install with BiocManager::install('anndata'))
suppressPackageStartupMessages({
  if (!requireNamespace("anndata", quietly = TRUE)) {
    message("Package 'anndata' not found. Install with: BiocManager::install('anndata')")
  }
})

# Load vitessceR for embedded image viewer (GitHub v0.1.0 required)
VITESSCE_AVAILABLE <- FALSE
suppressPackageStartupMessages({
  if (!requireNamespace("vitessceR", quietly = TRUE)) {
    message("Package 'vitessceR' not found. Install with: remotes::install_github('vitessce/vitessceR')")
    message("Version required: >= 0.1.0")
  } else {
    # Verify we have the GitHub version with required features
    tryCatch({
      pkg_version <- packageVersion("vitessceR")
      if (pkg_version < "0.1.0") {
        warning("vitessceR version ", pkg_version, " detected. Please upgrade to >= 0.1.0 from GitHub")
      }
      test_vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.16", name = "Test")
      message("[VITESSCE] Package verified working (version ", pkg_version, ")")
      VITESSCE_AVAILABLE <<- TRUE
    }, error = function(e) {
      warning("[VITESSCE] Package installed but failed basic test: ", conditionMessage(e))
    })

    # CRITICAL FIX: Ensure all child routes from MultiImageWrapper are registered
    # This fixes segmentation layers not loading properly and appearing at wrong sizes
    if (VITESSCE_AVAILABLE) {
      tryCatch({
        vitessceR::MultiImageWrapper$set("public", "get_routes", function() {
          routes <- list()
          if (length(self$image_wrappers) == 0L) return(routes)
          for (w in self$image_wrappers) {
            wr <- w$get_routes()
            if (length(wr) > 0L) routes <- c(routes, wr)
          }
          routes
        }, overwrite = TRUE)
        message("[VITESSCE] Applied MultiImageWrapper get_routes fix for segmentation layer routing")
      }, error = function(e) {
        warning("[VITESSCE] Failed to apply get_routes fix: ", conditionMessage(e))
      })
    }
  }
})

# ---- Runtime diagnostics (temporary) ----
APP_VERSION <- "vitessce-trajectory-v1"
try({
  message("[APP LOAD] Version=", APP_VERSION, " time=", Sys.time())
  message("[APP LOAD] getwd=", tryCatch(getwd(), error=function(e) NA))
  message("[APP LOAD] list.files(.) contains app.R? ", any(grepl('^app.R$', list.files('.'))))
  message("[APP LOAD] file.info(app.R)$mtime=", tryCatch(file.info('app.R')$mtime, error=function(e) NA))
  # Create a flag file to confirm this specific app.R executed.
  flag_path <- file.path(dirname(sys.frame(1)$ofile %||% 'app.R'), 'APP_LOADED_FLAG')
  writeLines(paste0('loaded ', Sys.time(), ' version=', APP_VERSION), flag_path)
  # Disable tracing that was causing issues with trajectory data
  # if (!isTRUE(getOption('traj_left_join_traced'))) {
  #   try(trace(dplyr::left_join,
  #             quote({
  #               message(sprintf('[TRACE left_join] classes lhs=%s rhs=%s by=%s',
  #                                paste(class(x), collapse='/'),
  #                                paste(class(y), collapse='/'),
  #                                paste(names(by %||% list()), collapse=',')))
  #             }), print = FALSE), silent = TRUE)
  #   options(traj_left_join_traced = TRUE)
  # }
}, silent = TRUE)


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
## for base64 encoding channel_config payload to Avivator
## (installed by scripts/install_shiny_deps.R)
suppressPackageStartupMessages({
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Run scripts/install_shiny_deps.R to install dependencies.")
  }
})


# ===== Islet Segmentation Viewer Helper Functions =====
# Cache for loaded GeoJSON data (persists across sessions for performance)
geojson_cache <- new.env()

# Load islet spatial lookup table (centroids and coordinates)
islet_spatial_lookup <- tryCatch({
  lookup_path <- file.path("..", "..", "data", "islet_spatial_lookup.csv")
  if (file.exists(lookup_path)) {
    df <- read.csv(lookup_path, stringsAsFactors = FALSE)
    message("[SEGMENTATION] Loaded islet_spatial_lookup.csv with ", nrow(df), " entries")
    df
  } else {
    message("[SEGMENTATION] islet_spatial_lookup.csv not found at ", lookup_path)
    NULL
  }
}, error = function(e) {
  message("[SEGMENTATION] Error loading islet_spatial_lookup.csv: ", conditionMessage(e))
  NULL
})

# Load GeoJSON for a specific case_id with caching
load_case_geojson <- function(case_id) {
  if (!SF_AVAILABLE) return(NULL)

  cache_key <- as.character(case_id)
  if (exists(cache_key, envir = geojson_cache)) {
    return(get(cache_key, envir = geojson_cache))
  }

  # Try uncompressed JSON first
  json_path <- file.path("..", "..", "data", "json", paste0(case_id, ".geojson"))
  if (!file.exists(json_path)) {
    # Try compressed version
    json_path <- file.path("..", "..", "data", "gson", paste0(case_id, ".geojson.gz"))
  }

  if (!file.exists(json_path)) {
    message("[SEGMENTATION] GeoJSON not found for case ", case_id)
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
# GeoJSON polygons are in pixel coordinates, so we convert µm → pixels
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

master_path <- file.path("..", "..", "data", "master_results.xlsx")

project_root <- tryCatch(normalizePath(file.path("..", ".."), mustWork = FALSE), error = function(e) NULL)

# Ensure reticulate/anndata can discover a Python binary when RETICULATE_PYTHON
# isn't set (shiny-server environments often lack PATH). This happens before
# any anndata import so reticulate locks onto the path once.
if (!nzchar(Sys.getenv("RETICULATE_PYTHON", ""))) {
  candidate_python <- NULL
  venv_path <- file.path(path.expand("~"), ".local", "share", "islet-explorer-py", "bin", "python")
  if (file.exists(venv_path)) candidate_python <- venv_path
  if (is.null(candidate_python) || !nzchar(candidate_python)) {
    candidate_python <- Sys.getenv("PYTHON", Sys.which("python3"))
    if (!nzchar(candidate_python)) candidate_python <- Sys.which("python")
  }
  if (nzchar(candidate_python)) {
    Sys.setenv(RETICULATE_PYTHON = candidate_python)
    message("[PYTHON] RETICULATE_PYTHON=", candidate_python)
  } else {
    message("[PYTHON] No python executable found. Set RETICULATE_PYTHON to a valid interpreter.")
  }
}

AI_SYSTEM_PROMPT <- paste(
  "You are an AI assistant embedded in the Islet Explorer Shiny app.",
  "Help users interpret plots, statistics, and data preparation steps without hallucinating.",
  "Favor concise, actionable answers rooted in the app's current state and available controls."
)

# ===== Credential Loading =====
# Set DEBUG_CREDENTIALS=1 to enable verbose logging
DEBUG_CREDS <- isTRUE(as.logical(Sys.getenv("DEBUG_CREDENTIALS", "0")))
VIEWER_DEBUG_ENABLED <- identical(Sys.getenv("VIEWER_DEBUG", "0"), "1")

if (DEBUG_CREDS) {
  cat("\n=== CREDENTIAL LOADING DEBUG ===\n")
  cat("Current working directory:", getwd(), "\n")
  cat("Home directory:", Sys.getenv("HOME"), "\n")
  cat("User:", Sys.getenv("USER"), "\n")
}

# Check if .Renviron exists - prioritize app directory for deployment
renviron_paths <- c(
  ".Renviron",  # App directory (for deployment)
  file.path(dirname(sys.frame(1)$ofile %||% "."), ".Renviron"),  # Script directory
  Sys.getenv("R_ENVIRON_USER", unset = ""),
  "~/.Renviron",
  "~/.Renviron.local",
  file.path(Sys.getenv("HOME"), ".Renviron")
)
renviron_paths <- unique(path.expand(renviron_paths[nzchar(renviron_paths)]))

if (DEBUG_CREDS) {
  cat("\nSearching for .Renviron files:\n")
  for (path in renviron_paths) {
    exists <- file.exists(path)
    cat("  ", path, ":", if(exists) "EXISTS" else "NOT FOUND", "\n")
    if (exists) {
      cat("    File size:", file.size(path), "bytes\n")
      cat("    Readable:", file.access(path, 4) == 0, "\n")
    }
  }
}

# Credential loader with optional debug output
load_user_env_credentials <- function() {
  if (DEBUG_CREDS) cat("\n[CREDENTIAL LOADER] Starting...\n")

  for (env_path in renviron_paths) {
    if (!file.exists(env_path)) next

    if (DEBUG_CREDS) cat("[CREDENTIAL LOADER] Attempting to load:", env_path, "\n")

    result <- tryCatch({
      readRenviron(env_path)

      if (DEBUG_CREDS) {
        lines <- readLines(env_path, warn = FALSE)
        key_line <- grep("^KEY\\s*=", lines, value = TRUE)
        base_line <- grep("^BASE\\s*=", lines, value = TRUE)
        cat("[CREDENTIAL LOADER] File contains", length(lines), "lines\n")
        cat("[CREDENTIAL LOADER] KEY line found:", length(key_line) > 0, "\n")
        cat("[CREDENTIAL LOADER] BASE line found:", length(base_line) > 0, "\n")

        key_after <- Sys.getenv("KEY", unset = "")
        base_after <- Sys.getenv("BASE", unset = "")
        cat("[CREDENTIAL LOADER] After loading:\n")
        cat("  KEY present:", nzchar(key_after), "\n")
        cat("  KEY length:", nchar(key_after), "\n")
        cat("  BASE present:", nzchar(base_after), "\n")
        cat("  BASE value:", base_after, "\n")
      }

      TRUE
    }, error = function(e) {
      if (DEBUG_CREDS) cat("[CREDENTIAL LOADER] ERROR loading", env_path, ":", conditionMessage(e), "\n")
      FALSE
    })

    if (result) {
      if (DEBUG_CREDS) cat("[CREDENTIAL LOADER] Successfully loaded from:", env_path, "\n")
      return(invisible(TRUE))
    }
  }

  if (DEBUG_CREDS) cat("[CREDENTIAL LOADER] No valid .Renviron file found\n")
  invisible(FALSE)
}

# Load credentials immediately
load_user_env_credentials()

# Simplified credential getters
get_llm_api_key <- function() {
  key_val <- Sys.getenv("KEY", unset = "")
  key_val <- stringr::str_trim(key_val)
  key_val <- stringr::str_replace_all(key_val, "\\s+", "")

  if (nzchar(key_val)) {
    if (DEBUG_CREDS) cat("[API KEY] Found key of length", nchar(key_val), "\n")
    return(key_val)
  }

  if (DEBUG_CREDS) cat("[API KEY] No key found in environment\n")
  return("")
}

get_llm_api_base <- function() {
  base_val <- Sys.getenv("BASE", unset = "")
  base_val <- stringr::str_trim(base_val)

  if (!nzchar(base_val)) {
    if (DEBUG_CREDS) cat("[API BASE] No BASE found, using default OpenAI endpoint\n")
    return("https://api.openai.com/v1")
  }

  # Remove trailing slashes
  base_val <- sub("/+\\z", "", base_val, perl = TRUE)

  # Add /v1 suffix if not present
  if (!grepl("/v[0-9]+$", base_val) && !grepl("/v[0-9]+/", base_val)) {
    base_val <- paste0(base_val, "/v1")
  }

  if (DEBUG_CREDS) cat("[API BASE] Using custom endpoint:", base_val, "\n")
  Sys.setenv(BASE = base_val)
  return(base_val)
}

# Initialize and verify
LLM_API_BASE <- get_llm_api_base()
test_key <- get_llm_api_key()

if (DEBUG_CREDS) {
  cat("\n=== FINAL CREDENTIAL CHECK ===\n")
  cat("API Base:", LLM_API_BASE, "\n")
  cat("API Key configured:", nzchar(test_key), "\n")
  cat("================================\n\n")
}

sanitize_openai_content <- function(text, max_chars = 4000) {
  if (is.null(text)) return("")
  trimmed <- stringr::str_trim(as.character(text))
  if (!nzchar(trimmed)) return("")
  if (nchar(trimmed, type = "chars") > max_chars) {
    return(paste0(substr(trimmed, 1, max_chars - 1), "…"))
  }
  trimmed
}

build_openai_messages <- function(history, system_prompt = AI_SYSTEM_PROMPT, max_turns = 6) {
  msgs <- list(list(role = "system", content = sanitize_openai_content(system_prompt)))
  if (!length(history)) return(msgs)
  filtered <- Filter(function(entry) entry$role %in% c("user", "assistant"), history)
  if (length(filtered) > max_turns * 2) {
    filtered <- filtered[(length(filtered) - max_turns * 2 + 1):length(filtered)]
  }
  for (entry in filtered) {
    sanitized <- sanitize_openai_content(entry$content)
    if (!nzchar(sanitized)) next
    msgs[[length(msgs) + 1]] <- list(role = entry$role, content = sanitized)
  }
  msgs
}

select_openai_models <- function(requested_model = "auto", history) {
  default_model <- Sys.getenv("OPENAI_DEFAULT_MODEL", unset = "gpt-oss-120b")
  fast_model <- Sys.getenv("OPENAI_FAST_MODEL", unset = "gpt-oss-20b")

  if (!identical(requested_model, "auto")) {
    return(unique(c(requested_model, default_model)))
  }

  last_user <- NULL
  if (length(history)) {
    user_entries <- rev(Filter(function(entry) identical(entry$role, "user"), history))
    if (length(user_entries)) {
      last_user <- user_entries[[1]]$content %||% ""
    }
  }

  convo_len <- sum(vapply(history, function(entry) !identical(entry$role, "system"), logical(1)))
  use_fast <- !is.null(last_user) && nchar(last_user, type = "chars") <= 600 && convo_len <= 6
  if (use_fast && nzchar(fast_model)) {
    return(unique(c(fast_model, default_model)))
  }
  unique(c(default_model))
}

call_openai_chat <- function(history, system_prompt = AI_SYSTEM_PROMPT,
                             model = "auto", temperature = 0.3,
                             max_output_tokens = 512,
                             stream = FALSE, stream_callback = NULL) {
  api_key <- get_llm_api_key()
  if (!nzchar(api_key)) {
    stop("OpenAI API key not configured.", call. = FALSE)
  }
  if (!httr2_available) {
    stop("OpenAI assistant requires the 'httr2' package. Install with install.packages('httr2').", call. = FALSE)
  }

  max_output_tokens <- suppressWarnings(as.integer(max_output_tokens))
  if (is.na(max_output_tokens) || max_output_tokens <= 0L) {
    max_output_tokens <- 512L
  } else if (max_output_tokens < 16L) {
    max_output_tokens <- 16L
  }

  messages <- build_openai_messages(history, system_prompt)
  input_messages <- lapply(messages, function(entry) {
    list(
      role = entry$role,
      content = list(list(type = "input_text", text = entry$content %||% ""))
    )
  })
  chat_messages <- lapply(messages, function(entry) {
    list(
      role = entry$role,
      content = entry$content %||% ""
    )
  })

  model_candidates <- select_openai_models(model, messages)
  fallback_model <- utils::tail(model_candidates, 1)
  last_error_message <- NULL

  base_candidate <- tryCatch(LLM_API_BASE, error = function(e) NULL)
  if (is.null(base_candidate) || !nzchar(base_candidate)) {
    base_candidate <- "https://api.openai.com/v1"
  }
  base_candidate <- sub("/+\\z", "", base_candidate, perl = TRUE)
  
  # UF Navigator uses standard OpenAI-compatible chat API, not /responses
  use_chat_api_only <- grepl("api\\.ai\\.it\\.ufl\\.edu", base_candidate, ignore.case = TRUE)
  
  responses_url <- paste0(base_candidate, "/responses")
  chat_url <- paste0(base_candidate, "/chat/completions")

  parse_usage <- function(body) {
    usage <- body$usage %||% NULL
    if (is.null(usage)) return(NULL)
    usage
  }

  stream_error <- function(message, fallback = FALSE) {
    err <- simpleError(message)
    attr(err, "fallback_to_nonstream") <- fallback
    err
  }

  try_stream_chat <- function(mdl) {
    message("[OPENAI] Streaming model '", mdl, "' via ", chat_url)
    chat_payload <- list(
      model = mdl,
      messages = chat_messages,
      temperature = temperature,
      max_tokens = max_output_tokens,
      stream = TRUE
    )
  payload_json <- jsonlite::toJSON(chat_payload, auto_unbox = TRUE, digits = NA)
  payload_raw <- charToRaw(payload_json)
    chunk_buffer <- ""
    accumulated <- ""
    usage_capture <- NULL

    append_segment <- function(segment) {
      lines <- strsplit(segment, "\n", fixed = TRUE)[[1]]
      if (!length(lines)) return()
      data_lines <- lines[startsWith(lines, "data:")]
      if (!length(data_lines)) return()
      for (line in data_lines) {
        data <- sub("^data:\\s*", "", line)
        if (identical(data, "[DONE]")) next
        if (!nzchar(data)) next
        parsed <- tryCatch(jsonlite::fromJSON(data, simplifyVector = FALSE), error = function(e) NULL)
        if (is.null(parsed)) next
        if (!is.null(parsed$usage)) {
          usage_capture <<- parsed$usage
          if (is.function(stream_callback)) {
            stream_callback(accumulated, usage_capture)
          }
        }
        if (!is.null(parsed$choices) && length(parsed$choices) > 0) {
          delta <- parsed$choices[[1]]$delta %||% list()
          piece <- delta$content
          if (!is.null(piece)) {
            piece <- paste(piece, collapse = "")
            accumulated <<- paste0(accumulated, piece)
            if (is.function(stream_callback)) {
              stream_callback(accumulated, usage_capture)
            }
          } else if (!is.null(delta$role) && is.function(stream_callback)) {
            stream_callback(accumulated, usage_capture)
          }
        }
      }
    }

    req_chat_stream <- httr2::request(chat_url) |>
      httr2::req_headers(
        Authorization = paste("Bearer", api_key),
        `Content-Type` = "application/json"
      ) |>
  httr2::req_body_raw(payload_raw, type = "application/json") |>
      httr2::req_timeout(60) |>
      httr2::req_retry(
        max_tries = 3,
        is_transient = function(resp) {
          status <- httr2::resp_status(resp)
          status == 429 || (status >= 500 && status < 600)
        }
      )

    resp_stream <- tryCatch(
      httr2::req_stream(req_chat_stream, function(chunk, ...) {
        chunk_buffer <<- paste0(chunk_buffer, rawToChar(chunk))
        repeat {
          split_pos <- regexpr("\n\n", chunk_buffer, fixed = TRUE)
          if (split_pos[1] == -1) break
          segment <- substr(chunk_buffer, 1, split_pos[1] - 1)
          chunk_buffer <<- substr(chunk_buffer, split_pos[1] + 2, nchar(chunk_buffer))
          append_segment(segment)
        }
        invisible(NULL)
      }),
      error = identity
    )

    if (inherits(resp_stream, "error")) {
      attr(resp_stream, "fallback_to_nonstream") <- TRUE
      return(resp_stream)
    }

    if (nzchar(chunk_buffer)) {
      append_segment(chunk_buffer)
      chunk_buffer <<- ""
    }

    status_chat <- httr2::resp_status(resp_stream)
    if (status_chat >= 300) {
      err_body_chat <- tryCatch(httr2::resp_body_json(resp_stream, simplifyVector = TRUE), error = function(e) NULL)
      status_reason_chat <- tryCatch(httr2::resp_status_desc(resp_stream), error = function(e) NULL)
      detail_chat <- if (!is.null(err_body_chat) && !is.null(err_body_chat$error) && !is.null(err_body_chat$error$message)) {
        err_body_chat$error$message
      } else if (!is.null(status_reason_chat) && nzchar(status_reason_chat)) {
        paste(status_chat, status_reason_chat)
      } else {
        paste("HTTP", status_chat)
      }
      return(stream_error(detail_chat, fallback = status_chat %in% c(404, 405)))
    }

    if (!nzchar(accumulated)) {
      return(stream_error("Empty completion returned from chat API.", fallback = TRUE))
    }

    list(
      text = stringr::str_trim(accumulated),
      usage = usage_capture
    )
  }

  for (mdl in model_candidates) {
    if (isTRUE(stream)) {
      stream_result <- try_stream_chat(mdl)
      if (!inherits(stream_result, "error")) {
        return(stream_result)
      }
      last_error_message <- conditionMessage(stream_result)
      if (!isTRUE(attr(stream_result, "fallback_to_nonstream"))) {
        next
      }
    }

    # For UF Navigator, skip /responses and use /chat/completions directly
    use_chat_api <- use_chat_api_only

    if (!use_chat_api) {
      payload <- list(
        model = mdl,
        input = input_messages,
        temperature = temperature,
        max_output_tokens = max_output_tokens,
        metadata = list(app = "Islet Explorer AI Assistant")
      )

      message("[OPENAI] Attempting model '", mdl, "' via ", responses_url)
      req <- httr2::request(responses_url) |>
        httr2::req_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ) |>
        httr2::req_body_json(payload, auto_unbox = TRUE) |>
        httr2::req_timeout(30) |>
        httr2::req_retry(
          max_tries = 3,
          is_transient = function(resp) {
            status <- httr2::resp_status(resp)
            status == 429 || (status >= 500 && status < 600)
          }
        )

      resp <- tryCatch(httr2::req_perform(req), error = identity)

      if (inherits(resp, "httr2_http")) {
        last_error_message <- conditionMessage(resp)
        if (!is.null(resp$resp)) {
          resp <- resp$resp
        } else {
          next
        }
      } else if (inherits(resp, "error")) {
        last_error_message <- conditionMessage(resp)
        next
      }

      status_code <- httr2::resp_status(resp)
      if (status_code >= 300) {
        err_body <- tryCatch(httr2::resp_body_json(resp, simplifyVector = TRUE), error = function(e) NULL)
        status_reason <- tryCatch(httr2::resp_status_desc(resp), error = function(e) NULL)
        detail <- if (!is.null(err_body) && !is.null(err_body$error) && !is.null(err_body$error$message)) {
          err_body$error$message
        } else if (!is.null(status_reason) && nzchar(status_reason)) {
          paste(status_code, status_reason)
        } else {
          paste("HTTP", status_code)
        }

        fallback_due_to_model <- status_code %in% c(400, 404) &&
          grepl("model", detail, ignore.case = TRUE) && !identical(mdl, fallback_model)
        fallback_to_chat <- status_code %in% c(404, 405) ||
          grepl("ResponsesAPIResponse|/responses|no-default-models", detail, ignore.case = TRUE)

        if (fallback_due_to_model) {
          message("[OPENAI] Falling back from model '", mdl, "' after error: ", detail)
          last_error_message <- detail
          next
        }

        if (fallback_to_chat) {
          use_chat_api <- TRUE
        } else {
          stop(detail, call. = FALSE)
        }
      } else {
        body <- httr2::resp_body_json(resp, simplifyVector = FALSE)

        extract_output_text <- function(x) {
          if (is.null(x) || !length(x)) return(NULL)
          pieces <- unlist(lapply(x, function(item) {
            if (!is.list(item)) return(if (is.character(item)) item else NULL)
            content <- item$content %||% list()
            unlist(lapply(content, function(chunk) {
              if (is.list(chunk)) {
                chunk$text %||% NULL
              } else if (is.character(chunk)) {
                chunk
              } else {
                NULL
              }
            }), use.names = FALSE)
          }), use.names = FALSE)
          pieces <- pieces[nzchar(pieces)]
          if (!length(pieces)) return(NULL)
          paste(pieces, collapse = "\n")
        }

        text <- NULL
        if (!is.null(body$output_text)) {
          if (is.character(body$output_text)) {
            text <- paste(body$output_text[nzchar(body$output_text)], collapse = "\n")
          } else if (is.list(body$output_text)) {
            text <- paste(unlist(body$output_text, use.names = FALSE), collapse = "\n")
          }
        }
        if (is.null(text) || !nzchar(text)) {
          text <- extract_output_text(body$output)
        }

        if (is.null(text) || !nzchar(text)) {
          stop("Empty completion returned from API.", call. = FALSE)
        }

        usage <- parse_usage(body)
        return(list(
          text = stringr::str_trim(as.character(text)),
          usage = usage
        ))
      }
    } # end if (!use_chat_api)
    
    if (use_chat_api) {
      message("[OPENAI] Retrying model '", mdl, "' via ", chat_url)
      chat_payload <- list(
        model = mdl,
        messages = chat_messages,
        temperature = temperature,
        max_tokens = max_output_tokens
      )

      req_chat <- httr2::request(chat_url) |>
        httr2::req_headers(
          Authorization = paste("Bearer", api_key),
          `Content-Type` = "application/json"
        ) |>
        httr2::req_body_json(chat_payload, auto_unbox = TRUE) |>
        httr2::req_timeout(30) |>
        httr2::req_retry(
          max_tries = 3,
          is_transient = function(resp) {
            status <- httr2::resp_status(resp)
            status == 429 || (status >= 500 && status < 600)
          }
        )

      resp_chat <- tryCatch(httr2::req_perform(req_chat), error = identity)

      if (inherits(resp_chat, "httr2_http")) {
        last_error_message <- conditionMessage(resp_chat)
        if (!is.null(resp_chat$resp)) {
          resp_chat <- resp_chat$resp
        } else {
          next
        }
      } else if (inherits(resp_chat, "error")) {
        last_error_message <- conditionMessage(resp_chat)
        next
      }

      status_chat <- httr2::resp_status(resp_chat)
      if (status_chat >= 300) {
        err_body_chat <- tryCatch(httr2::resp_body_json(resp_chat, simplifyVector = TRUE), error = function(e) NULL)
        status_reason_chat <- tryCatch(httr2::resp_status_desc(resp_chat), error = function(e) NULL)
        detail_chat <- if (!is.null(err_body_chat) && !is.null(err_body_chat$error) && !is.null(err_body_chat$error$message)) {
          err_body_chat$error$message
        } else if (!is.null(status_reason_chat) && nzchar(status_reason_chat)) {
          paste(status_chat, status_reason_chat)
        } else {
          paste("HTTP", status_chat)
        }

        fallback_due_to_model_chat <- status_chat %in% c(400, 404) &&
          grepl("model", detail_chat, ignore.case = TRUE) && !identical(mdl, fallback_model)

        if (fallback_due_to_model_chat) {
          message("[OPENAI] Falling back from model '", mdl, "' after chat error: ", detail_chat)
          last_error_message <- detail_chat
          next
        }

        stop(detail_chat, call. = FALSE)
      }

      body_chat <- httr2::resp_body_json(resp_chat, simplifyVector = FALSE)
      text_chat <- NULL
      if (!is.null(body_chat$choices) && length(body_chat$choices) > 0) {
        choice <- body_chat$choices[[1]]
        if (!is.null(choice$message) && !is.null(choice$message$content)) {
          msg_content <- choice$message$content
          if (is.character(msg_content)) {
            text_chat <- paste(msg_content[nzchar(msg_content)], collapse = "\n")
          } else if (is.list(msg_content)) {
            text_chat <- paste(unlist(lapply(msg_content, function(x) {
              if (is.list(x)) {
                x$text %||% NULL
              } else if (is.character(x)) {
                x
              } else {
                NULL
              }
            }), use.names = FALSE), collapse = "\n")
          }
        }
      }

      if (is.null(text_chat) || !nzchar(text_chat)) {
        stop("Empty completion returned from chat API.", call. = FALSE)
      }

      usage_chat <- parse_usage(body_chat)
      return(list(
        text = stringr::str_trim(as.character(text_chat)),
        usage = usage_chat
      ))
    }
  }

  stop(last_error_message %||% "OpenAI request failed after retries.", call. = FALSE)
}

# Restore basic viewer components needed for Avivator
load_channel_names <- function() {
  candidates <- unique(c(
    "Channel_names.txt",
    "Channel_names",
    file.path("..", "Channel_names.txt"),
    file.path("..", "Channel_names")
  ))

  for (cand in candidates) {
    if (!file.exists(cand)) next

    lines <- readLines(cand, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (!length(lines)) next

    parsed <- stringr::str_match(lines, "^\\s*(.*?)\\s*(?:\\(C(\\d+)\\))?\\s*$")
    labels <- parsed[, 2]
    idx <- suppressWarnings(as.integer(parsed[, 3]))

    if (all(is.na(idx))) {
      idx <- seq_along(labels)
    }

    ord <- order(idx, na.last = TRUE)
    labels <- labels[ord]
    labels <- labels[nzchar(labels)]
    if (!length(labels)) next

    message(sprintf("[CHANNELS] Loaded %d channel names from %s", length(labels), cand))
    return(labels)
  }

  message("[CHANNELS] Channel names not found. Embed names with: tiffcomment -set \"ChannelNames=<comma-separated>\" <file.ome.tif>")
  NULL
}

channel_names_vec <- load_channel_names()

# Local images setup
# Prefer static www/local_images so shiny-server can serve large OME-TIFFs with HTTP Range support.
# Only add a dynamic resource path if www/local_images does not exist.
local_images_env <- Sys.getenv("LOCAL_IMAGE_ROOT", unset = "")
local_images_root <- NULL
if (nzchar(local_images_env)) {
  local_images_root <- tryCatch(normalizePath(local_images_env, mustWork = TRUE), error = function(e) NULL)
}
if (is.null(local_images_root)) {
  candidate <- if (!is.null(project_root)) file.path(project_root, "local_images") else NULL
  if (!is.null(candidate) && dir.exists(candidate)) {
    local_images_root <- candidate
  } else {
    local_images_root <- tryCatch(normalizePath(file.path("..", "..", "local_images"), mustWork = FALSE), error = function(e) NULL)
  }
}

www_local_images_dir <- file.path("www", "local_images")
has_www_local_images <- dir.exists(www_local_images_dir)

# Register /images resource path for backup access (avoid conflict with www/local_images)
# This provides redundancy even if www/local_images exists
if (!is.null(local_images_root) && dir.exists(local_images_root)) {
  try({ 
    shiny::addResourcePath("images", local_images_root) 
    cat("[SETUP] Registered resource path /images ->", local_images_root, "\n")
  }, silent = TRUE)
}

# Register resource path for viewer URLs  
if (has_www_local_images) {
  try({
    # Register local_images path for both local and remote compatibility
    # Use 'img_data' to avoid conflict warning with www/local_images subdirectory
    shiny::addResourcePath("img_data", file.path("www", "local_images"))
    cat("[SETUP] Registered /img_data -> www/local_images\n")
  }, silent = TRUE)
}

label_exports_dir <- file.path("www", "LabelExports")
if (dir.exists(label_exports_dir)) {
  try({
    shiny::addResourcePath("label_exports", label_exports_dir)
    cat("[SETUP] Registered /label_exports ->", label_exports_dir, "\n")
    seg_files <- list.files(label_exports_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
    cat("[SETUP] Found", length(seg_files), "segmentation OME-TIFF files\n")
  }, silent = TRUE)
} else {
  cat("[SETUP] LabelExports directory not found at", label_exports_dir, "\n")
}

# Log the setup for debugging
cat("[SETUP] www/local_images exists:", has_www_local_images, "\n")
if (has_www_local_images) {
  local_files <- list.files(www_local_images_dir, pattern = "\\.(ome\\.)?tiff?$", ignore.case = TRUE)
  cat("[SETUP] Found", length(local_files), "images in www/local_images\n")
}

# Load segmentation annotations from data directory
load_segmentation_data <- function() {
  seg_path <- file.path("..", "..", "data", "annotations.tsv")
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
  lookup_path <- file.path("..", "..", "data", "islet_spatial_lookup.csv")
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

# Load at startup
segmentation_data <- load_segmentation_data()
islet_spatial_lookup <- load_islet_spatial_lookup()

# Keep legacy annotations_data for backward compatibility
annotations_data <- segmentation_data

# ---- Viewer helpers/defaults (Avivator) ----
# Optional default image URL; if NULL, we'll auto-pick the first available under /local_images
default_image_url <- NULL

# Environment detection for viewer configuration
detect_environment <- function(session = NULL) {
  # Check various indicators to determine if we're in a reverse proxy setup
  
  # Method 1: Check environment variables commonly set in reverse proxy setups
  proxy_indicators <- c(
    "HTTP_X_FORWARDED_FOR",
    "HTTP_X_FORWARDED_HOST", 
    "HTTP_X_FORWARDED_PROTO",
    "HTTP_X_REAL_IP",
    "SHINY_SERVER_VERSION"
  )
  
  has_proxy_env <- any(sapply(proxy_indicators, function(x) nzchar(Sys.getenv(x, ""))))
  
  # Method 2: Check if we're running on a non-localhost host
  is_remote_host <- FALSE
  if (!is.null(session) && !is.null(session$clientData$url_hostname)) {
    hostname <- session$clientData$url_hostname
    is_remote_host <- !(hostname %in% c("localhost", "127.0.0.1", "0.0.0.0"))
  }
  
  # Method 3: Check URL path structure (reverse proxy often has subpaths)
  has_subpath <- FALSE
  if (!is.null(session) && !is.null(session$clientData$url_pathname)) {
    pathname <- session$clientData$url_pathname
    # If pathname is not just "/" or "/app/", likely in a reverse proxy
    has_subpath <- !pathname %in% c("/", "/app/")
  }
  
  # Return environment info
  list(
    is_reverse_proxy = has_proxy_env || is_remote_host || has_subpath,
    hostname = if (!is.null(session)) session$clientData$url_hostname else "unknown",
    pathname = if (!is.null(session)) session$clientData$url_pathname else "/",
    indicators = list(
      proxy_env = has_proxy_env,
      remote_host = is_remote_host,
      subpath = has_subpath
    )
  )
}

# Lightweight helper to get image dimensions without loading full pixel data.
# Uses the 'magick' package if available; otherwise returns NULL.
get_image_dims <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("magick", quietly = TRUE)) return(NULL)
  dims <- tryCatch({
    img <- magick::image_read(path)
    info <- magick::image_info(img)
    # 'info' is a data.frame; width/height in pixels
    list(width = as.integer(info$width[[1]]), height = as.integer(info$height[[1]]))
  }, error = function(e) NULL)
  dims
}

# Extract OME-XML PhysicalSizeX/Y (µm) from TIFF ImageDescription and compute µm per pixel.
get_pixel_size_um <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("magick", quietly = TRUE) || !requireNamespace("xml2", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    img <- magick::image_read(path)
    attrs <- magick::image_attributes(img)
    desc <- NULL
    # magick attributes may contain TIFF tags; look for ImageDescription
    if (!is.null(attrs) && "tiff:ImageDescription" %in% names(attrs)) {
      desc <- attrs[["tiff:ImageDescription"]]
    } else if (!is.null(attrs) && "exif:ImageDescription" %in% names(attrs)) {
      desc <- attrs[["exif:ImageDescription"]]
    }
    if (is.null(desc) || !nzchar(desc)) return(NULL)
    # Parse OME-XML
    doc <- xml2::read_xml(desc)
    px <- xml2::xml_find_first(doc, "//Pixels")
    if (is.na(px)) return(NULL)
    phys_x <- suppressWarnings(as.numeric(xml2::xml_attr(px, "PhysicalSizeX")))
    phys_y <- suppressWarnings(as.numeric(xml2::xml_attr(px, "PhysicalSizeY")))
    unit_x <- xml2::xml_attr(px, "PhysicalSizeXUnit") %||% "µm"
    unit_y <- xml2::xml_attr(px, "PhysicalSizeYUnit") %||% "µm"
    # Only support µm units here
    if (is.na(phys_x) || is.na(phys_y)) return(NULL)
    if (!identical(unit_x, "µm") || !identical(unit_y, "µm")) return(NULL)
    list(x_um_per_px = phys_x, y_um_per_px = phys_y)
  }, error = function(e) NULL)
  res
}

# Helper function to find image file for case_id
find_image_for_case <- function(case_id, www_dir, local_root) {
  pattern <- sprintf(".*%s.*\\.(ome\\.)?tiff?$", case_id)
  
  # Try www directory first (preferred for HTTP Range support)
  if (!is.null(www_dir) && dir.exists(www_dir)) {
    matches <- list.files(www_dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
    if (length(matches) > 0) {
      cat("[IMAGE] Found in www:", matches[1], "\n")
      return(matches[1])
    }
  }
  
  # Fallback to LOCAL_IMAGE_ROOT
  if (!is.null(local_root) && dir.exists(local_root)) {
    matches <- list.files(local_root, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
    if (length(matches) > 0) {
      cat("[IMAGE] Found in LOCAL_IMAGE_ROOT:", matches[1], "\n")
      return(matches[1])
    }
  }
  
  cat("[IMAGE] Not found for case:", case_id, "\n")
  return(NULL)
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
  
  # Fallback to segmentation_data
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

resolve_avivator_base <- function() {
  # Return the URL path to the local avivator build under www/, or NULL if missing
  local_index <- file.path("www", "avivator", "index.html")
  if (file.exists(local_index)) return("avivator/index.html")
  NULL
}

# Build a base64-encoded channel configuration understood by embedded Avivator.
# Expected shape (before base64):
#   {
#     channelNames: ["DAPI", "CD31", ...],
#     primaryChannels: [
#       { index: 0, name: "DAPI", color: "#7F7F7F", visible: true },
#       { index: 25, name: "INS",  color: "#E41A1C", visible: true },
#       { index: 19, name: "GCG",  color: "#377EB8", visible: true },
#       { index: 13, name: "SST",  color: "#FFCC00", visible: true }
#     ]
#   }
build_channel_config_b64 <- function(names_vec) {
  if (is.null(names_vec) || length(names_vec) == 0) return(NULL)

  # Helper: RGB ints (0..255) to hex string
  rgb_hex <- function(r, g, b) sprintf("#%02X%02X%02X", as.integer(r), as.integer(g), as.integer(b))
  # Palette for key channels
  hex_red    <- rgb_hex(228,  26,  28)  # INS
  hex_blue   <- rgb_hex( 55, 126, 184)  # GCG
  hex_yellow <- rgb_hex(255, 204,   0)  # SST
  hex_grey   <- rgb_hex(127, 127, 127)  # DAPI

  # Locate indices for the key channels (names_vec is 1-based by channel number)
  find_idx0 <- function(tag) {
    hit <- which(toupper(names_vec) == tag)
    if (length(hit) == 0) return(NA_integer_)
    (hit[[1]] - 1L) # zero-based for viewer
  }
  idx0_dapi <- find_idx0("DAPI")
  idx0_ins  <- find_idx0("INS")
  idx0_gcg  <- find_idx0("GCG")
  idx0_sst  <- find_idx0("SST")

  prim <- list()
  add_pc <- function(idx0, name, hex) {
    if (!is.na(idx0) && idx0 >= 0) prim[[length(prim)+1]] <<- list(index = idx0, name = name, color = hex, visible = TRUE)
  }
  # Order: INS (red), GCG (blue), SST (yellow), DAPI (grey)
  add_pc(idx0_ins,  "INS",  hex_red)
  add_pc(idx0_gcg,  "GCG",  hex_blue)
  add_pc(idx0_sst,  "SST",  hex_yellow)
  add_pc(idx0_dapi, "DAPI", hex_grey)

  # Build config
  cfg <- list(
    channelNames = as.list(as.character(names_vec)),
    primaryChannels = prim
  )
  js <- jsonlite::toJSON(cfg, auto_unbox = TRUE)
  base64enc::base64encode(charToRaw(js))
}

# Construct an app-absolute URL (always rooted at the current Shiny pathname)
# so embedded viewers resolve assets from the correct base whether or not the
# app is hosted under a reverse proxy sub-path.
build_app_absolute_url <- function(session, rel_path) {
  if (is.null(rel_path) || !nzchar(rel_path)) return(NULL)
  rel_clean <- gsub("\\\\", "/", rel_path)
  rel_clean <- sub("^/+", "", rel_clean)

  base_path <- "/"
  path_candidate <- NULL
  if (!is.null(session)) {
    path_candidate <- tryCatch(session$clientData$url_pathname, error = function(e) NULL)
  }
  if (!is.null(path_candidate) && nzchar(path_candidate)) {
    base_path <- trimws(path_candidate)
  }
  # Strip any query/hash fragments that occasionally show up in clientData
  if (grepl("\\?", base_path, fixed = TRUE)) {
    base_path <- strsplit(base_path, "?", fixed = TRUE)[[1]][1]
  }
  if (grepl("#", base_path, fixed = TRUE)) {
    base_path <- strsplit(base_path, "#", fixed = TRUE)[[1]][1]
  }
  if (!nzchar(base_path)) base_path <- "/"
  if (substr(base_path, 1, 1) != "/") {
    base_path <- paste0("/", base_path)
  }
  base_path <- gsub("/+$", "", base_path)
  if (!nzchar(base_path)) base_path <- "/"
  if (!grepl("/$", base_path)) {
    base_path <- paste0(base_path, "/")
  }

  paste0(base_path, rel_clean)
}

# Build a fully-qualified HTTP(S) URL using the current session origin so that
# embedded clients (like Avivator) always receive a resolvable absolute URL.
build_public_http_url <- function(session, rel_path) {
  if (is.null(rel_path) || !nzchar(rel_path)) return(NULL)
  app_path <- build_app_absolute_url(session, rel_path)
  if (is.null(app_path)) return(NULL)

  proto <- NULL
  host <- NULL
  port <- NULL
  if (!is.null(session)) {
    proto <- tryCatch(session$clientData$url_protocol, error = function(e) NULL)
    host <- tryCatch(session$clientData$url_hostname, error = function(e) NULL)
    port <- tryCatch(session$clientData$url_port, error = function(e) NULL)
  }

  proto <- trimws(proto %||% "http:")
  proto <- sub(":$", "", proto)
  if (!nzchar(proto)) proto <- "http"

  host <- trimws(host %||% "127.0.0.1")
  if (!nzchar(host)) host <- "127.0.0.1"

  port <- trimws(port %||% "")
  default_port <- if (proto == "https") "443" else "80"
  port_fragment <- ""
  if (nzchar(port) && port != default_port) {
    port_fragment <- paste0(":", port)
  }

  paste0(proto, "://", host, port_fragment, app_path)
}

# Probe whether a relative viewer asset exists on disk and is reachable over HTTP.
probe_viewer_asset <- function(session, rel_path) {
  diag <- list(
    rel_path = rel_path,
    local_path = NULL,
    file_exists = NA,
    file_size = NA,
    http_url = NULL,
    http_status = NA_integer_,
    http_error = NULL
  )

  if (is.null(rel_path) || !nzchar(rel_path)) return(diag)

  # Resolve local filesystem path (only handles www/local_images or subdirectories)
  rel_parts <- strsplit(rel_path, "/", fixed = TRUE)[[1]]
  local_path <- do.call(file.path, as.list(c("www", rel_parts)))
  local_path <- suppressWarnings(normalizePath(local_path, winslash = "/", mustWork = FALSE))
  diag$local_path <- local_path
  diag$file_exists <- file.exists(local_path)
  if (isTRUE(diag$file_exists)) {
    info <- file.info(local_path)
    diag$file_size <- info$size
  }

  diag$http_url <- build_public_http_url(session, rel_path)
  if (!is.null(diag$http_url) && requireNamespace("curl", quietly = TRUE)) {
    h <- curl::new_handle()
    curl::handle_setopt(h, nobody = TRUE, customrequest = "HEAD", timeout = 5)
    tryCatch({
      resp <- curl::curl_fetch_memory(diag$http_url, handle = h)
      diag$http_status <- resp$status_code %||% NA_integer_
    }, error = function(e) {
      diag$http_error <- conditionMessage(e)
    })
  }

  diag
}

# ---------- Data loading and wrangling ----------

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

add_islet_key <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)

  # Start with an existing islet_key column if present; otherwise create NA placeholder
  if (!"islet_key" %in% names(df)) df$islet_key <- NA_character_

  # 1) Derive from 'region' if available, e.g., "Islet_200_union" -> "Islet_200"
  if ("region" %in% names(df)) {
    key_from_region <- stringr::str_extract(df$region, "Islet_\\d+")
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_region)
  }

  # 2) Derive from 'name' if available, e.g., "Islet_200" or any text containing digits
  if ("name" %in% names(df)) {
    key_from_name <- stringr::str_extract(df$name, "Islet_\\d+")
    only_digits  <- stringr::str_extract(df$name, "\\d+")
    fallback_name <- ifelse(!is.na(only_digits), paste0("Islet_", only_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, key_from_name, fallback_name)
  }

  # 3) Fallback from numeric/string 'islet_id' if present
  if ("islet_id" %in% names(df)) {
    id_str <- suppressWarnings(as.character(df$islet_id))
    id_digits <- stringr::str_extract(id_str, "\\d+")
    fallback_id <- ifelse(!is.na(id_digits), paste0("Islet_", id_digits), NA_character_)
    df$islet_key <- dplyr::coalesce(df$islet_key, fallback_id)
  }

  df
}

compute_diameter_um <- function(area_um2) {
  area_um2 <- suppressWarnings(as.numeric(area_um2))
  ifelse(is.finite(area_um2) & area_um2 > 0, 2 * sqrt(area_um2 / pi), NA_real_)
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
  # Diameter is computed ONLY from core area; if no core exists, diameter is NA
  size_area <- core_area %>%
    dplyr::mutate(islet_diam_um = compute_diameter_um(core_region_um2)) %>%
    dplyr::select(`Case ID`, `Donor Status`, islet_key, islet_diam_um)

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
      # Density is already in µm², no conversion needed
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

    core_rows <- targets_all %>% dplyr::filter(type == "islet_core") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_core = count)
    band_rows <- targets_all %>% dplyr::filter(type == "islet_band") %>% dplyr::select(dplyr::all_of(c(keys, "count"))) %>% dplyr::rename(count_band = count)
    union_missing <- core_rows %>%
      dplyr::inner_join(band_rows, by = keys) %>%
  { safe_left_join(., have_union, by = keys, context = "targets_union_missing:have_union") } %>%
      dplyr::filter(is.na(.has_union))
    if (nrow(union_missing) > 0) {
      synth <- union_missing %>%
        dplyr::transmute(`Case ID`, `Donor Status`, islet_key, type = "islet_union", class,
                          area_um2 = NA_real_, region_um2 = NA_real_,
                          area_density = NA_real_,
                          count = suppressWarnings(as.numeric(count_core)) + suppressWarnings(as.numeric(count_band)))
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

  # Markers with fraction positive / mean intensity (include LGALS3 sheet)
  markers <- master$markers %>% add_islet_key() %>% dplyr::filter(!is.na(islet_key))
  markers_all <- markers %>%
    dplyr::select(dplyr::all_of(c("Case ID", "Donor Status")), dplyr::all_of(aab_cols_markers), islet_key, region_type, marker, n_cells, pos_count, pos_frac) %>%
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
  comp <- comp %>% { safe_left_join(., size_area, by = c("Case ID", "Donor Status", "islet_key"), context = "comp:size_area") }

  message("[prep-final] size_area=", paste(class(size_area), collapse='/'),
    " targets_all=", paste(class(targets_all), collapse='/'),
    " markers_all=", paste(class(markers_all), collapse='/'),
    " comp=", paste(class(comp), collapse='/'))
  list(core_area = size_area, targets_all = targets_all, markers_all = markers_all, comp = comp)
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

# Unbiased pseudotime from features (no group labels). Tries DiffusionMap -> principal curve -> PC1.
## pseudotime removed

#' Compute Lamian-based pseudotime from AnnData
#' Uses Lamian's infer_tree_structure which accounts for multi-sample design
#' @param ad AnnData object (from anndata package)
#' @param features Character vector of feature names to use (NULL = use all)
#' @return Numeric vector of pseudotime values (same length as cells)
compute_lamian_pseudotime <- function(ad, features = NULL) {
  tryCatch({
    require(Lamian, quietly = TRUE)
    
    cat("Computing Lamian pseudotime with infer_tree_structure...\n")
    
    # Extract expression matrix (cells × genes)
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
    # Compute PCA (Lamian expects cells × PCs)
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
    # Transpose expression for Lamian (genes × cells expected)
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
    
    # REVERSE pseudotime direction (biological progression ND→Aab+→T1D)
    # Lamian origin.marker='INS' finds high-INS clusters, but trajectory may run backward
    # Reversal ensures: 0=healthy (high INS), 1=diseased (low INS)
    pseudotime <- 1 - pseudotime
    cat("  ✓ Pseudotime reversed for biological progression (0=healthy, 1=diseased)\n")
    
    cat(sprintf("✓ Lamian pseudotime computed (range: %.3f - %.3f)\n", 
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

summary_stats <- function(df, group_cols, value_col, stat = c("mean_se","mean_sd","median_iqr")) {
  stat <- match.arg(stat)
  df %>% group_by(across(all_of(group_cols))) %>%
    summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      se = ifelse(n > 0, sd / sqrt(pmax(n, 1)), NA_real_),
      median = median(.data[[value_col]], na.rm = TRUE),
      q1 = quantile(.data[[value_col]], 0.25, na.rm = TRUE, type = 7),
      q3 = quantile(.data[[value_col]], 0.75, na.rm = TRUE, type = 7),
      .groups = "drop"
    ) %>%
    mutate(
      y = dplyr::case_when(
        stat == "mean_se" ~ mean,
        stat == "mean_sd" ~ mean,
        stat == "median_iqr" ~ median,
        TRUE ~ mean
      ),
      ymin = dplyr::case_when(
        stat == "mean_se" ~ mean - se,
        stat == "mean_sd" ~ mean - sd,
        stat == "median_iqr" ~ q1,
        TRUE ~ mean - se
      ),
      ymax = dplyr::case_when(
        stat == "mean_se" ~ mean + se,
        stat == "mean_sd" ~ mean + sd,
        stat == "median_iqr" ~ q3,
        TRUE ~ mean + se
      )
    )
}

per_bin_anova <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  # capture unique bins with their numeric mid
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]]) %>%
    summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    # Require at least 2 donor groups present and at least some finite values
    if (n_distinct(na.omit(sub[[group_col]])) < 2 || sum(is.finite(sub[[value_col]])) < 2) return(NULL)
    fit <- tryCatch(aov(reformulate(group_col, response = value_col), data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    # Extract one-way ANOVA p-value from ANOVA table (first term)
    pval <- tryCatch({
      at <- anova(fit)
      as.numeric(at[["Pr(>F)"]][1])
    }, error = function(e) NA_real_)
    tibble(
      bin = as.character(b),
      mid = bmeta$mid[i],
      p_anova = pval
    )
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(tibble(bin = character(), mid = numeric(), p_anova = numeric()))
  arrange(out, mid)
}

per_bin_kendall <- function(df, bin_col, group_col, value_col, mid_col = "diam_mid") {
  if (nrow(df) == 0) return(tibble(bin = character(), mid = numeric(), p_kendall = numeric(), tau = numeric()))
  # Map group to ordered numeric along pseudotime ND < Aab+ < T1D
  code_group <- function(x) {
    x <- as.character(x)
    ifelse(x == "ND", 0, ifelse(x == "Aab+", 1, ifelse(x == "T1D", 2, NA_real_)))
  }
  bmeta <- df %>% filter(!is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]]) %>% summarise(mid = suppressWarnings(as.numeric(first(na.omit(.data[[mid_col]])))), .groups = "drop")
  res <- lapply(seq_len(nrow(bmeta)), function(i) {
    b <- bmeta[[bin_col]][i]
    sub <- df %>% filter(.data[[bin_col]] == b)
    x <- code_group(sub[[group_col]])
    y <- suppressWarnings(as.numeric(sub[[value_col]]))
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    if (length(unique(x)) < 2 || length(y) < 3) return(NULL)
    ct <- tryCatch(cor.test(x, y, method = "kendall", exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) return(NULL)
    tibble(bin = as.character(b), mid = bmeta$mid[i], p_kendall = unname(ct$p.value), tau = unname(ct$estimate))
  })
  out <- bind_rows(res)
  if (nrow(out) == 0) return(out)
  arrange(out, mid)
}

# ---------- UI ----------

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML("\n    /* Viewer and trajectory mode styles - fix tab positioning */\n    body.viewer-mode .col-sm-2 {\n      display: none !important;\n    }\n    body.viewer-mode .col-sm-10 {\n      width: 100% !important;\n      max-width: 100% !important;\n      flex: 0 0 100%;\n    }\n    body.viewer-mode .nav-tabs {\n      position: static !important;\n      width: 100% !important;\n    }\n    body.trajectory-mode .container-fluid > .row > .col-sm-3 {\n      display: none !important;\n    }\n    body.trajectory-mode .container-fluid > .row > .col-sm-9 {\n      width: 100% !important;\n      max-width: 100% !important;\n      flex: 0 0 100%;\n    }\n    \n    /* Global biomedical theme styling */\n    body {\n      background-color: #f8f9fa;\n      font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;\n      padding-top: 0;\n      min-width: 1200px;\n      overflow-x: auto;\n    }\n    \n    .container-fluid {\n      background-color: #ffffff;\n      border-radius: 12px;\n      box-shadow: 0 4px 12px rgba(0,0,0,0.08);\n      margin: 10px;\n      padding: 20px;\n      min-width: 1180px;\n    }\n    \n    /* Logo header styling */\n    .logo-header {\n      display: flex;\n      justify-content: flex-end;\n      align-items: center;\n      padding: 10px 0;\n      margin-bottom: 10px;\n      background: linear-gradient(135deg, #f8f9fa 0%, #ffffff 100%);\n      border-bottom: 2px solid #e3f2fd;\n    }\n    \n    /* Enhanced card styling with biomedical color scheme */\n    .card {\n      background: linear-gradient(145deg, #ffffff 0%, #f8f9fa 100%);\n      border: 1px solid #e3f2fd;\n      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.08);\n      transition: all 0.3s ease;\n      border-radius: 12px;\n    }\n    .card:hover {\n      box-shadow: 0 8px 20px rgba(44, 90, 160, 0.15);\n      transform: translateY(-2px);\n    }\n    \n    .card h5 {\n      color: #2c5aa0;\n      font-weight: 600;\n      border-bottom: 2px solid #e3f2fd;\n      padding-bottom: 8px;\n      margin-bottom: 15px;\n    }\n    \n    /* Sidebar styling with scientific theme */\n    .sidebar {\n      background: linear-gradient(180deg, #2c5aa0 0%, #1e3a72 100%);\n      color: white;\n      border-radius: 12px;\n      padding: 20px;\n      font-size: 14px;\n      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.2);\n    }\n    \n    .sidebar h4, .sidebar h5 {\n      color: #ffffff;\n      font-weight: 600;\n      border-bottom: 1px solid rgba(255,255,255,0.2);\n      padding-bottom: 8px;\n    }\n    \n    .sidebar .form-group {\n      margin-bottom: 18px;\n    }\n    \n    .sidebar label {\n      color: #e3f2fd;\n      font-size: 13px;\n      font-weight: 600;\n    }\n    \n    .sidebar .form-control {\n      background-color: rgba(255,255,255,0.9);\n      border: 1px solid #b3d9ff;\n      border-radius: 6px;\n      color: #2c5aa0;\n    }\n    \n    .sidebar .form-control:focus {\n      background-color: #ffffff;\n      border-color: #66b3ff;\n      box-shadow: 0 0 0 0.2rem rgba(44, 90, 160, 0.25);\n    }\n    \n    .sidebar .btn {\n      background: linear-gradient(145deg, #66b3ff 0%, #4da6ff 100%);\n      border: none;\n      border-radius: 6px;\n      color: white;\n      font-weight: 500;\n    }\n    \n    .sidebar .btn:hover {\n      background: linear-gradient(145deg, #4da6ff 0%, #3399ff 100%);\n      transform: translateY(-1px);\n    }\n    \n    /* Tab styling */\n    .nav-tabs {\n      border-bottom: 2px solid #e3f2fd;\n    }\n    \n    .nav-tabs .nav-link {\n      color: #2c5aa0;\n      font-weight: 500;\n      border: none;\n      border-radius: 8px 8px 0 0;\n      margin-right: 4px;\n      background-color: #f8f9fa;\n    }\n    \n    .nav-tabs .nav-link.active {\n      background: linear-gradient(145deg, #2c5aa0 0%, #1e3a72 100%);\n      color: white;\n      border-bottom: 3px solid #66b3ff;\n    }\n    \n    .nav-tabs .nav-link:hover {\n      background-color: #e3f2fd;\n      color: #1e3a72;\n    }\n    \n    /* Form controls styling */\n    .form-control {\n      border: 2px solid #e3f2fd;\n      border-radius: 6px;\n      transition: all 0.2s ease;\n    }\n    \n    .form-control:focus {\n      border-color: #66b3ff;\n      box-shadow: 0 0 0 0.2rem rgba(44, 90, 160, 0.15);\n    }\n    \n    /* Button styling */\n    .btn-primary {\n      background: linear-gradient(145deg, #2c5aa0 0%, #1e3a72 100%);\n      border: none;\n      border-radius: 8px;\n      font-weight: 500;\n      padding: 8px 16px;\n      transition: all 0.2s ease;\n    }\n    \n    .btn-primary:hover {\n      background: linear-gradient(145deg, #1e3a72 0%, #0f1f3d 100%);\n      transform: translateY(-2px);\n      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.3);\n    }\n    \n    /* Checkbox and radio button styling */\n    input[type='checkbox'], input[type='radio'] {\n      accent-color: #2c5aa0;\n    }\n    \n    /* Slider styling */\n    .irs--shiny {\n      color: #2c5aa0;\n    }\n    \n    .irs--shiny .irs-bar {\n      background: linear-gradient(90deg, #66b3ff 0%, #2c5aa0 100%);\n    }\n    \n    .irs--shiny .irs-handle {\n      background: #2c5aa0;\n      border: 3px solid #ffffff;\n    }\n    \n    /* Help text styling */\n    .help-block {\n      color: #b3d9ff;\n      font-size: 12px;\n      font-style: italic;\n    }\n    \n    /* Well and panel styling */\n    .well {\n      background: linear-gradient(145deg, #f8f9fa 0%, #e3f2fd 100%);\n      border: 1px solid #b3d9ff;\n      border-radius: 8px;\n    }\n    .ai-chat-panel {
      background: linear-gradient(145deg, #ffffff 0%, #f8f9fa 100%);\n      border: 1px solid #e3f2fd;\n      border-radius: 12px;\n      box-shadow: 0 4px 12px rgba(44, 90, 160, 0.08);\n      padding: 15px;\n      display: flex;\n      flex-direction: column;\n      gap: 12px;\n    }\n    .ai-chat-logo {\n      display: flex;\n      justify-content: center;\n      align-items: center;\n      padding-bottom: 8px;\n      border-bottom: 1px solid rgba(44, 90, 160, 0.15);\n    }\n    .ai-chat-logo img {\n      max-width: 100%;\n      height: auto;\n      max-height: 90px;\n    }\n    .ai-chat-header {\n      display: flex;\n      justify-content: center;\n      align-items: center;\n      margin-top: 8px;\n      margin-bottom: 4px;\n    }\n    .ai-chat-history {\n      flex: 1;\n      overflow-y: auto;\n      padding: 8px;\n      display: flex;\n      flex-direction: column;\n      gap: 10px;\n    }\n    .ai-chat-message {\n      padding: 10px 12px;\n      border-radius: 8px;\n      margin: 4px 0;\n      max-width: 90%;\n      word-wrap: break-word;\n    }\n    .ai-chat-message.user {\n      background: linear-gradient(145deg, #e3f2fd 0%, #bbdefb 100%);\n      border-left: 4px solid #2c5aa0;\n      align-self: flex-end;\n      margin-left: auto;\n    }\n    .ai-chat-message.assistant {\n      background: linear-gradient(145deg, #f5f5f5 0%, #eeeeee 100%);\n      border-left: 4px solid #66b3ff;\n      align-self: flex-start;\n      margin-right: auto;\n    }\n    .ai-chat-meta {\n      font-weight: 600;\n      font-size: 12px;\n      display: block;\n      margin-bottom: 4px;\n      text-transform: uppercase;\n      letter-spacing: 0.5px;\n    }\n    .ai-chat-message.user .ai-chat-meta {\n      color: #1e3a72;\n    }\n    .ai-chat-message.assistant .ai-chat-meta {\n      color: #2c5aa0;\n    }\n  "))),
  tags$head(tags$style(HTML("\n    .ai-chat-logo img {\n      max-width: 100%;\n      height: auto;\n      max-height: 120px;\n    }\n    /* Equal height panels */
    .equal-height-row {\n      display: flex;\n      flex-wrap: nowrap;\n    }\n    .equal-height-panel {\n      height: calc(100vh - 40px);\n      overflow-y: auto;\n    }\n    /* AI chat panel - match card heights */\n    .ai-chat-panel-container {\n      height: calc(100vh - 100px);\n      display: flex;\n      flex-direction: column;\n      overflow-y: auto;\n      margin-top: 20px;\n    }\n  "))),
  uiOutput("theme_css"),
  tags$script(HTML("
    $(document).on('shiny:connected', function() {
      function adjustLayout() {
        var activeTab = $('#tabs li.active a').text().trim();
        var sidebar = $('.equal-height-panel').first();
        var mainPanel = $('.main-content-panel');
        
        // Hide sidebar and expand main panel for non-Plot tabs
        if (activeTab !== 'Plot') {
          sidebar.hide();
          mainPanel.css({
            'width': '100%',
            'max-width': '100%',
            'flex': '0 0 100%'
          });
        } else {
          sidebar.show();
          mainPanel.css({
            'width': '88.33333333%',
            'max-width': '88.33333333%',
            'flex': '0 0 88.33333333%'
          });
        }
      }
      
      // Adjust on initial load
      setTimeout(adjustLayout, 100);
      
      // Adjust when tabs change
      $('a[data-toggle=\"tab\"]').on('shown.bs.tab', function() {
        adjustLayout();
      });
    });
  ")),
  fluidRow(class = "equal-height-row",
    # Left Sidebar Panel
    column(width = 1.4, class = "equal-height-panel",
      conditionalPanel(
        condition = "input.tabs == 'Plot'",
        div(class = "sidebar", style = "height: 100%;",
          h4("Data & Filters"),
          selectInput("mode", "Focus", choices = c("Markers", "Targets", "Composition"), selected = "Composition"),
          uiOutput("region_selector"),
          uiOutput("dynamic_selector"),
          uiOutput("metric_selector"),
      # color selector removed (pseudotime removed)
          hr(),
          h5("Autoantibody filter (Aab+ donors only)"),
          uiOutput("aab_warning"),
          checkboxGroupInput("aab_flags", NULL,
                             choices = c("GADA" = "AAb_GADA",
                                         "IA2A" = "AAb_IA2A",
                                         "ZnT8A" = "AAb_ZnT8A",
                                         "IAA" = "AAb_IAA",
                                         "mIAA" = "AAb_mIAA"),
                             selected = c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")),
          checkboxGroupInput("groups", "Donor Status", choices = c("ND", "Aab+", "T1D"), selected = c("ND", "Aab+", "T1D"), inline = TRUE),
          radioButtons("stat", "Statistic",
                       choices = c("Mean±SE" = "mean_se",
                                   "Mean±SD" = "mean_sd",
                                   "Median + IQR" = "median_iqr"),
                       selected = "mean_se"),
          checkboxInput("show_plot_outlier_table", "Show outlier table (if outliers detected)", value = FALSE),
          # Outlier table output
          uiOutput("plot_outlier_info"),
          hr(),
          h5("Export"),
          downloadButton("dl_summary", "Download summary CSV")
        )
      )
    ),
    # Main Panel - uses full width when sidebar is hidden
    column(width = 10.6, class = "equal-height-panel main-content-panel",
          tabsetPanel(id = "tabs",
        tabPanel("Plot",
                 fluidRow(
                   # Left card: Main plot (smaller to align bottom with left panel)
                   column(5,
                     div(class = "card", style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: calc(100vh - 100px); display: flex; flex-direction: column; gap: 12px;",
                       h5("Islet Size Distribution", style = "margin-top: 0; color: #333;"),
                       div(style = "flex: 1; min-height: 0; display: flex;", 
                         plotlyOutput("plt", height = "100%")
                       ),
                       div(style = "flex: 0 0 auto; display: flex; flex-direction: column; gap: 10px;",
                         # Reorganized controls with better spacing
                         fluidRow(
                           column(4, 
                             sliderInput("binwidth", "Diameter bin width (µm)", min = 1, max = 75, value = 50, step = 1),
                             sliderInput("diam_range", "Islet diameter range (µm)", min = 0, max = 500, value = c(0, 350), step = 10)
                           ),
                           column(4,
                             sliderInput("pt_size", "Point size", min = 0.3, max = 4.0, value = 0.8, step = 0.1),
                             sliderInput("pt_alpha", "Point transparency", min = 0.05, max = 1.0, value = 0.25, step = 0.05)
                           ),
                           column(4,
                             selectInput("plot_color_by", "Color points by:",
                                        choices = c("Donor Status" = "donor_status", 
                                                  "Donor ID" = "donor_id"),
                                        selected = "donor_status"),
                             radioButtons("add_smooth", "Trend line", choices = c("None", "LOESS"), selected = "None", inline = TRUE)
                           )
                         ),
                         div(style = "display: flex; flex-wrap: wrap; gap: 15px; justify-content: flex-end;",
                           checkboxInput("exclude_zero_top", "Exclude zero values", value = FALSE),
                           checkboxInput("show_points", "Show individual points", value = FALSE),
                           checkboxInput("log_scale", "Log scale y-axis", value = FALSE),
                           checkboxInput("log_scale_x", "Log2 scale x-axis", value = FALSE)
                         )
                       )
                     )
                   ),
                   # Middle card: Distribution plot (same height as main plot)
                   column(5,
                     div(class = "card", style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: calc(100vh - 100px); display: flex; flex-direction: column; gap: 12px;",
                       h5("Distribution Comparison", style = "margin-top: 0; color: #333;"),
                       div(style = "flex: 1; overflow-y: auto;",
                         plotlyOutput("dist", height = 400),
                         br(),
                         # Distribution controls below the plot
                         uiOutput("dist_ui")
                       )
                     )
                   ),
                   # Right card: AI Assistant
                   column(2,
                     div(class = "card ai-chat-panel", style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: calc(100vh - 100px); display: flex; flex-direction: column; gap: 12px;",
                       div(
                         class = "ai-chat-header",
                         style = "display: flex; align-items: center; padding-bottom: 12px; border-bottom: 1px solid rgba(44, 90, 160, 0.15);",
                         tags$img(src = "logo.png", alt = "Islet Explorer logo", style = "height: 80px; margin-right: 15px;"),
                         div(
                           style = "flex: 1;",
                           span("Powered by Mab Lab and University of Florida Navigator AI Toolkit", 
                                style = "font-size: 14px; font-weight: bold; color: #333;")
                         )
                       ),
                       uiOutput("chat_status", container = div, class = "ai-chat-status"),
                       div(style = "flex: 1; overflow-y: auto;",
                         uiOutput("chat_history", container = div, class = "ai-chat-history")
                       ),
                       div(
                         class = "ai-chat-model-picker",
                         selectInput(
                           "chat_model",
                           label = "Model",
                           choices = c(
                             "Navigator Fast (gpt-oss-20b)" = "gpt-oss-20b",
                             "Navigator Large (gpt-oss-210b)" = "gpt-oss-210b"
                           ),
                           selected = "gpt-oss-20b"
                         )
                       ),
                       div(
                         class = "ai-chat-input",
                         textAreaInput(
                           "chat_message",
                           label = NULL,
                           value = "",
                           placeholder = "Ask about donor trends, plots, or troubleshooting…",
                           height = "110px"
                         )
                       ),
                       div(
                         class = "ai-chat-controls",
                         actionButton("chat_send", "Send", class = "btn btn-primary"),
                         actionButton("chat_reset", "New Conversation", class = "btn btn-outline-secondary")
                       ),
                       uiOutput("chat_feedback", container = div, class = "ai-chat-feedback"),
                       uiOutput("chat_usage", container = div, class = "ai-chat-usage")
                     )
                   )
                 )
        ),
        tabPanel("Trajectory",
          tagList(
            checkboxInput(
              "show_outlier_table",
              "Show outlier table (if outliers detected)",
              value = FALSE
            ),
            uiOutput("traj_status"),
            fluidRow(
              # Left card: Scatterplot and heatmap
              column(8,
                div(class = "card", style = "padding: 15px; margin-bottom: 20px;",
                  # Four columns of selectors - two rows
                  fluidRow(
                    # Column 1: Feature selector and Trend lines
                    column(3,
                      uiOutput("traj_feature_selector"),
                      selectInput("traj_show_trend", "Trend lines:",
                                 choices = c("None" = "none", 
                                           "Overall" = "overall",
                                           "By Donor Status" = "by_donor"),
                                 selected = "by_donor")
                    ),
                    # Column 2: Color points by and Point size by
                    column(3,
                      selectInput("traj_color_by", "Color points by:",
                                 choices = c("Donor Status" = "donor_status", 
                                           "Donor ID" = "donor_id"),
                                 selected = "donor_status"),
                      selectInput("traj_point_size", "Point size by:",
                                 choices = c("Uniform" = "uniform",
                                           "Islet Diameter" = "islet_diam_um"),
                                 selected = "islet_diam_um")
                    ),
                    # Column 3: Two sliders stacked
                    column(3,
                      sliderInput("traj_alpha", "Point transparency:",
                                 min = 0.1, max = 1.0, value = 0.3, step = 0.05),
                      sliderInput("traj_point_size_slider", "Point size:",
                                 min = 0.5, max = 5.0, value = 2.3, step = 0.1)
                    ),
                    # Column 4: Legend display
                    column(3,
                      div(style = "border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #f9f9f9; min-height: 100px;",
                        h6("Legend", style = "margin-top: 0; margin-bottom: 10px; font-weight: bold; color: #333;"),
                        uiOutput("traj_legend")
                      )
                    )
                  ),
                  plotlyOutput("traj_scatter", height = 600),
                  br(),
                  plotOutput("traj_heatmap", height = 120)
                )
              ),
              # Right card: UMAP plots (height matched to left card)
              column(4,
                div(class = "card", style = "padding: 15px; margin-bottom: 20px; height: 950px; overflow-y: auto;",
                  fluidRow(
                    column(12,
                      h5("UMAP: Donor Status"),
                      plotOutput("traj_umap_donor", height = 400)
                    )
                  ),
                  br(),
                  fluidRow(
                    column(12,
                      h5("UMAP: Selected Feature"),
                      plotOutput("traj_umap_feature", height = 400)
                    )
                  )
                )
              )
            ),
            # Outlier table (moved below viewer)
            fluidRow(
              column(12,
                uiOutput("traj_outlier_info")
              )
            ),

            # Islet Segmentation Viewer - appears when a point is clicked
            uiOutput("segmentation_viewer_panel")

          )
        ),

        tabPanel("Viewer",
          div(style = "width: 100%;",
              uiOutput("local_image_picker"),
              uiOutput("vit_view")
          )
        ),
        
        tabPanel("Statistics",
                 fluidRow(
                   column(12,
                     wellPanel(
                       h4("Statistical Analysis Controls"),
                       fluidRow(
                         column(3, selectInput("alpha", "Significance level (α):", choices = c("0.05","0.01","0.001"), selected = "0.05")),
                         column(3, checkboxInput("stats_remove_outliers", "Remove outliers (>3 SD)", value = TRUE)),
                         column(3, br(), actionButton("run_tests", "Run Statistics", class = "btn-primary")),
                         column(3, br(), downloadButton("download_stats", "Download CSV"))
                       )
                     )
                   )
                 ),
                 fluidRow(
                   column(6,
                     wellPanel(
                       h4("Global ANOVA Results"),
                       tableOutput("stats_tbl")
                     )
                   ),
                   column(6,
                     wellPanel(
                       h4("Distribution by Donor Status"),
                       plotlyOutput("stats_plot", height = 320)
                     )
                   )
                 ),
                 fluidRow(
                   column(6,
                     wellPanel(
                       h4("Pairwise Comparisons"),
                       plotlyOutput("pairwise_plot", height = 320)
                     )
                   ),
                   column(6,
                     wellPanel(
                       h4("Area Under Curve by Donor Group"),
                       plotlyOutput("auc_plot", height = 320),
                       br(),
                       tableOutput("auc_table")
                     )
                   )
                 ),
                 fluidRow(
                   column(12,
                     wellPanel(
                       h4("Statistical Test Information"),
                       uiOutput("stats_explanation")
                     )
                   )
                 )
        )
      )  # Close tabsetPanel
    )  # Close Main Panel column
  )  # Close fluidRow
)    # Close fluidPage

# ---------- Server ----------

server <- function(input, output, session) {
  # Shared donor ID palette so colors are consistent across main and distribution plots
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
  # Variable to store forced image selection
  forced_image <- reactiveVal(NULL)
  
  validate_file <- reactive({
    shiny::validate(shiny::need(file.exists(master_path), paste("Not found:", master_path)))
    master_path
  })

  master <- reactive({
    req(validate_file())
    load_master(master_path)
  })

  prepared <- reactive({
    pd <- prep_data(master())
    # NA audits (console)
    try({
      audit_na(pd$targets_all, "targets_all")
      audit_na(pd$markers_all, "markers_all")
      audit_na(pd$comp, "comp")
    }, silent = TRUE)
    pd
  })

  # No initial greeting - start with empty chat history
  chat_history <- reactiveVal(list())
  chat_feedback <- reactiveVal(NULL)
  token_budget_limit <- {
    option_limit <- getOption("openai_token_budget", NULL)
    env_limit_raw <- Sys.getenv("OPENAI_TOKEN_BUDGET", unset = "")
    limit_candidate <- NA_integer_
    if (!is.null(option_limit)) {
      limit_candidate <- suppressWarnings(as.integer(option_limit))
    }
    if (nzchar(env_limit_raw)) {
      limit_candidate <- suppressWarnings(as.integer(env_limit_raw))
    }
    if (is.na(limit_candidate) || limit_candidate <= 0L) {
      NA_integer_
    } else {
      limit_candidate
    }
  }
  total_tokens_used <- reactiveVal(0L)
  streaming_active <- reactiveVal(FALSE)

  trim_history <- function(history, max_entries = 15) {
    if (length(history) <= max_entries) return(history)
    first_entry <- history[[1]]
    rest <- tail(history, max_entries - 1)
    c(list(first_entry), rest)
  }

  format_chat_content <- function(text) {
    if (is.null(text)) return(htmltools::HTML(""))
    text <- stringr::str_replace_all(text, "\r\n", "\n")
    htmltools::HTML(gsub("\n", "<br/>", htmltools::htmlEscape(text)))
  }

  output$chat_status <- renderUI({
    key_val <- get_llm_api_key()
    key_available <- nzchar(key_val)
  base_val <- tryCatch(LLM_API_BASE, error = function(e) NULL)
    norm_base <- if (!is.null(base_val) && nzchar(base_val)) sub("/+\\z", "", base_val, perl = TRUE) else ""
    host_label <- if (nzchar(norm_base)) {
      host <- sub("^https?://", "", norm_base, ignore.case = TRUE)
      sub("/.*$", "", host)
    } else {
      ""
    }
    provider_name <- if (nzchar(norm_base) && grepl("navigator", norm_base, ignore.case = TRUE)) {
      "Navigator Toolkit (UF)"
    } else {
      "OpenAI"
    }

    if (!key_available) {
      span("LLM key not found. Set KEY (and optionally BASE) in ~/.Renviron to enable assistant responses.")
    } else {
      NULL  # Don't show duplicate status text - it's already in the header
    }
  })

  output$chat_history <- renderUI({
    history <- chat_history()
    if (!length(history)) {
      return(div(
        class = "ai-chat-message assistant",
        span(class = "ai-chat-meta", "PANC FLOYD"),
        format_chat_content("I'm ready when you are!")
      ))
    }
    message_nodes <- lapply(history, function(entry) {
      role_label <- if (identical(entry$role, "user")) "You" else "PANC FLOYD"
      div(
        class = paste("ai-chat-message", entry$role),
        span(class = "ai-chat-meta", role_label),
        format_chat_content(entry$content)
      )
    })
    htmltools::tagList(message_nodes)
  })

  output$chat_feedback <- renderUI({
    msg <- chat_feedback()
    if (is.null(msg) || !nzchar(msg)) return(NULL)
    span(msg)
  })

  output$chat_usage <- renderUI({
    used <- total_tokens_used()
    limit <- token_budget_limit
    msg <- if (is.na(limit)) {
      sprintf("Tokens used: %s", format(used, big.mark = ","))
    } else {
      sprintf("Tokens used: %s / %s", format(used, big.mark = ","), format(limit, big.mark = ","))
    }
    span(msg)
  })

  observe({
    if (!isTRUE(streaming_active())) return()
    invalidateLater(120, session)
    session$flushReact()
  })

  observeEvent(input$chat_reset, {
    chat_history(list())
    chat_feedback(NULL)
    updateTextAreaInput(session, "chat_message", value = "")
  })

  observeEvent(input$chat_send, {
    if (!is.na(token_budget_limit) && total_tokens_used() >= token_budget_limit) {
      chat_feedback(sprintf(
        "Token budget reached (%s tokens). Adjust OPENAI_TOKEN_BUDGET or restart the session to continue.",
        format(total_tokens_used(), big.mark = ",")
      ))
      return()
    }

    user_text <- input$chat_message %||% ""
    user_text <- stringr::str_trim(user_text)
    if (!nzchar(user_text)) {
      chat_feedback("Please enter a question before sending.")
      return()
    }

    chat_feedback(NULL)
    updateTextAreaInput(session, "chat_message", value = "")

    history <- chat_history()
    history <- append(history, list(list(role = "user", content = user_text)))
    history <- trim_history(history)
    chat_history(history)

    history_for_api <- history

    if (!nzchar(get_llm_api_key())) {
      chat_feedback("OpenAI assistant is not configured. Define KEY in your environment (e.g., ~/.Renviron).")
      return()
    }

    placeholder_entry <- list(role = "assistant", content = "⏳ …")
  history_with_placeholder <- append(chat_history(), list(placeholder_entry))
  chat_index <- length(history_with_placeholder)
  chat_history(history_with_placeholder)

    final_usage <- NULL
    final_text <- ""

    stream_callback <- function(partial_text, usage = NULL) {
      if (!is.null(partial_text)) {
        final_text <<- partial_text
      }
      if (!is.null(usage)) {
        final_usage <<- usage
      }
      isolate({
        hist <- chat_history()
        if (length(hist) >= chat_index) {
          replacement <- if (!is.null(partial_text) && nzchar(partial_text)) partial_text else "⏳ …"
          hist[[chat_index]]$content <- replacement
          chat_history(hist)
        }
      })
      session$flushReact()
      invisible(NULL)
    }

    shinyjs::disable("chat_send")
    shinyjs::disable("chat_reset")
    streaming_active(TRUE)
    on.exit({
      shinyjs::enable("chat_send")
      shinyjs::enable("chat_reset")
    }, add = TRUE)
    on.exit(streaming_active(FALSE), add = TRUE)

    response <- withProgress(message = "Assistant is thinking...", value = 0, {
      tryCatch(
        call_openai_chat(
          history_for_api,
          model = input$chat_model %||% "auto",
          stream = TRUE,
          stream_callback = stream_callback
        ),
        error = function(e) e
      )
    })
    if (inherits(response, "error")) {
      err_text <- conditionMessage(response)
      message("[OPENAI] Assistant request failed: ", err_text)
      friendly <- "The assistant couldn’t reach the OpenAI service. Please try again shortly."
  base_val <- tryCatch(LLM_API_BASE, error = function(e) "")
      norm_base <- if (nzchar(base_val)) sub("/+\\z", "", base_val, perl = TRUE) else ""
      if (grepl("exceeded your current quota", err_text, ignore.case = TRUE)) {
        friendly <- paste(
          "The assistant can’t reply right now because the OpenAI API quota was exceeded.",
          "Please review plan and billing details at https://platform.openai.com/account/usage before trying again.")
      } else if (grepl("429", err_text, ignore.case = TRUE)) {
        friendly <- paste(
          "The assistant hit the OpenAI rate/usage limit.",
          "Please wait a moment or adjust your plan before trying again.")
      } else if (grepl("401", err_text, ignore.case = TRUE) || grepl("invalid api key", err_text, ignore.case = TRUE)) {
        friendly <- paste(
          "The assistant can’t authenticate with the LLM provider.",
          "Double-check the KEY value in your environment (e.g., ~/.Renviron), then reload the app.")
      } else if (grepl("could not resolve host", err_text, ignore.case = TRUE) ||
                 grepl("name or service not known", err_text, ignore.case = TRUE)) {
        endpoint <- if (nzchar(norm_base)) norm_base else "the configured endpoint"
        friendly <- paste(
          sprintf("The assistant could not resolve %s.", endpoint),
          "Confirm the BASE setting in your environment (e.g., ~/.Renviron) and network/DNS access.")
      } else if (grepl("timed? out", err_text, ignore.case = TRUE) ||
                 grepl("connection refused", err_text, ignore.case = TRUE)) {
        endpoint <- if (nzchar(norm_base)) norm_base else "the configured endpoint"
        friendly <- paste(
          sprintf("The assistant couldn’t reach %s.", endpoint),
          "Ensure the service is reachable and not blocking outbound requests.")
      }

      isolate({
        hist <- chat_history()
        if (length(hist) >= chat_index) {
          hist[[chat_index]]$content <- friendly
          chat_history(hist)
        } else {
          chat_history(append(hist, list(list(role = "assistant", content = friendly))))
        }
      })
      chat_feedback(friendly)
      return()
    }

    response_text <- NULL
    usage_info <- NULL
    if (is.list(response)) {
      response_text <- response$text %||% NULL
      usage_info <- response$usage %||% NULL
    }
    if (is.null(response_text)) {
      response_text <- response
    }
    response_text <- stringr::str_trim(paste(response_text, collapse = "\n"))
    if (!nzchar(response_text)) {
      response_text <- stringr::str_trim(final_text)
    }
    if (!nzchar(response_text)) {
      chat_feedback("Assistant returned an empty response. Try again.")
      isolate({
        hist <- chat_history()
        if (length(hist) >= chat_index) {
          hist[[chat_index]]$content <- "(no response)"
          chat_history(hist)
        }
      })
      return()
    }

    isolate({
      hist <- chat_history()
      if (length(hist) >= chat_index) {
        hist[[chat_index]]$content <- response_text
        chat_history(hist)
      } else {
        chat_history(append(hist, list(list(role = "assistant", content = response_text))))
      }
    })
    chat_feedback(NULL)

    usage_tokens <- NA_integer_
    effective_usage <- final_usage
    if (is.null(effective_usage)) effective_usage <- usage_info
    if (!is.null(effective_usage) && !is.null(effective_usage$total_tokens)) {
      usage_tokens <- suppressWarnings(as.integer(effective_usage$total_tokens))
    }
    if (!is.na(usage_tokens) && usage_tokens > 0L) {
      new_total <- total_tokens_used() + usage_tokens
      total_tokens_used(new_total)
      if (!is.na(token_budget_limit) && new_total >= token_budget_limit) {
        chat_feedback(sprintf(
          "Token budget reached (%s / %s tokens). Start a new session or raise OPENAI_TOKEN_BUDGET to continue.",
          format(new_total, big.mark = ","), format(token_budget_limit, big.mark = ",")
        ))
      }
    }
    chat_history(trim_history(chat_history()))
  })

  # ---------- Spatial Lookup Data ----------
  spatial_lookup <- reactive({
    lookup_path <- file.path("..", "..", "data", "islet_spatial_lookup.csv")
    if (file.exists(lookup_path)) {
      tryCatch({
       
        df <- read.csv(lookup_path, stringsAsFactors = FALSE)
        cat("[SPATIAL] Loaded", nrow(df), "islet spatial coordinates from", lookup_path, "\n")
        df
      }, error = function(e) {
        cat("[SPATIAL] Error loading spatial lookup:", e$message, "\n")
        NULL
      })
    } else {
      cat("[SPATIAL] Lookup file not found:", lookup_path, "\n")
      NULL
    }
  })

  # ---------- Trajectory (AnnData) ----------
  # Try multiple locations for the H5AD; allow ADATA_PATH override
  resolve_traj_path <- function() {
    cand <- c(
      file.path("..", "..", "data", "adata_ins_root.h5ad"),
      file.path("..", "..", "scripts", "adata_ins_root.h5ad"),
      if (!is.null(project_root)) file.path(project_root, "data", "adata_ins_root.h5ad") else NULL,
      if (!is.null(project_root)) file.path(project_root, "scripts", "adata_ins_root.h5ad") else NULL,
      Sys.getenv("ADATA_PATH", unset = "")
    )
    cand <- cand[nzchar(cand)]
    for (p in cand) if (file.exists(p)) return(p)
    NULL
  }
  traj_path <- resolve_traj_path()
  # OPTIMIZATION: Cache trajectory data with better error handling
  traj <- reactiveVal(NULL)
  traj_load_attempted <- reactiveVal(FALSE)
  
  observe({
    if (traj_load_attempted()) return()
    
    cur <- traj_path
    if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
    if (is.null(cur) || !file.exists(cur)) return()
    
    traj_load_attempted(TRUE)
    
    if (!requireNamespace("anndata", quietly = TRUE)) {
      traj(list(error = "R package 'anndata' not installed.", error_detail = NULL))
      return()
    }
    
    ad <- NULL
    load_err <- NULL
    norm_path <- tryCatch(normalizePath(cur, mustWork = TRUE), error = function(e) cur)
    
    ad <- tryCatch(
      anndata::read_h5ad(norm_path), 
      error = function(e) { 
        load_err <<- conditionMessage(e)
        NULL 
      }
    )
    
    if (is.null(ad)) {
      traj(list(error = sprintf("Failed to read H5AD: %s", norm_path), error_detail = load_err))
      return()
    }
    
    # Extract obs, obsm, uns
    obs <- tryCatch(ad$obs, error = function(e) NULL)
    obsm <- tryCatch(ad$obsm, error = function(e) NULL)
    uns <- tryCatch(ad$uns, error = function(e) NULL)
    var_names <- tryCatch(rownames(ad$var), error = function(e) NULL)
    
    if (is.null(obs)) {
      traj(list(error = "Could not extract obs data from AnnData", error_detail = NULL))
      return()
    }
    
    # Ensure obs is a data.frame
    if (!inherits(obs, "data.frame")) {
      obs <- tryCatch(as.data.frame(obs, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(obs)) {
        traj(list(error = "Obs coercion failed", error_detail = NULL))
        return()
      }
    }
    
    # Extract required columns
    cols <- intersect(
      c("combined_islet_id", "base_islet_id", "donor_status", "Case ID", 
        "case_id", "imageid", "dpt_pseudotime", "age", "gender", "total_cells"),
      colnames(obs)
    )
    obs2 <- as.data.frame(obs[, cols, drop = FALSE])
    
    # Parse combined_islet_id for Case ID and islet_key
    if (!"Case ID" %in% names(obs2)) obs2$`Case ID` <- NA_character_
    
    if ("imageid" %in% names(obs2)) {
      obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, as.character(obs2$imageid))
    }
    
    if ("combined_islet_id" %in% names(obs2)) {
      cid <- stringr::str_extract(obs2$combined_islet_id, "^[0-9]{3,4}")
      islet_num <- stringr::str_extract(obs2$combined_islet_id, "(?<=_Islet_)[0-9]{1,3}")
      obs2$`Case ID` <- dplyr::coalesce(obs2$`Case ID`, cid)
      obs2$islet_key <- ifelse(!is.na(islet_num), paste0("Islet_", islet_num), NA_character_)
    }
    
    # Normalize donor_status
    if ("donor_status" %in% names(obs2)) {
      obs2$donor_status <- stringr::str_trim(obs2$donor_status)
      obs2$donor_status <- dplyr::case_when(
        tolower(obs2$donor_status) %in% c("nd", "non-diabetic", "nondiabetic") ~ "ND",
        tolower(obs2$donor_status) %in% c("aab+", "aab", "autoantibody") ~ "Aab+",
        tolower(obs2$donor_status) %in% c("t1d", "type1", "diabetic") ~ "T1D",
        TRUE ~ obs2$donor_status
      )
    }
    
    # Create generic pseudotime column based on method selection
    # Default to PAGA (dpt_pseudotime) if available
    if ("dpt_pseudotime" %in% names(obs2)) {
      obs2$pseudotime <- obs2$dpt_pseudotime
    } else {
      obs2$pseudotime <- NA_real_
    }
    
    traj(list(obs = obs2, obsm = obsm, uns = uns, var_names = var_names, adata = ad))
  })

  # Reactive to handle pseudotime method changes
  traj_with_pseudotime <- reactive({
    tr <- traj()
    if (is.null(tr) || !is.null(tr$error)) return(tr)
    
    obs <- tr$obs
    ad <- tr$adata
    
    # Always use PAGA pseudotime (default)
    if ("dpt_pseudotime" %in% names(obs)) {
      obs$pseudotime <- obs$dpt_pseudotime
    }
    
    # Return updated trajectory data
    list(obs = obs, obsm = tr$obsm, uns = tr$uns, var_names = tr$var_names, adata = ad)
  })

  # Trajectory feature selector UI
  output$traj_feature_selector <- renderUI({
    tr <- traj()
    if (is.null(tr) || !is.null(tr$error)) {
      return(selectInput("traj_feature", "Select Feature:", choices = c("No features available" = ""), selected = ""))
    }
    
    # Get available features from AnnData
    available_features <- NULL
    if (!is.null(tr$var_names)) {
      available_features <- tr$var_names
    } else {
      # Fallback - try to get from the loaded trajectory
      tryCatch({
        adata_path <- resolve_traj_path()
        if (!is.null(adata_path) && file.exists(adata_path)) {
          adata <- anndata::read_h5ad(adata_path)
          available_features <- rownames(adata$var)
        }
      }, error = function(e) NULL)
    }
    
    if (is.null(available_features) || length(available_features) == 0) {
      return(selectInput("traj_feature", "Select Feature:", choices = c("No features available" = ""), selected = ""))
    }
    
    # Organize features into categories
    hormone_markers <- intersect(c("INS", "GCG", "SST"), available_features)
    immune_markers <- intersect(c("CD3e", "CD4", "CD8a", "CD68", "CD163", "CD20", "CD45", "HLADR"), available_features)
    vascular_markers <- intersect(c("CD31", "CD34", "SMA", "ColIV"), available_features)
    neural_markers <- intersect(c("B3TUBB", "GAP43", "PGP9.5"), available_features)
    other_markers <- setdiff(available_features, c(hormone_markers, immune_markers, vascular_markers, neural_markers))
    
    # Build choices list with categories
    choices <- list()
    if (length(hormone_markers) > 0) {
      choices[["Hormone Markers"]] <- setNames(hormone_markers, hormone_markers)
    }
    if (length(immune_markers) > 0) {
      choices[["Immune Markers"]] <- setNames(immune_markers, immune_markers)
    }
    if (length(vascular_markers) > 0) {
      choices[["Vascular Markers"]] <- setNames(vascular_markers, vascular_markers)
    }
    if (length(neural_markers) > 0) {
      choices[["Neural Markers"]] <- setNames(neural_markers, neural_markers)
    }
    if (length(other_markers) > 0) {
      choices[["Other Features"]] <- setNames(other_markers, other_markers)
    }
    
    # Default selection
    default_feature <- if ("INS" %in% available_features) "INS" else available_features[1]
    
    selectInput("traj_feature", "Select Feature:", 
                choices = choices, 
                selected = default_feature)
  })
  
  # Trajectory region selector  
  output$traj_region_selector <- renderUI({
    tagList(
      selectInput("traj_region", "Region:", 
                  choices = c("Islet" = "core", 
                             "Peri-Islet" = "band", 
                             "Islet+20um" = "union"), 
                  selected = "band"),
      tags$div(style = "font-size: 11px; color: #666; margin-top: -10px;",
               "Note: Single-cell trajectory data")
    )
  })
  
  # Build trajectory measurement from AnnData features
  traj_validate <- reactive({
    selected_feature <- input$traj_feature
    if (is.null(selected_feature) || !nzchar(selected_feature)) return(NULL)
    
    # Load AnnData and validate feature exists
    adata_path <- resolve_traj_path()
    if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)
    
    adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
    if (is.null(adata)) return(NULL)
    
    var_names <- rownames(adata$var)
    if (!selected_feature %in% var_names) return(NULL)
    
    # Return success indicator with feature name
    return(list(feature = selected_feature, status = "ready"))
  })

  # FIXED: Optimized trajectory data with spatial coordinates
  traj_data_clean <- reactive({
    ms <- traj_validate()
    if (is.null(ms)) return(NULL)
    
    selected_feature <- ms$feature
    if (is.null(selected_feature)) return(NULL)
    
    adata_path <- resolve_traj_path()
    if (is.null(adata_path) || !file.exists(adata_path)) return(NULL)
    
    tr <- traj_with_pseudotime()
    if (is.null(tr) || !is.null(tr$error)) return(NULL)
    
    # Use cached adata object instead of re-reading
    adata <- tr$adata
    if (is.null(adata)) {
      adata <- tryCatch(anndata::read_h5ad(adata_path), error = function(e) NULL)
      if (is.null(adata)) return(NULL)
    }
    
    var_names <- rownames(adata$var)
    if (!selected_feature %in% var_names) return(NULL)
    
    obs_df <- tr$obs
    if (is.null(obs_df)) return(NULL)
    
    # Extract expression values
    expression_vals <- tryCatch({
      if (selected_feature %in% var_names) {
        idx <- which(var_names == selected_feature)
        as.numeric(adata$X[, idx])
      } else {
        rep(NA_real_, nrow(obs_df))
      }
    }, error = function(e) {
      message("[traj_data_clean] Error extracting feature: ", e$message)
      rep(NA_real_, nrow(obs_df))
    })
    
    # Get donor_status with proper normalization
    donor_clean <- obs_df$donor_status
    if (is.null(donor_clean)) donor_clean <- rep(NA_character_, nrow(obs_df))
    
    # Extract UMAP coordinates
    umap_coords <- NULL
    if (!is.null(tr$obsm) && "X_umap" %in% names(tr$obsm)) {
      umap_coords <- tr$obsm$X_umap
    }
    
    # Get donor ID for coloring
    donor_ids <- as.character(obs_df$imageid)
    
    # Build result data.frame
    islet_keys <- as.character(obs_df$base_islet_id)
    islet_keys <- gsub("^Islet_Islet_", "Islet_", islet_keys)
    
    result_df <- data.frame(
      case_id = as.character(obs_df$imageid),
      islet_key = islet_keys,
      combined_islet_id = as.character(obs_df$combined_islet_id),
      pt = as.numeric(obs_df$pseudotime),
      donor_status = donor_clean,
      donor_id = donor_ids,
      value = as.numeric(expression_vals),
      feature_name = selected_feature,
      stringsAsFactors = FALSE
    )
    
    # Add UMAP coordinates
    if (!is.null(umap_coords) && nrow(umap_coords) == nrow(result_df)) {
      result_df$umap_1 <- as.numeric(umap_coords[, 1])
      result_df$umap_2 <- as.numeric(umap_coords[, 2])
    } else {
      result_df$umap_1 <- NA_real_
      result_df$umap_2 <- NA_real_
    }
    
    # OPTIMIZATION: Use data.table-style merge for diameter (faster than dplyr)
    prep_data <- prepared()
    if (!is.null(prep_data) && !is.null(prep_data$comp)) {
      size_lookup <- prep_data$comp[c("Case ID", "Donor Status", "islet_key", "islet_diam_um")]
      result_df$`Case ID` <- sprintf("%04d", as.numeric(result_df$case_id))
      result_df$`Donor Status` <- result_df$donor_status
      
      # Base R merge (faster for this use case)
      merged <- merge(result_df, size_lookup, 
                     by = c("Case ID", "Donor Status", "islet_key"), 
                     all.x = TRUE, sort = FALSE)
      
      if ("islet_diam_um" %in% colnames(merged)) {
        result_df$islet_diam_um <- merged$islet_diam_um[match(
          paste(result_df$`Case ID`, result_df$`Donor Status`, result_df$islet_key),
          paste(merged$`Case ID`, merged$`Donor Status`, merged$islet_key)
        )]
      }
    }
    
    if (!"islet_diam_um" %in% colnames(result_df)) {
      result_df$islet_diam_um <- NA_real_
    }
    
    # ENHANCEMENT: Add spatial coordinates from segmentation data
    if (!is.null(segmentation_data)) {
      result_df$centroid_x_um <- NA_real_
      result_df$centroid_y_um <- NA_real_
      result_df$has_spatial <- FALSE
      
      for (i in seq_len(nrow(result_df))) {
        annot <- get_islet_annotations(result_df$case_id[i], result_df$islet_key[i])
        if (!is.null(annot)) {
          result_df$centroid_x_um[i] <- annot$centroid_x
          result_df$centroid_y_um[i] <- annot$centroid_y
          result_df$has_spatial[i] <- TRUE
        }
      }
      
      cat("[traj_data_clean] Added spatial coordinates for", 
          sum(result_df$has_spatial), "islets\n")
    }
    
    # Filter valid rows
    valid_idx <- which(
      is.finite(result_df$pt) & 
      !is.na(result_df$value) & 
      !is.na(result_df$case_id) & 
      !is.na(result_df$islet_key)
    )
    
    if (length(valid_idx) == 0) {
      message("[traj_data_clean] No valid rows after filtering")
      return(NULL)
    }
    
    result_df <- result_df[valid_idx, , drop = FALSE]
    
    message("[traj_data_clean] Successfully processed ", nrow(result_df), 
            " observations for ", selected_feature)
    return(result_df)
  })

  output$traj_status <- renderUI({
    cur <- traj_path; if (is.null(cur) || !file.exists(cur)) cur <- resolve_traj_path()
    if (is.null(cur) || !file.exists(cur)) return(tags$div(style = "color:#b00;", "AnnData not found. Place 'adata_ins_root.h5ad' under data/ or scripts/, or set ADATA_PATH."))
    if (!requireNamespace("anndata", quietly = TRUE)) return(tags$div(style = "color:#b00;", "R package 'anndata' not installed. Install with BiocManager::install('anndata')."))
    tr <- traj_with_pseudotime(); if (is.null(tr)) return(tags$div(style = "color:#b00;", "AnnData object not available or H5AD failed to load."))
    if (!is.null(tr$error)) {
      det <- if (!is.null(tr$error_detail) && nzchar(tr$error_detail)) paste0(" Details: ", tr$error_detail) else ""
      install_hint <- if (grepl("anndata", tr$error, ignore.case = TRUE)) " Install with: BiocManager::install('anndata')" else ""
      return(tags$div(style = "color:#b00;", sprintf("Trajectory load error: %s%s%s", tr$error, det, install_hint)))
    }
    ms <- traj_validate(); jn <- traj_data_clean()
    nm <- if (is.null(ms)) 0 else nrow(ms)
    nj <- if (is.null(jn)) 0 else nrow(jn)
    # Distinct joined islets (case_id + islet)
    ndist <- if (is.null(jn)) 0 else nrow(dplyr::distinct(jn, case_id, islet_key))
    src <- if (!is.null(traj_path) && file.exists(traj_path)) traj_path else cur
    
    # Show pseudotime trajectory info
    # Just show data info, not trajectory range details
    tags$div(style = "color:#555;", 
      sprintf("Loaded AnnData (%d obs). Joined %d rows (%d unique islets) of %d islet metrics. Source: %s", 
              nrow(tr$obs), nj, ndist, nm, basename(src))
    )
  })

  output$traj_scatter <- renderPlotly({
    df <- traj_data_clean() 
    if (is.null(df) || nrow(df) == 0) {
      return(plotly_empty() %>% layout(title = "No trajectory data available"))
    }
    
    # Remove outliers (>3 SD from mean) and store for reporting
    df_original <- df
    value_mean <- mean(df$value, na.rm = TRUE)
    value_sd <- sd(df$value, na.rm = TRUE)
    outlier_threshold <- 3
    
    df$is_outlier <- abs(df$value - value_mean) > (outlier_threshold * value_sd)
    outliers <- df[df$is_outlier & !is.na(df$is_outlier), ]
    
    # Store outliers for display in table
    traj_outliers <<- if (nrow(outliers) > 0) {
      data.frame(
        Feature = outliers$feature_name,
        Case_ID = outliers$case_id,
        Islet = outliers$islet_key,
        Donor_Status = outliers$donor_status,
        Pseudotime = round(outliers$pt, 3),
        Value = round(outliers$value, 3),
        Z_Score = round((outliers$value - value_mean) / value_sd, 2),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
    
    # Remove outliers from plotting data
    df <- df[!df$is_outlier | is.na(df$is_outlier), ]
    
    cat(sprintf("[TRAJECTORY] Removed %d outliers (>3 SD)\n", nrow(df_original) - nrow(df)))
    
    # Determine color mapping and aesthetics
    color_by <- input$traj_color_by %||% "donor_status"
    size_by <- input$traj_point_size %||% "uniform"
    trend_type <- input$traj_show_trend %||% "by_donor"
    alpha_val <- input$traj_alpha %||% 0.3
    point_size <- input$traj_point_size_slider %||% 1.0
    # Always apply jitter (jitter checkbox removed from UI)
    use_jitter <- TRUE
    
    # Apply jitter to reduce overlap if requested
    if (use_jitter) {
      # Add small random noise to both x and y
      set.seed(42)  # Reproducible jitter
      jitter_amount_x <- diff(range(df$pt, na.rm = TRUE)) * 0.01  # 1% of x range
      jitter_amount_y <- diff(range(df$value, na.rm = TRUE)) * 0.02  # 2% of y range
      df$pt <- df$pt + runif(nrow(df), -jitter_amount_x, jitter_amount_x)
      df$value <- df$value + runif(nrow(df), -jitter_amount_y, jitter_amount_y)
    }
    
    # Create base aesthetic mapping
    aes_mapping <- aes(x = pt, y = value)
    
    # Add color aesthetic
    if (color_by == "donor_status") {
      df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
      aes_mapping$colour <- as.name("donor_status")
      color_title <- "Donor Status"
    } else if (color_by == "donor_id") {
      # Convert to factor for consistent coloring
      df$donor_id <- factor(df$donor_id)
      aes_mapping$colour <- as.name("donor_id")
      color_title <- "Donor ID"
    }
    
    # Add size aesthetic if not uniform
    if (size_by != "uniform") {
      if (size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
        # Check if diameter data is available for sizing
        if (!all(is.na(df$islet_diam_um))) {
          aes_mapping$size <- as.name("islet_diam_um")
          size_title <- "Islet Diameter"
        }
      }
    }
    
    # Create plot with user-controlled transparency and size
    # Smaller default point size for less overlap
    g <- ggplot(df, aes_mapping)
    if (size_by == "uniform") {
      g <- g + geom_point(alpha = alpha_val, size = point_size, stroke = 0)
    } else {
      g <- g + geom_point(alpha = alpha_val, stroke = 0)
    }
    g <- g +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "none",  # Legend moved to separate UI output in column 4
        plot.margin = unit(c(5, 5, 5, 5), "pt"),
        plot.title = element_text(margin = margin(b = 15))  # Space below title
      )
    
    # Apply color scales
    if (color_by == "donor_status") {
      g <- g + scale_color_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd"),
                                  name = color_title, drop = FALSE)
    } else if (color_by == "donor_id") {
      # Convert to factor for consistent coloring
      df$donor_id <- factor(df$donor_id)
      
      # Use Paired palette for better discrimination of donors
      # This gives maximum contrast between adjacent colors
      n_donors <- length(levels(df$donor_id))
      if (n_donors <= 12) {
        # Use Paired palette (max 12 colors, highly distinctive)
        colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
        names(colors) <- levels(df$donor_id)[1:length(colors)]
      } else {
        # Fall back to rainbow for many donors
        colors <- rainbow(n_donors, s = 0.8, v = 0.8)
        names(colors) <- levels(df$donor_id)
      }
      g <- g + scale_color_manual(values = colors, name = color_title, drop = FALSE)
    }
    
    # Add size scale - hide legend for size (only show color legend for Donor Status or ID)
    # Make the size range proportional to the point_size slider
    if (size_by != "uniform" && size_by == "islet_diam_um" && "islet_diam_um" %in% colnames(df)) {
      if (!all(is.na(df$islet_diam_um))) {
        # Scale the range based on point_size slider (default 1.0)
        # Base range is [0.5, 2.5], multiply by point_size to allow proportional scaling
        size_range <- c(0.5 * point_size, 2.5 * point_size)
        g <- g + scale_size_continuous(name = "Islet Diameter (µm)", range = size_range, guide = "none")
      }
    }
    
    # Get data range for constraining trend lines
    pt_range <- range(df$pt, na.rm = TRUE)
    
    # Add trend lines - ALWAYS use donor_status colors regardless of point coloring
    if (trend_type == "overall") {
      g <- g + geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "black", 
                          size = 1, fullrange = FALSE, span = 0.75)
    } else if (trend_type == "by_donor") {
      # Ensure donor_status factor exists for trend lines
      if (!"donor_status" %in% names(df) || all(is.na(df$donor_status))) {
        # Skip trend lines if no donor status
      } else {
        df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
        
        # Define trend line colors
        trend_colors <- c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")
        
        # When coloring by donor_id, we need to add trend lines with explicit colors
        # to avoid scale conflicts
        if (color_by == "donor_id") {
          # Add separate geom_smooth for each donor_status with explicit color
          for (status in c("ND", "Aab+", "T1D")) {
            df_subset <- df[df$donor_status == status & !is.na(df$donor_status), ]
            if (nrow(df_subset) > 0) {
              g <- g + geom_smooth(data = df_subset,
                                  aes(x = pt, y = value),
                                  method = "loess", se = TRUE, alpha = 0.15,
                                  size = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
                                  color = trend_colors[status])
            }
          }
        } else {
          # When coloring by donor_status, use the scale approach
          g <- g + geom_smooth(
                              method = "loess", se = TRUE, alpha = 0.15, 
                              size = 0.8, show.legend = FALSE, fullrange = FALSE, span = 0.75,
                              inherit.aes = FALSE,
                              mapping = aes(x = pt, y = value, group = donor_status, color = donor_status)) +
            scale_color_manual(values = trend_colors,
                             name = if (color_by == "donor_status") color_title else "Trend Lines",
                             drop = FALSE,
                             guide = if (color_by == "donor_status") "legend" else "none")
        }
      }
    }
    
    # Get current selection context for labels
    selected_feature <- input$traj_feature %||% "Selected feature"
    metric_label <- "Expression Level"
    
    # Dynamic x-axis label based on pseudotime method
    x_label <- "Pseudotime (PAGA trajectory, INS-rooted)"
    
    g <- g + labs(
      x = x_label,
      y = paste(selected_feature, metric_label),
      title = NULL
    ) + 
    coord_cartesian(xlim = c(0, 1))  # Set fixed range from 0 to 1
    
    # Apply dark theme if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0"),
        legend.background = element_rect(fill = "#111111", color = NA)
      )
    }
    
    # Create plotly with custom data for reliable click handling
    p <- ggplotly(g, tooltip = c("x", "y", "colour"), source = "traj_scatter")

    # CRITICAL FIX: Add customdata to plotly traces for click handling
    if ("combined_islet_id" %in% colnames(df) && "case_id" %in% colnames(df)) {
      # Store the original dataframe in a global variable for coordinate-based lookup
      traj_coord_lookup <<- data.frame(
        pt = df$pt,
        value = df$value,
        case_id = df$case_id,
        islet_key = df$islet_key,
        combined_islet_id = df$combined_islet_id,
        donor_status = df$donor_status,
        stringsAsFactors = FALSE
      )

      # FIXED: Match customdata to actual trace data by coordinates
      # Plotly splits data into separate traces by color groups
      for (i in seq_along(p$x$data)) {
        if (!is.null(p$x$data[[i]]$type) && p$x$data[[i]]$type == "scatter") {
          trace_x <- p$x$data[[i]]$x
          trace_y <- p$x$data[[i]]$y

          if (length(trace_x) > 0 && length(trace_y) > 0) {
            # Match each point in this trace to the original dataframe
            trace_customdata <- matrix(NA, nrow = length(trace_x), ncol = 4)

            for (j in seq_along(trace_x)) {
              # Find matching point in original data (use tolerance for jittered coordinates)
              tolerance <- 0.001
              matches <- which(
                abs(df$pt - trace_x[j]) < tolerance &
                abs(df$value - trace_y[j]) < tolerance
              )

              if (length(matches) > 0) {
                idx <- matches[1]
                trace_customdata[j, ] <- c(
                  as.character(df$case_id[idx]),
                  as.character(df$islet_key[idx]),
                  as.character(df$combined_islet_id[idx]),
                  as.character(df$donor_status[idx])
                )
              }
            }

            p$x$data[[i]]$customdata <- trace_customdata
            cat("[PLOT] Added customdata to trace", i, "with", length(trace_x), "points\n")
          }
        }
      }

      cat("[PLOT] Stored coordinate lookup table with", nrow(traj_coord_lookup), "entries\n")
    }

    # Apply layout and register click events (only once)
    p %>%
      layout(
        showlegend = FALSE,  # Legend hidden from plot, displayed in column 4 UI element
        title = NULL,  # No title
        margin = list(l = 60, r = 20, t = 20, b = 50, pad = 5),  # Reduced top margin since no title or legend
        xaxis = list(
          fixedrange = FALSE,
          automargin = TRUE,
          range = c(0, 1)  # Fixed range from 0 to 1
        ),
        yaxis = list(fixedrange = FALSE, automargin = TRUE)
      ) %>%
      event_register("plotly_click")
  })

  # ===== Islet Segmentation Viewer =====
  # Reactive value to store selected islet information
  selected_islet <- reactiveVal(NULL)

  # Click handler for trajectory scatter plot
  observeEvent(event_data("plotly_click", source = "traj_scatter"), {
    click_data <- event_data("plotly_click", source = "traj_scatter")
    if (is.null(click_data)) return()

    # Check if segmentation viewer is available
    if (!SF_AVAILABLE) {
      showNotification("Segmentation viewer requires the 'sf' package. Install with: install.packages('sf')",
                       type = "warning", duration = 5)
      return()
    }

    if (is.null(islet_spatial_lookup)) {
      showNotification("Islet spatial lookup data not available", type = "warning", duration = 5)
      return()
    }

    cat("[SEGMENTATION CLICK] Received click event\n")

    # Get clicked coordinates
    click_x <- click_data$x
    click_y <- click_data$y

    # First try to get from customdata
    custom <- click_data$customdata
    case_id <- NULL
    islet_key <- NULL

    if (!is.null(custom) && length(custom) >= 2) {
      case_id <- custom[[1]]
      islet_key <- custom[[2]]
      cat("[SEGMENTATION CLICK] From customdata: case_id=", case_id, ", islet_key=", islet_key, "\n")
    }

    # If customdata not available, look up from coordinate table
    if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
      if (exists("traj_coord_lookup") && !is.null(traj_coord_lookup) && nrow(traj_coord_lookup) > 0) {
        # Find nearest point in lookup table
        tolerance <- 0.01
        matches <- which(
          abs(traj_coord_lookup$pt - click_x) < tolerance &
          abs(traj_coord_lookup$value - click_y) < tolerance
        )

        if (length(matches) > 0) {
          idx <- matches[1]
          case_id <- traj_coord_lookup$case_id[idx]
          islet_key <- traj_coord_lookup$islet_key[idx]
          cat("[SEGMENTATION CLICK] From coord lookup: case_id=", case_id, ", islet_key=", islet_key, "\n")
        }
      }
    }

    if (is.null(case_id) || is.null(islet_key) || is.na(case_id) || is.na(islet_key)) {
      showNotification("Could not identify clicked islet", type = "warning", duration = 3)
      return()
    }

    # Look up centroid coordinates from spatial lookup
    spatial_match <- islet_spatial_lookup[
      islet_spatial_lookup$case_id == case_id & islet_spatial_lookup$islet_key == islet_key,
    ]

    if (nrow(spatial_match) == 0) {
      showNotification(paste("Spatial data not found for", islet_key, "in case", case_id),
                       type = "warning", duration = 3)
      return()
    }

    centroid_x <- spatial_match$centroid_x_um[1]
    centroid_y <- spatial_match$centroid_y_um[1]

    cat("[SEGMENTATION CLICK] Centroid: x=", centroid_x, ", y=", centroid_y, "\n")

    # Store selected islet info - this triggers the embedded viewer to update
    selected_islet(list(
      case_id = case_id,
      islet_key = islet_key,
      centroid_x = centroid_x,
      centroid_y = centroid_y
    ))
  })

  # Click handler for Plot tab scatter plot (when "Show individual points" is enabled)
  observeEvent(event_data("plotly_click", source = "plot_scatter"), {
    click_data <- event_data("plotly_click", source = "plot_scatter")
    if (is.null(click_data)) return()

    # Check if segmentation viewer is available
    if (!SF_AVAILABLE) {
      showNotification("Segmentation viewer requires the 'sf' package. Install with: install.packages('sf')",
                       type = "warning", duration = 5)
      return()
    }

    if (is.null(islet_spatial_lookup)) {
      showNotification("Islet spatial lookup data not available", type = "warning", duration = 5)
      return()
    }

    # Get the key from click data (format: "case_id|islet_key")
    click_key <- click_data$key
    if (is.null(click_key) || is.na(click_key) || !nzchar(click_key)) {
      showNotification("Click on individual points to view segmentation (enable 'Show individual points')",
                       type = "message", duration = 3)
      return()
    }

    cat("[PLOT CLICK] Received click with key:", click_key, "\n")

    # Parse the key to extract case_id and islet_key
    key_parts <- strsplit(click_key, "\\|")[[1]]
    if (length(key_parts) < 2) {
      showNotification("Could not parse islet identifier", type = "warning", duration = 3)
      return()
    }

    case_id <- key_parts[1]
    islet_key <- key_parts[2]

    cat("[PLOT CLICK] Parsed: case_id=", case_id, ", islet_key=", islet_key, "\n")

    # Look up centroid coordinates from spatial lookup
    spatial_match <- islet_spatial_lookup[
      islet_spatial_lookup$case_id == case_id & islet_spatial_lookup$islet_key == islet_key,
    ]

    if (nrow(spatial_match) == 0) {
      # Try with numeric case_id conversion
      case_id_num <- suppressWarnings(as.numeric(case_id))
      if (!is.na(case_id_num)) {
        spatial_match <- islet_spatial_lookup[
          islet_spatial_lookup$case_id == case_id_num & islet_spatial_lookup$islet_key == islet_key,
        ]
      }
    }

    if (nrow(spatial_match) == 0) {
      showNotification(paste("Spatial data not found for", islet_key, "in case", case_id),
                       type = "warning", duration = 3)
      return()
    }

    centroid_x <- spatial_match$centroid_x_um[1]
    centroid_y <- spatial_match$centroid_y_um[1]

    cat("[PLOT CLICK] Centroid: x=", centroid_x, ", y=", centroid_y, "\n")

    # Store selected islet info - this triggers the segmentation viewer to update
    selected_islet(list(
      case_id = case_id,
      islet_key = islet_key,
      centroid_x = centroid_x,
      centroid_y = centroid_y
    ))

    # Show modal dialog with segmentation view (Plot tab uses modal instead of embedded panel)
    # Uses the same islet_segmentation_view output that works for the Trajectory tab
    showModal(modalDialog(
      title = paste("Islet Segmentation:", islet_key, "(Case", case_id, ")"),
      plotOutput("islet_segmentation_view", height = "450px"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # Render the embedded segmentation viewer panel (appears when islet is selected)
  output$segmentation_viewer_panel <- renderUI({
    info <- selected_islet()
    if (is.null(info)) return(NULL)

    fluidRow(
      column(12,
        div(class = "card", style = "padding: 15px; margin-top: 20px; border: 2px solid #0066CC;",
          fluidRow(
            column(8,
              h4(paste("Islet Segmentation:", info$islet_key, "(Case", info$case_id, ")"),
                 style = "margin-top: 0; color: #0066CC;")
            ),
            column(4, style = "text-align: right;",
              actionButton("clear_segmentation", "Close", class = "btn btn-sm btn-outline-secondary")
            )
          ),
          fluidRow(
            # Left: Segmentation plot
            column(8,
              plotOutput("islet_segmentation_view", height = "450px")
            ),
            # Right: Legend and info
            column(4,
              div(style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px; height: 100%;",
                h5("Legend", style = "margin-top: 0; margin-bottom: 15px;"),
                div(style = "margin-bottom: 10px;",
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 25px; height: 4px; background-color: #0066CC; margin-right: 8px;"),
                    span("Islet boundary")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 25px; height: 4px; background-color: #00CCCC; margin-right: 8px;"),
                    span("Expanded (+20\u00b5m)")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 25px; height: 4px; background-color: #CC00CC; margin-right: 8px;"),
                    span("Nerve")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 25px; height: 4px; background-color: #CC0000; margin-right: 8px;"),
                    span("Capillary")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 25px; height: 4px; background-color: #00AA00; margin-right: 8px;"),
                    span("Lymphatic")
                  ),
                  div(style = "display: flex; align-items: center; margin-bottom: 8px;",
                    div(style = "width: 25px; height: 4px; background-color: #FFD700; margin-right: 8px;"),
                    span("Selected islet", style = "font-weight: bold;")
                  )
                ),
                hr(),
                h6("Islet Info", style = "margin-bottom: 10px;"),
                tags$dl(style = "font-size: 13px;",
                  tags$dt("Case ID:"), tags$dd(info$case_id),
                  tags$dt("Islet:"), tags$dd(info$islet_key),
                  tags$dt("Centroid X:"), tags$dd(paste0(round(info$centroid_x, 1), " \u00b5m")),
                  tags$dt("Centroid Y:"), tags$dd(paste0(round(info$centroid_y, 1), " \u00b5m"))
                ),
                hr(),
                p(style = "font-size: 11px; color: #666; margin-bottom: 0;",
                  "Click another point to view a different islet, or click Close to hide this panel.")
              )
            )
          )
        )
      )
    )
  })

  # Clear segmentation viewer when close button is clicked
  observeEvent(input$clear_segmentation, {
    selected_islet(NULL)
  })

  # Render islet segmentation plot
  output$islet_segmentation_view <- renderPlot({
    req(selected_islet())
    info <- selected_islet()

    cat("[SEGMENTATION RENDER] Rendering for case=", info$case_id, ", islet=", info$islet_key, "\n")

    # Load GeoJSON for this case
    geojson <- load_case_geojson(info$case_id)
    if (is.null(geojson)) {
      # Return empty plot with message
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("GeoJSON data not available for case", info$case_id),
                   size = 5, color = "gray50") +
          theme_void() +
          xlim(0, 1) + ylim(0, 1)
      )
    }

    # Get polygons in region around islet centroid
    buffer_um <- 250  # Show 250um around centroid
    polygons <- get_islet_region_polygons(geojson, info$centroid_x, info$centroid_y, buffer_um)

    # Calculate bounding box in pixels for proper axis limits
    centroid_x_px <- info$centroid_x / PIXEL_SIZE_UM
    centroid_y_px <- info$centroid_y / PIXEL_SIZE_UM
    buffer_px <- buffer_um / PIXEL_SIZE_UM

    xlim_range <- c(centroid_x_px - buffer_px, centroid_x_px + buffer_px)
    ylim_range <- c(centroid_y_px - buffer_px, centroid_y_px + buffer_px)

    # Define colors for each class
    colors <- c(
      "Islet" = "#0066CC",
      "IsletExpanded" = "#00CCCC",
      "Nerve" = "#CC00CC",
      "Capillary" = "#CC0000",
      "Lymphatic" = "#00AA00"
    )

    # Build plot with explicit limits and Y-axis reversed to match image coordinates
    p <- ggplot() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      coord_sf(xlim = xlim_range, ylim = rev(ylim_range), expand = FALSE)

    # Track if any polygons were added
    has_polygons <- FALSE

    # Add each class with appropriate styling
    # Order: IsletExpanded first (background), then Islet, then vessels/nerves on top
    layer_order <- c("IsletExpanded", "Islet", "Lymphatic", "Capillary", "Nerve")

    for (cls in layer_order) {
      if (!is.null(polygons[[cls]]) && nrow(polygons[[cls]]) > 0) {
        p <- p + geom_sf(data = polygons[[cls]],
                         fill = NA, color = colors[cls], linewidth = 0.8,
                         show.legend = FALSE)
        has_polygons <- TRUE
      }
    }

    if (!has_polygons) {
      # No polygons found in region
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste("No segmentation data found near", info$islet_key),
                   size = 5, color = "gray50") +
          theme_void() +
          xlim(0, 1) + ylim(0, 1)
      )
    }

    # Highlight the clicked islet by matching name
    if (!is.null(polygons$Islet) && nrow(polygons$Islet) > 0) {
      # Check if there's a "name" column or "id" that matches islet_key
      islet_names <- if ("name" %in% names(polygons$Islet)) {
        polygons$Islet$name
      } else if ("id" %in% names(polygons$Islet)) {
        polygons$Islet$id
      } else {
        NULL
      }

      if (!is.null(islet_names)) {
        # Try exact match first
        clicked_idx <- which(islet_names == info$islet_key)
        # If no exact match, try pattern match
        if (length(clicked_idx) == 0) {
          clicked_idx <- grep(info$islet_key, islet_names, fixed = TRUE)
        }

        if (length(clicked_idx) > 0) {
          clicked_islet <- polygons$Islet[clicked_idx[1], ]
          p <- p + geom_sf(data = clicked_islet, fill = NA,
                           color = "#FFD700", linewidth = 2.5,
                           show.legend = FALSE)
        }
      }
    }

    # Add crosshairs at centroid (centroid_x_px and centroid_y_px calculated above)
    p <- p +
      geom_vline(xintercept = centroid_x_px, color = "gray50", linetype = "dashed", linewidth = 0.3) +
      geom_hline(yintercept = centroid_y_px, color = "gray50", linetype = "dashed", linewidth = 0.3)

    # Add labels (coordinates displayed in pixels to match QuPath)
    p <- p +
      labs(
        title = paste(info$islet_key, "- Case", info$case_id),
        x = "X position (pixels)",
        y = "Y position (pixels)"
      )

    p
  })

  # Outlier table display - hidden by default, shown only when checkbox is checked
  output$traj_outlier_info <- renderUI({
    # Only show if checkbox is checked AND outliers exist
    show_table <- isTRUE(input$show_outlier_table)
    
    if (show_table && exists("traj_outliers") && !is.null(traj_outliers) && nrow(traj_outliers) > 0) {
      tagList(
        div(style = "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; padding: 10px; margin-top: 10px;",
          h5(style = "color: #856404; margin-top: 0;", 
             sprintf("⚠️ %d Outlier%s Removed (>3 SD)", nrow(traj_outliers), ifelse(nrow(traj_outliers) > 1, "s", ""))),
          p(style = "color: #856404; font-size: 12px; margin-bottom: 10px;",
            "The following data points were excluded from the plot because they exceed 3 standard deviations from the mean:"),
          div(style = "max-height: 200px; overflow-y: auto;",
            renderTable({
              traj_outliers
            }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = 'xs', width = "100%")
          )
        )
      )
    } else {
      NULL
    }
  })

  # Legend output for trajectory scatterplot
  output$traj_legend <- renderUI({
    color_by <- input$traj_color_by %||% "donor_status"
    
    if (color_by == "donor_status") {
      # Donor Status legend
      tagList(
        div(style = "display: flex; align-items: center; margin-bottom: 8px;",
          div(style = "width: 20px; height: 20px; background-color: #2ca02c; border-radius: 3px; margin-right: 8px;"),
          span("ND", style = "font-size: 13px;")
        ),
        div(style = "display: flex; align-items: center; margin-bottom: 8px;",
          div(style = "width: 20px; height: 20px; background-color: #ffcc00; border-radius: 3px; margin-right: 8px;"),
          span("Aab+", style = "font-size: 13px;")
        ),
        div(style = "display: flex; align-items: center;",
          div(style = "width: 20px; height: 20px; background-color: #9467bd; border-radius: 3px; margin-right: 8px;"),
          span("T1D", style = "font-size: 13px;")
        )
      )
    } else if (color_by == "donor_id") {
      # Get unique donor IDs from the data
      df <- traj_data_clean()
      if (!is.null(df) && nrow(df) > 0) {
        donor_ids <- sort(unique(df$donor_id))
        n_donors <- length(donor_ids)
        
        # Generate colors matching the plot
        if (n_donors <= 12) {
          colors <- RColorBrewer::brewer.pal(min(12, max(3, n_donors)), "Paired")
        } else {
          colors <- rainbow(n_donors, s = 0.8, v = 0.8)
        }
        
        # Split into two columns
        mid_point <- ceiling(n_donors / 2)
        col1_ids <- donor_ids[1:mid_point]
        col2_ids <- if (n_donors > mid_point) donor_ids[(mid_point + 1):n_donors] else NULL
        
        # Create legend items for column 1
        col1_items <- lapply(seq_along(col1_ids), function(i) {
          div(style = "display: flex; align-items: center; margin-bottom: 6px;",
            div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[i]),
            ),
            span(col1_ids[i], style = "font-size: 12px;")
          )
        })
        
        # Create legend items for column 2
        col2_items <- if (!is.null(col2_ids)) {
          lapply(seq_along(col2_ids), function(i) {
            idx <- mid_point + i
            div(style = "display: flex; align-items: center; margin-bottom: 6px;",
              div(style = sprintf("width: 18px; height: 18px; background-color: %s; border-radius: 3px; margin-right: 6px;", colors[idx]),
              ),
              span(col2_ids[i], style = "font-size: 12px;")
            )
          })
        } else {
          NULL
        }
        
        # Create two-column layout
        tagList(
          div(style = "display: flex; gap: 10px;",
            div(style = "flex: 1;", col1_items),
            if (!is.null(col2_items)) div(style = "flex: 1;", col2_items)
          )
        )
      } else {
        p("No data available", style = "font-size: 12px; color: #999;")
      }
    }
  })

  # UMAP plot colored by donor status
  output$traj_umap_donor <- renderPlot({
    df <- traj_data_clean()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Check if UMAP coordinates are available
    if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
      return(ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
        theme_void())
    }
    
    # Create UMAP plot colored by donor status
    df$donor_status <- factor(df$donor_status, levels = c("ND", "Aab+", "T1D"))
    
    g <- ggplot(df, aes(x = umap_1, y = umap_2, color = donor_status)) +
      geom_point(alpha = 0.7, size = 1.2) +
      scale_color_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")) +
      labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP: Donor Status") +
      theme_minimal(base_size = 12) +
      theme(aspect.ratio = 1)
    
    # Apply dark theme if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0"),
        legend.background = element_rect(fill = "#111111", color = NA)
      )
    }
    
    return(g)
  })
  
  # UMAP plot colored by selected feature
  output$traj_umap_feature <- renderPlot({
    df <- traj_data_clean()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Check if UMAP coordinates are available
    if (all(is.na(df$umap_1)) || all(is.na(df$umap_2))) {
      return(ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "UMAP coordinates not available", size = 5) +
        theme_void())
    }
    
    selected_feature <- input$traj_feature %||% "Selected feature"
    
    # Create UMAP plot colored by feature expression
    g <- ggplot(df, aes(x = umap_1, y = umap_2, color = value)) +
      geom_point(alpha = 0.7, size = 1.2) +
      labs(x = "UMAP 1", y = "UMAP 2", 
           title = paste("UMAP:", selected_feature),
           color = "Expression") +
      theme_minimal(base_size = 12) +
      theme(aspect.ratio = 1)
    
    # Use viridis or gradient color scale
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      g <- g + scale_color_viridis_c(option = "plasma", na.value = "#bbbbbb")
    } else {
      g <- g + scale_color_gradient(low = "#132B43", high = "#56B1F7", na.value = "#bbbbbb")
    }
    
    # Apply dark theme if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0"),
        legend.background = element_rect(fill = "#111111", color = NA)
      )
    }
    
    return(g)
  })

  output$traj_heatmap <- renderPlot({
    df <- traj_data_clean() 
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # Always show binned average (checkbox removed)
    
    # Create binned analysis
    enc <- function(s) ifelse(s == "ND", 0, ifelse(s == "Aab+", 1, ifelse(s == "T1D", 2, NA_real_)))
    df$ds_code <- enc(df$donor_status)
    
    # Normalize pseudotime to 0-1 range for binning
    pt_range <- range(df$pt, na.rm = TRUE)
    df$pt_norm <- (df$pt - pt_range[1]) / (pt_range[2] - pt_range[1])
    
    nb <- 25  # Number of bins
    brks <- seq(0, 1, length.out = nb + 1)
    df$pt_bin <- cut(pmax(0, pmin(1, df$pt_norm)), breaks = brks, include.lowest = TRUE, right = FALSE)
    
    # Calculate averages per bin using base R
    bin_levels <- levels(df$pt_bin)
    hm_list <- list()
    
    for (i in seq_along(bin_levels)) {
        bin_name <- bin_levels[i]
        bin_data <- df[!is.na(df$pt_bin) & df$pt_bin == bin_name, ]
        
        if (nrow(bin_data) >= 3) {  # Only bins with at least 3 observations
          hm_list[[length(hm_list) + 1]] <- data.frame(
            pt_bin = factor(bin_name, levels = bin_levels),
            avg_donor_status = mean(bin_data$ds_code, na.rm = TRUE),
            count = nrow(bin_data),
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(hm_list) == 0) return(NULL)
      hm <- do.call(rbind, hm_list)
      
      if (nrow(hm) == 0) return(NULL)
      
      # Calculate bin midpoints in actual pseudotime range
      mids <- head(brks, -1) + diff(brks)[1]/2
      hm$x <- mids * (pt_range[2] - pt_range[1]) + pt_range[1]
      
      # Create heatmap
      g <- ggplot(hm, aes(x = x, y = 1, fill = avg_donor_status)) +
        geom_tile(height = 1) +
        scale_fill_gradientn(
          colors = c("#2ca02c", "#ffcc00", "#9467bd"),
          limits = c(0, 2), 
          na.value = "#dddddd",
          name = "Average\nDonor Type",
          breaks = c(0, 1, 2),
          labels = c("ND", "Aab+", "T1D")
        ) +
        scale_x_continuous(limits = c(pt_range[1], pt_range[2]), expand = c(0, 0)) +
        labs(x = "Pseudotime", y = "", title = "Donor Status Progression Along Pseudotime") +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 12, hjust = 0.5)
        )
      
      # Apply dark theme if selected
      if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
        g <- g + theme(
          plot.background = element_rect(fill = "#000000", colour = NA),
          panel.background = element_rect(fill = "#000000", colour = NA),
          axis.text = element_text(color = "#e6e6e6"),
          axis.title = element_text(color = "#f0f0f0"),
          plot.title = element_text(color = "#f0f0f0"),
          legend.text = element_text(color = "#e6e6e6"),
          legend.title = element_text(color = "#f0f0f0"),
          legend.background = element_rect(fill = "#111111", color = NA)
        )
      }
      
      return(g)
  })

  # Selected islet viewer removed
 
  # Selected islet vitessce viewer removed
 
  viewer_info <- reactive({
    # Get forced image value if set (for trajectory integration)
    forced_img <- forced_image()

    # Detect environment to determine URL construction strategy
    env_info <- detect_environment(session)

    base <- resolve_avivator_base()
    info <- list(
      base = base,
      selection = NULL,
      mode = "picker",
      ok = FALSE,
      iframe_src = NULL,
      image_url = NULL,
      image_public_url = NULL,
      image_rel = NULL,
      image_app_url = NULL,
      asset_diag = NULL,
      env = env_info
    )
    if (is.null(base)) return(info)
    cat("[VIEWER] Environment detection - reverse_proxy:", env_info$is_reverse_proxy, 
        "hostname:", env_info$hostname, "pathname:", env_info$pathname, "\n")
    
    params <- list()
    # Select image: prefer forced image first, then user selection
    rel_url <- NULL
    
    # Check for forced image first (use stored value to ensure consistency)
    sel_basename <- forced_img
    if (is.null(sel_basename)) {
      # Fall back to input selection
      sel_basename <- input$selected_image
    }
    
    cat("[VIEWER] === Image Selection Debug ===\n")
    cat("[VIEWER] Forced image:", if(is.null(forced_img)) "NULL" else paste0("'", forced_img, "'"), "\n")
    cat("[VIEWER] Input image:", if(is.null(input$selected_image)) "NULL" else paste0("'", input$selected_image, "'"), "\n")
    cat("[VIEWER] Selected basename:", if(is.null(sel_basename)) "NULL" else paste0("'", sel_basename, "'"), "\n")
    cat("[VIEWER] === End Selection Debug ===\n")
    
    if (!is.null(sel_basename) && nzchar(sel_basename)) {
      rel_url <- NULL
      www_candidate <- file.path("www", "local_images", sel_basename)
      resource_candidate <- NULL
      if (!is.null(local_images_root) && dir.exists(local_images_root)) {
        resource_candidate <- file.path(local_images_root, sel_basename)
      }

      if (file.exists(www_candidate)) {
        rel_url <- paste("local_images", sel_basename, sep = "/")
        cat("[VIEWER] Using www/local_images asset:", www_candidate, "\n")
      } else if (!is.null(resource_candidate) && file.exists(resource_candidate)) {
        rel_url <- paste("images", sel_basename, sep = "/")
        cat("[VIEWER] Using /images resource asset:", resource_candidate, "\n")
      }

      if (!is.null(rel_url)) {
        info$selection <- sel_basename
        cat("[VIEWER] Using rel URL:", rel_url, "\n")
      } else {
        fallback_rel <- paste("local_images", sel_basename, sep = "/")
        rel_url <- fallback_rel
        info$selection <- sel_basename
        cat("[VIEWER] WARNING: image not found on disk, falling back to rel:", rel_url, "\n")
      }
    }
    
    asset_diag <- probe_viewer_asset(session, rel_url)
    info$asset_diag <- asset_diag
    if (!is.null(asset_diag$rel_path)) {
      cat("[VIEWER] Asset diag: exists=", asset_diag$file_exists,
          "size=", asset_diag$file_size,
          "http=", asset_diag$http_status %||% "NA",
          if (!is.null(asset_diag$http_error)) paste("err=", asset_diag$http_error) else "",
          "\n")
    }
    
    resolved_app_url <- NULL
    resolved_public_url <- NULL
    if (!is.null(rel_url) && nzchar(rel_url)) {
      info$image_rel <- rel_url
      resolved_app_url <- build_app_absolute_url(session, rel_url)
      resolved_public_url <- build_public_http_url(session, rel_url)
      info$image_app_url <- resolved_app_url
      info$image_public_url <- resolved_public_url
      cat("[VIEWER] App URL:", resolved_app_url %||% "NULL", "\n")
      cat("[VIEWER] Public URL:", resolved_public_url %||% "NULL", "\n")
    }

    sel_url <- NULL
    if (!is.null(resolved_public_url) && nzchar(resolved_public_url)) {
      sel_url <- resolved_public_url
    } else if (!is.null(resolved_app_url) && nzchar(resolved_app_url)) {
      sel_url <- resolved_app_url
    } else if (!is.null(default_image_url) && nzchar(default_image_url)) {
      sel_url <- default_image_url
    }

    if (!is.null(sel_url)) {
      # Keep slashes unencoded (reserved=FALSE) so relative paths work correctly
      params[["image_url"]] <- utils::URLencode(sel_url, reserved = FALSE)
      info$image_param <- sel_url
      if (is.null(info$image_public_url)) {
        info$image_public_url <- if (!identical(sel_url, resolved_app_url)) sel_url else NULL
      }
      info$image_url <- resolved_app_url %||% sel_url
    }
    # Channel config based on Channel_names mapping (base64-encoded per viewer expectation)
    ch_b64 <- tryCatch(build_channel_config_b64(channel_names_vec), error = function(e) NULL)
    if (!is.null(ch_b64) && nzchar(ch_b64)) {
      # Encode base64 for transport; viewer decodes via decodeURIComponent + atob
      params[["channel_config"]] <- utils::URLencode(ch_b64, reserved = TRUE)
    }
    
    query <- NULL
    if (length(params)) {
      parts <- vapply(names(params), function(nm) sprintf("%s=%s", nm, params[[nm]]), character(1))
      query <- paste(parts, collapse = "&")
    }
    info$iframe_src <- if (!is.null(query)) paste0(base, "?", query) else base
    info$ok <- TRUE
    info
  })

  observe({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) {
      shinyjs::runjs("document.body.classList.add('viewer-mode');")
    } else {
      shinyjs::runjs("document.body.classList.remove('viewer-mode');")
    }
    
    # Hide sidebar for trajectory tab
    if (!is.null(input$tabs) && identical(input$tabs, "Trajectory")) {
      shinyjs::runjs("document.body.classList.add('trajectory-mode');")
    } else {
      shinyjs::runjs("document.body.classList.remove('trajectory-mode');")
    }
  })

  # Base per-islet dataset (no binning), aligned with current selections, without zero filtering
  raw_df_base <- reactive({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) {
      return(prepared()$targets_all[FALSE, , drop = FALSE])
    }
    pd <- prepared()
    mode <- input$mode %||% "Composition"
    groups <- input$groups
    if (is.null(groups) || !length(groups)) {
      groups <- c("ND", "Aab+", "T1D")
    }
    w <- input$which
    if (is.null(w) || length(w) == 0) w <- NA_character_
    # Resolve defaults before branching so initial renders don't error while inputs initialize
    if (identical(mode, "Composition")) {
      if (!is.character(w) || length(w) != 1 || !nzchar(w) || !(w %in% c("Ins_any","Glu_any","Stt_any"))) {
        w <- "Ins_any"
      }
    } else if (identical(mode, "Targets")) {
      if (!is.character(w) || length(w) != 1 || !nzchar(w)) {
        candidate_classes <- pd$targets_all %>% dplyr::pull(class) %>% unique() %>% na.omit() %>% sort()
        if (length(candidate_classes) > 0) w <- candidate_classes[1]
      }
    } else if (identical(mode, "Markers")) {
      if (!is.character(w) || length(w) != 1 || !nzchar(w)) {
        candidate_markers <- pd$markers_all %>% dplyr::pull(marker) %>% unique() %>% na.omit() %>% sort()
        if (length(candidate_markers) > 0) w <- candidate_markers[1]
      }
    }
    region_val <- input$region %||% "band"
    # Resolve composition default early
    if (identical(mode, "Targets")) {
      req(!is.null(w) && nzchar(w))
      region_tag <- paste0("islet_", tolower(region_val))
      df <- pd$targets_all %>% dplyr::filter(`Donor Status` %in% groups, class == w, tolower(type) == region_tag)
      df <- df %>% dplyr::filter(is.finite(islet_diam_um))
      df <- df %>% dplyr::mutate(value = as.numeric(if (!is.null(input$target_metric) && input$target_metric == "Counts") count else area_density))
    } else if (identical(mode, "Markers")) {
      req(!is.null(w) && nzchar(w))
      region_tag <- paste0("islet_", tolower(region_val))
      df <- pd$markers_all %>% dplyr::filter(`Donor Status` %in% groups, marker == w, tolower(region_type) == region_tag)
      df <- df %>% dplyr::filter(is.finite(islet_diam_um))
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") {
        df <- df %>% dplyr::mutate(value = suppressWarnings(as.numeric(pos_count)))
      } else {
        num <- suppressWarnings(as.numeric(df$pos_count))
        den <- suppressWarnings(as.numeric(df$n_cells))
        df <- df %>% dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
      }
    } else {
      df <- pd$comp %>% dplyr::filter(`Donor Status` %in% groups) %>% dplyr::filter(is.finite(islet_diam_um))
      num <- suppressWarnings(as.numeric(df[[w]]))
      den <- suppressWarnings(as.numeric(df$cells_total))
      df <- df %>% dplyr::mutate(value = ifelse(is.finite(num) & is.finite(den) & den > 0, 100.0 * num / den, NA_real_))
    }
    out <- df %>% dplyr::mutate(donor_status = `Donor Status`) %>% dplyr::filter(!is.na(value))
    # Apply diameter range filter
    if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
      diam_min <- as.numeric(input$diam_range[1])
      diam_max <- as.numeric(input$diam_range[2])
      if (is.finite(diam_min) && is.finite(diam_max)) {
        out <- out %>% dplyr::filter(is.finite(islet_diam_um) & islet_diam_um >= diam_min & islet_diam_um <= diam_max)
      }
    }
    # Apply AAb filters (only within the Aab+ donor group) if selected
    # Reverse logic: EXCLUDE donors with UNCHECKED antibodies
    flags <- input$aab_flags
    all_aab_cols <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
    
    # Keep ND and T1D unchanged; filter Aab+ by EXCLUDING unchecked AAbs
    others <- out %>% dplyr::filter(donor_status != "Aab+")
    aabp   <- out %>% dplyr::filter(donor_status == "Aab+")
    
    if (nrow(aabp) > 0) {
      if (is.null(flags) || length(flags) == 0) {
        # If no flags selected, exclude all Aab+ donors
        aabp <- aabp[0, , drop = FALSE]  # Empty dataframe
      } else {
        # Get unchecked antibodies
        unchecked <- setdiff(all_aab_cols, flags)
        unchecked_avail <- intersect(unchecked, colnames(aabp))
        
        if (length(unchecked_avail) > 0) {
          # Exclude donors who have any of the unchecked antibodies
          exclude_mat <- as.data.frame(aabp[, unchecked_avail, drop = FALSE])
          for (cc in colnames(exclude_mat)) exclude_mat[[cc]] <- as.logical(exclude_mat[[cc]])
          has_unchecked <- rowSums(exclude_mat, na.rm = TRUE) > 0
          aabp <- aabp[!has_unchecked, , drop = FALSE]
        }
      }
    }
    out <- dplyr::bind_rows(others, aabp)
    attr(out, "selection_used") <- w
    attr(out, "mode_used") <- mode
    out
  })

  # Raw per-islet dataset, with top-plot zero filtering applied if selected
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

  # App-wide CSS for theme background
  output$theme_css <- renderUI({
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      tags$style(HTML("body { background-color: #000000; color: #e6e6e6; } .well { background-color: #111111; }"))
    } else {
      tags$style(HTML("body { background-color: #ffffff; color: #111111; }"))
    }
  })

  output$dynamic_selector <- renderUI({
    if (identical(input$mode, "Targets")) {
      # Choices should be independent of Region to preserve selection when switching regions
      classes <- prepared()$targets_all %>% pull(class) %>% unique() %>% na.omit() %>% sort()
      if (length(classes) == 0) classes <- character(0)
      # Preserve current selection if still valid
      sel <- if (!is.null(input$which) && input$which %in% classes) input$which else if (length(classes)>0) classes[1] else NULL
      selectInput("which", "Target class", choices = classes, selected = sel)
    } else if (identical(input$mode, "Markers")) {
      # Choices should be independent of Region to preserve selection when switching regions
      markers <- prepared()$markers_all %>% pull(marker) %>% unique() %>% na.omit() %>% sort()
      if (length(markers) == 0) markers <- character(0)
      sel <- if (!is.null(input$which) && input$which %in% markers) input$which else if (length(markers)>0) markers[1] else NULL
      selectInput("which", "Marker", choices = markers, selected = sel)
    } else {
      selectInput("which", "Composition measure", choices = c("Ins_frac" = "Ins_any", "Glu_frac" = "Glu_any", "Stt_frac" = "Stt_any"), selected = "Ins_any")
    }
  })

  output$aab_warning <- renderUI({
    # Check if current AAb filter results in no Aab+ donors
    df <- raw_df_base()
    if (!is.null(df) && nrow(df) > 0) {
      aabp_donors <- df %>% dplyr::filter(donor_status == "Aab+")
      if (nrow(aabp_donors) == 0 && "Aab+" %in% input$groups) {
        div(style = "color: #d9534f; font-weight: bold; margin-bottom: 10px;", 
            "⚠ No matching donors.")
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  output$region_selector <- renderUI({
    if (identical(input$mode, "Composition")) return(NULL)
    # Display labels while preserving underlying values used in code
    selectInput(
      "region", "Region",
      choices = c("Islet" = "core", "Peri-Islet" = "band", "Islet+20um" = "union"),
      selected = "band"
    )
  })

  output$metric_selector <- renderUI({
    if (identical(input$mode, "Targets")) {
      radioButtons("target_metric", "Targets metric", choices = c("Counts", "Density"), selected = "Density", inline = TRUE)
    } else if (identical(input$mode, "Markers")) {
      radioButtons("marker_metric", "Markers metric", choices = c("Counts", "% positive"), selected = "% positive", inline = TRUE)
    } else {
      NULL
    }
  })

  # color selector removed

  # Dataset for plotting: build from raw_df then bin and normalize
  plot_df <- reactive({
    df <- raw_df()
    if (is.null(df) || !nrow(df)) {
      return(df %||% data.frame())
    }
    mode <- attr(df, "mode_used") %||% (input$mode %||% "Composition")
    selection_label <- attr(df, "selection_used") %||% input$which
    if (is.null(selection_label) || length(selection_label) == 0 || is.na(selection_label[1]) || !nzchar(selection_label[1])) {
      selection_label <- if (identical(mode, "Composition")) "Ins_any" else ""
    }
    selection_label <- selection_label[1]
    bw <- suppressWarnings(as.numeric(input$binwidth))
    if (length(bw) != 1 || !is.finite(bw) || bw <= 0) {
      bw <- 50
    }
    
    # Identify outliers (>3 SD from mean) but keep them in data, just mark them
    df_original <- df
    value_mean <- mean(df$value, na.rm = TRUE)
    value_sd <- sd(df$value, na.rm = TRUE)
    outlier_threshold <- 3
    
    df$is_outlier <- abs(df$value - value_mean) > (outlier_threshold * value_sd)
    outliers <- df[df$is_outlier & !is.na(df$is_outlier), ]
    
    # Store outliers for display in table
    plot_outliers <<- if (nrow(outliers) > 0) {
      data.frame(
        Mode = mode,
        Selection = selection_label,
        Case_ID = outliers$`Case ID`,
        Islet = outliers$islet_key,
        Donor_Status = outliers$donor_status,
        Diameter_um = round(outliers$islet_diam_um, 1),
        Value = round(outliers$value, 3),
        Z_Score = round((outliers$value - value_mean) / value_sd, 2),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
    
    # DO NOT remove outliers - keep them for coloring red in the plot
    cat(sprintf("[PLOT] Identified %d outliers (>3 SD) - will be colored red\n", nrow(outliers)))
    
    # attach bins on the fly
    df <- bin_islet_sizes(df, "islet_diam_um", bw)
    
    # Identify the metric column name to group by for normalization
    metric_col <- if (identical(mode, "Targets")) {
      paste0(input$target_region, "__", input$target_class)
    } else if (identical(mode, "Markers")) {
      paste0(input$marker_region, "__", input$marker_name)
    } else {
      paste0("comp__", input$comp_choice)
    }
    
    # Apply normalization for main plot if requested (default to "none")
    norm_mode <- "none"  # Normalization selector removed
    if (identical(norm_mode, "robust") && all(c("Case ID") %in% colnames(df))) {
      # Group by BOTH donor AND metric to normalize per-metric per-donor
      df <- df %>% dplyr::mutate(metric_id = metric_col) %>%
        dplyr::group_by(`Case ID`, metric_id) %>% dplyr::mutate(
          .med = suppressWarnings(stats::median(value, na.rm = TRUE)),
          .mad = suppressWarnings(stats::mad(value, center = .med, constant = 1, na.rm = TRUE)),
          .r_sd = ifelse(is.finite(.mad) & .mad > 0, .mad * 1.4826, NA_real_),
          value = ifelse(is.finite(.r_sd) & .r_sd > 0, (value - .med) / .r_sd, value)
        ) %>% dplyr::ungroup() %>% dplyr::select(-.med, -.mad, -.r_sd, -metric_id)
    } else if (identical(norm_mode, "global")) {
      mu <- mean(df$value, na.rm = TRUE)
      sdv <- sd(df$value, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      df$value <- (df$value - mu) / sdv
    }
    df
  })

  # Summary with mean and SE per bin per group
  summary_df <- reactive({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) {
      return(data.frame())
    }
    df <- plot_df()
    if (is.null(df) || !nrow(df)) return(data.frame())
    stat_choice <- input$stat %||% "mean_se"
    if (!stat_choice %in% c("mean_se","mean_sd","median_iqr")) {
      stat_choice <- "mean_se"
    }
    summary_stats(df, group_cols = c("donor_status", "diam_bin", "diam_mid"), value_col = "value", stat = stat_choice)
  })

  output$plt <- renderPlotly({
    sm <- summary_df()
    if (is.null(sm) || !nrow(sm)) {
      return(plotly_empty() %>% layout(title = "Loading trajectory data..."))
    }

    # Order groups ND, Aab+, T1D and set colors
    grp_levels <- c("ND","Aab+","T1D")
    sm$donor_status <- factor(sm$donor_status, levels = grp_levels)
    
    # Determine color scheme based on selector
    color_by <- input$plot_color_by %||% "donor_status"
    
    # Build color map based on selection
    if (color_by == "donor_status") {
      color_map <- c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")
    } else {
      # For donor_id, we'll apply colors to raw points only
      color_map <- c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")
    }

    # Build y-label reflecting normalization choice
    ylab_base <- if (identical(input$mode, "Targets")) {
      if (!is.null(input$target_metric) && input$target_metric == "Counts") "Target count" else "Target density (per µm²)"
    } else if (identical(input$mode, "Markers")) {
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") "Positive cell count" else "n positive / total (%)"
    } else {
      "% composition"
    }
  cnorm <- "none"  # Normalization selector removed, always use "none"
    ylab <- switch(cnorm,
                   none = ylab_base,
                   global = paste0(ylab_base, " (scaled z)"),
                   robust = paste0(ylab_base, " (robust z)"),
                   ylab_base)

    # Build title reflecting selection; for Composition, indicate fraction (%)
    title_text <- if (identical(input$mode, "Targets")) {
      paste0(input$which, " vs Islet Size")
    } else if (identical(input$mode, "Markers")) {
      paste0(input$which, " vs Islet Size")
    } else {
      wsel <- input$which
      if (is.null(wsel) || length(wsel) != 1 || !nzchar(wsel)) wsel <- "Ins_any"
      nm <- switch(wsel,
                   Ins_any = "Insulin+ fraction",
                   Glu_any = "Glucagon+ fraction",
                   Stt_any = "Somatostatin+ fraction",
                   wsel)
      paste0(nm, " vs Islet Size")
    }

    # Create base plot with summary lines/points
    # When coloring by donor_id, keep summary lines by donor_status
    p <- ggplot(sm, aes(x = diam_mid, y = y, color = donor_status, group = donor_status)) +
      geom_line(alpha = ifelse(!is.null(input$add_smooth) && input$add_smooth == "LOESS", 0, 1)) +
      geom_point() +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) +
   labs(x = "Islet diameter (µm)",
     y = ylab,
     color = NULL,
     title = title_text) +
      scale_color_manual(values = color_map, breaks = grp_levels, drop = FALSE) +
      theme_minimal(base_size = 14) +
      theme(legend.position = c(1, 1), 
            legend.justification = c(1, 1),
            legend.direction = "vertical")

    # Apply axis scales BEFORE adding points so transformations are applied correctly
    # Determine x-axis max for breaks
    xmax <- suppressWarnings(max(sm$diam_mid, na.rm = TRUE))
    if (!is.finite(xmax)) xmax <- 300
    
    # Apply log2 scale to x-axis if selected
    if (!is.null(input$log_scale_x) && input$log_scale_x) {
      # For log2 scale, use powers of 2 as breaks
      max_pow <- ceiling(log2(xmax))
      major_breaks <- 2^(0:max_pow)
      p <- p + scale_x_continuous(trans = "log2", breaks = major_breaks)
    } else {
      # Linear scale with standard breaks
      major_breaks <- seq(0, ceiling(xmax/50)*50, by = 50)
      minor_breaks <- seq(0, ceiling(xmax/10)*10, by = 10)
      p <- p + scale_x_continuous(breaks = major_breaks, minor_breaks = minor_breaks)
    }
    
    # Apply log scale to y-axis if selected
    if (!is.null(input$log_scale) && input$log_scale) {
      p <- p + scale_y_log10()
    }

    if (isTRUE(input$show_points)) {
      raw <- plot_df()
      raw$donor_status <- factor(raw$donor_status, levels = grp_levels)

      # Add islet_click_key for segmentation viewer click handling
      # Combines Case ID and islet_key with a separator for parsing in click handler
      if ("Case ID" %in% names(raw) && "islet_key" %in% names(raw)) {
        raw$islet_click_key <- paste(raw$`Case ID`, raw$islet_key, sep = "|")
      } else {
        raw$islet_click_key <- NA_character_
      }

      # Separate normal points from outliers
      raw_normal <- raw[!raw$is_outlier | is.na(raw$is_outlier), ]
      raw_outliers <- raw[raw$is_outlier & !is.na(raw$is_outlier), ]
      
      # Adjust jitter width based on whether log2 scale is enabled
      # In log space, jitter should be multiplicative rather than additive
      if (!is.null(input$log_scale_x) && input$log_scale_x) {
        # For log2 scale, use a fraction of the value (multiplicative jitter)
        # This creates proportional jitter that looks consistent across the log scale
        jw <- 0.15  # 15% multiplicative jitter in log space
      } else {
        # For linear scale, use standard additive jitter
        jw <- max(1, as.numeric(input$binwidth) * 0.35)
      }
      
      # Add slight vertical jitter to avoid rows when plotting counts
      yr <- suppressWarnings(range(raw$value, na.rm = TRUE))
      ydiff <- if (all(is.finite(yr))) diff(yr) else 0
      is_counts <- (identical(input$mode, "Targets") && !is.null(input$target_metric) && input$target_metric == "Counts") ||
                   (identical(input$mode, "Markers") && !is.null(input$marker_metric) && input$marker_metric == "Counts")
      jh <- if (is_counts) {
        # 2% of range or at least 0.5 in data units
        mx <- max(0.02 * ydiff, 0.5)
        if (!is.finite(mx) || mx <= 0) 0.5 else mx
      } else {
        # 1% of range
        mx <- 0.01 * ydiff
        if (!is.finite(mx) || mx < 0) 0 else mx
      }
      
      # Apply color scheme for NORMAL points
      if (color_by == "donor_id") {
        # Color by donor ID
        if ("Case ID" %in% colnames(raw_normal)) {
          raw_normal$donor_id <- factor(raw_normal$`Case ID`)
          donor_colors <- get_donor_palette(levels(raw_normal$donor_id))
          # Combine donor_status colors (for lines) with donor_id colors (for points)
          combined_colors <- c(
            "ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd",
            donor_colors
          )
          p <- p +
            geom_point(data = raw_normal, aes(x = diam_mid, y = value, color = donor_id, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh), size = input$pt_size, alpha = input$pt_alpha, inherit.aes = FALSE) +
            scale_color_manual(values = combined_colors,
                             breaks = levels(raw_normal$donor_id),
                             name = "Donor ID")
        } else {
          # Fallback to donor_status if Case ID not available
          p <- p +
            geom_point(data = raw_normal, aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                       position = position_jitter(width = jw, height = jh), size = input$pt_size, alpha = input$pt_alpha, inherit.aes = FALSE)
        }
      } else {
        # Color by donor_status (default)
        p <- p +
          geom_point(data = raw_normal, aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                     position = position_jitter(width = jw, height = jh), size = input$pt_size, alpha = input$pt_alpha, inherit.aes = FALSE)
      }
      
      # Add OUTLIERS on top (colored by donor status, slightly larger/brighter)
      if (nrow(raw_outliers) > 0) {
        p <- p +
          geom_point(data = raw_outliers, aes(x = diam_mid, y = value, color = donor_status, key = islet_click_key),
                     position = position_jitter(width = jw, height = jh),
                     size = input$pt_size * 1.2, # Slightly larger
                     alpha = min(1.0, input$pt_alpha * 1.5), # More visible
                     inherit.aes = FALSE)
      }
    }

    # Add theme for grid lines
    p <- p + theme(panel.grid.minor = element_line(size = 0.2, colour = if (!is.null(input$theme_bg) && input$theme_bg == "Dark") "#222222" else "#eeeeee"))

    # Highlight significant bins if stats were run
    st <- NULL
    if (!is.null(st) && nrow(st) > 0) {
      alpha_num <- as.numeric(input$alpha)
      st <- st %>% mutate(p_adj = p.adjust(p_value, method = "BH"), sig = p_adj <= alpha_num)
      sig <- st %>% filter(sig)
      if (nrow(sig) > 0) {
        ytop <- suppressWarnings(max(sm$ymax, na.rm = TRUE))
        ytop <- ifelse(is.finite(ytop), ytop, suppressWarnings(max(sm$y, na.rm = TRUE)))
        sig$y <- ytop * 1.05
        p <- p + geom_text(data = sig, aes(x = mid, y = y), label = ifelse(unique(st$test) == 'kendall', '816', '*'), inherit.aes = FALSE, color = "#444444")
      }
    }

    # Optional smoothing overlay (kept subtle)
    if (!is.null(input$add_smooth) && input$add_smooth == "LOESS") {
      p <- p + geom_smooth(se = FALSE, method = "loess", span = 0.6, size = 0.9)
    }

    # Apply dark theme inside the plot if selected
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      p <- p + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.text = element_text(color = "#e6e6e6"),
        legend.title = element_text(color = "#f0f0f0")
      )
    }

    gg <- ggplotly(p, source = "plot_scatter")
    gg <- gg %>% layout(legend = list(orientation = "h", x = 0, y = -0.15))
    gg
  })

  # Stats: continuous tests (no binning). Global LM + post-hoc pairwise on residuals
  stats_run <- eventReactive(input$run_tests, {
    rdf <- raw_df()
    if (is.null(rdf) || !nrow(rdf)) return(NULL)
    
    # Remove outliers if requested (same as plot_df)
    if (isTRUE(input$stats_remove_outliers)) {
      value_mean <- mean(rdf$value, na.rm = TRUE)
      value_sd <- sd(rdf$value, na.rm = TRUE)
      outlier_threshold <- 3
      
      rdf$is_outlier <- abs(rdf$value - value_mean) > (outlier_threshold * value_sd)
      n_outliers <- sum(rdf$is_outlier, na.rm = TRUE)
      
      # Remove outliers
      rdf <- rdf[!rdf$is_outlier | is.na(rdf$is_outlier), ]
      
      cat(sprintf("[STATISTICS] Removed %d outliers (>3 SD)\n", n_outliers))
    }
    
    # Ensure donor_status factor ordering
    rdf$donor_status <- factor(rdf$donor_status, levels = c("ND","Aab+","T1D"))
    # Global: value ~ donor_status + islet_diam_um (Type I with donor_status first)

    fit <- tryCatch(lm(value ~ donor_status + islet_diam_um, data = rdf), error = function(e) NULL)
    p_global <- NA_real_
    if (!is.null(fit)) {
      at <- tryCatch(anova(fit), error = function(e) NULL)
      if (!is.null(at)) p_global <- suppressWarnings(as.numeric(at[["Pr(>F)"]][1]))
    }
    # Post-hoc: residualize value ~ islet_diam_um, then pairwise t-tests across groups
    res <- tryCatch({
      fit_res <- lm(value ~ islet_diam_um, data = rdf)
      r <- resid(fit_res)
      pt <- pairwise.t.test(r, rdf$donor_status, p.adjust.method = "BH")
      mat <- as.data.frame(as.table(pt$p.value))
      colnames(mat) <- c("group1","group2","p_value")
      mat <- mat[!is.na(mat$p_value), , drop = FALSE]
      mat
    }, error = function(e) NULL)
    # Assemble results
    global <- data.frame(type = "global", contrast = "donor_status", p_value = p_global, stringsAsFactors = FALSE)
    if (!is.null(res) && nrow(res)) {
      pairs <- data.frame(type = "pairwise",
                          contrast = paste(res$group1, "vs", res$group2),
                          p_value = res$p_value,
                          stringsAsFactors = FALSE)
      out <- rbind(global, pairs)
    } else {
      out <- global
    }
    out
  }, ignoreInit = TRUE)

  stats_data <- reactive({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_run()
    if (is.null(st) || nrow(st) == 0) return(NULL)
    alpha_num <- as.numeric(input$alpha)
    st$p_adj <- NA_real_
    if (any(st$type == "pairwise")) {
      idx <- which(st$type == "pairwise")
      st$p_adj[idx] <- p.adjust(st$p_value[idx], method = "BH")
    }
    st$sig <- ifelse(!is.na(st$p_adj), st$p_adj <= alpha_num, st$p_value <= alpha_num)
    st
  })

  # Statistical test explanation
  output$stats_explanation <- renderUI({
    st <- stats_data()
    if (is.null(st) || nrow(st) == 0) {
      return(tags$div(
        style = "padding: 15px; color: #666;",
        tags$p("Click 'Run Statistics' to perform statistical analysis.")
      ))
    }
    
    alpha_num <- as.numeric(input$alpha)
    n_total <- nrow(st)
    n_sig <- sum(st$sig, na.rm = TRUE)
    
    tags$div(
      style = "padding: 15px;",
      tags$h5("Analysis Method:", style = "margin-top: 0;"),
      tags$p(
        style = "margin-bottom: 10px;",
        tags$strong("Global ANOVA:"), " Tests if there are significant differences in the selected metric across the three donor status groups (ND, Aab+, T1D) while controlling for islet diameter.",
        tags$br(),
        tags$em("Model: value ~ donor_status + islet_diam_um")
      ),
      tags$p(
        style = "margin-bottom: 10px;",
        tags$strong("Pairwise t-tests:"), " Compares each pair of donor status groups (ND vs Aab+, ND vs T1D, Aab+ vs T1D) on residuals from the diameter-adjusted model.",
        tags$br(),
        tags$em("Multiple testing correction: Benjamini-Hochberg (FDR) method")
      ),
      tags$hr(),
      tags$h5("Results Summary:"),
      tags$p(
        sprintf("Significance threshold: α = %s", input$alpha),
        tags$br(),
        sprintf("Significant tests: %d / %d (%.1f%%)", n_sig, n_total, 100 * n_sig / n_total)
      ),
      tags$p(
        style = "font-size: 90%; color: #666;",
        tags$strong("Interpretation:"),
        " A significant result (p < α) suggests the selected metric differs between donor status groups after accounting for islet size. ",
        "Non-significant results may indicate no biological difference or insufficient statistical power."
      )
    )
  })

  output$stats_tbl <- renderTable({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_data()
    if (is.null(st) || nrow(st) == 0) return(NULL)
    fmt <- function(x) ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 2))
    display <- st
    display$p_value <- fmt(display$p_value)
    display$p_adj <- fmt(display$p_adj)
    display
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$stats_plot <- renderPlotly({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_data()
    req(st, nrow(st) > 0)
    # Force dependency on alpha to trigger re-render
    alpha_num <- as.numeric(input$alpha)
    req(alpha_num)
    
    df <- st %>%
      mutate(neglog = ifelse(!is.na(p_value) & p_value > 0, -log10(p_value), NA_real_)) %>%
      filter(!is.na(neglog))
    if (nrow(df) == 0) return(NULL)
    thresh <- -log10(alpha_num)
    g <- ggplot(df, aes(x = reorder(contrast, neglog), y = neglog, fill = type)) +
      geom_col(show.legend = TRUE) +
      geom_hline(yintercept = thresh, linetype = "dashed", color = "#444444") +
      coord_flip() +
      labs(x = "Contrast", y = expression(-log[10](p)), fill = "Test type",
           title = "Test significance overview") +
      theme_minimal(base_size = 14)
    ggplotly(g)
  })

  output$pairwise_plot <- renderPlotly({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    st <- stats_data()
    req(st, nrow(st) > 0)
    # Force dependency on alpha to trigger re-render
    alpha_num <- as.numeric(input$alpha)
    req(alpha_num)
    
    df <- st %>% filter(type == "pairwise" & !is.na(p_value))
    if (nrow(df) == 0) return(NULL)
    df <- df %>% mutate(sig = ifelse(is.na(p_adj), p_value <= alpha_num, p_adj <= alpha_num))
    g <- ggplot(df, aes(x = reorder(contrast, p_value), y = p_value, color = sig)) +
      geom_point(size = 3) +
      geom_hline(yintercept = alpha_num, linetype = "dashed", color = "#444444") +
      scale_y_log10() +
      scale_color_manual(values = c(`FALSE` = "#1f77b4", `TRUE` = "#d62728"),
                         labels = c(`FALSE` = "NS", `TRUE` = "Significant")) +
      coord_flip() +
      labs(x = "Pairwise contrast", y = "p-value (log scale)", color = "", title = "Pairwise comparison p-values") +
      theme_minimal(base_size = 14)
    ggplotly(g)
  })

  # Area Under Curve by Donor Group
  output$auc_plot <- renderPlotly({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    
    # Get the summary data (binned by diameter)
    sdf <- summary_df()
    if (is.null(sdf) || nrow(sdf) == 0) return(NULL)
    
    # Calculate area under each line using trapezoidal rule
    # summary_df returns columns: donor_status, diam_bin, diam_mid, y, ymin, ymax
    auc_by_group <- tryCatch({
      sdf %>%
        dplyr::group_by(donor_status) %>%
        dplyr::arrange(diam_mid) %>%
        dplyr::summarise(
          auc = if (n() < 2) 0 else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
          .groups = "drop"
        )
    }, error = function(e) {
      cat("[AUC ERROR]", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(auc_by_group) || nrow(auc_by_group) == 0) return(NULL)
    
    # Create bar plot
    g <- ggplot(auc_by_group, aes(x = donor_status, y = auc, fill = donor_status)) +
      geom_col(width = 0.7, color = "black", alpha = 0.8) +
      scale_fill_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd")) +
      labs(
        x = "Donor Status",
        y = "Area Under Curve (AUC)",
        title = "Integrated Area by Donor Group"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
    
    ggplotly(g)
  })
  
  # AUC summary table with statistical comparisons
  output$auc_table <- renderTable({
    if (!is.null(input$tabs) && identical(input$tabs, "Viewer")) return(NULL)
    
    sdf <- summary_df()
    if (is.null(sdf) || nrow(sdf) == 0) return(NULL)
    
    # Calculate AUC for each donor group
    # summary_df returns columns: donor_status, diam_bin, diam_mid, y, ymin, ymax
    auc_by_group <- tryCatch({
      result <- sdf %>%
        dplyr::group_by(donor_status) %>%
        dplyr::arrange(diam_mid) %>%
        dplyr::summarise(
          auc = if (n() < 2) 0 else sum(diff(diam_mid) * (head(y, -1) + tail(y, -1))) / 2,
          n_bins = n(),
          .groups = "drop"
        )
      # Only return result if it has valid data
      if (nrow(result) > 0) result else NULL
    }, error = function(e) {
      cat("[AUC TABLE ERROR]", conditionMessage(e), "\n")
      NULL
    })
    
    # Only show error if truly no data
    if (is.null(auc_by_group)) {
      return(data.frame(Message = "Insufficient data for AUC calculation"))
    }
    
    # Calculate pairwise percent differences
    nd_auc <- auc_by_group$auc[auc_by_group$donor_status == "ND"]
    aab_auc <- auc_by_group$auc[auc_by_group$donor_status == "Aab+"]
    t1d_auc <- auc_by_group$auc[auc_by_group$donor_status == "T1D"]
    
    result <- data.frame(
      "Donor Group" = auc_by_group$donor_status,
      "AUC" = sprintf("%.2f", auc_by_group$auc),
      "N Bins" = auc_by_group$n_bins,
      check.names = FALSE
    )
    
    # Add comparison row
    if (length(nd_auc) > 0 && length(t1d_auc) > 0) {
      pct_change_nd_t1d <- ((t1d_auc - nd_auc) / nd_auc) * 100
      result <- rbind(result, data.frame(
        "Donor Group" = "T1D vs ND",
        "AUC" = sprintf("%.1f%%", pct_change_nd_t1d),
        "N Bins" = "",
        check.names = FALSE
      ))
    }
    
    result
  }, striped = TRUE, bordered = TRUE)

  # Vitessce embed
  # Image picker for local OME-TIFFs under www/local_images (preferred) or LOCAL_IMAGE_ROOT fallback
  output$local_image_picker <- renderUI({
    # Only render when on Viewer tab
    if (is.null(input$tabs) || input$tabs != "Viewer") {
      return(NULL)
    }
    
    vi <- viewer_info()
    cat("[VIEWER RENDER] Base:", vi$base, "Image URL (app):", vi$image_url, 
        "Public:", vi$image_public_url %||% "NULL", "\n")
    cat("[VIEWER RENDER] Full iframe src:", vi$iframe_src, "\n")
    
    if (is.null(vi$base)) {
      return(tagList(
        tags$div(style = "color:#b00;", "Local Avivator static build not found under shiny/www/avivator."),
        tags$div(style = "color:#666; font-size:90%;", "To install: run scripts/install_avivator.sh (requires Node ≥ 18) or place a prebuilt bundle under shiny/www/avivator.")
      ))
    }
    
    # Build list of available images
    available_images <- character(0)
    images_dir <- file.path("www", "local_images")
    if (dir.exists(images_dir)) {
      image_files <- list.files(images_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
      available_images <- c(available_images, image_files)
    }
    
    # Also check LOCAL_IMAGE_ROOT if set
    local_root <- Sys.getenv("LOCAL_IMAGE_ROOT", unset = "")
    if (nzchar(local_root) && dir.exists(local_root)) {
      root_files <- list.files(local_root, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
      available_images <- unique(c(available_images, root_files))
    }
    
    # Check if image is accessible and not too large for web serving
    image_basename <- NULL
    if (!is.null(vi$image_rel)) {
      image_basename <- basename(vi$image_rel)
    } else if (!is.null(vi$image_url)) {
      image_basename <- basename(vi$image_url)
    } else if (!is.null(vi$image_public_url)) {
      image_basename <- basename(vi$image_public_url)
    }

    if (!is.null(image_basename) && nzchar(image_basename)) {
      local_image_path <- file.path("www", "local_images", image_basename)
      file_exists <- file.exists(local_image_path)
      
      cat("[VIEWER DEBUG] Looking for image file:", local_image_path, "Exists:", file_exists, "\n")
      
      if (file_exists) {
        # Log file size for debugging
        file_size <- file.info(local_image_path)$size
        file_size_gb <- round(file_size / (1024^3), 2)
        cat("[VIEWER DEBUG] Image size:", file_size_gb, "GB\n")
      }
    }
    
    # Determine current selection
    current_selection <- if (!is.null(forced_image())) {
      forced_image()
    } else if (!is.null(input$selected_image)) {
      input$selected_image
    } else if (length(available_images) > 0) {
      available_images[1]
    } else {
      NULL
    }
    
    # Build UI with dropdown if images are available
    tagList(
      if (length(available_images) > 0) {
        fluidRow(
          column(12,
            selectInput(
              "selected_image",
              "Select Image:",
              choices = available_images,
              selected = current_selection,
              width = "100%"
            )
          )
        )
      } else {
        tags$div(
          style = "padding: 10px; background-color: #fff3cd; border: 1px solid #ffc107; margin-bottom: 10px;",
          tags$strong("No images found."),
          " Place OME-TIFF files in ",
          tags$code("app/shiny_app/www/local_images/"),
          " or set ",
          tags$code("LOCAL_IMAGE_ROOT"),
          " environment variable."
        )
      },
      if (VIEWER_DEBUG_ENABLED) {
        tags$details(
          style = "margin: 10px 0;",
          tags$summary("Viewer debug info"),
          tags$pre(
            style = "max-height:200px; overflow:auto;",
            jsonlite::toJSON(list(
              env = vi$env,
              image_rel = vi$image_rel,
              image_app = vi$image_app_url %||% vi$image_url,
              image_public = vi$image_public_url,
              iframe_src = vi$iframe_src,
              asset_diag = vi$asset_diag
            ), auto_unbox = TRUE, pretty = TRUE)
          )
        )
      },
      # Only render iframe if we have an image selected
      if (!is.null(vi$image_url) && nzchar(vi$image_url)) {
        tags$iframe(
          id = "avivator_viewer",
          src = vi$iframe_src,
          width = "100%",
          frameBorder = 0,
          allowfullscreen = NA,
          style = "width:100%; height:calc(100vh - 180px); border:0;"
        )
      } else {
        tags$div(
          style = "padding: 20px; text-align: center; color: #666;",
          "Select an image from the dropdown above to view it in AVIVATOR."
        )
      }
    )
  })

  
  # Pseudotime scatter when counts are selected (scaled counts vs pseudotime)
  # Distribution plot UI (violin/box)
  output$dist_ui <- renderUI({
    tagList(
      fluidRow(
        column(4,
          sliderInput("dist_pt_size", "Point size", min = 0.3, max = 4.0, value = 0.7, step = 0.1),
          sliderInput("dist_pt_alpha", "Point transparency", min = 0.05, max = 1.0, value = 0.25, step = 0.05)
        ),
        column(4,
          selectInput("dist_color_by", "Color points by:",
                     choices = c("Donor Status" = "donor_status", 
                               "Donor ID" = "donor_id"),
                     selected = "donor_status"),
          radioButtons("dist_type", "Plot type", choices = c("Violin", "Box"), selected = "Violin", inline = TRUE)
        ),
        column(4,
          checkboxInput("exclude_zero_dist", "Exclude zero values", value = FALSE),
          checkboxInput("dist_show_points", "Show individual points", value = TRUE),
          checkboxInput("dist_log_scale", "Log scale y-axis", value = FALSE)
        )
      )
    )
  })

  # Distribution plot: violin/box per donor group for current raw_df()
  output$dist <- renderPlotly({
    # Use the SAME data source as the main plot (respects Statistic selector)
    rdf <- raw_df()
    if (is.null(rdf) || nrow(rdf) == 0) return(NULL)
    
    # Explicit dependencies to ensure reactivity
    req(input$mode, input$which)

    # Force reactivity on distribution-specific controls
    dist_color_by_val <- input$dist_color_by
    dist_show_points_val <- input$dist_show_points

    cat(sprintf("[DIST] Rendering with mode=%s, which=%s, region=%s, color_by=%s\n",
                input$mode %||% "NULL",
                input$which %||% "NULL",
                input$region %||% "NULL",
                dist_color_by_val %||% "NULL"))
    
    # Apply exclude_zero_dist separately (since raw_df uses exclude_zero_top)
    if (isTRUE(input$exclude_zero_dist) && !isTRUE(input$exclude_zero_top)) {
      rdf <- rdf %>% dplyr::filter(value != 0)
    }
    
    # Identify the metric column name to group by for normalization
    metric_col <- if (identical(input$mode, "Targets")) {
      paste0(input$region, "__", input$which)
    } else if (identical(input$mode, "Markers")) {
      paste0(input$region, "__", input$which)
    } else {
      paste0("comp__", input$which)
    }
    
    # Apply normalization consistent with main plot (default to "none")
    norm_mode <- "none"  # Normalization selector removed
    if (identical(norm_mode, "robust") && all(c("Case ID") %in% colnames(rdf))) {
      # Group by BOTH donor AND metric to normalize per-metric per-donor
      rdf <- rdf %>% dplyr::mutate(metric_id = metric_col) %>%
        dplyr::group_by(`Case ID`, metric_id) %>% dplyr::mutate(
          .med = suppressWarnings(stats::median(value, na.rm = TRUE)),
          .mad = suppressWarnings(stats::mad(value, center = .med, constant = 1, na.rm = TRUE)),
          .r_sd = ifelse(is.finite(.mad) & .mad > 0, .mad * 1.4826, NA_real_),
          value = ifelse(is.finite(.r_sd) & .r_sd > 0, (value - .med) / .r_sd, value)
        ) %>% dplyr::ungroup() %>% dplyr::select(-.med, -.mad, -.r_sd, -metric_id)
    } else if (identical(norm_mode, "global")) {
      mu <- mean(rdf$value, na.rm = TRUE)
      sdv <- sd(rdf$value, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      rdf$value <- (rdf$value - mu) / sdv
    }
    
    # Determine color scheme based on selector
    color_by <- input$dist_color_by %||% "donor_status"
    
    # Ensure factor order for donor_status
    rdf$donor_status <- factor(rdf$donor_status, levels = c("ND","Aab+","T1D"))
    
    # Create donor_id column BEFORE using it in aes() to avoid errors
    if (color_by == "donor_id" && "Case ID" %in% colnames(rdf)) {
      rdf$donor_id <- factor(rdf$`Case ID`)
      cat(sprintf("[DIST] Created donor_id with %d levels: %s\n", 
                  length(levels(rdf$donor_id)), 
                  paste(levels(rdf$donor_id), collapse=", ")))
    }
    
    ylab_base <- if (identical(input$mode, "Targets")) {
      if (!is.null(input$target_metric) && input$target_metric == "Counts") "Target count" else "Target density (per µm²)"
    } else if (identical(input$mode, "Markers")) {
      if (!is.null(input$marker_metric) && input$marker_metric == "Counts") "Positive cell count" else "n positive / total (%)"
    } else {
      "% composition"
    }
    ylab <- ylab_base  # Normalization selector removed, always use base label
    
    # Standard ggplot-based violin/box with jitter; keep styling
    g <- ggplot(rdf, aes(x = donor_status, y = value, fill = donor_status))
    
    if (!is.null(input$dist_type) && input$dist_type == "Box") {
      g <- g + geom_boxplot(
        outlier.shape = NA,
        outlier.size = 0,
        outlier.colour = NA,
        outlier.fill = NA,
        outlier.alpha = 0,
        outlier.stroke = 0,
        alpha = 0.8
      )
    } else {
      g <- g + geom_violin(trim = FALSE, alpha = 0.8)
    }
    
    g <- g + scale_fill_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd"), guide = "none")
    
    if (isTRUE(input$dist_show_points)) {
      if (color_by == "donor_id" && "donor_id" %in% colnames(rdf)) {
        donor_colors <- get_donor_palette(levels(rdf$donor_id))
        cat(sprintf("[DIST] Applying donor_id colors for %d donors\n", length(donor_colors)))
        g <- g + geom_jitter(aes(x = donor_status, y = value, color = donor_id), width = 0.18,
                             alpha = ifelse(is.null(input$dist_pt_alpha), 0.5, input$dist_pt_alpha),
                             size = ifelse(is.null(input$dist_pt_size), 1.5, input$dist_pt_size),
                             stroke = 0, inherit.aes = FALSE) +
          scale_color_manual(values = donor_colors,
                             breaks = levels(rdf$donor_id),
                             name = "Donor ID",
                             guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
      } else {
        g <- g + geom_jitter(aes(x = donor_status, y = value, color = donor_status), width = 0.18,
                             alpha = ifelse(is.null(input$dist_pt_alpha), 0.25, input$dist_pt_alpha),
                             size = ifelse(is.null(input$dist_pt_size), 0.7, input$dist_pt_size),
                             stroke = 0, inherit.aes = FALSE) +
          scale_color_manual(values = c("ND" = "#2ca02c", "Aab+" = "#ffcc00", "T1D" = "#9467bd"), guide = "none")
      }
    }
    
    g <- g + labs(x = "Donor Status", y = ylab, title = "Distribution across donor groups") +
      theme_minimal(base_size = 14) +
      theme(legend.position = if (color_by == "donor_id") "right" else "none")
    
    if (!is.null(input$dist_log_scale) && input$dist_log_scale) {
      g <- g + scale_y_log10()
    }
    
    if (!is.null(input$theme_bg) && input$theme_bg == "Dark") {
      g <- g + theme(
        plot.background = element_rect(fill = "#000000", colour = NA),
        panel.background = element_rect(fill = "#000000", colour = NA),
        panel.grid.major = element_line(color = "#333333"),
        axis.text = element_text(color = "#e6e6e6"),
        axis.title = element_text(color = "#f0f0f0"),
        plot.title = element_text(color = "#f0f0f0"),
        legend.position = if (color_by == "donor_id") "right" else "none"
      )
    }
    
    p <- ggplotly(g, tooltip = c("x", "y", "colour"))
    
    if (!is.null(input$dist_type) && input$dist_type == "Box") {
      for (i in seq_along(p$x$data)) {
        if (!is.null(p$x$data[[i]]$type) && p$x$data[[i]]$type == "box") {
          p$x$data[[i]]$marker$opacity <- 0
          p$x$data[[i]]$marker$size <- 0
          p$x$data[[i]]$marker$color <- "rgba(0,0,0,0)"
          p$x$data[[i]]$marker$outliercolor <- "rgba(0,0,0,0)"
          p$x$data[[i]]$boxpoints <- FALSE
        }
      }
    }
    
    p <- p %>% layout(showlegend = isTRUE(color_by == "donor_id"))
    p
  })

  # Helper function to create selection descriptions for CSV export
  get_selection_description <- function() {
    desc <- c()
    desc <- c(desc, paste("Generated:", Sys.time()))
    desc <- c(desc, paste("Focus:", input$mode %||% "Unknown"))
    desc <- c(desc, paste("Region:", input$region %||% "Unknown"))
    
    if (!is.null(input$dynamic_metric)) {
      desc <- c(desc, paste("Metric:", input$dynamic_metric))
    }
    
    desc <- c(desc, paste("Donor Groups:", paste(input$groups %||% c(), collapse = ", ")))
    
    if (length(input$aab_flags) > 0) {
      # Reverse logic: showing which antibodies are EXCLUDED (unchecked)
      all_aab <- c("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA")
      excluded <- setdiff(all_aab, input$aab_flags)
      if (length(excluded) > 0) {
        aab_desc <- paste("AAb Filter: Excluding donors with", paste(excluded, collapse = ", "))
        desc <- c(desc, aab_desc)
      } else {
        desc <- c(desc, "AAb Filter: All Aab+ donors included")
      }
    }
    
    desc <- c(desc, paste("Statistic:", input$stat %||% "mean_se"))
    desc <- c(desc, paste("Bin Width:", input$binwidth %||% "50", "µm"))
    if (!is.null(input$diam_range) && length(input$diam_range) == 2) {
      desc <- c(desc, paste("Diameter Range:", input$diam_range[1], "-", input$diam_range[2], "µm"))
    }
    
    return(desc)
  }

  # pseudotime stats removed

  output$dl_summary <- downloadHandler(
    filename = function() {
      paste0("summary_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      df <- summary_df()
      if (is.null(df) || nrow(df) == 0) df <- data.frame()
      
      # Add selection descriptions as comments at the top
      descriptions <- get_selection_description()
      
      # Write descriptions as comments first
      writeLines(paste("#", descriptions), file)
      writeLines("", file, sep = "\n")
      
      # Append the actual data
      write.table(df, file, append = TRUE, sep = ",", row.names = FALSE)
    }
  )

  output$dl_stats <- downloadHandler(
    filename = function() {
      paste0("stats_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      st <- stats_data()
      if (is.null(st) || nrow(st) == 0) st <- data.frame()
      
      # Add selection descriptions as comments at the top
      descriptions <- get_selection_description()
      
      # Write descriptions as comments first
      writeLines(paste("#", descriptions), file)
      writeLines("", file, sep = "\n")
      
      # Append the actual data
      write.table(st, file, append = TRUE, sep = ",", row.names = FALSE)
    }
  )
  
  # Download handler for the new Statistics tab download button
  output$download_stats <- downloadHandler(
    filename = function() {
      paste0("statistics_", gsub("[^0-9A-Za-z]+","_", Sys.time()), ".csv")
    },
    content = function(file) {
      st <- stats_data()
      if (is.null(st) || nrow(st) == 0) st <- data.frame()
      
      # Add selection descriptions as comments at the top
      descriptions <- get_selection_description()
      
      # Write descriptions as comments first
      writeLines(paste("#", descriptions), file)
      writeLines("", file, sep = "\n")
      
      # Append the actual data
      write.table(st, file, append = TRUE, sep = ",", row.names = FALSE)
    }
  )
  
  # Outlier table display for Plot tab - hidden by default, shown only when checkbox is checked
  output$plot_outlier_info <- renderUI({
    show_table <- input$show_plot_outlier_table %||% FALSE
    
    if (show_table && exists("plot_outliers") && !is.null(plot_outliers) && nrow(plot_outliers) > 0) {
      tags$div(
        style = "margin-top: 15px; padding: 10px; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 5px;",
        tags$h6(
          style = "margin-top: 0; color: #856404;",
          sprintf("⚠️ %d Outlier%s Removed (>3 SD)", nrow(plot_outliers), ifelse(nrow(plot_outliers) > 1, "s", ""))),
        tags$div(
          style = "max-height: 200px; overflow-y: auto;",
          renderTable(
            plot_outliers
          )
        )
      )
    } else {
      NULL
    }
  })
}

shinyApp(ui = ui, server = server)
