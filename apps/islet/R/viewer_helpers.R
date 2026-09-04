# viewer_helpers.R
# Avivator viewer helper functions: channel loading, local-image setup,
# environment detection, URL construction, and asset probing.
#
# Globals referenced (must be defined before sourcing):
#   project_root            – top-level project directory (or NULL)
#   VIEWER_DEBUG_ENABLED    – logical flag for extra logging
#
# Globals *created* by this file:
#   channel_names_vec       – character vector of channel names (or NULL)
#   local_images_root       – path to local_images directory (or NULL)
#   has_www_local_images    – logical: does www/local_images/ exist?
#   default_image_url       – NULL; set later by server logic

# ---- Channel names -------------------------------------------------------

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

    return(labels)
  }

  message("[CHANNELS] Channel names not found. Embed names with: tiffcomment -set \"ChannelNames=<comma-separated>\" <file.ome.tif>")
  NULL
}

channel_names_vec <- load_channel_names()

# ---- Local images setup ---------------------------------------------------

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
  }, silent = TRUE)
}

# Register resource path for viewer URLs
if (has_www_local_images) {
  try({
    # Register local_images path for both local and remote compatibility
    # Use 'img_data' to avoid conflict warning with www/local_images subdirectory
    shiny::addResourcePath("img_data", file.path("www", "local_images"))
  }, silent = TRUE)
}

label_exports_dir <- file.path("www", "LabelExports")
if (dir.exists(label_exports_dir)) {
  try({
    shiny::addResourcePath("label_exports", label_exports_dir)
  }, silent = TRUE)
}

# ---- Viewer helpers/defaults (Avivator) -----------------------------------

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

# Extract OME-XML PhysicalSizeX/Y (um) from TIFF ImageDescription and compute um per pixel.
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
    unit_x <- xml2::xml_attr(px, "PhysicalSizeXUnit") %||% "\u00b5m"
    unit_y <- xml2::xml_attr(px, "PhysicalSizeYUnit") %||% "\u00b5m"
    # Only support um units here
    if (is.na(phys_x) || is.na(phys_y)) return(NULL)
    if (!identical(unit_x, "\u00b5m") || !identical(unit_y, "\u00b5m")) return(NULL)
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
      return(matches[1])
    }
  }

  # Fallback to LOCAL_IMAGE_ROOT
  if (!is.null(local_root) && dir.exists(local_root)) {
    matches <- list.files(local_root, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    }
  }

  return(NULL)
}

# ---- Viewer backend: Avivator (default) or TissUUmaps ---------------------
# TissUUmaps is a separate Flask service (see deploy/islet-tissuumaps.container)
# that tiles pre-generated per-channel pyramids and adds region drawing, GeoJSON
# import/export and per-cell marker overlays. Its projects are `.tmap` files
# built by scripts/senior/build_tissuumaps_project.py — one per donor, with the
# channel layers named after their markers.
#
# Set ISLET_VIEWER_BACKEND=tissuumaps to switch the Viewer tab over. Anything
# else (including unset) keeps Avivator, so existing deployments are unaffected.

resolve_viewer_backend <- function() {
  backend <- tolower(trimws(Sys.getenv("ISLET_VIEWER_BACKEND", "avivator")))
  if (!nzchar(backend)) backend <- "avivator"
  if (!backend %in% c("avivator", "tissuumaps")) {
    warning(sprintf("Unknown ISLET_VIEWER_BACKEND '%s'; falling back to avivator", backend))
    backend <- "avivator"
  }
  backend
}

VIEWER_BACKEND <- resolve_viewer_backend()

#' Base URL of the TissUUmaps service.
#'
#' Defaults to a same-origin `/tissuumaps/` route (RC proxies it to the
#' container, so no CORS and no mixed content). Override with TISSUUMAPS_URL,
#' e.g. `http://127.0.0.1:5100` for local development.
resolve_tissuumaps_base <- function() {
  env_url <- Sys.getenv("TISSUUMAPS_URL", "")
  if (nzchar(env_url)) return(sub("/+$", "", env_url))
  "/tissuumaps"
}

#' Filesystem location of the `.tmap` projects, if the app can see it.
#'
#' Only used to populate the donor dropdown; the tiles themselves are always
#' fetched from the TissUUmaps service, never from here.
resolve_tissuumaps_project_dir <- function() {
  env_dir <- Sys.getenv("TISSUUMAPS_PROJECT_DIR", "")
  if (nzchar(env_dir) && dir.exists(env_dir)) return(env_dir)
  candidates <- c(
    if (!is.null(project_root)) file.path(project_root, "data", "tissuumaps"),
    file.path("..", "..", "data", "tissuumaps")
  )
  for (cand in candidates) {
    if (!is.null(cand) && dir.exists(cand)) return(cand)
  }
  NULL
}

#' List available TissUUmaps projects as bare case ids.
#'
#' Falls back to the OME-TIFF stems when the project directory isn't visible to
#' the app (a split deployment where only the container mounts the volume) —
#' build_tissuumaps_project.py names each project after the same case id.
list_available_tmap_projects <- function() {
  dir <- resolve_tissuumaps_project_dir()
  if (!is.null(dir)) {
    projects <- list.files(dir, pattern = "\\.tmap$", ignore.case = TRUE)
    if (length(projects)) return(sort(sub("\\.tmap$", "", projects, ignore.case = TRUE)))
  }

  images_dir <- file.path("www", "local_images")
  if (dir.exists(images_dir)) {
    files <- list.files(images_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
    if (length(files)) return(sort(sub("\\.ome\\.tiff?$", "", files, ignore.case = TRUE)))
  }

  character(0)
}

#' Full iframe URL for one TissUUmaps project.
#'
#' TissUUmaps switches itself into embedded mode when `window != window.top`
#' (collapses the navbar, shrinks the side menu, disables passive wheel events),
#' so no extra query parameter is needed to make it iframe-friendly.
build_tissuumaps_iframe_url <- function(project) {
  if (is.null(project) || !nzchar(project)) return(NULL)
  base <- resolve_tissuumaps_base()
  stem <- sub("\\.tmap$", "", project, ignore.case = TRUE)
  url <- paste0(base, "/", utils::URLencode(stem, reserved = TRUE), ".tmap")

  # Projects stored in a subdirectory of the server's slide root need ?path=;
  # TissUUmaps resolves every relative asset in the .tmap against it.
  sub_path <- Sys.getenv("TISSUUMAPS_PROJECT_PATH", "")
  if (nzchar(sub_path)) {
    url <- paste0(url, "?path=", utils::URLencode(sub_path, reserved = TRUE))
  }
  url
}

resolve_avivator_base <- function() {
  # Phase 2: if the operator points us at a hosted Avivator instance (self-
  # hosted or the public demo), use that. Falls back to the bundled static
  # build under www/ for legacy single-box deployments.
  env_url <- Sys.getenv("AVIVATOR_URL", "")
  if (nzchar(env_url)) return(env_url)

  local_index <- file.path("www", "avivator", "index.html")
  if (file.exists(local_index)) return("avivator/index.html")
  NULL
}

# ---- Phase 2: OME-Zarr on MinIO ------------------------------------------
# When images are hosted as OME-Zarr stores in an S3-compatible bucket we no
# longer need nginx to serve large OME-TIFFs with Range headers. The Shiny
# app just hands Avivator a fully-qualified URL pointing at the `.zarr`
# directory.

#' Base URL for the OME-Zarr bucket (MinIO or cloud S3).
#'
#' Defaults to the production MinIO deployment. Override per-environment via
#' the `MINIO_IMAGE_URL` env var. Example values:
#'
#'   * MinIO on prod: http://10.15.152.7:9000/islet-images
#'   * nginx-proxied: https://islet-explorer.example.org/images
#'   * S3: https://islet-explorer.s3.amazonaws.com
resolve_minio_base <- function() {
  env_url <- Sys.getenv("MINIO_IMAGE_URL", "")
  if (nzchar(env_url)) return(sub("/+$", "", env_url))
  # Reasonable default that matches the production nginx location block.
  "http://10.15.152.7:9000/islet-images"
}

#' Build the Zarr store URL for a given image name.
#'
#' `image_name` may include or omit the `.zarr` suffix. The returned URL is
#' always fully qualified so Avivator/Viv can pass it straight to `source=`.
build_zarr_image_url <- function(image_name) {
  if (is.null(image_name) || !nzchar(image_name)) return(NULL)
  base <- resolve_minio_base()
  stem <- sub("\\.ome\\.tiff?$", "", image_name, ignore.case = TRUE)
  stem <- sub("\\.zarr$", "", stem, ignore.case = TRUE)
  paste0(base, "/", stem, ".zarr")
}

#' Construct the full iframe URL for Avivator, pointing at a Zarr store.
#'
#' Keeps the channel configuration payload (base64-encoded JSON) that the
#' bundled Avivator build expects, so disease-specific channel assignments
#' (INS/GCG/SST/DAPI) still apply.
build_avivator_iframe_url <- function(image_name, channel_config_b64 = NULL) {
  base <- resolve_avivator_base()
  if (is.null(base)) return(NULL)
  zarr_url <- build_zarr_image_url(image_name)
  if (is.null(zarr_url)) return(base)

  params <- list(
    source = utils::URLencode(zarr_url, reserved = TRUE),
    # image_url is retained so the legacy Avivator bundle (which reads this
    # query key) still binds without a code change on the JS side.
    image_url = utils::URLencode(zarr_url, reserved = FALSE)
  )
  if (!is.null(channel_config_b64) && nzchar(channel_config_b64)) {
    params$channel_config <- utils::URLencode(channel_config_b64, reserved = TRUE)
  }

  query <- paste(
    vapply(names(params), function(k) sprintf("%s=%s", k, params[[k]]), character(1)),
    collapse = "&"
  )
  paste0(base, if (grepl("\\?", base)) "&" else "?", query)
}

#' List available OME-Zarr images.
#'
#' Priority order:
#'   1. `ISLET_ZARR_MANIFEST` env var pointing at a newline-delimited file of
#'      image stems (one per line) — most reliable for MinIO deployments.
#'   2. `data/zarr_images/*.zarr` on the local filesystem (dev / prod where
#'      the scaling conversion has been run).
#'   3. Legacy `www/local_images/*.ome.tif[f]` (fallback for unconverted
#'      deployments — the stems are still valid MinIO keys if the operator
#'      has uploaded them).
list_available_zarr_images <- function() {
  manifest <- Sys.getenv("ISLET_ZARR_MANIFEST", "")
  if (nzchar(manifest) && file.exists(manifest)) {
    lines <- readLines(manifest, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines) & !grepl("^#", lines)]
    if (length(lines) > 0) return(lines)
  }

  zarr_root <- file.path("..", "..", "data", "app_data", "zarr_images")
  if (!is.null(project_root)) {
    alt <- file.path(project_root, "data", "app_data", "zarr_images")
    if (dir.exists(alt)) zarr_root <- alt
  }
  if (dir.exists(zarr_root)) {
    dirs <- list.dirs(zarr_root, full.names = FALSE, recursive = FALSE)
    dirs <- sub("\\.zarr$", "", dirs, ignore.case = TRUE)
    if (length(dirs) > 0) return(dirs)
  }

  # Last-resort legacy fallback: expose OME-TIFF stems.
  legacy_dir <- file.path("www", "local_images")
  if (dir.exists(legacy_dir)) {
    files <- list.files(legacy_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
    return(sub("\\.ome\\.tiff?$", "", files, ignore.case = TRUE))
  }

  character(0)
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
build_channel_config_b64 <- function(names_vec = NULL) {
  # Emit the default-on marker set BY NAME only (no positional `channelNames`
  # override, no hardcoded per-channel indices). The embedded Avivator resolves
  # each primaryChannel against the SELECTED image's own OME-XML channel names
  # (its `findIndex(name)` fallback), and — with `channelNames` omitted — keeps
  # each image's own channel labels. This is required because the cohort mixes
  # panels: batch-1 slides are 59-ch with SST@35, batch-2 slides put SST@1 and
  # some are 58-ch (GCG@45). A single positional sidecar mislabels/mis-defaults
  # batch-2; name resolution is correct for every donor regardless of panel.
  # `names_vec` is accepted for backward compatibility but no longer used.
  rgb_hex <- function(r, g, b) sprintf("#%02X%02X%02X", as.integer(r), as.integer(g), as.integer(b))

  # Order: INS (red), GCG (blue), SST (yellow), DAPI (grey)
  prim <- list(
    list(name = "INS",  color = rgb_hex(228,  26,  28), visible = TRUE),
    list(name = "GCG",  color = rgb_hex( 55, 126, 184), visible = TRUE),
    list(name = "SST",  color = rgb_hex(255, 204,   0), visible = TRUE),
    list(name = "DAPI", color = rgb_hex(127, 127, 127), visible = TRUE)
  )

  cfg <- list(primaryChannels = prim)
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
