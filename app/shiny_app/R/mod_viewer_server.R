viewer_server <- function(id, forced_image, current_tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    viewer_info <- reactive({
      forced_img <- forced_image()
      env_info <- detect_environment(session)
      base <- resolve_avivator_base()

      info <- list(
        base = base, selection = NULL, mode = "picker", ok = FALSE,
        iframe_src = NULL, image_url = NULL, image_public_url = NULL,
        image_rel = NULL, image_app_url = NULL, asset_diag = NULL, env = env_info
      )
      if (is.null(base)) return(info)

      cat("[VIEWER] Environment detection - reverse_proxy:", env_info$is_reverse_proxy,
          "hostname:", env_info$hostname, "pathname:", env_info$pathname, "\n")

      params <- list()
      rel_url <- NULL

      sel_basename <- forced_img
      if (is.null(sel_basename)) sel_basename <- input$selected_image

      cat("[VIEWER] === Image Selection Debug ===\n")
      cat("[VIEWER] Forced image:", if (is.null(forced_img)) "NULL" else paste0("'", forced_img, "'"), "\n")
      cat("[VIEWER] Input image:", if (is.null(input$selected_image)) "NULL" else paste0("'", input$selected_image, "'"), "\n")
      cat("[VIEWER] Selected basename:", if (is.null(sel_basename)) "NULL" else paste0("'", sel_basename, "'"), "\n")

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
        } else {
          rel_url <- paste("local_images", sel_basename, sep = "/")
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
            if (!is.null(asset_diag$http_error)) paste("err=", asset_diag$http_error) else "", "\n")
      }

      resolved_app_url <- NULL
      resolved_public_url <- NULL
      if (!is.null(rel_url) && nzchar(rel_url)) {
        info$image_rel <- rel_url
        resolved_app_url <- build_app_absolute_url(session, rel_url)
        resolved_public_url <- build_public_http_url(session, rel_url)
        info$image_app_url <- resolved_app_url
        info$image_public_url <- resolved_public_url
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
        params[["image_url"]] <- utils::URLencode(sel_url, reserved = FALSE)
        info$image_param <- sel_url
        if (is.null(info$image_public_url)) {
          info$image_public_url <- if (!identical(sel_url, resolved_app_url)) sel_url else NULL
        }
        info$image_url <- resolved_app_url %||% sel_url
      }

      ch_b64 <- tryCatch(build_channel_config_b64(channel_names_vec), error = function(e) NULL)
      if (!is.null(ch_b64) && nzchar(ch_b64)) {
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

    output$local_image_picker <- renderUI({
      tab <- current_tab()
      if (is.null(tab) || tab != "Viewer") return(NULL)

      vi <- viewer_info()
      cat("[VIEWER RENDER] Base:", vi$base, "Image URL (app):", vi$image_url,
          "Public:", vi$image_public_url %||% "NULL", "\n")

      if (is.null(vi$base)) {
        return(tagList(
          tags$div(style = "color:#b00;", "Local Avivator static build not found under shiny/www/avivator."),
          tags$div(style = "color:#666; font-size:90%;",
                   "To install: run scripts/install_avivator.sh (requires Node >= 18) or place a prebuilt bundle under shiny/www/avivator.")
        ))
      }

      available_images <- character(0)
      images_dir <- file.path("www", "local_images")
      if (dir.exists(images_dir)) {
        image_files <- list.files(images_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
        available_images <- c(available_images, image_files)
      }
      local_root <- Sys.getenv("LOCAL_IMAGE_ROOT", unset = "")
      if (nzchar(local_root) && dir.exists(local_root)) {
        root_files <- list.files(local_root, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE)
        available_images <- unique(c(available_images, root_files))
      }

      image_basename <- NULL
      if (!is.null(vi$image_rel)) {
        image_basename <- basename(vi$image_rel)
      } else if (!is.null(vi$image_url)) {
        image_basename <- basename(vi$image_url)
      }

      if (!is.null(image_basename) && nzchar(image_basename)) {
        local_image_path <- file.path("www", "local_images", image_basename)
        file_exists <- file.exists(local_image_path)
        cat("[VIEWER DEBUG] Looking for image file:", local_image_path, "Exists:", file_exists, "\n")
        if (file_exists) {
          file_size_gb <- round(file.info(local_image_path)$size / (1024^3), 2)
          cat("[VIEWER DEBUG] Image size:", file_size_gb, "GB\n")
        }
      }

      current_selection <- if (!is.null(forced_image())) {
        forced_image()
      } else if (!is.null(input$selected_image)) {
        input$selected_image
      } else if (length(available_images) > 0) {
        available_images[1]
      } else {
        NULL
      }

      tagList(
        if (length(available_images) > 0) {
          fluidRow(
            column(12,
              selectInput(ns("selected_image"), "Select Image:",
                          choices = available_images, selected = current_selection, width = "100%")
            )
          )
        } else {
          tags$div(
            style = "padding: 10px; background-color: #fff3cd; border: 1px solid #ffc107; margin-bottom: 10px;",
            tags$strong("No images found."),
            " Place OME-TIFF files in ", tags$code("app/shiny_app/www/local_images/"),
            " or set ", tags$code("LOCAL_IMAGE_ROOT"), " environment variable."
          )
        },
        if (VIEWER_DEBUG_ENABLED) {
          tags$details(
            style = "margin: 10px 0;",
            tags$summary("Viewer debug info"),
            tags$pre(
              style = "max-height:200px; overflow:auto;",
              jsonlite::toJSON(list(
                env = vi$env, image_rel = vi$image_rel,
                image_app = vi$image_app_url %||% vi$image_url,
                image_public = vi$image_public_url,
                iframe_src = vi$iframe_src,
                asset_diag = vi$asset_diag
              ), auto_unbox = TRUE, pretty = TRUE)
            )
          )
        },
        if (!is.null(vi$image_url) && nzchar(vi$image_url)) {
          tags$div(
            id = ns("viewer_wrap"),
            style = "position: relative;",
            tags$div(
              style = "padding: 6px 10px; background-color: #f6f8fa; color: #666; font-size: 90%; margin-bottom: 6px; border: 1px solid #e1e4e8; border-radius: 3px;",
              "Click inside the viewer to enable pan/zoom and wheel zoom."
            ),
            tags$iframe(
              id = "avivator_viewer", src = vi$iframe_src,
              width = "100%", frameBorder = 0, allowfullscreen = NA,
              tabindex = "0", allow = "fullscreen",
              style = "width:100%; height:calc(100vh - 180px); border:0;"
            ),
            tags$script(HTML(sprintf("(function(){
  var wrap = document.getElementById('%s');
  if (!wrap) return;
  var iframe = wrap.querySelector('#avivator_viewer');
  if (!iframe) return;
  function focusViewer(){ try { if (iframe.contentWindow) iframe.contentWindow.focus(); } catch (e) {} }
  iframe.addEventListener('load', focusViewer);
  wrap.addEventListener('mouseenter', focusViewer);
  wrap.addEventListener('wheel', function(){
    try {
      if (document.activeElement !== iframe) focusViewer();
    } catch (e) {}
  }, { passive: true });
  setTimeout(focusViewer, 200);
})();", ns("viewer_wrap"))))
          )
        } else {
          tags$div(
            style = "padding: 20px; text-align: center; color: #666;",
            "Select an image from the dropdown above to view it in AVIVATOR."
          )
        }
      )
    })
  })
}
