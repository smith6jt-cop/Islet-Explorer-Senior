viewer_server <- function(id, forced_image, current_tab) {
  # Auto-launch, as asked for 2026-09-14. Returns FALSE and moves on if a server
  # cannot be reached or started -- the tab then shows its ordinary "no projects"
  # message rather than the app failing to start. See R/tissuumaps_service.R for
  # why this holds a lock until the server ANSWERS: shiny-server runs one R
  # worker per session, and a naive check-then-spawn was measured starting six
  # servers from six concurrent workers.
  if (identical(tolower(Sys.getenv("ISLET_VIEWER_BACKEND", "avivator")), "tissuumaps")) {
    isTRUE(try(ensure_tissuumaps(), silent = TRUE))
  }

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Which viewer renders the pixels. Fixed for the session (env var), so read
    # it once — every branch below is non-reactive on this value.
    use_tissuumaps <- identical(VIEWER_BACKEND, "tissuumaps")
    viewer_label <- if (use_tissuumaps) "TissUUmaps" else "AVIVATOR"

    # Viewer base: the TissUUmaps service URL, or the Avivator static build
    # location (a file existence check) — both fixed for the session, so
    # compute once, non-reactively.
    viewer_base <- if (use_tissuumaps) resolve_tissuumaps_base() else resolve_avivator_base()
    iframe_id <- ns("avivator_viewer")
    wrap_id   <- ns("viewer_wrap")

    # Resolve the active build once per session. Missing membership must never
    # turn a filesystem listing into permission to display every donor.
    viewer_cohort_r <- reactive(read_viewer_cohort())
    cohort_note_r <- reactive(ie_require_measured(
      viewer_cohort_r()$availability, guard = "viewer_cohort"))

    # Images available to pick from. TissUUmaps picks `.tmap` projects (bare
    # case ids); Avivator picks OME-TIFF basenames from www/local_images +
    # LOCAL_IMAGE_ROOT.
    available_images_r <- reactive({
      if (!is.null(cohort_note_r())) return(character(0))
      cohort <- viewer_cohort_r()
      if (use_tissuumaps) return(list_available_tmap_projects(cohort))

      imgs <- character(0)
      images_dir <- file.path("www", "local_images")
      if (dir.exists(images_dir)) {
        imgs <- c(imgs, filter_viewer_images(
          list.files(images_dir, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE), cohort))
      }
      local_root <- Sys.getenv("LOCAL_IMAGE_ROOT", unset = "")
      if (nzchar(local_root) && dir.exists(local_root)) {
        imgs <- unique(c(imgs, filter_viewer_images(
          list.files(local_root, pattern = "\\.ome\\.tiff?$", ignore.case = TRUE), cohort)))
      }
      imgs
    })

    # The basename currently being viewed (forced image wins, then the user's
    # dropdown choice, then the first available image).
    selected_basename_r <- reactive({
      available <- available_images_r()
      valid_choice <- function(value) {
        if (use_tissuumaps) {
          valid <- is.character(value) && length(value) == 1L && !is.na(value) &&
            islet_safe_path_segment(value) && value %in% available
          if (!valid && !is.null(value)) {
            islet_log_validation_failure("viewer_image_invalid")
          }
          return(valid)
        }
        islet_valid_asset_basename(value, available)
      }
      fi <- forced_image()
      if (!is.null(fi) && valid_choice(fi)) return(fi)
      si <- input$selected_image
      if (!is.null(si) && valid_choice(si)) return(si)
      if (length(available) > 0) return(available[[1]])
      NULL
    })

    # Coarse "is there something to show" flag. Computed reactively, then
    # mirrored into a reactiveVal so it only invalidates dependents when the
    # boolean actually FLIPS. Switching between images keeps it TRUE, so
    # output$viewer_frame does NOT re-render — which would destroy + recreate
    # the WebGL <iframe> and cause the macOS WebRender flicker. Image-URL
    # changes are instead pushed into the existing iframe (see observeEvent).
    has_image_now <- reactive({
      sb <- selected_basename_r()
      !is.null(viewer_base) && !is.null(sb) && nzchar(sb)
    })
    has_image_stable <- reactiveVal(FALSE)
    observe({ has_image_stable(isTRUE(has_image_now())) })

    viewer_info <- reactive({
      info <- list(
        base = viewer_base, selection = NULL, ok = FALSE,
        iframe_src = NULL, image_url = NULL, image_public_url = NULL,
        image_rel = NULL, image_app_url = NULL, asset_diag = NULL, env = NULL
      )
      if (is.null(viewer_base)) return(info)

      sel_basename <- selected_basename_r()

      # TissUUmaps serves its own tiles, so there is no local asset to resolve:
      # the project name is all we need to build the URL.
      if (use_tissuumaps) {
        if (is.null(sel_basename) || !nzchar(sel_basename)) return(info)
        info$selection  <- sel_basename
        info$iframe_src <- build_tissuumaps_iframe_url(sel_basename)
        info$image_url  <- info$iframe_src
        info$ok         <- !is.null(info$iframe_src)
        return(info)
      }

      rel_url <- NULL
      if (!is.null(sel_basename) && nzchar(sel_basename)) {
        www_candidate <- file.path("www", "local_images", sel_basename)
        resource_candidate <- NULL
        if (!is.null(local_images_root) && dir.exists(local_images_root)) {
          resource_candidate <- file.path(local_images_root, sel_basename)
        }
        if (file.exists(www_candidate)) {
          rel_url <- paste("local_images", sel_basename, sep = "/")
        } else if (!is.null(resource_candidate) && file.exists(resource_candidate)) {
          rel_url <- paste("images", sel_basename, sep = "/")
        } else {
          rel_url <- paste("local_images", sel_basename, sep = "/")
        }
        info$selection <- sel_basename
      }

      resolved_app_url <- NULL
      resolved_public_url <- NULL
      if (!is.null(rel_url) && nzchar(rel_url)) {
        info$image_rel <- rel_url
        # The session URL origin (protocol/host/port/pathname) is fixed for the
        # lifetime of the session. Isolate the clientData reads so any churn in
        # session$clientData$url_* cannot re-trigger this reactive — which would
        # otherwise rebuild the viewer iframe and flicker the WebGL canvas.
        resolved_app_url    <- isolate(build_app_absolute_url(session, rel_url))
        resolved_public_url <- isolate(build_public_http_url(session, rel_url))
        info$image_app_url    <- resolved_app_url
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

      params <- list()
      if (!is.null(sel_url)) {
        params[["image_url"]] <- utils::URLencode(sel_url, reserved = FALSE)
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
      info$iframe_src <- if (!is.null(query)) paste0(viewer_base, "?", query) else viewer_base
      info$ok <- TRUE

      # Diagnostics only — detect_environment() and probe_viewer_asset() both
      # read session$clientData / do a blocking curl HEAD, so keep them off the
      # hot path unless explicitly debugging.
      if (VIEWER_DEBUG_ENABLED) {
        info$env <- isolate(detect_environment(session))
        info$asset_diag <- isolate(probe_viewer_asset(session, rel_url))
      }
      info
    })

    # ---- Image picker (mounts once; independent of the iframe) --------------
    output$image_selector <- renderUI({
      tab <- current_tab()
      if (is.null(tab) || tab != "Viewer") return(NULL)
      if (is.null(viewer_base)) return(NULL)  # missing-build message lives in viewer_frame
      if (!is.null(cohort_note_r())) return(NULL)  # named note lives in viewer_frame

      available_images <- available_images_r()
      # isolate(): the selector reflects the current choice on first mount but
      # must not take a reactive dependency on it, or picking an image would
      # rebuild the dropdown.
      current_selection <- isolate(selected_basename_r())

      tagList(
        if (length(available_images) > 0) {
          fluidRow(
            column(12,
              selectInput(ns("selected_image"),
                          if (use_tissuumaps) "Select Donor:" else "Select Image:",
                          choices = available_images, selected = current_selection, width = "100%")
            )
          )
        } else if (use_tissuumaps) {
          tags$div(
            style = "padding: 10px; margin-bottom: 10px",
            tags$strong("No TissUUmaps projects found."),
            " Build one per donor with ", tags$code("scripts/senior/build_tissuumaps_project.py"),
            " and point ", tags$code("TISSUUMAPS_PROJECT_DIR"), " at the output directory."
          )
        } else {
          tags$div(
            style = "padding: 10px; margin-bottom: 10px",
            tags$strong("No images found."),
            " Place OME-TIFF files in ", tags$code("apps/islet/www/local_images/"),
            " or set ", tags$code("LOCAL_IMAGE_ROOT"), " environment variable."
          )
        },
        if (VIEWER_DEBUG_ENABLED) {
          vi <- isolate(viewer_info())
          tags$details(
            style = "margin: 10px 0",
            tags$summary("Viewer debug info"),
            tags$pre(
              style = "max-height:200px; overflow:auto",
              jsonlite::toJSON(list(
                env = vi$env, image_rel = vi$image_rel,
                image_app = vi$image_app_url %||% vi$image_url,
                image_public = vi$image_public_url,
                iframe_src = vi$iframe_src,
                asset_diag = vi$asset_diag
              ), auto_unbox = TRUE, pretty = TRUE)
            )
          )
        }
      )
    })

    # ---- Avivator iframe (mounts once when an image first becomes available) -
    output$viewer_frame <- renderUI({
      tab <- current_tab()
      if (is.null(tab) || tab != "Viewer") return(NULL)

      if (is.null(viewer_base)) {
        return(tagList(
          tags$div("Local Avivator static build not found under shiny/www/avivator."),
          tags$div(style = "font-size:90%",
                   "To install: run scripts/install_avivator.sh (requires Node >= 18) or place a prebuilt bundle under shiny/www/avivator.")
        ))
      }

      cohort_note <- cohort_note_r()
      if (!is.null(cohort_note)) return(tags$div(format(cohort_note)))

      # Depend on the FLIP-only boolean (mounts once); read the URL via isolate
      # so this renderUI does NOT re-fire — and recreate the iframe — when the
      # src changes. URL changes are pushed in place via the observeEvent below.
      ready <- has_image_stable()
      if (!isTRUE(ready)) {
        return(tags$div(
          style = "padding: 20px; text-align: center",
          sprintf("Select an image from the dropdown above to view it in %s.", viewer_label)
        ))
      }

      src0 <- isolate(viewer_info()$iframe_src)
      if (is.null(src0) || !nzchar(src0)) {
        return(tags$div(
          style = "padding: 20px; text-align: center",
          sprintf("Select an image from the dropdown above to view it in %s.", viewer_label)
        ))
      }

      tags$div(
        id = wrap_id,
        style = "position: relative",
        tags$div(
          style = "padding: 6px 10px; font-size: 90%; margin-bottom: 6px",
          if (use_tissuumaps) {
            paste("Click inside the viewer to enable pan/zoom and wheel zoom.",
                  "Islet outlines load automatically; use the Regions tab to draw or export annotations,",
                  "and the Markers tab to overlay single cells.")
          } else {
            "Click inside the viewer to enable pan/zoom and wheel zoom."
          }
        ),
        tags$iframe(
          id = iframe_id, src = src0, `data-relsrc` = src0,
          title = "Interactive pancreatic tissue image viewer",
          width = "100%", frameBorder = 0, allowfullscreen = NA,
          tabindex = "0", allow = "fullscreen",
          style = "width:100%; height:calc(100vh - 180px)"
        ),
        tags$script(HTML(sprintf("(function(){
  var wrap = document.getElementById('%s');
  if (!wrap) return;
  var iframe = wrap.querySelector('#%s');
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
})();", wrap_id, iframe_id)))
      )
    })

    # ---- Push src changes into the existing iframe (no DOM recreation) -------
    observeEvent(viewer_info()$iframe_src, {
      src <- viewer_info()$iframe_src
      if (is.null(src) || !nzchar(src)) return()
      session$sendCustomMessage("islet_update_viewer_src",
                                list(id = iframe_id, src = src))
    }, ignoreNULL = TRUE)
  })
}
