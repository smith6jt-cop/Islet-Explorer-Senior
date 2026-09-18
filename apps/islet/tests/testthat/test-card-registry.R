# Exercise the real registry dispatch and card reactives without scientific data.
# Missing artifacts must produce a named note, including with no activation args.
skip_if_not_installed("shiny")

.card_registry_sources <- function() {
  env <- new.env(parent = globalenv())
  for (name in c("NS", "tagList", "div", "tags", "moduleServer", "reactive")) {
    env[[name]] <- getExportedValue("shiny", name)
  }
  files <- c(paste0("card_", names(.card_registry_outputs), ".R"),
             "mod_crossmodal_ui.R", "mod_crossmodal_server.R",
             "mod_lymphnode_ui.R", "mod_lymphnode_server.R")
  for (file in files) sys.source(.find_source(file.path("R", file)), envir = env)

  # Read the authoritative keys without running the data-loading globals.
  env$SHOW_INFLAMMATION_TAB <- FALSE
  for (expr in parse(.find_source("R/00_globals.R"))) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        as.character(expr[[2L]]) %in% c("IE_PIPELINE", "IE_REFERENCE", "IE_TAB_KEYS")) {
      eval(expr, envir = env)
    }
  }
  env
}

.card_registry_outputs <- c(
  s1_size_composition = "s1_body",
  s3_compartments = "s3_compartments_body", x1_gradient = "x1_gradient_body",
  x2_overlay = "x2_overlay_content", x3_identity = "x3_identity_content",
  g1_grid = "g1_body"
)
.card_registry_probes <- c(
  s1_size_composition = "layer", s3_compartments = "selected",
  x1_gradient = "selected", x2_overlay = "selection", x3_identity = "selection",
  g1_grid = "layer"
)
.card_registries <- function(env) {
  list(
    list(module = "crossmodal", tab = env$IE_CROSSMODAL_TAB, cards = env$IE_CROSSMODAL_CARDS),
    list(module = "lymphnode", tab = env$IE_LYMPHNODE_TAB, cards = env$IE_LYMPHNODE_CARDS)
  )
}

.expect_card_note <- function(output, id) {
  html <- output[[.card_registry_outputs[[id]]]]$html
  expect_true(is.character(html) && nzchar(html), info = paste(id, "rendered a note"))
  expect_match(html, paste0("[", if (id %in% c("x1_gradient", "s3_compartments")) "guard_", id, "]"),
               fixed = TRUE, info = paste(id, "names its guard"))
  expect_match(html, "Artifact:", fixed = TRUE, info = paste(id, "names its missing input"))
}

test_that("every registered card is registered once and is actually implemented", {
  env <- .card_registry_sources()
  ids <- character()
  for (registry in .card_registries(env)) {
    expect_true(registry$tab %in% env$IE_TAB_KEYS)
    html <- htmltools::renderTags(env[[paste0(registry$module, "_ui")]](registry$module))$html
    for (card in registry$cards) {
      ids <- c(ids, card$id)
      ui <- paste0("card_", card$id, "_ui")
      server <- paste0("card_", card$id, "_server")
      expect_match(html, paste0('data-card-id="', card$id, '"'), fixed = TRUE)
      if (card$id %in% c("ledger", "n1_sections")) {
        expect_false(exists(ui, envir = env, mode = "function"))
        expect_false(exists(server, envir = env, mode = "function"))
        expect_match(html, card$title, fixed = TRUE)
        expect_match(html, paste0("unavailable(not_built): no ", ui, "() is defined yet."),
                     fixed = TRUE)
      } else {
        expect_true(is.function(env[[ui]]), info = ui)
        expect_true(is.function(env[[server]]), info = server)
        args <- formals(env[[server]])
        expect_identical(tail(names(args), 3L), c("active_tab", "tab_key", "..."))
        expect_null(args$active_tab)
        expect_null(args$tab_key)
      }
    }
  }
  expect_false(anyDuplicated(ids) > 0L)

  # EVERY REGISTERED CARD MUST BE BUILT. This test used to allow two placeholders by name --
  # "ledger" and "n1_sections" -- and so blessed the thing it should have caught: both were the
  # FIRST card on their tab, neither had a card_*.R file or a server anywhere, and both rendered
  # unavailable on every donor. A registry entry is a promise that a card exists. Unregistered is
  # honest; "unavailable" over a card nobody wrote is not. Derived, so a new placeholder fails
  # here rather than in front of a reader.
  for (id in ids) {
    expect_true(is.function(env[[paste0("card_", id, "_ui")]]),
                info = paste0("card_", id, "_ui is registered but not implemented"))
    expect_true(is.function(env[[paste0("card_", id, "_server")]]),
                info = paste0("card_", id, "_server is registered but not implemented"))
  }
  expect_setequal(ids, names(.card_registry_outputs))
})

test_that("both registries reach every implemented server and gate its actual output on their tab", {
  withr::local_envvar(c(ISLET_LAYERS_DIR = tempfile("missing-card-layers-")))
  for (registry_index in 1:2) {
    env <- .card_registry_sources()
    registry <- .card_registries(env)[[registry_index]]
    calls <- new.env(parent = emptyenv())
    # Observe dispatch, but run the original server with the original arguments.
    for (id in names(.card_registry_outputs)) {
      local({
        card_id <- id
        name <- paste0("card_", card_id, "_server")
        original <- env[[name]]
        env[[name]] <- function(...) {
          value <- original(...)
          calls[[card_id]] <- list(tab_key = list(...)$tab_key, value = value)
          value
        }
      })
    }
    current_tab <- shiny::reactiveVal(NULL)
    expected <- setdiff(vapply(registry$cards, `[[`, "", "id"), c("ledger", "n1_sections"))
    shiny::testServer(env[[paste0(registry$module, "_server")]],
      args = list(id = registry$module, prepared = shiny::reactive(NULL), active_tab = current_tab), {
        expect_setequal(ls(calls), expected)
        for (id in expected) {
          expect_identical(calls[[id]]$tab_key, registry$tab, info = id)
          expect_false(calls[[id]]$value$active(), info = paste(id, "before a tab is selected"))
        }
        # Include every other real tab, particularly the former Plot/Spatial
        # guesses and the other registry. Return home to check reactivation.
        for (tab in c(registry$tab, setdiff(env$IE_TAB_KEYS, registry$tab), registry$tab)) {
          current_tab(tab)
          session$flushReact()
          for (id in expected) {
            card <- calls[[id]]$value
            expect_identical(card$active(), identical(tab, registry$tab),
                             info = paste(registry$module, id, "on", tab))
            probe <- card[[.card_registry_probes[[id]]]]
            if (identical(tab, registry$tab)) {
              expect_type(probe(), "list")
              .expect_card_note(output, id)
            } else {
              expect_error(probe(), class = "shiny.silent.error",
                           info = paste(id, "must not load on", tab))
            }
          }
        }
      })
  }
})

test_that("standalone cards render with optional activation args and obey an arbitrary mount key", {
  withr::local_envvar(c(ISLET_LAYERS_DIR = tempfile("missing-card-layers-")))
  env <- .card_registry_sources()
  scenarios <- list(
    list(), list(active_tab = function() "Elsewhere"), list(tab_key = "TestMount"),
    list(active_tab = function() "TestMount", tab_key = "TestMount")
  )
  for (id in names(.card_registry_outputs)) {
    for (args in scenarios) {
      run_card <- function(input, output, session) {
        card <- do.call(env[[paste0("card_", id, "_server")]],
                        c(list(input = input, output = output, session = session), args))
      }
      shiny::testServer(run_card, {
        session$flushReact()
        expect_true(card$active(), info = id)
        .expect_card_note(output, id)
      })
    }
  }
})

test_that("the within-donor trajectory card is mounted on the Trajectory tab, not a registry", {
  env <- .card_registry_sources()
  # Moved 2026-09-18: it is a within-donor islet ordering, not a cross-modality comparison, and it
  # was on Cross-Modal only because that was the registry that first carried it. It is now wired
  # directly by trajectory_server, so it appears in NEITHER registry -- and this test is what stops
  # it being quietly re-added to one, or dropped from the tab it moved to.
  for (registry in .card_registries(env)) {
    expect_false("s2_trajectory" %in% vapply(registry, `[[`, "", "id"))
  }
  expect_true(is.function(env$card_s2_trajectory_ui))
  expect_true(is.function(env$card_s2_trajectory_server))

  ui <- paste(readLines("../../R/mod_trajectory_ui.R"), collapse = " ")
  expect_true(grepl("card_s2_trajectory_ui(ns)", ui, fixed = TRUE))
  server <- paste(readLines("../../R/mod_trajectory_server.R"), collapse = " ")
  expect_true(grepl("card_s2_trajectory_server(", server, fixed = TRUE))
  # It must take its tab key rather than guess one: card_s1 shipped with tab_key = "spatial",
  # which matches nothing in IE_TAB_KEYS, so its gate was never true and it could not draw.
  expect_true(grepl('tab_key = "Trajectory"', server, fixed = TRUE))
})
