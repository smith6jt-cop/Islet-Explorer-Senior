# Exercise the real registry dispatch and card reactives without scientific data.
# Missing artifacts must produce a named note, including with no activation args.
#
# NOTHING HERE NAMES THE CARDS. Every card this file checks is read from the two
# registries, its sources from R/card_*.R, its body output from its own UI and
# its gated reactive from what its server returns. The three hardcoded lists that
# used to sit here -- the file list, an id -> output map and an id -> probe map --
# were each a second copy of a fact the card already carried, and a card added to
# a registry by another hand was simply absent from all three.
skip_if_not_installed("shiny")

.card_registry_sources <- function() {
  env <- new.env(parent = globalenv())
  for (name in c("NS", "tagList", "div", "tags", "moduleServer", "reactive")) {
    env[[name]] <- getExportedValue("shiny", name)
  }
  app_r <- dirname(.find_source("R/00_globals.R"))
  files <- c(sort(list.files(app_r, pattern = "^card_.*[.][Rr]$")),
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

.card_registries <- function(env) {
  list(
    list(module = "crossmodal", tab = env$IE_CROSSMODAL_TAB, cards = env$IE_CROSSMODAL_CARDS),
    list(module = "lymphnode", tab = env$IE_LYMPHNODE_TAB, cards = env$IE_LYMPHNODE_CARDS)
  )
}

#' Every card id both registries mount, in registry order.
.card_ids <- function(env) {
  unlist(lapply(.card_registries(env), function(r) vapply(r$cards, `[[`, "", "id")))
}

#' The output a card renders its own body into -- read from the card's UI.
#'
#' The convention every card follows is one output named `*_body` or `*_content`
#' carrying the whole card, with any others being controls or sub-panels. A card
#' that declares neither, or two, fails here rather than being read on an output
#' that was never the card's answer.
.card_body_output <- function(id) {
  ui <- NULL
  for (expr in parse(.find_source(file.path("R", paste0("card_", id, ".R"))),
                     keep.source = FALSE)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name(paste0("card_", id, "_ui")))) ui <- expr[[3L]]
  }
  outputs <- character()
  rec <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    parts <- as.list(e)
    head <- parts[[1L]]
    name <- if (is.name(head)) as.character(head)
            else if (is.call(head) && identical(head[[1L]], as.name("::"))) as.character(head[[3L]])
            else ""
    if (grepl("Output$", name) && length(parts) >= 2L) {
      arg <- parts[[2L]]
      if (is.call(arg) && identical(arg[[1L]], as.name("ns")) && is.character(arg[[2L]])) {
        outputs <<- c(outputs, arg[[2L]])
      }
    }
    for (i in seq_along(parts)) {
      if (!identical(unname(parts[i]), list(quote(expr = )))) rec(parts[[i]])
    }
    invisible(NULL)
  }
  rec(ui)
  body <- grep("_(body|content)$", outputs, value = TRUE)
  expect_length(body, 1L)
  body[[1L]]
}

#' The returned reactive the tab gate must close over.
#'
#' Read from the card's own return value, not from a map beside it: a card hands
#' back `layer`, `selected` or `selection` as the reactive that loads its
#' artifact, and that is the one that must refuse to run off its tab.
.card_gate <- function(card, id) {
  gate <- intersect(c("layer", "selected", "selection"), names(card))
  expect_gt(length(gate), 0L)
  card[[gate[[1L]]]]
}

.expect_card_note <- function(output, id) {
  html <- output[[.card_body_output(id)]]$html
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
      # NO PLACEHOLDER EXEMPTION. "ledger" and "n1_sections" were named here and
      # allowed to have no card_*.R at all; both registries dropped them on
      # 2026-09-18 and a real n1_sections card is being written, so the exemption
      # is dead code that would re-bless the next unbuilt entry. Every registered
      # card is checked the same way.
      expect_true(is.function(env[[ui]]), info = ui)
      expect_true(is.function(env[[server]]), info = server)
      args <- formals(env[[server]])
      expect_identical(tail(names(args), 3L), c("active_tab", "tab_key", "..."))
      expect_null(args$active_tab)
      expect_null(args$tab_key)
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
  # Each card names one body output of its own, and no two cards share one.
  bodies <- vapply(ids, .card_body_output, "")
  expect_false(anyDuplicated(bodies) > 0L)
})

test_that("both registries reach every implemented server and gate its actual output on their tab", {
  withr::local_envvar(c(ISLET_LAYERS_DIR = tempfile("missing-card-layers-")))
  for (registry_index in 1:2) {
    env <- .card_registry_sources()
    registry <- .card_registries(env)[[registry_index]]
    calls <- new.env(parent = emptyenv())
    # Observe dispatch, but run the original server with the original arguments.
    for (id in .card_ids(env)) {
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
    expected <- vapply(registry$cards, `[[`, "", "id")
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
            probe <- .card_gate(card, id)
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
  for (id in .card_ids(env)) {
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
    # registry$cards, NOT registry: a registry is list(module=, tab=, cards=), so
    # vapply over it subscripted the character scalar "crossmodal" and the test
    # ERRORED instead of asserting. Line 116 above had it right all along.
    expect_false("s2_trajectory" %in% vapply(registry$cards, `[[`, "", "id"))
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

# ===========================================================================
# The app's own wiring: what it CALLS must exist, and the rail's Help entries
# must open rather than throw.
#
# tests/parse_check.R parses every source file, which proves each one is
# syntactically valid and NOTHING about whether a call has a callee. Two
# functions were called and defined nowhere in the repo at once --
# `ie_tab_tips_modal()` from app.R's rail Help observer, on every tab, and
# `load_cell_distance()` from mod_inflammation_server.R, latent only because
# SHOW_INFLAMMATION_TAB defaults to 0 -- and both parsed cleanly for days.
# ===========================================================================

# An empty argument (`df[, 1]`, a formal with no default) is the missing symbol.
# Touching its VALUE raises "argument is missing", so it is tested on the
# one-element sublist, which never forces it.
.arg_is_empty <- function(parts, i) identical(unname(parts[i]), list(quote(expr = )))

#' Every name an expression could BIND, at any depth: assignment targets, for
#' variables, assign("x", ...) and the formals of any function inside it.
#' Deliberately a superset -- this test hunts names that are bound NOWHERE, and
#' over-collecting here can only make it quieter, never wrong.
.names_bound_in <- function(expr) {
  out <- character()
  rec <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    parts <- as.list(e)
    head <- parts[[1L]]
    if (is.name(head)) {
      nm <- as.character(head)
      if (nm %in% c("<-", "=", "<<-") && length(parts) >= 3L && !.arg_is_empty(parts, 2L) &&
          (is.name(parts[[2L]]) || is.character(parts[[2L]]))) {
        out <<- c(out, as.character(parts[[2L]]))
      } else if (nm == "for" && length(parts) >= 3L && !.arg_is_empty(parts, 2L) &&
                 is.name(parts[[2L]])) {
        out <<- c(out, as.character(parts[[2L]]))
      } else if (nm == "assign" && length(parts) >= 2L && !.arg_is_empty(parts, 2L) &&
                 is.character(parts[[2L]])) {
        out <<- c(out, parts[[2L]])
      } else if (nm == "function" && length(parts) >= 2L) {
        out <<- c(out, names(as.list(parts[[2L]])))
      }
    }
    for (i in seq_along(parts)) if (!.arg_is_empty(parts, i)) rec(parts[[i]])
    invisible(NULL)
  }
  rec(expr)
  unique(out)
}

#' Every function NAME called in `expr`, minus the ones bound in the scope that
#' reaches it. A `pkg::fn()` call is skipped: its head is a call, not a name, and
#' R resolves it itself.
.unbound_calls_in <- function(expr, scope) {
  out <- character()
  walk <- function(e, scope) {
    if (!is.call(e)) return(invisible(NULL))
    parts <- as.list(e)
    head <- parts[[1L]]
    if (is.name(head)) {
      nm <- as.character(head)
      if (!nm %in% scope) out <<- c(out, nm)
      if (nm == "function") {
        formals_list <- as.list(parts[[2L]])
        inner <- c(scope, names(formals_list), .names_bound_in(parts[[3L]]))
        for (i in seq_along(formals_list)) {
          if (!.arg_is_empty(formals_list, i)) walk(formals_list[[i]], scope)
        }
        walk(parts[[3L]], inner)
        return(invisible(NULL))
      }
    } else {
      walk(head, scope)
    }
    for (i in seq_along(parts)[-1L]) if (!.arg_is_empty(parts, i)) walk(parts[[i]], scope)
    invisible(NULL)
  }
  walk(expr, scope)
  unique(out)
}

#' ieCore, then every app source, then app.R -- the files Shiny itself loads.
.app_source_files <- function() {
  app_dir <- dirname(.find_source("app.R"))
  c(sort(list.files(.iecore_dir, pattern = "[.][Rr]$", full.names = TRUE)),
    sort(list.files(file.path(app_dir, "R"), pattern = "[.][Rr]$", full.names = TRUE)),
    file.path(app_dir, "app.R"))
}

#' The packages the app ATTACHES, read out of R/00_globals.R rather than listed
#' here, so this check follows the app instead of a copy of it.
.app_attached_packages <- function() {
  out <- character()
  rec <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    parts <- as.list(e)
    if (is.name(parts[[1L]]) && as.character(parts[[1L]]) %in% c("library", "require") &&
        length(parts) >= 2L && !.arg_is_empty(parts, 2L)) {
      arg <- parts[[2L]]
      if (is.name(arg) || is.character(arg)) out <<- c(out, as.character(arg))
    }
    for (i in seq_along(parts)) if (!.arg_is_empty(parts, i)) rec(parts[[i]])
    invisible(NULL)
  }
  for (expr in parse(.find_source("R/00_globals.R"))) rec(expr)
  unique(out)
}

test_that("every function the app calls has a definition somewhere", {
  packages <- .app_attached_packages()
  # Never pass by skipping quietly on ONE package: a missing package would make
  # its functions look undefined, so the check is all of them or none of them.
  for (pkg in packages) skip_if_not_installed(pkg)

  files <- .app_source_files()
  expect_gt(length(files), 20L)

  # Top-level bindings from every file, because Shiny sources them all into one
  # environment and order between them is not guaranteed.
  defined <- unlist(lapply(files, function(f) {
    unlist(lapply(parse(f, keep.source = FALSE), function(e) {
      if (is.call(e) && length(e) >= 3L && is.name(e[[1L]]) &&
          as.character(e[[1L]]) %in% c("<-", "=", "<<-") &&
          (is.name(e[[2L]]) || is.character(e[[2L]]))) as.character(e[[2L]]) else NULL
    }))
  }))
  exported <- unlist(lapply(packages, getNamespaceExports))

  missing <- do.call(rbind, lapply(files, function(f) {
    called <- unlist(lapply(parse(f, keep.source = FALSE), function(e) {
      # A top-level `local({...})` or `if (...) {...}` binds names for the block
      # it encloses, so the expression's own bindings join the scope it sees.
      .unbound_calls_in(e, c(defined, .names_bound_in(e)))
    }))
    called <- setdiff(unique(called), c(defined, exported))
    called <- called[!vapply(called, exists, logical(1), mode = "function")]
    if (!length(called)) NULL else data.frame(name = called, file = basename(f))
  }))

  expect_null(missing, info = if (is.null(missing)) "" else
    paste0("called but defined nowhere: ",
           paste(sprintf("%s() in %s", missing$name, missing$file), collapse = "; ")))
})

#' The `observeEvent(input$help_*, ...)` expressions, lifted out of app.R's own
#' server body. Copying them into the test would test the copy.
.app_help_observers <- function() {
  server <- NULL
  for (expr in parse(.find_source("app.R"), keep.source = FALSE)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("server"))) server <- expr[[3L]]
  }
  if (is.null(server)) stop("app.R defines no `server`")
  found <- list()
  rec <- function(e) {
    if (!is.call(e)) return(invisible(NULL))
    parts <- as.list(e)
    if (identical(parts[[1L]], as.name("observeEvent")) && length(parts) >= 2L &&
        !.arg_is_empty(parts, 2L)) {
      trigger <- paste(deparse(parts[[2L]]), collapse = "")
      if (grepl("^input\\$help_", trigger)) found[[trigger]] <<- e
    }
    for (i in seq_along(parts)) if (!.arg_is_empty(parts, i)) rec(parts[[i]])
    invisible(NULL)
  }
  rec(server[[3L]])
  found
}

test_that("the rail's Help entries open a modal on every tab instead of throwing", {
  env <- .card_registry_sources()
  for (name in c("observeEvent", "modalDialog", "modalButton", "icon", "req")) {
    env[[name]] <- getExportedValue("shiny", name)
  }
  # The Help modal and its two helpers, evaluated from app.R itself.
  for (expr in parse(.find_source("app.R"), keep.source = FALSE)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        as.character(expr[[2L]]) %in% c("IE_TAB_TIPS", ".ie_tab_card_titles",
                                        ".ie_tab_states_ui", "ie_tab_tips_modal")) {
      eval(expr, envir = env)
    }
  }
  expect_true(is.function(env$ie_tab_tips_modal))

  observers <- .app_help_observers()
  expect_setequal(names(observers), c("input$help_rules", "input$help_tips"))

  seen <- new.env(parent = emptyenv())
  seen$modals <- list()
  env$showModal <- function(ui, session = NULL) {
    seen$modals[[length(seen$modals) + 1L]] <- ui
    invisible(ui)
  }
  # The rules modal reads the domain's rule file; this test is about the
  # observers' wiring, and the call-graph test above is what proves the real
  # builder exists.
  env$phenotype_rules_modal_ui <- function() {
    shiny::modalDialog(title = "Cell Phenotype Rules", "stub")
  }

  run <- function(input, output, session) {
    active_tab <- shiny::reactive(input$tabs)
    for (observer in observers) eval(observer)
  }
  # The observers are evaluated with app.R's own definitions in scope, and with
  # the test's `observers` reachable, so the code under test is the source file's
  # and not a copy of it.
  run_env <- new.env(parent = env)
  run_env$observers <- observers
  environment(run) <- run_env

  shiny::testServer(run, {
    # One flush BEFORE the first input: an observeEvent registered inside
    # testServer consumes its ignoreInit pass on the first flush, so without this
    # every assertion below would read one step behind the input it just set.
    session$flushReact()
    session$setInputs(help_rules = 1L)
    session$flushReact()
    expect_length(seen$modals, 1L)
    expect_match(htmltools::renderTags(seen$modals[[1L]])$html,
                 "Cell Phenotype Rules", fixed = TRUE)

    tabs <- c(env$IE_TAB_KEYS, "NotATab")
    for (i in seq_along(tabs)) {
      session$setInputs(tabs = tabs[i])
      session$setInputs(help_tips = i)   # observeEvent fires on a CHANGE
      session$flushReact()
      expect_length(seen$modals, i + 1L)
      html <- htmltools::renderTags(seen$modals[[i + 1L]])$html
      label <- Filter(function(d) identical(d$key, tabs[i]),
                      c(env$IE_PIPELINE, env$IE_REFERENCE))
      label <- if (length(label)) label[[1L]]$label else "This tab"
      expect_match(html, paste0(label, " — how to read it"), fixed = TRUE,
                   info = tabs[i])
      # Every tab's modal says what an absence means, in ieCore's own words.
      expect_match(html, ie_state_label("unavailable"), fixed = TRUE, info = tabs[i])
    }
  })
})

test_that("the tips modal reads its content back from the registries and ieCore", {
  env <- .card_registry_sources()
  for (expr in parse(.find_source("app.R"), keep.source = FALSE)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        as.character(expr[[2L]]) %in% c("IE_TAB_TIPS", ".ie_tab_card_titles",
                                        ".ie_tab_states_ui", "ie_tab_tips_modal")) {
      eval(expr, envir = env)
    }
  }
  # Card titles come from the registry, not from a second list that could drift.
  for (registry in .card_registries(env)) {
    titles <- env$.ie_tab_card_titles(registry$tab)
    expect_identical(titles, vapply(registry$cards, `[[`, "", "title"))
    html <- htmltools::renderTags(env$ie_tab_tips_modal(registry$tab))$html
    for (title in titles) expect_match(html, title, fixed = TRUE)
  }
  expect_identical(env$.ie_tab_card_titles("Plot"), character(0))

  # Every state ieCore declares is named, so a new state cannot appear in a card
  # note while the modal still explains four of five.
  html <- htmltools::renderTags(env$ie_tab_tips_modal("Plot"))$html
  for (state in IE_STATES) expect_match(html, ie_state_label(state), fixed = TRUE)

  # A key that is not a tab still opens, and says nothing it cannot support.
  html <- htmltools::renderTags(env$ie_tab_tips_modal(NULL))$html
  expect_match(html, "This tab — how to read it", fixed = TRUE)
  expect_false(grepl("Cards on this tab", html, fixed = TRUE))
})

#' Every character literal in an expression, with no evaluation.
#'
#' NA_character_ is a character literal too, and an NA in here turns every
#' comparison below into NA rather than FALSE -- which is how the first draft of
#' this test reported three violations that did not exist.
.string_literals_in <- function(expr) {
  out <- character()
  rec <- function(e) {
    if (is.character(e)) { out <<- c(out, e[!is.na(e)]); return(invisible(NULL)) }
    if (!is.call(e)) return(invisible(NULL))
    parts <- as.list(e)
    for (i in seq_along(parts)) {
      if (!identical(unname(parts[i]), list(quote(expr = )))) rec(parts[[i]])
    }
    invisible(NULL)
  }
  rec(expr)
  out
}

test_that("the .ie-card class marks a registry card and nothing else", {
  # tests/tab_smoke.R counts the cards a tab mounted with
  # `.ie-card-stack > .ie-card`, so that selector is a contract, not a style
  # choice: card_g1_grid.R drew its four tissue x modality cells as .ie-card
  # inside an .ie-card-stack of its own, and Cross-Modal reported ELEVEN cards
  # against a registry of seven and failed its own harness with every card
  # drawing correctly. One .ie-card per card, in its UI; one .ie-card-stack per
  # registry, in the registry.
  app_r <- dirname(.find_source("R/00_globals.R"))
  cards <- list.files(app_r, pattern = "^card_.*[.][Rr]$")
  expect_gt(length(cards), 0L)

  for (file in cards) {
    id <- sub("^card_(.*)[.][Rr]$", "\\1", file)
    roots <- character()
    others <- character()
    for (expr in parse(file.path(app_r, file), keep.source = FALSE)) {
      literals <- .string_literals_in(expr)
      target <- if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
                    is.name(expr[[2L]])) as.character(expr[[2L]]) else ""
      if (identical(target, paste0("card_", id, "_ui"))) {
        roots <- c(roots, literals)
      } else {
        others <- c(others, literals)
      }
    }
    expect_identical(sum(roots == "ie-card"), 1L, info = file)
    expect_identical(sum(others == "ie-card"), 0L, info = file)
    expect_identical(sum(c(roots, others) == "ie-card-stack"), 0L, info = file)
  }

  for (file in c("mod_crossmodal_ui.R", "mod_lymphnode_ui.R")) {
    literals <- unlist(lapply(parse(file.path(app_r, file), keep.source = FALSE),
                              .string_literals_in))
    expect_identical(sum(literals == "ie-card-stack"), 1L, info = file)
    expect_identical(sum(literals == "ie-card"), 0L, info = file)
  }
})

# ===========================================================================
# load_cell_distance(), the other function that was called and never defined.
#
# It belongs in a test-spatial-helpers.R; it is here because this is the file
# this lane owns, and an untested definition is how the first one was lost.
# ===========================================================================

test_that("load_cell_distance returns the frame the Inflammation tab reads, or NULL", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  env <- new.env(parent = globalenv())
  sys.source(.find_source("R/spatial_helpers.R"), envir = env)
  expect_true(is.function(env$load_cell_distance))

  # No connection at all: a missing view is NULL, which is what the module's
  # need() reports. Never an empty frame -- "no view" and "no cells" are
  # different statements and a reader must not see one as the other.
  expect_null(env$load_cell_distance())

  con <- DBI::dbConnect(duckdb::duckdb())
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
  assign("con", con, envir = globalenv())
  withr::defer(rm("con", envir = globalenv()))

  # A view without the columns this frame promises is NULL and names them,
  # rather than a frame of NA that cannot be told from unmeasured cells.
  DBI::dbExecute(con, "CREATE VIEW cell_distance AS SELECT 1 AS case_id")
  expect_message(expect_null(env$load_cell_distance()), "cell_distance_columns_absent")
  DBI::dbExecute(con, "DROP VIEW cell_distance")

  DBI::dbExecute(con, "CREATE TABLE cd (case_id VARCHAR, donor_status VARCHAR,
                                        cell_region VARCHAR, phenotype VARCHAR,
                                        dist_islet DOUBLE, ins_pos BOOLEAN)")
  DBI::dbExecute(con, "INSERT INTO cd VALUES
    ('6505','ND','core','Beta cell',-12.5,true),
    ('6505','ND','peri','CD8a Tcell',8.25,false),
    ('6520','T1D','tissue','Macrophage',NULL,false)")
  DBI::dbExecute(con, "CREATE VIEW cell_distance AS SELECT * FROM cd")

  df <- env$load_cell_distance()
  expect_identical(colnames(df),
                   c("case_id", "donor_status", "cell_region", "phenotype", "signed_dist"))
  expect_identical(nrow(df), 3L)
  # Character, not factor: table() over a factor emits a zero for every level the
  # subset lacks, and a phenotype that is absent would be indistinguishable from
  # one that was counted and found zero.
  for (col in c("case_id", "donor_status", "cell_region", "phenotype")) {
    expect_type(df[[col]], "character")
  }
  # SIGNED, and a distance that was never measured stays NA. Zero is the islet
  # edge, so a coerced zero would place an unmeasured cell exactly on it.
  expect_identical(df$signed_dist[df$phenotype == "Beta cell"], -12.5)
  expect_true(is.na(df$signed_dist[df$phenotype == "Macrophage"]))
  expect_identical(sum(is.na(df$signed_dist)), 1L)

  # The same object on the second call: the module holds it in one reactive and
  # the frame is the whole view.
  expect_identical(env$load_cell_distance(), df)
})
