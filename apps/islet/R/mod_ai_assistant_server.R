ai_assistant_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

    strip_markdown <- function(text) {
      if (is.null(text) || !nzchar(text)) return(text)
      # Remove fenced code block markers (```lang ... ```)
      text <- gsub("```[a-zA-Z]*\n?", "", text)
      # Remove inline code backticks
      text <- gsub("`([^`]+)`", "\\1", text)
      # Remove bold (**...**) and (__...__)
      text <- gsub("\\*\\*([^*]+)\\*\\*", "\\1", text)
      text <- gsub("__([^_]+)__", "\\1", text)
      # Remove italic (*...*) and (_..._)
      text <- gsub("\\*([^*]+)\\*", "\\1", text)
      text <- gsub("(?<![a-zA-Z])_([^_]+)_(?![a-zA-Z])", "\\1", text, perl = TRUE)
      # Remove headers
      text <- gsub("(?m)^#{1,6}\\s+", "", text, perl = TRUE)
      # Remove bullet markers (- or * at start of line), keep the text
      text <- gsub("(?m)^[*\\-]\\s+", "", text, perl = TRUE)
      # Remove numbered list markers
      text <- gsub("(?m)^\\d+\\.\\s+", "", text, perl = TRUE)
      # Convert markdown links [text](url) to just text
      text <- gsub("\\[([^]]+)\\]\\([^)]+\\)", "\\1", text)
      # Remove horizontal rules
      text <- gsub("(?m)^[-*]{3,}$", "", text, perl = TRUE)
      text
    }

    format_chat_content <- function(text) {
      if (is.null(text)) return(htmltools::HTML(""))
      text <- strip_markdown(text)
      text <- stringr::str_replace_all(text, "\r\n", "\n")
      htmltools::HTML(gsub("\n", "<br/>", htmltools::htmlEscape(text)))
    }

    output$chat_status <- renderUI({
      key_val <- get_llm_api_key()
      key_available <- nzchar(key_val)
      if (!key_available) {
        span("LLM key not found. Set KEY (and optionally BASE) in ~/.Renviron to enable assistant responses.")
      } else {
        NULL
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

      placeholder_entry <- list(role = "assistant", content = "\u23f3 \u2026")
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
            replacement <- if (!is.null(partial_text) && nzchar(partial_text)) partial_text else "\u23f3 \u2026"
            hist[[chat_index]]$content <- replacement
            chat_history(hist)
          }
        })
        session$flushReact()
        invisible(NULL)
      }

      shinyjs::disable(ns("chat_send"))
      shinyjs::disable(ns("chat_reset"))
      streaming_active(TRUE)
      on.exit({
        shinyjs::enable(ns("chat_send"))
        shinyjs::enable(ns("chat_reset"))
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
        friendly <- paste0("Assistant error: ", err_text)
        base_val <- tryCatch(LLM_API_BASE, error = function(e) "")
        norm_base <- if (nzchar(base_val)) sub("/+\\z", "", base_val, perl = TRUE) else ""
        if (grepl("exceeded your current quota", err_text, ignore.case = TRUE)) {
          friendly <- paste(
            "The assistant can't reply right now because the OpenAI API quota was exceeded.",
            "Please review plan and billing details at https://platform.openai.com/account/usage before trying again.")
        } else if (grepl("429", err_text, ignore.case = TRUE)) {
          friendly <- paste(
            "The assistant hit the OpenAI rate/usage limit.",
            "Please wait a moment or adjust your plan before trying again.")
        } else if (grepl("401", err_text, ignore.case = TRUE) || grepl("invalid api key", err_text, ignore.case = TRUE)) {
          friendly <- paste(
            "The assistant can't authenticate with the LLM provider.",
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
            sprintf("The assistant couldn't reach %s.", endpoint),
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
  })
}
