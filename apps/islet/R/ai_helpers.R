# ai_helpers.R
# AI assistant helper functions for the Islet Explorer Shiny app.
# Extracted from app.R lines 335-920.
#
# Dependencies from global.R:
#   - DEBUG_CREDS   (logical flag, set via DEBUG_CREDENTIALS env var)
#   - httr2_available (logical, TRUE if httr2 namespace can be loaded)

AI_SYSTEM_PROMPT <- paste(
  "You are an AI assistant embedded in the Islet Explorer Shiny app.",
  "Help users interpret plots, statistics, and data preparation steps without hallucinating.",
  "Favor concise, actionable answers rooted in the app's current state and available controls."
)

# ===== Credential Loading =====
# Set DEBUG_CREDENTIALS=1 to enable verbose logging

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
    return(paste0(substr(trimmed, 1, max_chars - 1), "\u2026"))
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
                             max_output_tokens = 4096,
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
    is_reasoning <- FALSE
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
            # If transitioning from reasoning to content, clear the thinking indicator
            if (is_reasoning) {
              accumulated <<- ""
              is_reasoning <<- FALSE
            }
            accumulated <<- paste0(accumulated, piece)
            if (is.function(stream_callback)) {
              stream_callback(accumulated, usage_capture)
            }
          } else if (!is.null(delta$reasoning_content) && is.function(stream_callback)) {
            # Reasoning model: show thinking indicator while reasoning
            is_reasoning <<- TRUE
            stream_callback("\u2728 Thinking...", usage_capture)
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
      httr2::req_timeout(120) |>
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
      return(stream_error(detail_chat, fallback = status_chat %in% c(401, 404, 405)))
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
        httr2::req_timeout(120) |>
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

        fallback_due_to_model <- status_code %in% c(400, 401, 404) &&
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
        httr2::req_timeout(120) |>
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

        fallback_due_to_model_chat <- status_chat %in% c(400, 401, 404) &&
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
