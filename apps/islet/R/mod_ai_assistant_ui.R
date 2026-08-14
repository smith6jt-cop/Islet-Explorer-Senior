ai_assistant_ui <- function(id) {
  ns <- NS(id)
  div(class = "card ai-chat-panel",
    style = "margin-bottom: 20px; padding: 15px; border: 1px solid #ddd; border-radius: 8px; height: clamp(715px, 84.5vh - 97px, 1560px); display: flex; flex-direction: column; gap: 12px;",
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
    uiOutput(ns("chat_status"), container = div, class = "ai-chat-status"),
    div(style = "flex: 1; overflow-y: auto;",
      uiOutput(ns("chat_history"), container = div, class = "ai-chat-history")
    ),
    div(
      class = "ai-chat-model-picker",
      selectInput(
        ns("chat_model"),
        label = "Model",
        choices = c(
          "Navigator Fast (gpt-oss-20b)" = "gpt-oss-20b",
          "Navigator Large (gpt-oss-120b)" = "gpt-oss-120b"
        ),
        selected = "gpt-oss-20b"
      )
    ),
    div(
      class = "ai-chat-input",
      textAreaInput(
        ns("chat_message"),
        label = NULL,
        value = "",
        placeholder = "Ask about donor trends, plots, or troubleshooting\u2026",
        height = "110px"
      )
    ),
    div(
      class = "ai-chat-controls",
      actionButton(ns("chat_send"), "Send", class = "btn btn-primary"),
      actionButton(ns("chat_reset"), "New Conversation", class = "btn btn-outline-secondary")
    ),
    uiOutput(ns("chat_feedback"), container = div, class = "ai-chat-feedback"),
    uiOutput(ns("chat_usage"), container = div, class = "ai-chat-usage")
  )
}
