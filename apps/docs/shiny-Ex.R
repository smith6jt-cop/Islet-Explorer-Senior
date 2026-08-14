pkgname <- "shiny"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('shiny')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ExtendedTask")
### * ExtendedTask

flush(stderr()); flush(stdout())

### Name: ExtendedTask
### Title: Task or computation that proceeds in the background
### Aliases: ExtendedTask

### ** Examples

## Don't show: 
if (rlang::is_interactive() && rlang::is_installed("mirai")) (if (getRversion() >= "3.4") withAutoprint else force)({ # examplesIf
## End(Don't show)
library(shiny)
library(bslib)
library(mirai)

# Set background processes for running tasks
daemons(1)
# Reset when the app is stopped
onStop(function() daemons(0))

ui <- page_fluid(
  titlePanel("Extended Task Demo"),
  p(
    'Click the button below to perform a "calculation"',
    "that takes a while to perform."
  ),
  input_task_button("recalculate", "Recalculate"),
  p(textOutput("result"))
)

server <- function(input, output) {
  rand_task <- ExtendedTask$new(function() {
    mirai(
      {
        # Slow operation goes here
        Sys.sleep(2)
        sample(1:100, 1)
      }
    )
  })

  # Make button state reflect task.
  # If using R >=4.1, you can do this instead:
  # rand_task <- ExtendedTask$new(...) |> bind_task_button("recalculate")
  bind_task_button(rand_task, "recalculate")

  observeEvent(input$recalculate, {
    # Invoke the extended in an observer
    rand_task$invoke()
  })

  output$result <- renderText({
    # React to updated results when the task completes
    number <- rand_task$result()
    paste0("Your number is ", number, ".")
  })
}

shinyApp(ui, server)
## Don't show: 
}) # examplesIf
## End(Don't show)



cleanEx()
nameEx("MockShinySession")
### * MockShinySession

flush(stderr()); flush(stdout())

### Name: MockShinySession
### Title: Mock Shiny Session
### Aliases: MockShinySession

### ** Examples


## ------------------------------------------------
## Method `MockShinySession$setInputs`
## ------------------------------------------------

## Not run: 
##D session$setInputs(x=1, y=2)
## End(Not run)



cleanEx()
nameEx("Progress")
### * Progress

flush(stderr()); flush(stdout())

### Name: Progress
### Title: Reporting progress (object-oriented API)
### Aliases: Progress

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  plotOutput("plot")
)

server <- function(input, output, session) {
  output$plot <- renderPlot({
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())

    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')

    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
    plot(cars)
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("actionButton")
### * actionButton

flush(stderr()); flush(stdout())

### Name: actionButton
### Title: Action button/link
### Aliases: actionButton actionLink

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("obs", "Number of observations", 0, 1000, 500),
  actionButton("goButton", "Go!", class = "btn-success"),
  plotOutput("distPlot")
)

server <- function(input, output) {
  output$distPlot <- renderPlot({
    # Take a dependency on input$goButton. This will run once initially,
    # because the value changes from NULL to 0.
    input$goButton

    # Use isolate() to avoid dependency on input$obs
    dist <- isolate(rnorm(input$obs))
    hist(dist)
  })
}

shinyApp(ui, server)

}

## Example of adding extra class values
actionButton("largeButton", "Large Primary Button", class = "btn-primary btn-lg")
actionLink("infoLink", "Information Link", class = "btn-info")




cleanEx()
nameEx("bindCache")
### * bindCache

flush(stderr()); flush(stdout())

### Name: bindCache
### Title: Add caching with reactivity to an object
### Aliases: bindCache

### ** Examples

## Not run: 
##D rc <- bindCache(
##D   x = reactive({
##D     Sys.sleep(2)   # Pretend this is expensive
##D     input$x * 100
##D   }),
##D   input$x
##D )
##D 
##D # Can make it prettier with the %>% operator
##D library(magrittr)
##D 
##D rc <- reactive({
##D   Sys.sleep(2)
##D   input$x * 100
##D }) %>%
##D   bindCache(input$x)
##D 
## End(Not run)

## Only run app examples in interactive R sessions
if (interactive()) {

# Basic example
shinyApp(
  ui = fluidPage(
    sliderInput("x", "x", 1, 10, 5),
    sliderInput("y", "y", 1, 10, 5),
    div("x * y: "),
    verbatimTextOutput("txt")
  ),
  server = function(input, output) {
    r <- reactive({
      # The value expression is an _expensive_ computation
      message("Doing expensive computation...")
      Sys.sleep(2)
      input$x * input$y
    }) %>%
      bindCache(input$x, input$y)

    output$txt <- renderText(r())
  }
)


# Caching renderText
shinyApp(
  ui = fluidPage(
    sliderInput("x", "x", 1, 10, 5),
    sliderInput("y", "y", 1, 10, 5),
    div("x * y: "),
    verbatimTextOutput("txt")
  ),
  server = function(input, output) {
    output$txt <- renderText({
      message("Doing expensive computation...")
      Sys.sleep(2)
      input$x * input$y
    }) %>%
      bindCache(input$x, input$y)
  }
)


# Demo of using events and caching with an actionButton
shinyApp(
  ui = fluidPage(
    sliderInput("x", "x", 1, 10, 5),
    sliderInput("y", "y", 1, 10, 5),
    actionButton("go", "Go"),
    div("x * y: "),
    verbatimTextOutput("txt")
  ),
  server = function(input, output) {
    r <- reactive({
      message("Doing expensive computation...")
      Sys.sleep(2)
      input$x * input$y
    }) %>%
      bindCache(input$x, input$y) %>%
      bindEvent(input$go)
      # The cached, eventified reactive takes a reactive dependency on
      # input$go, but doesn't use it for the cache key. It uses input$x and
      # input$y for the cache key, but doesn't take a reactive dependency on
      # them, because the reactive dependency is superseded by addEvent().

    output$txt <- renderText(r())
  }
)

}




cleanEx()
nameEx("bookmarkButton")
### * bookmarkButton

flush(stderr()); flush(stdout())

### Name: bookmarkButton
### Title: Create a button for bookmarking/sharing
### Aliases: bookmarkButton

### ** Examples

## Only run these examples in interactive sessions
if (interactive()) {

# This example shows how to use multiple bookmark buttons. If you only need
# a single bookmark button, see examples in ?enableBookmarking.
ui <- function(request) {
  fluidPage(
    tabsetPanel(id = "tabs",
      tabPanel("One",
        checkboxInput("chk1", "Checkbox 1"),
        bookmarkButton(id = "bookmark1")
      ),
      tabPanel("Two",
        checkboxInput("chk2", "Checkbox 2"),
        bookmarkButton(id = "bookmark2")
      )
    )
  )
}
server <- function(input, output, session) {
  # Need to exclude the buttons from themselves being bookmarked
  setBookmarkExclude(c("bookmark1", "bookmark2"))

  # Trigger bookmarking with either button
  observeEvent(input$bookmark1, {
    session$doBookmark()
  })
  observeEvent(input$bookmark2, {
    session$doBookmark()
  })
}
enableBookmarking(store = "url")
shinyApp(ui, server)
}



cleanEx()
nameEx("brushedPoints")
### * brushedPoints

flush(stderr()); flush(stdout())

### Name: brushedPoints
### Title: Find rows of data selected on an interactive plot.
### Aliases: brushedPoints nearPoints

### ** Examples

## Not run: 
##D # Note that in practice, these examples would need to go in reactives
##D # or observers.
##D 
##D # This would select all points within 5 pixels of the click
##D nearPoints(mtcars, input$plot_click)
##D 
##D # Select just the nearest point within 10 pixels of the click
##D nearPoints(mtcars, input$plot_click, threshold = 10, maxpoints = 1)
##D 
## End(Not run)



cleanEx()
nameEx("busyIndicatorOptions")
### * busyIndicatorOptions

flush(stderr()); flush(stdout())

### Name: busyIndicatorOptions
### Title: Customize busy indicator options
### Aliases: busyIndicatorOptions

### ** Examples

## Don't show: 
if (rlang::is_interactive()) (if (getRversion() >= "3.4") withAutoprint else force)({ # examplesIf
## End(Don't show)

library(bslib)

card_ui <- function(id, spinner_type = id) {
  card(
    busyIndicatorOptions(spinner_type = spinner_type),
    card_header(paste("Spinner:", spinner_type)),
    plotOutput(shiny::NS(id, "plot"))
  )
}

card_server <- function(id, simulate = reactive()) {
  moduleServer(
    id = id,
    function(input, output, session) {
      output$plot <- renderPlot({
        Sys.sleep(1)
        simulate()
        plot(x = rnorm(100), y = rnorm(100))
      })
    }
  )
}

ui <- page_fillable(
  useBusyIndicators(),
  input_task_button("simulate", "Simulate", icon = icon("refresh")),
  layout_columns(
    card_ui("ring"),
    card_ui("bars"),
    card_ui("dots"),
    card_ui("pulse"),
    col_widths = 6
  )
)

server <- function(input, output, session) {
  simulate <- reactive(input$simulate)
  card_server("ring", simulate)
  card_server("bars", simulate)
  card_server("dots", simulate)
  card_server("pulse", simulate)
}

shinyApp(ui, server)
## Don't show: 
}) # examplesIf
## End(Don't show)



cleanEx()
nameEx("checkboxGroupInput")
### * checkboxGroupInput

flush(stderr()); flush(stdout())

### Name: checkboxGroupInput
### Title: Checkbox Group Input Control
### Aliases: checkboxGroupInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  checkboxGroupInput("variable", "Variables to show:",
                     c("Cylinders" = "cyl",
                       "Transmission" = "am",
                       "Gears" = "gear")),
  tableOutput("data")
)

server <- function(input, output, session) {
  output$data <- renderTable({
    mtcars[, c("mpg", input$variable), drop = FALSE]
  }, rownames = TRUE)
}

shinyApp(ui, server)

ui <- fluidPage(
  checkboxGroupInput("icons", "Choose icons:",
    choiceNames =
      list(icon("calendar"), icon("bed"),
           icon("cog"), icon("bug")),
    choiceValues =
      list("calendar", "bed", "cog", "bug")
  ),
  textOutput("txt")
)

server <- function(input, output, session) {
  output$txt <- renderText({
    icons <- paste(input$icons, collapse = ", ")
    paste("You chose", icons)
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("checkboxInput")
### * checkboxInput

flush(stderr()); flush(stdout())

### Name: checkboxInput
### Title: Checkbox Input Control
### Aliases: checkboxInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  checkboxInput("somevalue", "Some value", FALSE),
  verbatimTextOutput("value")
)
server <- function(input, output) {
  output$value <- renderText({ input$somevalue })
}
shinyApp(ui, server)
}




cleanEx()
nameEx("column")
### * column

flush(stderr()); flush(stdout())

### Name: column
### Title: Create a column within a UI definition
### Aliases: column

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  fluidRow(
    column(4,
      sliderInput("obs", "Number of observations:",
                  min = 1, max = 1000, value = 500)
    ),
    column(8,
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs))
  })
}

shinyApp(ui, server)



ui <- fluidPage(
  fluidRow(
    column(width = 4,
      "4"
    ),
    column(width = 3, offset = 2,
      "3 offset 2"
    )
  )
)
shinyApp(ui, server = function(input, output) { })
}



cleanEx()
nameEx("conditionalPanel")
### * conditionalPanel

flush(stderr()); flush(stdout())

### Name: conditionalPanel
### Title: Conditional Panel
### Aliases: conditionalPanel

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  ui <- fluidPage(
    sidebarPanel(
      selectInput("plotType", "Plot Type",
        c(Scatter = "scatter", Histogram = "hist")
      ),
      # Only show this panel if the plot type is a histogram
      conditionalPanel(
        condition = "input.plotType == 'hist'",
        selectInput(
          "breaks", "Breaks",
          c("Sturges", "Scott", "Freedman-Diaconis", "[Custom]" = "custom")
        ),
        # Only show this panel if Custom is selected
        conditionalPanel(
          condition = "input.breaks == 'custom'",
          sliderInput("breakCount", "Break Count", min = 1, max = 50, value = 10)
        )
      )
    ),
    mainPanel(
      plotOutput("plot")
    )
  )

  server <- function(input, output) {
    x <- rnorm(100)
    y <- rnorm(100)

    output$plot <- renderPlot({
      if (input$plotType == "scatter") {
        plot(x, y)
      } else {
        breaks <- input$breaks
        if (breaks == "custom") {
          breaks <- input$breakCount
        }

        hist(x, breaks = breaks)
      }
    })
  }

  shinyApp(ui, server)
}



cleanEx()
nameEx("createRenderFunction")
### * createRenderFunction

flush(stderr()); flush(stdout())

### Name: createRenderFunction
### Title: Implement custom render functions
### Aliases: createRenderFunction quoToFunction installExprFunction

### ** Examples

# A custom render function that repeats the supplied value 3 times
renderTriple <- function(expr) {
  # Wrap user-supplied reactive expression into a function
  func <- quoToFunction(rlang::enquo0(expr))

  createRenderFunction(
    func,
    transform = function(value, session, name, ...) {
      paste(rep(value, 3), collapse=", ")
    },
    outputFunc = textOutput
  )
}

# For better legacy support, consider using installExprFunction() over quoToFunction()
renderTripleLegacy <- function(expr, env = parent.frame(), quoted = FALSE) {
  func <- installExprFunction(expr, "func", env, quoted)

  createRenderFunction(
    func,
    transform = function(value, session, name, ...) {
      paste(rep(value, 3), collapse=", ")
    },
    outputFunc = textOutput
  )
}

# Test render function from the console
reactiveConsole(TRUE)

v <- reactiveVal("basic")
r <- renderTriple({ v() })
r()
#> [1] "basic, basic, basic"

# User can supply quoted code via rlang::quo(). Note that evaluation of the
# expression happens when r2() is invoked, not when r2 is created.
q <- rlang::quo({ v() })
r2 <- rlang::inject(renderTriple(!!q))
v("rlang")
r2()
#> [1] "rlang, rlang, rlang"

# Supplying quoted code without rlang::quo() requires installExprFunction()
expr <- quote({ v() })
r3 <- renderTripleLegacy(expr, quoted = TRUE)
v("legacy")
r3()
#> [1] "legacy, legacy, legacy"

# The legacy approach also supports with quosures (env is ignored in this case)
q <- rlang::quo({ v() })
r4 <- renderTripleLegacy(q, quoted = TRUE)
v("legacy-rlang")
r4()
#> [1] "legacy-rlang, legacy-rlang, legacy-rlang"

# Turn off reactivity in the console
reactiveConsole(FALSE)




cleanEx()
nameEx("dateInput")
### * dateInput

flush(stderr()); flush(stdout())

### Name: dateInput
### Title: Create date input
### Aliases: dateInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  dateInput("date1", "Date:", value = "2012-02-29"),

  # Default value is the date in client's time zone
  dateInput("date2", "Date:"),

  # value is always yyyy-mm-dd, even if the display format is different
  dateInput("date3", "Date:", value = "2012-02-29", format = "mm/dd/yy"),

  # Pass in a Date object
  dateInput("date4", "Date:", value = Sys.Date()-10),

  # Use different language and different first day of week
  dateInput("date5", "Date:",
          language = "ru",
          weekstart = 1),

  # Start with decade view instead of default month view
  dateInput("date6", "Date:",
            startview = "decade"),

  # Disable Mondays and Tuesdays.
  dateInput("date7", "Date:", daysofweekdisabled = c(1,2)),

  # Disable specific dates.
  dateInput("date8", "Date:", value = "2012-02-29",
            datesdisabled = c("2012-03-01", "2012-03-02"))
)

shinyApp(ui, server = function(input, output) { })
}




cleanEx()
nameEx("dateRangeInput")
### * dateRangeInput

flush(stderr()); flush(stdout())

### Name: dateRangeInput
### Title: Create date range input
### Aliases: dateRangeInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  dateRangeInput("daterange1", "Date range:",
                 start = "2001-01-01",
                 end   = "2010-12-31"),

  # Default start and end is the current date in the client's time zone
  dateRangeInput("daterange2", "Date range:"),

  # start and end are always specified in yyyy-mm-dd, even if the display
  # format is different
  dateRangeInput("daterange3", "Date range:",
                 start  = "2001-01-01",
                 end    = "2010-12-31",
                 min    = "2001-01-01",
                 max    = "2012-12-21",
                 format = "mm/dd/yy",
                 separator = " - "),

  # Pass in Date objects
  dateRangeInput("daterange4", "Date range:",
                 start = Sys.Date()-10,
                 end = Sys.Date()+10),

  # Use different language and different first day of week
  dateRangeInput("daterange5", "Date range:",
                 language = "de",
                 weekstart = 1),

  # Start with decade view instead of default month view
  dateRangeInput("daterange6", "Date range:",
                 startview = "decade")
)

shinyApp(ui, server = function(input, output) { })
}




cleanEx()
nameEx("debounce")
### * debounce

flush(stderr()); flush(stdout())

### Name: debounce
### Title: Slow down a reactive expression with debounce/throttle
### Aliases: debounce throttle

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

library(shiny)
library(magrittr)

ui <- fluidPage(
  plotOutput("plot", click = clickOpts("hover")),
  helpText("Quickly click on the plot above, while watching the result table below:"),
  tableOutput("result")
)

server <- function(input, output, session) {
  hover <- reactive({
    if (is.null(input$hover))
      list(x = NA, y = NA)
    else
      input$hover
  })
  hover_d <- hover %>% debounce(1000)
  hover_t <- hover %>% throttle(1000)

  output$plot <- renderPlot({
    plot(cars)
  })

  output$result <- renderTable({
    data.frame(
      mode = c("raw", "throttle", "debounce"),
      x = c(hover()$x, hover_t()$x, hover_d()$x),
      y = c(hover()$y, hover_t()$y, hover_d()$y)
    )
  })
}

shinyApp(ui, server)
}




cleanEx()
nameEx("devmode")
### * devmode

flush(stderr()); flush(stdout())

### Name: devmode
### Title: Shiny Developer Mode
### Aliases: devmode in_devmode with_devmode devmode_inform
###   register_devmode_option get_devmode_option

### ** Examples

# Enable Shiny Developer mode
devmode()

in_devmode() # TRUE/FALSE?

# Execute code in a temporary shiny dev mode
with_devmode(TRUE, in_devmode()) # TRUE

# Ex: Within shiny, we register the option "shiny.minified"
#   to default to `FALSE` when in Dev Mode
## Not run: 
##D register_devmode_option(
##D   "shiny.minified",
##D   devmode_message = paste0(
##D     "Using full shiny javascript file. ",
##D     "To use the minified version, call `options(shiny.minified = TRUE)`"
##D   ),
##D   devmode_default = FALSE
##D )
## End(Not run)

# Used within `shiny::runApp(launch.browser)`
get_devmode_option("shiny.minified", TRUE) # TRUE if Dev mode is off
is_minified <- with_devmode(TRUE, {
  get_devmode_option("shiny.minified", TRUE)
})
is_minified # FALSE




cleanEx()
nameEx("downloadButton")
### * downloadButton

flush(stderr()); flush(stdout())

### Name: downloadButton
### Title: Create a download button or link
### Aliases: downloadButton downloadLink

### ** Examples

## Not run: 
##D ui <- fluidPage(
##D   p("Choose a dataset to download."),
##D   selectInput("dataset", "Dataset", choices = c("mtcars", "airquality")),
##D   downloadButton("downloadData", "Download")
##D )
##D 
##D server <- function(input, output) {
##D   # The requested dataset
##D   data <- reactive({
##D     get(input$dataset)
##D   })
##D 
##D   output$downloadData <- downloadHandler(
##D     filename = function() {
##D       # Use the selected dataset as the suggested file name
##D       paste0(input$dataset, ".csv")
##D     },
##D     content = function(file) {
##D       # Write the dataset to the `file` that will be downloaded
##D       write.csv(data(), file)
##D     }
##D   )
##D }
##D 
##D shinyApp(ui, server)
## End(Not run)




cleanEx()
nameEx("downloadHandler")
### * downloadHandler

flush(stderr()); flush(stdout())

### Name: downloadHandler
### Title: File Downloads
### Aliases: downloadHandler

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  downloadButton("downloadData", "Download")
)

server <- function(input, output) {
  # Our dataset
  data <- mtcars

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data, file)
    }
  )
}

shinyApp(ui, server)
}




cleanEx()
nameEx("enableBookmarking")
### * enableBookmarking

flush(stderr()); flush(stdout())

### Name: enableBookmarking
### Title: Enable bookmarking for a Shiny application
### Aliases: enableBookmarking

### ** Examples

## Only run these examples in interactive R sessions
if (interactive()) {

# Basic example with state encoded in URL
ui <- function(request) {
  fluidPage(
    textInput("txt", "Text"),
    checkboxInput("chk", "Checkbox"),
    bookmarkButton()
  )
}
server <- function(input, output, session) { }
enableBookmarking("url")
shinyApp(ui, server)


# An alternative to calling enableBookmarking(): use shinyApp's
# enableBookmarking argument
shinyApp(ui, server, enableBookmarking = "url")


# Same basic example with state saved to disk
enableBookmarking("server")
shinyApp(ui, server)


# Save/restore arbitrary values
ui <- function(req) {
  fluidPage(
    textInput("txt", "Text"),
    checkboxInput("chk", "Checkbox"),
    bookmarkButton(),
    br(),
    textOutput("lastSaved")
  )
}
server <- function(input, output, session) {
  vals <- reactiveValues(savedTime = NULL)
  output$lastSaved <- renderText({
    if (!is.null(vals$savedTime))
      paste("Last saved at", vals$savedTime)
    else
      ""
  })

  onBookmark(function(state) {
    vals$savedTime <- Sys.time()
    # state is a mutable reference object, and we can add arbitrary values
    # to it.
    state$values$time <- vals$savedTime
  })
  onRestore(function(state) {
    vals$savedTime <- state$values$time
  })
}
enableBookmarking(store = "url")
shinyApp(ui, server)


# Usable with dynamic UI (set the slider, then change the text input,
# click the bookmark button)
ui <- function(request) {
  fluidPage(
    sliderInput("slider", "Slider", 1, 100, 50),
    uiOutput("ui"),
    bookmarkButton()
  )
}
server <- function(input, output, session) {
  output$ui <- renderUI({
    textInput("txt", "Text", input$slider)
  })
}
enableBookmarking("url")
shinyApp(ui, server)


# Exclude specific inputs (The only input that will be saved in this
# example is chk)
ui <- function(request) {
  fluidPage(
    passwordInput("pw", "Password"), # Passwords are never saved
    sliderInput("slider", "Slider", 1, 100, 50), # Manually excluded below
    checkboxInput("chk", "Checkbox"),
    bookmarkButton()
  )
}
server <- function(input, output, session) {
  setBookmarkExclude("slider")
}
enableBookmarking("url")
shinyApp(ui, server)


# Update the browser's location bar every time an input changes. This should
# not be used with enableBookmarking("server"), because that would create a
# new saved state on disk every time the user changes an input.
ui <- function(req) {
  fluidPage(
    textInput("txt", "Text"),
    checkboxInput("chk", "Checkbox")
  )
}
server <- function(input, output, session) {
  observe({
    # Trigger this observer every time an input changes
    reactiveValuesToList(input)
    session$doBookmark()
  })
  onBookmarked(function(url) {
    updateQueryString(url)
  })
}
enableBookmarking("url")
shinyApp(ui, server)


# Save/restore uploaded files
ui <- function(request) {
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose CSV File", multiple = TRUE,
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        tags$hr(),
        checkboxInput("header", "Header", TRUE),
        bookmarkButton()
      ),
      mainPanel(
        tableOutput("contents")
      )
    )
  )
}
server <- function(input, output) {
  output$contents <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)

    if (nrow(inFile) == 1) {
      read.csv(inFile$datapath, header = input$header)
    } else {
      data.frame(x = "multiple files")
    }
  })
}
enableBookmarking("server")
shinyApp(ui, server)

}



cleanEx()
nameEx("exportTestValues")
### * exportTestValues

flush(stderr()); flush(stdout())

### Name: exportTestValues
### Title: Register expressions for export in test mode
### Aliases: exportTestValues

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {

options(shiny.testmode = TRUE)

# This application shows the test snapshot URL; clicking on it will
# fetch the input, output, and exported values in JSON format.
shinyApp(
  ui = basicPage(
    h4("Snapshot URL: "),
    uiOutput("url"),
    h4("Current values:"),
    verbatimTextOutput("values"),
    actionButton("inc", "Increment x")
  ),

  server = function(input, output, session) {
    vals <- reactiveValues(x = 1)
    y <- reactive({ vals$x + 1 })

    observeEvent(input$inc, {
      vals$x <<- vals$x + 1
    })

    exportTestValues(
      x = vals$x,
      y = y()
    )

    output$url <- renderUI({
      url <- session$getTestSnapshotUrl(format="json")
      a(href = url, url)
    })

    output$values <- renderText({
      paste0("vals$x: ", vals$x, "\ny: ", y())
    })
  }
)
}



cleanEx()
nameEx("fileInput")
### * fileInput

flush(stderr()); flush(stdout())

### Name: fileInput
### Title: File Upload Control
### Aliases: fileInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File", accept = ".csv"),
      checkboxInput("header", "Header", TRUE)
    ),
    mainPanel(
      tableOutput("contents")
    )
  )
)

server <- function(input, output) {
  output$contents <- renderTable({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)

    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))

    read.csv(file$datapath, header = input$header)
  })
}

shinyApp(ui, server)
}




cleanEx()
nameEx("fillPage")
### * fillPage

flush(stderr()); flush(stdout())

### Name: fillPage
### Title: Create a page that fills the window
### Aliases: fillPage

### ** Examples

fillPage(
  tags$style(type = "text/css",
    ".half-fill { width: 50%; height: 100%; }",
    "#one { float: left; background-color: #ddddff; }",
    "#two { float: right; background-color: #ccffcc; }"
  ),
  div(id = "one", class = "half-fill",
    "Left half"
  ),
  div(id = "two", class = "half-fill",
    "Right half"
  ),
  padding = 10
)

fillPage(
  fillRow(
    div(style = "background-color: red; width: 100%; height: 100%;"),
    div(style = "background-color: blue; width: 100%; height: 100%;")
  )
)



cleanEx()
nameEx("fillRow")
### * fillRow

flush(stderr()); flush(stdout())

### Name: fillRow
### Title: Flex Box-based row/column layouts
### Aliases: fillRow fillCol

### ** Examples

# Only run this example in interactive R sessions.
if (interactive()) {

ui <- fillPage(fillRow(
  plotOutput("plotLeft", height = "100%"),
  fillCol(
    plotOutput("plotTopRight", height = "100%"),
    plotOutput("plotBottomRight", height = "100%")
  )
))

server <- function(input, output, session) {
  output$plotLeft <- renderPlot(plot(cars))
  output$plotTopRight <- renderPlot(plot(pressure))
  output$plotBottomRight <- renderPlot(plot(AirPassengers))
}

shinyApp(ui, server)

}



cleanEx()
nameEx("fixedPage")
### * fixedPage

flush(stderr()); flush(stdout())

### Name: fixedPage
### Title: Create a page with a fixed layout
### Aliases: fixedPage fixedRow

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fixedPage(
  title = "Hello, Shiny!",
  fixedRow(
    column(width = 4,
      "4"
    ),
    column(width = 3, offset = 2,
      "3 offset 2"
    )
  )
)

shinyApp(ui, server = function(input, output) { })
}




cleanEx()
nameEx("flowLayout")
### * flowLayout

flush(stderr()); flush(stdout())

### Name: flowLayout
### Title: Flow layout
### Aliases: flowLayout

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- flowLayout(
  numericInput("rows", "How many rows?", 5),
  selectInput("letter", "Which letter?", LETTERS),
  sliderInput("value", "What value?", 0, 100, 50)
)
shinyApp(ui, server = function(input, output) { })
}



cleanEx()
nameEx("fluidPage")
### * fluidPage

flush(stderr()); flush(stdout())

### Name: fluidPage
### Title: Create a page with fluid layout
### Aliases: fluidPage fluidRow

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

# Example of UI with fluidPage
ui <- fluidPage(

  # Application title
  titlePanel("Hello Shiny!"),

  sidebarLayout(

    # Sidebar with a slider input
    sidebarPanel(
      sliderInput("obs",
                  "Number of observations:",
                  min = 0,
                  max = 1000,
                  value = 500)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs))
  })
}

# Complete app with UI and server components
shinyApp(ui, server)


# UI demonstrating column layouts
ui <- fluidPage(
  title = "Hello Shiny!",
  fluidRow(
    column(width = 4,
      "4"
    ),
    column(width = 3, offset = 2,
      "3 offset 2"
    )
  )
)

shinyApp(ui, server = function(input, output) { })
}



cleanEx()
nameEx("freezeReactiveValue")
### * freezeReactiveValue

flush(stderr()); flush(stdout())

### Name: freezeReactiveVal
### Title: Freeze a reactive value
### Aliases: freezeReactiveVal freezeReactiveValue

### ** Examples

## Only run this examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  selectInput("data", "Data Set", c("mtcars", "pressure")),
  checkboxGroupInput("cols", "Columns (select 2)", character(0)),
  plotOutput("plot")
)

server <- function(input, output, session) {
  observe({
    data <- get(input$data)
    # Sets a flag on input$cols to essentially do req(FALSE) if input$cols
    # is accessed. Without this, an error will momentarily show whenever a
    # new data set is selected.
    freezeReactiveValue(input, "cols")
    updateCheckboxGroupInput(session, "cols", choices = names(data))
  })

  output$plot <- renderPlot({
    # When a new data set is selected, input$cols will have been invalidated
    # above, and this will essentially do the same as req(FALSE), causing
    # this observer to stop and raise a silent exception.
    cols <- input$cols
    data <- get(input$data)

    if (length(cols) == 2) {
      plot(data[[ cols[1] ]], data[[ cols[2] ]])
    }
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("getCurrentOutputInfo")
### * getCurrentOutputInfo

flush(stderr()); flush(stdout())

### Name: getCurrentOutputInfo
### Title: Get output information
### Aliases: getCurrentOutputInfo

### ** Examples


if (interactive()) {
  shinyApp(
    fluidPage(
      tags$style(HTML("body {background-color: black; color: white; }")),
      tags$style(HTML("body a {color: purple}")),
      tags$style(HTML("#info {background-color: teal; color: orange; }")),
      plotOutput("p"),
      "Computed CSS styles for the output named info:",
      tagAppendAttributes(
        textOutput("info"),
        class = "shiny-report-theme"
      )
    ),
    function(input, output) {
      output$p <- renderPlot({
        info <- getCurrentOutputInfo()
        par(bg = info$bg(), fg = info$fg(), col.axis = info$fg(), col.main = info$fg())
        plot(1:10, col = info$accent(), pch = 19)
        title("A simple R plot that uses its CSS styling")
      })
      output$info <- renderText({
        info <- getCurrentOutputInfo()
        jsonlite::toJSON(
          list(
            bg = info$bg(),
            fg = info$fg(),
            accent = info$accent(),
            font = info$font()
          ),
          auto_unbox = TRUE
        )
      })
    }
  )
}





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("getQueryString")
### * getQueryString

flush(stderr()); flush(stdout())

### Name: getQueryString
### Title: Get the query string / hash component from the URL
### Aliases: getQueryString getUrlHash

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {

  ## App 1: getQueryString
  ## Printing the value of the query string
  ## (Use the back and forward buttons to see how the browser
  ## keeps a record of each state)
  shinyApp(
    ui = fluidPage(
      textInput("txt", "Enter new query string"),
      helpText("Format: ?param1=val1&param2=val2"),
      actionButton("go", "Update"),
      hr(),
      verbatimTextOutput("query")
    ),
    server = function(input, output, session) {
      observeEvent(input$go, {
        updateQueryString(input$txt, mode = "push")
      })
      output$query <- renderText({
        query <- getQueryString()
        queryText <- paste(names(query), query,
                       sep = "=", collapse=", ")
        paste("Your query string is:\n", queryText)
      })
    }
  )

  ## App 2: getUrlHash
  ## Printing the value of the URL hash
  ## (Use the back and forward buttons to see how the browser
  ## keeps a record of each state)
  shinyApp(
    ui = fluidPage(
      textInput("txt", "Enter new hash"),
      helpText("Format: #hash"),
      actionButton("go", "Update"),
      hr(),
      verbatimTextOutput("hash")
    ),
    server = function(input, output, session) {
      observeEvent(input$go, {
        updateQueryString(input$txt, mode = "push")
      })
      output$hash <- renderText({
        hash <- getUrlHash()
        paste("Your hash is:\n", hash)
      })
    }
  )
}



cleanEx()
nameEx("helpText")
### * helpText

flush(stderr()); flush(stdout())

### Name: helpText
### Title: Create a help text element
### Aliases: helpText

### ** Examples

helpText("Note: while the data view will show only",
         "the specified number of observations, the",
         "summary will be based on the full dataset.")



cleanEx()
nameEx("htmlOutput")
### * htmlOutput

flush(stderr()); flush(stdout())

### Name: htmlOutput
### Title: Create an HTML output element
### Aliases: htmlOutput uiOutput

### ** Examples

htmlOutput("summary")

# Using a custom container and class
tags$ul(
  htmlOutput("summary", container = tags$li, class = "custom-li-output")
)



cleanEx()
nameEx("httpResponse")
### * httpResponse

flush(stderr()); flush(stdout())

### Name: httpResponse
### Title: Create an HTTP response object
### Aliases: httpResponse
### Keywords: internal

### ** Examples

httpResponse(status = 405L,
  content_type = "text/plain",
  content = "The requested method was not allowed"
)




cleanEx()
nameEx("icon")
### * icon

flush(stderr()); flush(stdout())

### Name: icon
### Title: Create an icon
### Aliases: icon

### ** Examples

# add an icon to a submit button
submitButton("Update View", icon = icon("redo"))

navbarPage("App Title",
  tabPanel("Plot", icon = icon("bar-chart-o")),
  tabPanel("Summary", icon = icon("list-alt")),
  tabPanel("Table", icon = icon("table"))
)



cleanEx()
nameEx("insertTab")
### * insertTab

flush(stderr()); flush(stdout())

### Name: insertTab
### Title: Dynamically insert/remove a tabPanel
### Aliases: insertTab prependTab appendTab removeTab

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {

# example app for inserting/removing a tab
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      actionButton("add", "Add 'Dynamic' tab"),
      actionButton("remove", "Remove 'Foo' tab")
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
        tabPanel("Hello", "This is the hello tab"),
        tabPanel("Foo", "This is the foo tab"),
        tabPanel("Bar", "This is the bar tab")
      )
    )
  )
)
server <- function(input, output, session) {
  observeEvent(input$add, {
    insertTab(inputId = "tabs",
      tabPanel("Dynamic", "This a dynamically-added tab"),
      target = "Bar"
    )
  })
  observeEvent(input$remove, {
    removeTab(inputId = "tabs", target = "Foo")
  })
}

shinyApp(ui, server)


# example app for prepending/appending a navbarMenu
ui <- navbarPage("Navbar page", id = "tabs",
  tabPanel("Home",
    actionButton("prepend", "Prepend a navbarMenu"),
    actionButton("append", "Append a navbarMenu")
  )
)
server <- function(input, output, session) {
  observeEvent(input$prepend, {
    id <- paste0("Dropdown", input$prepend, "p")
    prependTab(inputId = "tabs",
      navbarMenu(id,
        tabPanel("Drop1", paste("Drop1 page from", id)),
        tabPanel("Drop2", paste("Drop2 page from", id)),
        "------",
        "Header",
        tabPanel("Drop3", paste("Drop3 page from", id))
      )
    )
  })
  observeEvent(input$append, {
    id <- paste0("Dropdown", input$append, "a")
    appendTab(inputId = "tabs",
      navbarMenu(id,
        tabPanel("Drop1", paste("Drop1 page from", id)),
        tabPanel("Drop2", paste("Drop2 page from", id)),
        "------",
        "Header",
        tabPanel("Drop3", paste("Drop3 page from", id))
      )
    )
  })
}

shinyApp(ui, server)

}



cleanEx()
nameEx("insertUI")
### * insertUI

flush(stderr()); flush(stdout())

### Name: insertUI
### Title: Insert and remove UI objects
### Aliases: insertUI removeUI

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
# Define UI
ui <- fluidPage(
  actionButton("add", "Add UI")
)

# Server logic
server <- function(input, output, session) {
  observeEvent(input$add, {
    insertUI(
      selector = "#add",
      where = "afterEnd",
      ui = textInput(paste0("txt", input$add),
                     "Insert some text")
    )
  })
}

# Complete app with UI and server components
shinyApp(ui, server)
}

if (interactive()) {
# Define UI
ui <- fluidPage(
  actionButton("rmv", "Remove UI"),
  textInput("txt", "This is no longer useful")
)

# Server logic
server <- function(input, output, session) {
  observeEvent(input$rmv, {
    removeUI(
      selector = "div:has(> #txt)"
    )
  })
}

# Complete app with UI and server components
shinyApp(ui, server)
}



cleanEx()
nameEx("invalidateLater")
### * invalidateLater

flush(stderr()); flush(stdout())

### Name: invalidateLater
### Title: Scheduled Invalidation
### Aliases: invalidateLater

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("n", "Number of observations", 2, 1000, 500),
  plotOutput("plot")
)

server <- function(input, output, session) {

  observe({
    # Re-execute this reactive expression after 1000 milliseconds
    invalidateLater(1000, session)

    # Do something each time this is invalidated.
    # The isolate() makes this observer _not_ get invalidated and re-executed
    # when input$n changes.
    print(paste("The value of input$n is", isolate(input$n)))
  })

  # Generate a new histogram at timed intervals, but not when
  # input$n changes.
  output$plot <- renderPlot({
    # Re-execute this reactive expression after 2000 milliseconds
    invalidateLater(2000)
    hist(rnorm(isolate(input$n)))
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("isolate")
### * isolate

flush(stderr()); flush(stdout())

### Name: isolate
### Title: Create a non-reactive scope for an expression
### Aliases: isolate

### ** Examples

## Not run: 
##D observe({
##D   input$saveButton  # Do take a dependency on input$saveButton
##D 
##D   # isolate a simple expression
##D   data <- get(isolate(input$dataset))  # No dependency on input$dataset
##D   writeToDatabase(data)
##D })
##D 
##D observe({
##D   input$saveButton  # Do take a dependency on input$saveButton
##D 
##D   # isolate a whole block
##D   data <- isolate({
##D     a <- input$valueA   # No dependency on input$valueA or input$valueB
##D     b <- input$valueB
##D     c(a=a, b=b)
##D   })
##D   writeToDatabase(data)
##D })
##D 
##D observe({
##D   x <- 1
##D   # x outside of isolate() is affected
##D   isolate(x <- 2)
##D   print(x) # 2
##D 
##D   y <- 1
##D   # Use local() to avoid affecting calling environment
##D   isolate(local(y <- 2))
##D   print(y) # 1
##D })
##D 
## End(Not run)

# Can also use isolate to call reactive expressions from the R console
values <- reactiveValues(A=1)
fun <- reactive({ as.character(values$A) })
isolate(fun())
# "1"

# isolate also works if the reactive expression accesses values from the
# input object, like input$x



cleanEx()
nameEx("makeReactiveBinding")
### * makeReactiveBinding

flush(stderr()); flush(stdout())

### Name: makeReactiveBinding
### Title: Make a reactive variable
### Aliases: makeReactiveBinding
### Keywords: internal

### ** Examples

reactiveConsole(TRUE)

a <- 10
makeReactiveBinding("a")

b <- reactive(a * -1)
observe(print(b()))

a <- 20
a <- 30

reactiveConsole(FALSE)



cleanEx()
nameEx("markdown")
### * markdown

flush(stderr()); flush(stdout())

### Name: markdown
### Title: Insert inline Markdown
### Aliases: markdown

### ** Examples

ui <- fluidPage(
  markdown("
    # Markdown Example

    This is a markdown paragraph, and will be contained within a `<p>` tag
    in the UI.

    The following is an unordered list, which will be represented in the UI as
    a `<ul>` with `<li>` children:

    * a bullet
    * another

    [Links](https://developer.mozilla.org/en-US/docs/Web/HTML/Element/a) work;
    so does *emphasis*.

    To see more of what's possible, check out [commonmark.org/help](https://commonmark.org/help).
    ")
)



cleanEx()
nameEx("modalDialog")
### * modalDialog

flush(stderr()); flush(stdout())

### Name: modalDialog
### Title: Create a modal dialog UI
### Aliases: modalDialog modalButton

### ** Examples

if (interactive()) {
# Display an important message that can be dismissed only by clicking the
# dismiss button.
shinyApp(
  ui = basicPage(
    actionButton("show", "Show modal dialog")
  ),
  server = function(input, output) {
    observeEvent(input$show, {
      showModal(modalDialog(
        title = "Important message",
        "This is an important message!"
      ))
    })
  }
)


# Display a message that can be dismissed by clicking outside the modal dialog,
# or by pressing Esc.
shinyApp(
  ui = basicPage(
    actionButton("show", "Show modal dialog")
  ),
  server = function(input, output) {
    observeEvent(input$show, {
      showModal(modalDialog(
        title = "Somewhat important message",
        "This is a somewhat important message.",
        easyClose = TRUE,
        footer = NULL
      ))
    })
  }
)


# Display a modal that requires valid input before continuing.
shinyApp(
  ui = basicPage(
    actionButton("show", "Show modal dialog"),
    verbatimTextOutput("dataInfo")
  ),

  server = function(input, output) {
    # reactiveValues object for storing current data set.
    vals <- reactiveValues(data = NULL)

    # Return the UI for a modal dialog with data selection input. If 'failed' is
    # TRUE, then display a message that the previous value was invalid.
    dataModal <- function(failed = FALSE) {
      modalDialog(
        textInput("dataset", "Choose data set",
          placeholder = 'Try "mtcars" or "abc"'
        ),
        span('(Try the name of a valid data object like "mtcars", ',
             'then a name of a non-existent object like "abc")'),
        if (failed)
          div(tags$b("Invalid name of data object", style = "color: red;")),

        footer = tagList(
          modalButton("Cancel"),
          actionButton("ok", "OK")
        )
      )
    }

    # Show modal when button is clicked.
    observeEvent(input$show, {
      showModal(dataModal())
    })

    # When OK button is pressed, attempt to load the data set. If successful,
    # remove the modal. If not show another modal, but this time with a failure
    # message.
    observeEvent(input$ok, {
      # Check that data object exists and is data frame.
      if (!is.null(input$dataset) && nzchar(input$dataset) &&
          exists(input$dataset) && is.data.frame(get(input$dataset))) {
        vals$data <- get(input$dataset)
        removeModal()
      } else {
        showModal(dataModal(failed = TRUE))
      }
    })

    # Display information about selected data
    output$dataInfo <- renderPrint({
      if (is.null(vals$data))
        "No data selected"
      else
        summary(vals$data)
    })
  }
)
}



cleanEx()
nameEx("moduleServer")
### * moduleServer

flush(stderr()); flush(stdout())

### Name: moduleServer
### Title: Shiny modules
### Aliases: moduleServer

### ** Examples

# Define the UI for a module
counterUI <- function(id, label = "Counter") {
  ns <- NS(id)
  tagList(
    actionButton(ns("button"), label = label),
    verbatimTextOutput(ns("out"))
  )
}

# Define the server logic for a module
counterServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      count <- reactiveVal(0)
      observeEvent(input$button, {
        count(count() + 1)
      })
      output$out <- renderText({
        count()
      })
      count
    }
  )
}

# Use the module in an app
ui <- fluidPage(
  counterUI("counter1", "Counter #1"),
  counterUI("counter2", "Counter #2")
)
server <- function(input, output, session) {
  counterServer("counter1")
  counterServer("counter2")
}
if (interactive()) {
  shinyApp(ui, server)
}



# If you want to pass extra parameters to the module's server logic, you can
# add them to your function. In this case `prefix` is text that will be
# printed before the count.
counterServer2 <- function(id, prefix = NULL) {
  moduleServer(
    id,
    function(input, output, session) {
      count <- reactiveVal(0)
      observeEvent(input$button, {
        count(count() + 1)
      })
      output$out <- renderText({
        paste0(prefix, count())
      })
      count
    }
  )
}

ui <- fluidPage(
  counterUI("counter", "Counter"),
)
server <- function(input, output, session) {
  counterServer2("counter", "The current count is: ")
}
if (interactive()) {
  shinyApp(ui, server)
}




cleanEx()
nameEx("navbarPage")
### * navbarPage

flush(stderr()); flush(stdout())

### Name: navbarPage
### Title: Create a page with a top level navigation bar
### Aliases: navbarPage navbarMenu

### ** Examples

navbarPage("App Title",
  tabPanel("Plot"),
  tabPanel("Summary"),
  tabPanel("Table")
)

navbarPage("App Title",
  tabPanel("Plot"),
  navbarMenu("More",
    tabPanel("Summary"),
    "----",
    "Section header",
    tabPanel("Table")
  )
)



cleanEx()
nameEx("navlistPanel")
### * navlistPanel

flush(stderr()); flush(stdout())

### Name: navlistPanel
### Title: Create a navigation list panel
### Aliases: navlistPanel

### ** Examples

fluidPage(

  titlePanel("Application Title"),

  navlistPanel(
    "Header",
    tabPanel("First"),
    tabPanel("Second"),
    tabPanel("Third")
  )
)



cleanEx()
nameEx("numericInput")
### * numericInput

flush(stderr()); flush(stdout())

### Name: numericInput
### Title: Create a numeric input control
### Aliases: numericInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  numericInput("obs", "Observations:", 10, min = 1, max = 100),
  verbatimTextOutput("value")
)
server <- function(input, output) {
  output$value <- renderText({ input$obs })
}
shinyApp(ui, server)
}




cleanEx()
nameEx("observe")
### * observe

flush(stderr()); flush(stdout())

### Name: observe
### Title: Create a reactive observer
### Aliases: observe

### ** Examples

values <- reactiveValues(A=1)

obsB <- observe({
  print(values$A + 1)
})

# To store expressions for later conversion to observe, use rlang::quo()
myquo <- rlang::quo({ print(values$A + 3) })
obsC <- rlang::inject(observe(!!myquo))

# (Legacy) Can use quoted expressions
obsD <- observe(quote({ print(values$A + 2) }), quoted = TRUE)

# In a normal Shiny app, the web client will trigger flush events. If you
# are at the console, you can force a flush with flushReact()
shiny:::flushReact()



cleanEx()
nameEx("observeEvent")
### * observeEvent

flush(stderr()); flush(stdout())

### Name: observeEvent
### Title: Event handler
### Aliases: observeEvent eventReactive

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

  ## App 1: Sample usage
  shinyApp(
    ui = fluidPage(
      column(4,
        numericInput("x", "Value", 5),
        br(),
        actionButton("button", "Show")
      ),
      column(8, tableOutput("table"))
    ),
    server = function(input, output) {
      # Take an action every time button is pressed;
      # here, we just print a message to the console
      observeEvent(input$button, {
        cat("Showing", input$x, "rows\n")
      })
      # The observeEvent() above is equivalent to:
      # observe({
      #    cat("Showing", input$x, "rows\n")
      #   }) %>%
      #   bindEvent(input$button)

      # Take a reactive dependency on input$button, but
      # not on any of the stuff inside the function
      df <- eventReactive(input$button, {
        head(cars, input$x)
      })
      output$table <- renderTable({
        df()
      })
    }
  )

  ## App 2: Using `once`
  shinyApp(
    ui = basicPage( actionButton("go", "Go")),
    server = function(input, output, session) {
      observeEvent(input$go, {
        print(paste("This will only be printed once; all",
              "subsequent button clicks won't do anything"))
      }, once = TRUE)
      # The observeEvent() above is equivalent to:
      # observe({
      #   print(paste("This will only be printed once; all",
      #         "subsequent button clicks won't do anything"))
      #   }) %>%
      #   bindEvent(input$go, once = TRUE)
    }
  )

  ## App 3: Using `ignoreInit` and `once`
  shinyApp(
    ui = basicPage(actionButton("go", "Go")),
    server = function(input, output, session) {
      observeEvent(input$go, {
        insertUI("#go", "afterEnd",
                 actionButton("dynamic", "click to remove"))

        # set up an observer that depends on the dynamic
        # input, so that it doesn't run when the input is
        # created, and only runs once after that (since
        # the side effect is remove the input from the DOM)
        observeEvent(input$dynamic, {
          removeUI("#dynamic")
        }, ignoreInit = TRUE, once = TRUE)
      })
    }
  )
}



cleanEx()
nameEx("onBookmark")
### * onBookmark

flush(stderr()); flush(stdout())

### Name: onBookmark
### Title: Add callbacks for Shiny session bookmarking events
### Aliases: onBookmark onBookmarked onRestore onRestored

### ** Examples

## Only run these examples in interactive sessions
if (interactive()) {

# Basic use of onBookmark and onRestore: This app saves the time in its
# arbitrary values, and restores that time when the app is restored.
ui <- function(req) {
  fluidPage(
    textInput("txt", "Input text"),
    bookmarkButton()
  )
}
server <- function(input, output) {
  onBookmark(function(state) {
    savedTime <- as.character(Sys.time())
    cat("Last saved at", savedTime, "\n")
    # state is a mutable reference object, and we can add arbitrary values to
    # it.
    state$values$time <- savedTime
  })

  onRestore(function(state) {
    cat("Restoring from state bookmarked at", state$values$time, "\n")
  })
}
enableBookmarking("url")
shinyApp(ui, server)



ui <- function(req) {
  fluidPage(
    textInput("txt", "Input text"),
    bookmarkButton()
  )
}
server <- function(input, output, session) {
  lastUpdateTime <- NULL

  observeEvent(input$txt, {
    updateTextInput(session, "txt",
      label = paste0("Input text (Changed ", as.character(Sys.time()), ")")
    )
  })

  onBookmark(function(state) {
    # Save content to a file
    messageFile <- file.path(state$dir, "message.txt")
    cat(as.character(Sys.time()), file = messageFile)
  })

  onRestored(function(state) {
    # Read the file
    messageFile <- file.path(state$dir, "message.txt")
    timeText <- readChar(messageFile, 1000)

    # updateTextInput must be called in onRestored, as opposed to onRestore,
    # because onRestored happens after the client browser is ready.
    updateTextInput(session, "txt",
      label = paste0("Input text (Changed ", timeText, ")")
    )
  })
}
# "server" bookmarking is needed for writing to disk.
enableBookmarking("server")
shinyApp(ui, server)


# This app has a module, and both the module and the main app code have
# onBookmark and onRestore functions which write and read state$values$hash. The
# module's version of state$values$hash does not conflict with the app's version
# of state$values$hash.
#
# A basic module that captializes text.
capitalizerUI <- function(id) {
  ns <- NS(id)
  wellPanel(
    h4("Text captializer module"),
    textInput(ns("text"), "Enter text:"),
    verbatimTextOutput(ns("out"))
  )
}
capitalizerServer <- function(input, output, session) {
  output$out <- renderText({
    toupper(input$text)
  })
  onBookmark(function(state) {
    state$values$hash <- rlang::hash(input$text)
  })
  onRestore(function(state) {
    if (identical(rlang::hash(input$text), state$values$hash)) {
      message("Module's input text matches hash ", state$values$hash)
    } else {
      message("Module's input text does not match hash ", state$values$hash)
    }
  })
}
# Main app code
ui <- function(request) {
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        capitalizerUI("tc"),
        textInput("text", "Enter text (not in module):"),
        bookmarkButton()
      ),
      mainPanel()
    )
  )
}
server <- function(input, output, session) {
  callModule(capitalizerServer, "tc")
  onBookmark(function(state) {
    state$values$hash <- rlang::hash(input$text)
  })
  onRestore(function(state) {
    if (identical(rlang::hash(input$text), state$values$hash)) {
      message("App's input text matches hash ", state$values$hash)
    } else {
      message("App's input text does not match hash ", state$values$hash)
    }
  })
}
enableBookmarking(store = "url")
shinyApp(ui, server)
}



cleanEx()
nameEx("onFlush")
### * onFlush

flush(stderr()); flush(stdout())

### Name: onFlush
### Title: Add callbacks for Shiny session events
### Aliases: onFlush onFlushed onSessionEnded onUnhandledError

### ** Examples

## Don't show: 
if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)({ # examplesIf
## End(Don't show)
library(shiny)

ui <- fixedPage(
  markdown(c(
    "Set the number to 8 or higher to cause an error",
    "in the `renderText()` output."
  )),
  sliderInput("number", "Number", 0, 10, 4),
  textOutput("text"),
  hr(),
  markdown(c(
    "Click the button below to crash the app with an unhandled error",
    "in an `observe()` block."
  )),
  actionButton("crash", "Crash the app!")
)

log_event <- function(level, ...) {
  ts <- strftime(Sys.time(), " [%F %T] ")
  message(level, ts, ...)
}

server <- function(input, output, session) {
  log_event("INFO", "Session started")

  onUnhandledError(function(err) {
    # log the unhandled error
    level <- if (inherits(err, "shiny.error.fatal")) "FATAL" else "ERROR"
    log_event(level, conditionMessage(err))
  })

  onStop(function() {
    log_event("INFO", "Session ended")
  })

  observeEvent(input$crash, stop("Oops, an unhandled error happened!"))

  output$text <- renderText({
    if (input$number > 7) {
      stop("that's too high!")
    }
    sprintf("You picked number %d.", input$number)
  })
}

shinyApp(ui, server)
## Don't show: 
}) # examplesIf
## End(Don't show)



cleanEx()
nameEx("onStop")
### * onStop

flush(stderr()); flush(stdout())

### Name: onStop
### Title: Run code after an application or session ends
### Aliases: onStop

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  # Open this application in multiple browsers, then close the browsers.
  shinyApp(
    ui = basicPage("onStop demo"),

    server = function(input, output, session) {
      onStop(function() cat("Session stopped\n"))
    },

    onStart = function() {
      cat("Doing application setup\n")

      onStop(function() {
        cat("Doing application cleanup\n")
      })
    }
  )
}
# In the example above, onStop() is called inside of onStart(). This is
# the pattern that should be used when creating a shinyApp() object from
# a function, or at the console. If instead you are writing an app.R which
# will be invoked with runApp(), you can do it that way, or put the onStop()
# before the shinyApp() call, as shown below.

## Not run: 
##D # ==== app.R ====
##D cat("Doing application setup\n")
##D onStop(function() {
##D   cat("Doing application cleanup\n")
##D })
##D 
##D shinyApp(
##D   ui = basicPage("onStop demo"),
##D 
##D   server = function(input, output, session) {
##D     onStop(function() cat("Session stopped\n"))
##D   }
##D )
##D # ==== end app.R ====
##D 
##D 
##D # Similarly, if you have a global.R, you can call onStop() from there.
##D # ==== global.R ====
##D cat("Doing application setup\n")
##D onStop(function() {
##D   cat("Doing application cleanup\n")
##D })
##D # ==== end global.R ====
## End(Not run)



cleanEx()
nameEx("outputOptions")
### * outputOptions

flush(stderr()); flush(stdout())

### Name: outputOptions
### Title: Set options for an output object.
### Aliases: outputOptions

### ** Examples

## Not run: 
##D # Get the list of options for all observers within output
##D outputOptions(output)
##D 
##D # Disable suspend for output$myplot
##D outputOptions(output, "myplot", suspendWhenHidden = FALSE)
##D 
##D # Change priority for output$myplot
##D outputOptions(output, "myplot", priority = 10)
##D 
##D # Get the list of options for output$myplot
##D outputOptions(output, "myplot")
## End(Not run)




cleanEx()
nameEx("parseQueryString")
### * parseQueryString

flush(stderr()); flush(stdout())

### Name: parseQueryString
### Title: Parse a GET query string from a URL
### Aliases: parseQueryString

### ** Examples

parseQueryString("?foo=1&bar=b%20a%20r")

## Not run: 
##D # Example of usage within a Shiny app
##D function(input, output, session) {
##D 
##D   output$queryText <- renderText({
##D     query <- parseQueryString(session$clientData$url_search)
##D 
##D     # Ways of accessing the values
##D     if (as.numeric(query$foo) == 1) {
##D       # Do something
##D     }
##D     if (query[["bar"]] == "targetstring") {
##D       # Do something else
##D     }
##D 
##D     # Return a string with key-value pairs
##D     paste(names(query), query, sep = "=", collapse=", ")
##D   })
##D }
## End(Not run)




cleanEx()
nameEx("passwordInput")
### * passwordInput

flush(stderr()); flush(stdout())

### Name: passwordInput
### Title: Create a password input control
### Aliases: passwordInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  passwordInput("password", "Password:"),
  actionButton("go", "Go"),
  verbatimTextOutput("value")
)
server <- function(input, output) {
  output$value <- renderText({
    req(input$go)
    isolate(input$password)
  })
}
shinyApp(ui, server)
}



cleanEx()
nameEx("plotOutput")
### * plotOutput

flush(stderr()); flush(stdout())

### Name: plotOutput
### Title: Create an plot or image output element
### Aliases: plotOutput imageOutput

### ** Examples

# Only run these examples in interactive R sessions
if (interactive()) {

# A basic shiny app with a plotOutput
shinyApp(
  ui = fluidPage(
    sidebarLayout(
      sidebarPanel(
        actionButton("newplot", "New plot")
      ),
      mainPanel(
        plotOutput("plot")
      )
    )
  ),
  server = function(input, output) {
    output$plot <- renderPlot({
      input$newplot
      # Add a little noise to the cars data
      cars2 <- cars + rnorm(nrow(cars))
      plot(cars2)
    })
  }
)


# A demonstration of clicking, hovering, and brushing
shinyApp(
  ui = basicPage(
    fluidRow(
      column(width = 4,
        plotOutput("plot", height=300,
          click = "plot_click",  # Equiv, to click=clickOpts(id="plot_click")
          hover = hoverOpts(id = "plot_hover", delayType = "throttle"),
          brush = brushOpts(id = "plot_brush")
        ),
        h4("Clicked points"),
        tableOutput("plot_clickedpoints"),
        h4("Brushed points"),
        tableOutput("plot_brushedpoints")
      ),
      column(width = 4,
        verbatimTextOutput("plot_clickinfo"),
        verbatimTextOutput("plot_hoverinfo")
      ),
      column(width = 4,
        wellPanel(actionButton("newplot", "New plot")),
        verbatimTextOutput("plot_brushinfo")
      )
    )
  ),
  server = function(input, output, session) {
    data <- reactive({
      input$newplot
      # Add a little noise to the cars data so the points move
      cars + rnorm(nrow(cars))
    })
    output$plot <- renderPlot({
      d <- data()
      plot(d$speed, d$dist)
    })
    output$plot_clickinfo <- renderPrint({
      cat("Click:\n")
      str(input$plot_click)
    })
    output$plot_hoverinfo <- renderPrint({
      cat("Hover (throttled):\n")
      str(input$plot_hover)
    })
    output$plot_brushinfo <- renderPrint({
      cat("Brush (debounced):\n")
      str(input$plot_brush)
    })
    output$plot_clickedpoints <- renderTable({
      # For base graphics, we need to specify columns, though for ggplot2,
      # it's usually not necessary.
      res <- nearPoints(data(), input$plot_click, "speed", "dist")
      if (nrow(res) == 0)
        return()
      res
    })
    output$plot_brushedpoints <- renderTable({
      res <- brushedPoints(data(), input$plot_brush, "speed", "dist")
      if (nrow(res) == 0)
        return()
      res
    })
  }
)


# Demo of clicking, hovering, brushing with imageOutput
# Note that coordinates are in pixels
shinyApp(
  ui = basicPage(
    fluidRow(
      column(width = 4,
        imageOutput("image", height=300,
          click = "image_click",
          hover = hoverOpts(
            id = "image_hover",
            delay = 500,
            delayType = "throttle"
          ),
          brush = brushOpts(id = "image_brush")
        )
      ),
      column(width = 4,
        verbatimTextOutput("image_clickinfo"),
        verbatimTextOutput("image_hoverinfo")
      ),
      column(width = 4,
        wellPanel(actionButton("newimage", "New image")),
        verbatimTextOutput("image_brushinfo")
      )
    )
  ),
  server = function(input, output, session) {
    output$image <- renderImage({
      input$newimage

      # Get width and height of image output
      width  <- session$clientData$output_image_width
      height <- session$clientData$output_image_height

      # Write to a temporary PNG file
      outfile <- tempfile(fileext = ".png")

      png(outfile, width=width, height=height)
      plot(rnorm(200), rnorm(200))
      dev.off()

      # Return a list containing information about the image
      list(
        src = outfile,
        contentType = "image/png",
        width = width,
        height = height,
        alt = "This is alternate text"
      )
    })
    output$image_clickinfo <- renderPrint({
      cat("Click:\n")
      str(input$image_click)
    })
    output$image_hoverinfo <- renderPrint({
      cat("Hover (throttled):\n")
      str(input$image_hover)
    })
    output$image_brushinfo <- renderPrint({
      cat("Brush (debounced):\n")
      str(input$image_brush)
    })
  }
)

}



cleanEx()
nameEx("radioButtons")
### * radioButtons

flush(stderr()); flush(stdout())

### Name: radioButtons
### Title: Create radio buttons
### Aliases: radioButtons

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  radioButtons("dist", "Distribution type:",
               c("Normal" = "norm",
                 "Uniform" = "unif",
                 "Log-normal" = "lnorm",
                 "Exponential" = "exp")),
  plotOutput("distPlot")
)

server <- function(input, output) {
  output$distPlot <- renderPlot({
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)

    hist(dist(500))
  })
}

shinyApp(ui, server)

ui <- fluidPage(
  radioButtons("rb", "Choose one:",
               choiceNames = list(
                 icon("calendar"),
                 HTML("<p style='color:red;'>Red Text</p>"),
                 "Normal text"
               ),
               choiceValues = list(
                 "icon", "html", "text"
               )),
  textOutput("txt")
)

server <- function(input, output) {
  output$txt <- renderText({
    paste("You chose", input$rb)
  })
}

shinyApp(ui, server)
}




cleanEx()
nameEx("reactive")
### * reactive

flush(stderr()); flush(stdout())

### Name: reactive
### Title: Create a reactive expression
### Aliases: reactive is.reactive

### ** Examples

library(rlang)
values <- reactiveValues(A=1)

reactiveB <- reactive({
  values$A + 1
})
# View the values from the R console with isolate()
isolate(reactiveB())
# 2

# To store expressions for later conversion to reactive, use quote()
myquo <- rlang::quo(values$A + 2)
# Unexpected value! Sending a quosure directly will not work as expected.
reactiveC <- reactive(myquo)
# We'd hope for `3`, but instead we get the quosure that was supplied.
isolate(reactiveC())

# Instead, the quosure should be `rlang::inject()`ed
reactiveD <- rlang::inject(reactive(!!myquo))
isolate(reactiveD())
# 3

# (Legacy) Can use quoted expressions
expr <- quote({ values$A + 3 })
reactiveE <- reactive(expr, quoted = TRUE)
isolate(reactiveE())
# 4




cleanEx()
nameEx("reactiveConsole")
### * reactiveConsole

flush(stderr()); flush(stdout())

### Name: reactiveConsole
### Title: Activate reactivity in the console
### Aliases: reactiveConsole
### Keywords: internal

### ** Examples

reactiveConsole(TRUE)
x <- reactiveVal(10)
y <- observe({
  message("The value of x is ", x())
})
x(20)
x(30)
reactiveConsole(FALSE)



cleanEx()
nameEx("reactiveFileReader")
### * reactiveFileReader

flush(stderr()); flush(stdout())

### Name: reactiveFileReader
### Title: Reactive file reader
### Aliases: reactiveFileReader

### ** Examples

## Not run: 
##D # Per-session reactive file reader
##D function(input, output, session) {
##D   fileData <- reactiveFileReader(1000, session, 'data.csv', read.csv)
##D 
##D   output$data <- renderTable({
##D     fileData()
##D   })
##D }
##D 
##D # Cross-session reactive file reader. In this example, all sessions share
##D # the same reader, so read.csv only gets executed once no matter how many
##D # user sessions are connected.
##D fileData <- reactiveFileReader(1000, NULL, 'data.csv', read.csv)
##D function(input, output, session) {
##D   output$data <- renderTable({
##D     fileData()
##D   })
##D }
## End(Not run)



cleanEx()
nameEx("reactivePoll")
### * reactivePoll

flush(stderr()); flush(stdout())

### Name: reactivePoll
### Title: Reactive polling
### Aliases: reactivePoll

### ** Examples

function(input, output, session) {

  data <- reactivePoll(1000, session,
    # This function returns the time that log_file was last modified
    checkFunc = function() {
      if (file.exists(log_file))
        file.info(log_file)$mtime[1]
      else
        ""
    },
    # This function returns the content of log_file
    valueFunc = function() {
      read.csv(log_file)
    }
  )

  output$dataTable <- renderTable({
    data()
  })
}



cleanEx()
nameEx("reactiveTimer")
### * reactiveTimer

flush(stderr()); flush(stdout())

### Name: reactiveTimer
### Title: Timer
### Aliases: reactiveTimer

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("n", "Number of observations", 2, 1000, 500),
  plotOutput("plot")
)

server <- function(input, output) {

  # Anything that calls autoInvalidate will automatically invalidate
  # every 2 seconds.
  autoInvalidate <- reactiveTimer(2000)

  observe({
    # Invalidate and re-execute this reactive expression every time the
    # timer fires.
    autoInvalidate()

    # Do something each time this is invalidated.
    # The isolate() makes this observer _not_ get invalidated and re-executed
    # when input$n changes.
    print(paste("The value of input$n is", isolate(input$n)))
  })

  # Generate a new histogram each time the timer fires, but not when
  # input$n changes.
  output$plot <- renderPlot({
    autoInvalidate()
    hist(rnorm(isolate(input$n)))
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("reactiveVal")
### * reactiveVal

flush(stderr()); flush(stdout())

### Name: reactiveVal
### Title: Create a (single) reactive value
### Aliases: reactiveVal

### ** Examples


## Not run: 
##D 
##D # Create the object by calling reactiveVal
##D r <- reactiveVal()
##D 
##D # Set the value by calling with an argument
##D r(10)
##D 
##D # Read the value by calling without arguments
##D r()
##D 
## End(Not run)

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  actionButton("minus", "-1"),
  actionButton("plus", "+1"),
  br(),
  textOutput("value")
)

# The comments below show the equivalent logic using reactiveValues()
server <- function(input, output, session) {
  value <- reactiveVal(0)       # rv <- reactiveValues(value = 0)

  observeEvent(input$minus, {
    newValue <- value() - 1     # newValue <- rv$value - 1
    value(newValue)             # rv$value <- newValue
  })

  observeEvent(input$plus, {
    newValue <- value() + 1     # newValue <- rv$value + 1
    value(newValue)             # rv$value <- newValue
  })

  output$value <- renderText({
    value()                     # rv$value
  })
}

shinyApp(ui, server)

}




cleanEx()
nameEx("reactiveValues")
### * reactiveValues

flush(stderr()); flush(stdout())

### Name: reactiveValues
### Title: Create an object for storing reactive values
### Aliases: reactiveValues

### ** Examples

# Create the object with no values
values <- reactiveValues()

# Assign values to 'a' and 'b'
values$a <- 3
values[['b']] <- 4

## Not run: 
##D # From within a reactive context, you can access values with:
##D values$a
##D values[['a']]
## End(Not run)

# If not in a reactive context (e.g., at the console), you can use isolate()
# to retrieve the value:
isolate(values$a)
isolate(values[['a']])

# Set values upon creation
values <- reactiveValues(a = 1, b = 2)
isolate(values$a)




cleanEx()
nameEx("reactiveValuesToList")
### * reactiveValuesToList

flush(stderr()); flush(stdout())

### Name: reactiveValuesToList
### Title: Convert a reactivevalues object to a list
### Aliases: reactiveValuesToList

### ** Examples

values <- reactiveValues(a = 1)
## Not run: 
##D reactiveValuesToList(values)
## End(Not run)

# To get the objects without taking dependencies on them, use isolate().
# isolate() can also be used when calling from outside a reactive context (e.g.
# at the console)
isolate(reactiveValuesToList(values))



cleanEx()
nameEx("registerInputHandler")
### * registerInputHandler

flush(stderr()); flush(stdout())

### Name: registerInputHandler
### Title: Register an Input Handler
### Aliases: registerInputHandler

### ** Examples

## Not run: 
##D # Register an input handler which rounds a input number to the nearest integer
##D # In a package, this should be called from the .onLoad function.
##D registerInputHandler("mypackage.validint", function(x, shinysession, name) {
##D   if (is.null(x)) return(NA)
##D   round(x)
##D })
##D 
##D ## On the Javascript side, the associated input binding must have a corresponding getType method:
##D # getType: function(el) {
##D #   return "mypackage.validint";
##D # }
##D 
## End(Not run)



cleanEx()
nameEx("renderCachedPlot")
### * renderCachedPlot

flush(stderr()); flush(stdout())

### Name: renderCachedPlot
### Title: Plot output with cached images
### Aliases: renderCachedPlot

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

# A basic example that uses the default app-scoped memory cache.
# The cache will be shared among all simultaneous users of the application.
shinyApp(
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        sliderInput("n", "Number of points", 4, 32, value = 8, step = 4)
      ),
      mainPanel(plotOutput("plot"))
    )
  ),
  function(input, output, session) {
    output$plot <- renderCachedPlot({
        Sys.sleep(2)  # Add an artificial delay
        seqn <- seq_len(input$n)
        plot(mtcars$wt[seqn], mtcars$mpg[seqn],
             xlim = range(mtcars$wt), ylim = range(mtcars$mpg))
      },
      cacheKeyExpr = { list(input$n) }
    )
  }
)



# An example uses a data object shared across sessions. mydata() is part of
# the cache key, so when its value changes, plots that were previously
# stored in the cache will no longer be used (unless mydata() changes back
# to its previous value).
mydata <- reactiveVal(data.frame(x = rnorm(400), y = rnorm(400)))

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Number of points", 50, 400, 100, step = 50),
      actionButton("newdata", "New data")
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$newdata, {
    mydata(data.frame(x = rnorm(400), y = rnorm(400)))
  })

  output$plot <- renderCachedPlot(
    {
      Sys.sleep(2)
      d <- mydata()
      seqn <- seq_len(input$n)
      plot(d$x[seqn], d$y[seqn], xlim = range(d$x), ylim = range(d$y))
    },
    cacheKeyExpr = { list(input$n, mydata()) },
  )
}

shinyApp(ui, server)


# A basic application with two plots, where each plot in each session has
# a separate cache.
shinyApp(
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        sliderInput("n", "Number of points", 4, 32, value = 8, step = 4)
      ),
      mainPanel(
        plotOutput("plot1"),
        plotOutput("plot2")
      )
    )
  ),
  function(input, output, session) {
    output$plot1 <- renderCachedPlot({
        Sys.sleep(2)  # Add an artificial delay
        seqn <- seq_len(input$n)
        plot(mtcars$wt[seqn], mtcars$mpg[seqn],
             xlim = range(mtcars$wt), ylim = range(mtcars$mpg))
      },
      cacheKeyExpr = { list(input$n) },
      cache = cachem::cache_mem()
    )
    output$plot2 <- renderCachedPlot({
        Sys.sleep(2)  # Add an artificial delay
        seqn <- seq_len(input$n)
        plot(mtcars$wt[seqn], mtcars$mpg[seqn],
             xlim = range(mtcars$wt), ylim = range(mtcars$mpg))
      },
      cacheKeyExpr = { list(input$n) },
      cache = cachem::cache_mem()
    )
  }
)

}

## Not run: 
##D # At the top of app.R, this set the application-scoped cache to be a memory
##D # cache that is 20 MB in size, and where cached objects expire after one
##D # hour.
##D shinyOptions(cache = cachem::cache_mem(max_size = 20e6, max_age = 3600))
##D 
##D # At the top of app.R, this set the application-scoped cache to be a disk
##D # cache that can be shared among multiple concurrent R processes, and is
##D # deleted when the system reboots.
##D shinyOptions(cache = cachem::cache_disk(file.path(dirname(tempdir()), "myapp-cache")))
##D 
##D # At the top of app.R, this set the application-scoped cache to be a disk
##D # cache that can be shared among multiple concurrent R processes, and
##D # persists on disk across reboots.
##D shinyOptions(cache = cachem::cache_disk("./myapp-cache"))
##D 
##D # At the top of the server function, this set the session-scoped cache to be
##D # a memory cache that is 5 MB in size.
##D server <- function(input, output, session) {
##D   shinyOptions(cache = cachem::cache_mem(max_size = 5e6))
##D 
##D   output$plot <- renderCachedPlot(
##D     ...,
##D     cache = "session"
##D   )
##D }
##D 
## End(Not run)



cleanEx()
nameEx("renderDataTable")
### * renderDataTable

flush(stderr()); flush(stdout())

### Name: dataTableOutput
### Title: Table output with the JavaScript DataTables library
### Aliases: dataTableOutput renderDataTable
### Keywords: internal

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  # pass a callback function to DataTables using I()
  shinyApp(
    ui = fluidPage(
      fluidRow(
        column(12,
          dataTableOutput('table')
        )
      )
    ),
    server = function(input, output) {
      output$table <- renderDataTable(iris,
        options = list(
          pageLength = 5,
          initComplete = I("function(settings, json) {alert('Done.');}")
        )
      )
    }
  )
}



cleanEx()
nameEx("renderImage")
### * renderImage

flush(stderr()); flush(stdout())

### Name: renderImage
### Title: Image file output
### Aliases: renderImage

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

ui <- fluidPage(
  sliderInput("n", "Number of observations", 2, 1000, 500),
  plotOutput("plot1"),
  plotOutput("plot2"),
  plotOutput("plot3")
)

server <- function(input, output, session) {

  # A plot of fixed size
  output$plot1 <- renderImage({
    # A temp file to save the output. It will be deleted after renderImage
    # sends it, because deleteFile=TRUE.
    outfile <- tempfile(fileext='.png')

    # Generate a png
    png(outfile, width=400, height=400)
    hist(rnorm(input$n))
    dev.off()

    # Return a list
    list(src = outfile,
         alt = "This is alternate text")
  }, deleteFile = TRUE)

  # A dynamically-sized plot
  output$plot2 <- renderImage({
    # Read plot2's width and height. These are reactive values, so this
    # expression will re-run whenever these values change.
    width  <- session$clientData$output_plot2_width
    height <- session$clientData$output_plot2_height

    # A temp file to save the output.
    outfile <- tempfile(fileext='.png')

    png(outfile, width=width, height=height)
    hist(rnorm(input$n))
    dev.off()

    # Return a list containing the filename
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")
  }, deleteFile = TRUE)

  # Send a pre-rendered image, and don't delete the image after sending it
  # NOTE: For this example to work, it would require files in a subdirectory
  # named images/
  output$plot3 <- renderImage({
    # When input$n is 1, filename is ./images/image1.jpeg
    filename <- normalizePath(file.path('./images',
                              paste('image', input$n, '.jpeg', sep='')))

    # Return a list containing the filename
    list(src = filename)
  }, deleteFile = FALSE)
}

shinyApp(ui, server)
}



cleanEx()
nameEx("renderPrint")
### * renderPrint

flush(stderr()); flush(stdout())

### Name: renderPrint
### Title: Text Output
### Aliases: renderPrint renderText

### ** Examples

isolate({

# renderPrint captures any print output, converts it to a string, and
# returns it
visFun <- renderPrint({ "foo" })
visFun()
# '[1] "foo"'

invisFun <- renderPrint({ invisible("foo") })
invisFun()
# ''

multiprintFun <- renderPrint({
  print("foo");
  "bar"
})
multiprintFun()
# '[1] "foo"\n[1] "bar"'

nullFun <- renderPrint({ NULL })
nullFun()
# 'NULL'

invisNullFun <- renderPrint({ invisible(NULL) })
invisNullFun()
# ''

vecFun <- renderPrint({ 1:5 })
vecFun()
# '[1] 1 2 3 4 5'


# Contrast with renderText, which takes the value returned from the function
# and uses cat() to convert it to a string
visFun <- renderText({ "foo" })
visFun()
# 'foo'

invisFun <- renderText({ invisible("foo") })
invisFun()
# 'foo'

multiprintFun <- renderText({
  print("foo");
  "bar"
})
multiprintFun()
# 'bar'

nullFun <- renderText({ NULL })
nullFun()
# ''

invisNullFun <- renderText({ invisible(NULL) })
invisNullFun()
# ''

vecFun <- renderText({ 1:5 })
vecFun()
# '1 2 3 4 5'

})



cleanEx()
nameEx("renderTable")
### * renderTable

flush(stderr()); flush(stdout())

### Name: tableOutput
### Title: Table Output
### Aliases: tableOutput renderTable

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  # table example
  shinyApp(
    ui = fluidPage(
      fluidRow(
        column(12,
          tableOutput('table')
        )
      )
    ),
    server = function(input, output) {
      output$table <- renderTable(iris)
    }
  )
}



cleanEx()
nameEx("renderUI")
### * renderUI

flush(stderr()); flush(stdout())

### Name: renderUI
### Title: UI Output
### Aliases: renderUI

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  uiOutput("moreControls")
)

server <- function(input, output) {
  output$moreControls <- renderUI({
    tagList(
      sliderInput("n", "N", 1, 1000, 500),
      textInput("label", "Label")
    )
  })
}
shinyApp(ui, server)
}




cleanEx()
nameEx("repeatable")
### * repeatable

flush(stderr()); flush(stdout())

### Name: repeatable
### Title: Make a random number generator repeatable
### Aliases: repeatable

### ** Examples

rnormA <- repeatable(rnorm)
rnormB <- repeatable(rnorm)
rnormA(3)  # [1]  1.8285879 -0.7468041 -0.4639111
rnormA(3)  # [1]  1.8285879 -0.7468041 -0.4639111
rnormA(5)  # [1]  1.8285879 -0.7468041 -0.4639111 -1.6510126 -1.4686924
rnormB(5)  # [1] -0.7946034  0.2568374 -0.6567597  1.2451387 -0.8375699



cleanEx()
nameEx("req")
### * req

flush(stderr()); flush(stdout())

### Name: req
### Title: Check for required values
### Aliases: req

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
  ui <- fluidPage(
    textInput('data', 'Enter a dataset from the "datasets" package', 'cars'),
    p('(E.g. "cars", "mtcars", "pressure", "faithful")'), hr(),
    tableOutput('tbl')
  )

  server <- function(input, output) {
    output$tbl <- renderTable({

      ## to require that the user types something, use: `req(input$data)`
      ## but better: require that input$data is valid and leave the last
      ## valid table up
      req(exists(input$data, "package:datasets", inherits = FALSE),
          cancelOutput = TRUE)

      head(get(input$data, "package:datasets", inherits = FALSE))
    })
  }

  shinyApp(ui, server)
}



cleanEx()
nameEx("resourcePaths")
### * resourcePaths

flush(stderr()); flush(stdout())

### Name: addResourcePath
### Title: Resource Publishing
### Aliases: addResourcePath resourcePaths removeResourcePath

### ** Examples

addResourcePath('datasets', system.file('data', package='datasets'))
resourcePaths()
removeResourcePath('datasets')
resourcePaths()

# make sure all resources are removed
lapply(names(resourcePaths()), removeResourcePath)



cleanEx()
nameEx("runApp")
### * runApp

flush(stderr()); flush(stdout())

### Name: runApp
### Title: Run Shiny Application
### Aliases: runApp

### ** Examples

## Not run: 
##D # Start app in the current working directory
##D runApp()
##D 
##D # Start app in a subdirectory called myapp
##D runApp("myapp")
## End(Not run)

## Only run this example in interactive R sessions
if (interactive()) {
  options(device.ask.default = FALSE)

  # Apps can be run without a server.r and ui.r file
  runApp(list(
    ui = bootstrapPage(
      numericInput('n', 'Number of obs', 100),
      plotOutput('plot')
    ),
    server = function(input, output) {
      output$plot <- renderPlot({ hist(runif(input$n)) })
    }
  ))


  # Running a Shiny app object
  app <- shinyApp(
    ui = bootstrapPage(
      numericInput('n', 'Number of obs', 100),
      plotOutput('plot')
    ),
    server = function(input, output) {
      output$plot <- renderPlot({ hist(runif(input$n)) })
    }
  )
  runApp(app)
}



cleanEx()
nameEx("runExample")
### * runExample

flush(stderr()); flush(stdout())

### Name: runExample
### Title: Run Shiny Example Applications
### Aliases: runExample

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  # List all available examples
  runExample()

  # Run one of the examples
  runExample("01_hello")

  # Print the directory containing the code for all examples
  system.file("examples", package="shiny")
}



cleanEx()
nameEx("runGadget")
### * runGadget

flush(stderr()); flush(stdout())

### Name: runGadget
### Title: Run a gadget
### Aliases: runGadget

### ** Examples

## Not run: 
##D library(shiny)
##D 
##D ui <- fillPage(...)
##D 
##D server <- function(input, output, session) {
##D   ...
##D }
##D 
##D # Either pass ui/server as separate arguments...
##D runGadget(ui, server)
##D 
##D # ...or as a single app object
##D runGadget(shinyApp(ui, server))
## End(Not run)



cleanEx()
nameEx("runUrl")
### * runUrl

flush(stderr()); flush(stdout())

### Name: runUrl
### Title: Run a Shiny application from a URL
### Aliases: runUrl runGist runGitHub

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  runUrl('https://github.com/rstudio/shiny_example/archive/main.tar.gz')

  # Can run an app from a subdirectory in the archive
  runUrl("https://github.com/rstudio/shiny_example/archive/main.zip",
    subdir = "inst/shinyapp/")
}
## Only run this example in interactive R sessions
if (interactive()) {
  runGist(3239667)
  runGist("https://gist.github.com/jcheng5/3239667")

  # Old URL format without username
  runGist("https://gist.github.com/3239667")
}

## Only run this example in interactive R sessions
if (interactive()) {
  runGitHub("shiny_example", "rstudio")
  # or runGitHub("rstudio/shiny_example")

  # Can run an app from a subdirectory in the repo
  runGitHub("shiny_example", "rstudio", subdir = "inst/shinyapp/")
}



cleanEx()
nameEx("safeError")
### * safeError

flush(stderr()); flush(stdout())

### Name: safeError
### Title: Declare an error safe for the user to see
### Aliases: safeError

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

# uncomment the desired line to experiment with shiny.sanitize.errors
# options(shiny.sanitize.errors = TRUE)
# options(shiny.sanitize.errors = FALSE)

# Define UI
ui <- fluidPage(
  textInput('number', 'Enter your favorite number from 1 to 10', '5'),
  textOutput('normalError'),
  textOutput('safeError')
)

# Server logic
server <- function(input, output) {
  output$normalError <- renderText({
    number <- input$number
    if (number %in% 1:10) {
      return(paste('You chose', number, '!'))
    } else {
      stop(
        paste(number, 'is not a number between 1 and 10')
      )
    }
  })
  output$safeError <- renderText({
    number <- input$number
    if (number %in% 1:10) {
      return(paste('You chose', number, '!'))
    } else {
      stop(safeError(
        paste(number, 'is not a number between 1 and 10')
      ))
    }
  })
}

# Complete app with UI and server components
shinyApp(ui, server)
}



cleanEx()
nameEx("selectInput")
### * selectInput

flush(stderr()); flush(stdout())

### Name: selectInput
### Title: Create a select list input control
### Aliases: selectInput selectizeInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

# basic example
shinyApp(
  ui = fluidPage(
    selectInput("variable", "Variable:",
                c("Cylinders" = "cyl",
                  "Transmission" = "am",
                  "Gears" = "gear")),
    tableOutput("data")
  ),
  server = function(input, output) {
    output$data <- renderTable({
      mtcars[, c("mpg", input$variable), drop = FALSE]
    }, rownames = TRUE)
  }
)

# demoing group support in the `choices` arg
shinyApp(
  ui = fluidPage(
    selectInput("state", "Choose a state:",
      list(`East Coast` = list("NY", "NJ", "CT"),
           `West Coast` = list("WA", "OR", "CA"),
           `Midwest` = list("MN", "WI", "IA"))
    ),
    textOutput("result")
  ),
  server = function(input, output) {
    output$result <- renderText({
      paste("You chose", input$state)
    })
  }
)
}




cleanEx()
nameEx("shinyApp")
### * shinyApp

flush(stderr()); flush(stdout())

### Name: shinyApp
### Title: Create a Shiny app object
### Aliases: shinyApp shinyAppDir shinyAppFile

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  options(device.ask.default = FALSE)

  shinyApp(
    ui = fluidPage(
      numericInput("n", "n", 1),
      plotOutput("plot")
    ),
    server = function(input, output) {
      output$plot <- renderPlot( plot(head(cars, input$n)) )
    }
  )

  shinyAppDir(system.file("examples/01_hello", package="shiny"))


  # The object can be passed to runApp()
  app <- shinyApp(
    ui = fluidPage(
      numericInput("n", "n", 1),
      plotOutput("plot")
    ),
    server = function(input, output) {
      output$plot <- renderPlot( plot(head(cars, input$n)) )
    }
  )

  runApp(app)
}



cleanEx()
nameEx("shinyServer")
### * shinyServer

flush(stderr()); flush(stdout())

### Name: shinyServer
### Title: Define Server Functionality
### Aliases: shinyServer
### Keywords: internal

### ** Examples

## Not run: 
##D # A very simple Shiny app that takes a message from the user
##D # and outputs an uppercase version of it.
##D shinyServer(function(input, output, session) {
##D   output$uppercase <- renderText({
##D     toupper(input$message)
##D   })
##D })
##D 
##D 
##D # It is also possible for a server.R file to simply return the function,
##D # without calling shinyServer().
##D # For example, the server.R file could contain just the following:
##D function(input, output, session) {
##D   output$uppercase <- renderText({
##D     toupper(input$message)
##D   })
##D }
## End(Not run)



cleanEx()
nameEx("showNotification")
### * showNotification

flush(stderr()); flush(stdout())

### Name: showNotification
### Title: Show or remove a notification
### Aliases: showNotification removeNotification

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
# Show a message when button is clicked
shinyApp(
  ui = fluidPage(
    actionButton("show", "Show")
  ),
  server = function(input, output) {
    observeEvent(input$show, {
      showNotification("Message text",
        action = a(href = "javascript:location.reload();", "Reload page")
      )
    })
  }
)

# App with show and remove buttons
shinyApp(
  ui = fluidPage(
    actionButton("show", "Show"),
    actionButton("remove", "Remove")
  ),
  server = function(input, output) {
    # A queue of notification IDs
    ids <- character(0)
    # A counter
    n <- 0

    observeEvent(input$show, {
      # Save the ID for removal later
      id <- showNotification(paste("Message", n), duration = NULL)
      ids <<- c(ids, id)
      n <<- n + 1
    })

    observeEvent(input$remove, {
      if (length(ids) > 0)
        removeNotification(ids[1])
      ids <<- ids[-1]
    })
  }
)
}



cleanEx()
nameEx("showTab")
### * showTab

flush(stderr()); flush(stdout())

### Name: showTab
### Title: Dynamically hide/show a tabPanel
### Aliases: showTab hideTab

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {

ui <- navbarPage("Navbar page", id = "tabs",
  tabPanel("Home",
    actionButton("hideTab", "Hide 'Foo' tab"),
    actionButton("showTab", "Show 'Foo' tab"),
    actionButton("hideMenu", "Hide 'More' navbarMenu"),
    actionButton("showMenu", "Show 'More' navbarMenu")
  ),
  tabPanel("Foo", "This is the foo tab"),
  tabPanel("Bar", "This is the bar tab"),
  navbarMenu("More",
    tabPanel("Table", "Table page"),
    tabPanel("About", "About page"),
    "------",
    "Even more!",
    tabPanel("Email", "Email page")
  )
)

server <- function(input, output, session) {
  observeEvent(input$hideTab, {
    hideTab(inputId = "tabs", target = "Foo")
  })

  observeEvent(input$showTab, {
    showTab(inputId = "tabs", target = "Foo")
  })

  observeEvent(input$hideMenu, {
    hideTab(inputId = "tabs", target = "More")
  })

  observeEvent(input$showMenu, {
    showTab(inputId = "tabs", target = "More")
  })
}

shinyApp(ui, server)
}




cleanEx()
nameEx("sidebarLayout")
### * sidebarLayout

flush(stderr()); flush(stdout())

### Name: sidebarLayout
### Title: Layout a sidebar and main area
### Aliases: sidebarLayout sidebarPanel mainPanel

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

# Define UI
ui <- fluidPage(

  # Application title
  titlePanel("Hello Shiny!"),

  sidebarLayout(

    # Sidebar with a slider input
    sidebarPanel(
      sliderInput("obs",
                  "Number of observations:",
                  min = 0,
                  max = 1000,
                  value = 500)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs))
  })
}

# Complete app with UI and server components
shinyApp(ui, server)
}



cleanEx()
nameEx("sizeGrowthRatio")
### * sizeGrowthRatio

flush(stderr()); flush(stdout())

### Name: sizeGrowthRatio
### Title: Create a sizing function that grows at a given ratio
### Aliases: sizeGrowthRatio

### ** Examples

f <- sizeGrowthRatio(500, 500, 1.25)
f(c(400, 400))
f(c(500, 500))
f(c(530, 550))
f(c(625, 700))




cleanEx()
nameEx("sliderInput")
### * sliderInput

flush(stderr()); flush(stdout())

### Name: sliderInput
### Title: Slider Input Widget
### Aliases: sliderInput animationOptions

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

ui <- fluidPage(
  sliderInput("obs", "Number of observations:",
    min = 0, max = 1000, value = 500
  ),
  plotOutput("distPlot")
)

# Server logic
server <- function(input, output) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs))
  })
}

# Complete app with UI and server components
shinyApp(ui, server)
}




cleanEx()
nameEx("splitLayout")
### * splitLayout

flush(stderr()); flush(stdout())

### Name: splitLayout
### Title: Split layout
### Aliases: splitLayout

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

# Server code used for all examples
server <- function(input, output) {
  output$plot1 <- renderPlot(plot(cars))
  output$plot2 <- renderPlot(plot(pressure))
  output$plot3 <- renderPlot(plot(AirPassengers))
}

# Equal sizing
ui <- splitLayout(
  plotOutput("plot1"),
  plotOutput("plot2")
)
shinyApp(ui, server)

# Custom widths
ui <- splitLayout(cellWidths = c("25%", "75%"),
  plotOutput("plot1"),
  plotOutput("plot2")
)
shinyApp(ui, server)

# All cells at 300 pixels wide, with cell padding
# and a border around everything
ui <- splitLayout(
  style = "border: 1px solid silver;",
  cellWidths = 300,
  cellArgs = list(style = "padding: 6px"),
  plotOutput("plot1"),
  plotOutput("plot2"),
  plotOutput("plot3")
)
shinyApp(ui, server)
}



cleanEx()
nameEx("stacktrace")
### * stacktrace

flush(stderr()); flush(stdout())

### Name: stacktrace
### Title: Stack trace manipulation functions
### Aliases: stacktrace captureStackTraces withLogErrors printError
###   printStackTrace conditionStackTrace conditionStackTrace<-
###   ..stacktraceon.. ..stacktraceoff..
### Keywords: internal

### ** Examples

# Keeps tryCatch and withVisible related calls off the
# pretty-printed stack trace

visibleFunction1 <- function() {
  stop("Kaboom!")
}

visibleFunction2 <- function() {
  visibleFunction1()
}

hiddenFunction <- function(expr) {
  expr
}

# An example without ..stacktraceon/off.. manipulation.
# The outer "try" is just to prevent example() from stopping.
try({
  # The withLogErrors call ensures that stack traces are captured
  # and that errors that bubble up are logged using warning().
  withLogErrors({
    # tryCatch and withVisible are just here to add some noise to
    # the stack trace.
    tryCatch(
      withVisible({
        hiddenFunction(visibleFunction2())
      })
    )
  })
})

# Now the same example, but with ..stacktraceon/off.. to hide some
# of the less-interesting bits (tryCatch and withVisible).
..stacktraceoff..({
  try({
    withLogErrors({
      tryCatch(
        withVisible(
          hiddenFunction(
            ..stacktraceon..(visibleFunction2())
          )
        )
      )
    })
  })
})





cleanEx()
nameEx("submitButton")
### * submitButton

flush(stderr()); flush(stdout())

### Name: submitButton
### Title: Create a submit button
### Aliases: submitButton

### ** Examples

if (interactive()) {

shinyApp(
  ui = basicPage(
    numericInput("num", label = "Make changes", value = 1),
    submitButton("Update View", icon("refresh")),
    helpText("When you click the button above, you should see",
             "the output below update to reflect the value you",
             "entered at the top:"),
    verbatimTextOutput("value")
  ),
  server = function(input, output) {

    # submit buttons do not have a value of their own,
    # they control when the app accesses values of other widgets.
    # input$num is the value of the number widget.
    output$value <- renderPrint({ input$num })
  }
)
}



cleanEx()
nameEx("tabPanel")
### * tabPanel

flush(stderr()); flush(stdout())

### Name: tabPanel
### Title: Create a tab panel
### Aliases: tabPanel tabPanelBody

### ** Examples

# Show a tabset that includes a plot, summary, and
# table view of the generated distribution
mainPanel(
  tabsetPanel(
    tabPanel("Plot", plotOutput("plot")),
    tabPanel("Summary", verbatimTextOutput("summary")),
    tabPanel("Table", tableOutput("table"))
  )
)



cleanEx()
nameEx("tabsetPanel")
### * tabsetPanel

flush(stderr()); flush(stdout())

### Name: tabsetPanel
### Title: Create a tabset panel
### Aliases: tabsetPanel

### ** Examples

# Show a tabset that includes a plot, summary, and
# table view of the generated distribution
mainPanel(
  tabsetPanel(
    tabPanel("Plot", plotOutput("plot")),
    tabPanel("Summary", verbatimTextOutput("summary")),
    tabPanel("Table", tableOutput("table"))
  )
)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      radioButtons("controller", "Controller", 1:3, 1)
    ),
    mainPanel(
      tabsetPanel(
        id = "hidden_tabs",
        # Hide the tab values.
        # Can only switch tabs by using `updateTabsetPanel()`
        type = "hidden",
        tabPanelBody("panel1", "Panel 1 content"),
        tabPanelBody("panel2", "Panel 2 content"),
        tabPanelBody("panel3", "Panel 3 content")
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$controller, {
    updateTabsetPanel(session, "hidden_tabs", selected = paste0("panel", input$controller))
  })
}

if (interactive()) {
  shinyApp(ui, server)
}



cleanEx()
nameEx("testServer")
### * testServer

flush(stderr()); flush(stdout())

### Name: testServer
### Title: Reactive testing for Shiny server functions and modules
### Aliases: testServer

### ** Examples

# Testing a server function  ----------------------------------------------
server <- function(input, output, session) {
  x <- reactive(input$a * input$b)
}

testServer(server, {
  session$setInputs(a = 2, b = 3)
  stopifnot(x() == 6)
})


# Testing a module --------------------------------------------------------
myModuleServer <- function(id, multiplier = 2, prefix = "I am ") {
  moduleServer(id, function(input, output, session) {
    myreactive <- reactive({
      input$x * multiplier
    })
    output$txt <- renderText({
      paste0(prefix, myreactive())
    })
  })
}

testServer(myModuleServer, args = list(multiplier = 2), {
  session$setInputs(x = 1)
  # You're also free to use third-party
  # testing packages like testthat:
  #   expect_equal(myreactive(), 2)
  stopifnot(myreactive() == 2)
  stopifnot(output$txt == "I am 2")

  session$setInputs(x = 2)
  stopifnot(myreactive() == 4)
  stopifnot(output$txt == "I am 4")
  # Any additional arguments, below, are passed along to the module.
})



cleanEx()
nameEx("textAreaInput")
### * textAreaInput

flush(stderr()); flush(stdout())

### Name: textAreaInput
### Title: Create a textarea input control
### Aliases: textAreaInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  textAreaInput("caption", "Caption", "Data Summary", width = "1000px"),
  verbatimTextOutput("value")
)
server <- function(input, output) {
  output$value <- renderText({ input$caption })
}
shinyApp(ui, server)

}




cleanEx()
nameEx("textInput")
### * textInput

flush(stderr()); flush(stdout())

### Name: textInput
### Title: Create a text input control
### Aliases: textInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  textInput("caption", "Caption", "Data Summary"),
  verbatimTextOutput("value")
)
server <- function(input, output) {
  output$value <- renderText({ input$caption })
}
shinyApp(ui, server)
}




cleanEx()
nameEx("textOutput")
### * textOutput

flush(stderr()); flush(stdout())

### Name: textOutput
### Title: Create a text output element
### Aliases: textOutput verbatimTextOutput

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  shinyApp(
    ui = basicPage(
      textInput("txt", "Enter the text to display below:"),
      textOutput("text"),
      verbatimTextOutput("verb")
    ),
    server = function(input, output) {
      output$text <- renderText({ input$txt })
      output$verb <- renderText({ input$txt })
    }
  )
}



cleanEx()
nameEx("titlePanel")
### * titlePanel

flush(stderr()); flush(stdout())

### Name: titlePanel
### Title: Create a panel containing an application title.
### Aliases: titlePanel

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  titlePanel("Hello Shiny!")
)
shinyApp(ui, server = function(input, output) { })
}



cleanEx()
nameEx("updateActionButton")
### * updateActionButton

flush(stderr()); flush(stdout())

### Name: updateActionButton
### Title: Change the label or icon of an action button on the client
### Aliases: updateActionButton updateActionLink

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  actionButton("update", "Update other buttons and link"),
  br(),
  actionButton("goButton", "Go"),
  br(),
  actionButton("goButton2", "Go 2", icon = icon("area-chart")),
  br(),
  actionButton("goButton3", "Go 3"),
  br(),
  actionLink("goLink", "Go Link")
)

server <- function(input, output, session) {
  observe({
    req(input$update)

    # Updates goButton's label and icon
    updateActionButton(session, "goButton",
      label = "New label",
      icon = icon("calendar"))

    # Leaves goButton2's label unchanged and
    # removes its icon
    updateActionButton(session, "goButton2",
      icon = character(0))

    # Leaves goButton3's icon, if it exists,
    # unchanged and changes its label
    updateActionButton(session, "goButton3",
      label = "New label 3")

    # Updates goLink's label and icon
    updateActionButton(session, "goLink",
      label = "New link label",
      icon = icon("link"))
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateCheckboxGroupInput")
### * updateCheckboxGroupInput

flush(stderr()); flush(stdout())

### Name: updateCheckboxGroupInput
### Title: Change the value of a checkbox group input on the client
### Aliases: updateCheckboxGroupInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  p("The first checkbox group controls the second"),
  checkboxGroupInput("inCheckboxGroup", "Input checkbox",
    c("Item A", "Item B", "Item C")),
  checkboxGroupInput("inCheckboxGroup2", "Input checkbox 2",
    c("Item A", "Item B", "Item C"))
)

server <- function(input, output, session) {
  observe({
    x <- input$inCheckboxGroup

    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)

    # Can also set the label and select items
    updateCheckboxGroupInput(session, "inCheckboxGroup2",
      label = paste("Checkboxgroup label", length(x)),
      choices = x,
      selected = x
    )
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateCheckboxInput")
### * updateCheckboxInput

flush(stderr()); flush(stdout())

### Name: updateCheckboxInput
### Title: Change the value of a checkbox input on the client
### Aliases: updateCheckboxInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("controller", "Controller", 0, 1, 0, step = 1),
  checkboxInput("inCheckbox", "Input checkbox")
)

server <- function(input, output, session) {
  observe({
    # TRUE if input$controller is odd, FALSE if even.
    x_even <- input$controller %% 2 == 1

    updateCheckboxInput(session, "inCheckbox", value = x_even)
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateDateInput")
### * updateDateInput

flush(stderr()); flush(stdout())

### Name: updateDateInput
### Title: Change the value of a date input on the client
### Aliases: updateDateInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("n", "Day of month", 1, 30, 10),
  dateInput("inDate", "Input date")
)

server <- function(input, output, session) {
  observe({
    date <- as.Date(paste0("2013-04-", input$n))
    updateDateInput(session, "inDate",
      label = paste("Date label", input$n),
      value = date,
      min   = date - 3,
      max   = date + 3
    )
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateDateRangeInput")
### * updateDateRangeInput

flush(stderr()); flush(stdout())

### Name: updateDateRangeInput
### Title: Change the start and end values of a date range input on the
###   client
### Aliases: updateDateRangeInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("n", "Day of month", 1, 30, 10),
  dateRangeInput("inDateRange", "Input date range")
)

server <- function(input, output, session) {
  observe({
    date <- as.Date(paste0("2013-04-", input$n))

    updateDateRangeInput(session, "inDateRange",
      label = paste("Date range label", input$n),
      start = date - 1,
      end = date + 1,
      min = date - 5,
      max = date + 5
    )
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateNumericInput")
### * updateNumericInput

flush(stderr()); flush(stdout())

### Name: updateNumericInput
### Title: Change the value of a number input on the client
### Aliases: updateNumericInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("controller", "Controller", 0, 20, 10),
  numericInput("inNumber", "Input number", 0),
  numericInput("inNumber2", "Input number 2", 0)
)

server <- function(input, output, session) {

  observeEvent(input$controller, {
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    x <- input$controller

    updateNumericInput(session, "inNumber", value = x)

    updateNumericInput(session, "inNumber2",
      label = paste("Number label ", x),
      value = x, min = x-10, max = x+10, step = 5)
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateQueryString")
### * updateQueryString

flush(stderr()); flush(stdout())

### Name: updateQueryString
### Title: Update URL in browser's location bar
### Aliases: updateQueryString

### ** Examples

## Only run these examples in interactive sessions
if (interactive()) {

  ## App 1: Doing "live" bookmarking
  ## Update the browser's location bar every time an input changes.
  ## This should not be used with enableBookmarking("server"),
  ## because that would create a new saved state on disk every time
  ## the user changes an input.
  enableBookmarking("url")
  shinyApp(
    ui = function(req) {
      fluidPage(
        textInput("txt", "Text"),
        checkboxInput("chk", "Checkbox")
      )
    },
    server = function(input, output, session) {
      observe({
        # Trigger this observer every time an input changes
        reactiveValuesToList(input)
        session$doBookmark()
      })
      onBookmarked(function(url) {
        updateQueryString(url)
      })
    }
  )

  ## App 2: Printing the value of the query string
  ## (Use the back and forward buttons to see how the browser
  ## keeps a record of each state)
  shinyApp(
    ui = fluidPage(
      textInput("txt", "Enter new query string"),
      helpText("Format: ?param1=val1&param2=val2"),
      actionButton("go", "Update"),
      hr(),
      verbatimTextOutput("query")
    ),
    server = function(input, output, session) {
      observeEvent(input$go, {
        updateQueryString(input$txt, mode = "push")
      })
      output$query <- renderText({
        query <- getQueryString()
        queryText <- paste(names(query), query,
                       sep = "=", collapse=", ")
        paste("Your query string is:\n", queryText)
      })
    }
  )
}



cleanEx()
nameEx("updateRadioButtons")
### * updateRadioButtons

flush(stderr()); flush(stdout())

### Name: updateRadioButtons
### Title: Change the value of a radio input on the client
### Aliases: updateRadioButtons

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  p("The first radio button group controls the second"),
  radioButtons("inRadioButtons", "Input radio buttons",
    c("Item A", "Item B", "Item C")),
  radioButtons("inRadioButtons2", "Input radio buttons 2",
    c("Item A", "Item B", "Item C"))
)

server <- function(input, output, session) {
  observe({
    x <- input$inRadioButtons

    # Can also set the label and select items
    updateRadioButtons(session, "inRadioButtons2",
      label = paste("radioButtons label", x),
      choices = x,
      selected = x
    )
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateSelectInput")
### * updateSelectInput

flush(stderr()); flush(stdout())

### Name: updateSelectInput
### Title: Change the value of a select input on the client
### Aliases: updateSelectInput updateSelectizeInput updateVarSelectInput
###   updateVarSelectizeInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  p("The checkbox group controls the select input"),
  checkboxGroupInput("inCheckboxGroup", "Input checkbox",
    c("Item A", "Item B", "Item C")),
  selectInput("inSelect", "Select input",
    c("Item A", "Item B", "Item C"))
)

server <- function(input, output, session) {
  observe({
    x <- input$inCheckboxGroup

    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)

    # Can also set the label and select items
    updateSelectInput(session, "inSelect",
      label = paste("Select input label", length(x)),
      choices = x,
      selected = tail(x, 1)
    )
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateSliderInput")
### * updateSliderInput

flush(stderr()); flush(stdout())

### Name: updateSliderInput
### Title: Update Slider Input Widget
### Aliases: updateSliderInput

### ** Examples

## Only run this example in interactive R sessions
if (interactive()) {
  shinyApp(
    ui = fluidPage(
      sidebarLayout(
        sidebarPanel(
          p("The first slider controls the second"),
          sliderInput("control", "Controller:", min=0, max=20, value=10,
                       step=1),
          sliderInput("receive", "Receiver:", min=0, max=20, value=10,
                       step=1)
        ),
        mainPanel()
      )
    ),
    server = function(input, output, session) {
      observe({
        val <- input$control
        # Control the value, min, max, and step.
        # Step size is 2 when input value is even; 1 when value is odd.
        updateSliderInput(session, "receive", value = val,
          min = floor(val/2), max = val+4, step = (val+1)%%2 + 1)
      })
    }
  )
}



cleanEx()
nameEx("updateTabsetPanel")
### * updateTabsetPanel

flush(stderr()); flush(stdout())

### Name: updateTabsetPanel
### Title: Change the selected tab on the client
### Aliases: updateTabsetPanel updateNavbarPage updateNavlistPanel

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    sliderInput("controller", "Controller", 1, 3, 1)
  ),
  mainPanel(
    tabsetPanel(id = "inTabset",
      tabPanel(title = "Panel 1", value = "panel1", "Panel 1 content"),
      tabPanel(title = "Panel 2", value = "panel2", "Panel 2 content"),
      tabPanel(title = "Panel 3", value = "panel3", "Panel 3 content")
    )
  )
))

server <- function(input, output, session) {
  observeEvent(input$controller, {
    updateTabsetPanel(session, "inTabset",
      selected = paste0("panel", input$controller)
    )
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateTextAreaInput")
### * updateTextAreaInput

flush(stderr()); flush(stdout())

### Name: updateTextAreaInput
### Title: Change the value of a textarea input on the client
### Aliases: updateTextAreaInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("controller", "Controller", 0, 20, 10),
  textAreaInput("inText", "Input textarea"),
  textAreaInput("inText2", "Input textarea 2")
)

server <- function(input, output, session) {
  observe({
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    x <- input$controller

    # This will change the value of input$inText, based on x
    updateTextAreaInput(session, "inText", value = paste("New text", x))

    # Can also set the label, this time for input$inText2
    updateTextAreaInput(session, "inText2",
      label = paste("New label", x),
      value = paste("New text", x))
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("updateTextInput")
### * updateTextInput

flush(stderr()); flush(stdout())

### Name: updateTextInput
### Title: Change the value of a text input on the client
### Aliases: updateTextInput

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  sliderInput("controller", "Controller", 0, 20, 10),
  textInput("inText", "Input text"),
  textInput("inText2", "Input text 2")
)

server <- function(input, output, session) {
  observe({
    # We'll use the input$controller variable multiple times, so save it as x
    # for convenience.
    x <- input$controller

    # This will change the value of input$inText, based on x
    updateTextInput(session, "inText", value = paste("New text", x))

    # Can also set the label, this time for input$inText2
    updateTextInput(session, "inText2",
      label = paste("New label", x),
      value = paste("New text", x))
  })
}

shinyApp(ui, server)
}



cleanEx()
nameEx("useBusyIndicators")
### * useBusyIndicators

flush(stderr()); flush(stdout())

### Name: useBusyIndicators
### Title: Enable/disable busy indication
### Aliases: useBusyIndicators

### ** Examples

## Don't show: 
if (rlang::is_interactive()) (if (getRversion() >= "3.4") withAutoprint else force)({ # examplesIf
## End(Don't show)

library(bslib)

ui <- page_fillable(
  useBusyIndicators(),
  card(
    card_header(
      "A plot",
      input_task_button("simulate", "Simulate"),
      class = "d-flex justify-content-between align-items-center"
    ),
    plotOutput("p"),
  )
)

server <- function(input, output) {
  output$p <- renderPlot({
    input$simulate
    Sys.sleep(4)
    plot(x = rnorm(100), y = rnorm(100))
  })
}

shinyApp(ui, server)
## Don't show: 
}) # examplesIf
## End(Don't show)



cleanEx()
nameEx("validate")
### * validate

flush(stderr()); flush(stdout())

### Name: validate
### Title: Validate input values and other conditions
### Aliases: validate need

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

ui <- fluidPage(
  checkboxGroupInput('in1', 'Check some letters', choices = head(LETTERS)),
  selectizeInput('in2', 'Select a state', choices = c("", state.name)),
  plotOutput('plot')
)

server <- function(input, output) {
  output$plot <- renderPlot({
    validate(
      need(input$in1, 'Check at least one letter!'),
      need(input$in2 != '', 'Please choose a state.')
    )
    plot(1:10, main = paste(c(input$in1, input$in2), collapse = ', '))
  })
}

shinyApp(ui, server)

}



cleanEx()
nameEx("varSelectInput")
### * varSelectInput

flush(stderr()); flush(stdout())

### Name: varSelectInput
### Title: Select variables from a data frame
### Aliases: varSelectInput varSelectizeInput

### ** Examples


## Only run examples in interactive R sessions
if (interactive()) {

library(ggplot2)

# single selection
shinyApp(
  ui = fluidPage(
    varSelectInput("variable", "Variable:", mtcars),
    plotOutput("data")
  ),
  server = function(input, output) {
    output$data <- renderPlot({
      ggplot(mtcars, aes(!!input$variable)) + geom_histogram()
    })
  }
)


# multiple selections
## Not run: 
##D shinyApp(
##D  ui = fluidPage(
##D    varSelectInput("variables", "Variable:", mtcars, multiple = TRUE),
##D    tableOutput("data")
##D  ),
##D  server = function(input, output) {
##D    output$data <- renderTable({
##D       if (length(input$variables) == 0) return(mtcars)
##D       mtcars %>% dplyr::select(!!!input$variables)
##D    }, rownames = TRUE)
##D  }
##D )
## End(Not run)

}



cleanEx()
nameEx("verticalLayout")
### * verticalLayout

flush(stderr()); flush(stdout())

### Name: verticalLayout
### Title: Lay out UI elements vertically
### Aliases: verticalLayout

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {

ui <- fluidPage(
  verticalLayout(
    a(href="http://example.com/link1", "Link One"),
    a(href="http://example.com/link2", "Link Two"),
    a(href="http://example.com/link3", "Link Three")
  )
)
shinyApp(ui, server = function(input, output) { })
}



cleanEx()
nameEx("withMathJax")
### * withMathJax

flush(stderr()); flush(stdout())

### Name: withMathJax
### Title: Load the MathJax library and typeset math expressions
### Aliases: withMathJax

### ** Examples

withMathJax(helpText("Some math here $$\\alpha+\\beta$$"))
# now we can just write "static" content without withMathJax()
div("more math here $$\\sqrt{2}$$")



cleanEx()
nameEx("withProgress")
### * withProgress

flush(stderr()); flush(stdout())

### Name: withProgress
### Title: Reporting progress (functional API)
### Aliases: withProgress setProgress incProgress

### ** Examples

## Only run examples in interactive R sessions
if (interactive()) {
options(device.ask.default = FALSE)

ui <- fluidPage(
  plotOutput("plot")
)

server <- function(input, output) {
  output$plot <- renderPlot({
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
      for (i in 1:15) {
        incProgress(1/15)
        Sys.sleep(0.25)
      }
    })
    plot(cars)
  })
}

shinyApp(ui, server)
}



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
