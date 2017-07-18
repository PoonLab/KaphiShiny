library(shiny)
library(Kaphi)
library(phylocanvas)

models <- list(
  "Coalescent" = list(
    Constant = "const.coalescent"
  ),
  "Compartmental" = list(
    SIRD = "sir.dynamic",
    SIRND = "sir.nondynamic",
    SEIR = "seir",
    SIS = "sis"
  ),
  "Networks" = list(
  ),
  "Speciation" = list(
    Yule = "yule", 
    BirthDeath = "bd",
    BiSSE = "bisse",
    MuSSE = "musse",
    QuaSSE = "quasse",
    GeoSSE = "geosse",
    BiSSness = "bisseness",
    ClaSSE  = "classe"
  )
)

ui <- fluidPage(
  
  # Page title
  titlePanel(
    title = h1(strong("Kaphi"), "- Kernel-embedded ABC-SMC for phylodynamic inference"),
    windowTitle = "Kaphi - Kernel-embedded ABC-SMC for phylodynamic inference"
  ),
  
  sidebarLayout(
    
    sidebarPanel( 
      # Allowing independent scrolling in the sidebar
      id = "sidebarPanel",
      style = "overflow-y:scroll; max-height:90vh",
      # Row for newick text/file input 
      fluidRow(
        h3(strong(em("Newick Input"))),
        textInput(inputId = "newickString", label = "Enter a Newick String"), 
        fileInput(inputId = "newickFile", label = "Choose a Newick File"),
        actionButton(inputId = "processString", label = "Process String"),
        actionButton(inputId = "processFile", label = "Process File")
      ),
      # Row for model selection 
      fluidRow(
        h3(strong(em("Model Selection"))),
        selectInput("generalModel", "General Model", names(models)),
        selectInput("specificModel", "Specific Model", models[[1]])
      ),
      # Row for configuration download and upload
      fluidRow(
        h3(strong(em("SMC Configuration Download and Upload"))),
        downloadButton(outputId = "downloadDefaultConfigurationFile", label = "Download Default Configuration File"),
        fileInput(inputId = "configFile", label = "Upload Configuration File")
      ),
      # Row for running simulation
      fluidRow(
        h3(strong(em("Run Kaphi"))),
        actionButton(inputId = "runKaphi", label = "Run Kaphi")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        # Tab for tree plot
        tabPanel(
          title = "Tree Plot",
          selectInput(
            inputId = "treeFormat", 
            label = "Tree Format", 
            choices = c(
              Rectangular = "rectangular",
              Circular = "circular",
              Radial = "radial",
              Hierarchical = "hierarchical"
            )
          ),
          fluidRow(
            column(
              6,
              sliderInput("width", "Panel Width (px)", min = 1, max = 10000, value = 1000)
            ),
            column(
              6,
              sliderInput("height", "Panel Height (px)", min = 1, max = 10000, value = 1000)
            )
          ),
          uiOutput("treeVisualization")
        ), 
        # Tab for prior distributions
        tabPanel(
          title = "Priors Distributions",
          uiOutput(outputId = "priorsDistributionsPlots")
        ), 
        # Tab for feedback, diagnosis, and results 
        tabPanel(
          title = "SMC-ABC Run",
          verbatimTextOutput("consoleHeading"),
          verbatimTextOutput("console"),
          uiOutput("downloadTraceFileButton")
        )
      )
    )
    
  )
  
)  

server <- function(input, output, session) {
  
  newickInput <- reactiveValues(data = NULL)
  
  # Reading tree from newick string
  observeEvent(
    input$processString,
    {
      newickInput$data <- read.tree(text = input$newickString)
    }
  )
  
  # Reading tree from newick file
  observeEvent(
    input$processFile,
    {
      newickFile <- input$newickFile
      newickInput$data <- read.tree(newickFile$datapath)
    }
  )
  
  # Plotting newick input
  output$tree <- renderPhylocanvas({
    if (is.null(newickInput$data)) return()
    phylocanvas(newickInput$data, treetype = input$treeFormat)
  })
  
  # Rendering newick input
  output$treeVisualization <- renderUI({
    phylocanvasOutput("tree", width = input$width, height = input$height)
  })
  
  # Handling general and specific model selection
  observe({
    updateSelectInput(session, "specificModel", choices = models[[input$generalModel]])
  })
  
  # Downloading the config file corresponding to the choosen specific model
  output$downloadDefaultConfigurationFile <- downloadHandler(
    filename = function() {
      paste0(input$specificModel, ".yaml")
    },
    content = function(file) {
      file.copy(paste0("configs/", input$specificModel, ".yaml"), file)
    }
  )
  
  # Handling config file upload by the user
  observeEvent(
    input$configFile,
    {
      # Loading configuration file
      configFile <- input$configFile
      config <- load.config(configFile$datapath)
      config <- set.model(config, input$specificModel)
      # Plotting prior distributions (heavily inspired by plot.smc.config)
      y <- rbind(sapply(1:1000, function(x) sample.priors(config)))
      h <- apply(y, 1, density) 
      output$priorsDistributionsPlots <- renderUI({
        nTabs = length(names(config$priors))
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(names(config$priors)[[i]]),
            plotOutput(outputId = paste0(names(config$priors)[[i]], "Plot"))
          )
        })
        do.call(tabsetPanel, tabs)
      })
      observe(
        lapply(seq_len(length(names(config$priors))), function(i) {
          q <- quantile(y[i,], c(0.05, 0.95))
          output[[paste0(names(config$priors)[[i]], "Plot")]] <- renderPlot(
            plot(h[[i]], xlab=names(h)[i], main=paste0("Sample prior distribution of ", names(config$priors)[[i]]), xlim=q)
          )
        })
      )
    }
  )
  
  uniqueTraceFileName <- Sys.time()
  trace <- reactiveValues()
  # Running Kaphi
  observeEvent(
    input$runKaphi,
    {
      # Loading configuration file
      configFile <- input$configFile
      config <- load.config(configFile$datapath)
      config <- set.model(config, input$specificModel)
      # Loading tree input
      if (is.null(newickInput$data)) return()
      obs.tree <- newickInput$data
      obs.tree <- parse.input.tree(obs.tree, config)
      # Initializing workspace
      ws <- init.workspace(obs.tree, config)
      # Running ABC-SMC and outputing the console output to the user
      output$console <- renderPrint({
        logText()
        return(print(trace[["log"]]))
      })
      logText <- reactive({
        trace[["log"]] <- capture.output(res <- run.smc(ws, trace.file = sprintf("tmp/%s.tsv", uniqueTraceFileName), model=input$specificModel))
      })
      output$consoleHeading <- renderText("Trace of SMC-ABC Run:")
      # Rendering download button to download trace file
      output$downloadTraceFileButton <- renderUI({
        downloadButton(outputId = "downloadTraceFile", label = "Download Trace File")
      })
    }
  )
  
  # Downloading the generated trace file
  output$downloadTraceFile <- downloadHandler(
    filename = function() {
      sprintf("%s.tsv", uniqueTraceFileName)
    },
    content = function(file) {
      file.copy(sprintf("tmp/%s.tsv", uniqueTraceFileName), file)
    }
  )
  
  # Deleting user trace files after the user ends their session
  session$onSessionEnded(function() {
    if (file.exists(sprintf("tmp/%s.tsv", uniqueTraceFileName))) {
      file.remove(sprintf("tmp/%s.tsv", uniqueTraceFileName))
    }
  })
  
}

shinyApp(ui = ui, server = server)
