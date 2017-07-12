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
          title = "Prior Distributions",
          verbatimTextOutput(outputId = "test1"),
          verbatimTextOutput(outputId = "test2"),
          verbatimTextOutput(outputId = "test3")
        ), 
        # Tab for feedback/diagnosis
        tabPanel(title = "Feedback/Diagnosis"),
        # Tab for simulation results
        tabPanel(title = "Simulation Results")
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
      # Load configuration file
      configFile <- input$configFile
      config <- load.config(configFile$datapath)
      config <- set.model(config, input$specificModel)
      # Load tree input
      if (is.null(newickInput$data)) return()
      obs.tree <- newickInput$data
      obs.tree <- parse.input.tree(obs.tree, config)
      # Initialize workspace
      ws <- init.workspace(obs.tree, config)
      # Run ABC-SMC
      # res <- run.smc(ws, model=input$specificModel, verbose=F)
    }
  )
  
}

shinyApp(ui = ui, server = server)
