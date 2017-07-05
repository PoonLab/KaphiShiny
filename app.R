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
      # Row for model selection and configuration form download
      fluidRow(
        h3(strong(em("Model Selection and Configuration Form Download"))),
        selectInput("generalModel", "General Model", names(models)),
        selectInput("specificModel", "Specific Model", models[[1]])
      ),
      # Row for configuration upload
      fluidRow(
        h3(strong(em("Configuration Upload")))
      ),
      # Row for running simulation
      fluidRow(
        h3(strong(em("Running Simulation")))
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
        tabPanel(title = "Prior Distributions"), 
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
      inFile <- input$newickFile
      newickInput$data <- read.tree(inFile$datapath)
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
  
}

shinyApp(ui = ui, server = server)
