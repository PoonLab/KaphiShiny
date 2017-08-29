library(shiny)
library(Kaphi)
library(phylocanvas)

distributions <- list(
  exp = list(
    rate = list(Lower = 0, Upper = Inf, Default = 1)
  ),
  gamma = list(
    rate = list(Lower = 0, Upper = Inf, Default = 1),
    shape = list(Lower = 0, Upper = Inf, Default = 1)
  ),
  lnorm = list(
    mean = list(Lower = 0, Upper = Inf, Default = 1),
    sd = list(Lower = 0, Upper = Inf, Default = 1)
  ),
  norm = list(
    mean = list(Lower = 0, Upper = Inf, Default = 1),
    sd = list(Lower = 0, Upper = Inf, Default = 1)
  )
)

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

parameters <- list(
  const.coalescent = list(
    "Ne.tau"
  ),
  sir.dynamic = list(
    "beta",
    "gamma",
    "mu"
  ),
  sir.nondynamic = list(
    "beta",
    "gamma"
  ),
  seir= list(
    "beta",
    "gamma",
    "mu", 
    "alpha" 
  ),
  sis = list(
    "beta",
    "gamma",
    "mu"
  ),
  yule = list(
    "lambda"
  ), 
  bd = list(
    "lambda",
    "mu" 
  ),
  bisse = list(
    "lambda0",
    "lambda1",
    "mu0",
    "mu1",
    "q01",
    "q10"
  ),
  musse = list(
    "lambda1", 
    "lambda2",
    "lambda3", 
    "mu1",
    "mu2",
    "mu3",
    "q12",
    "q13",
    "q21",
    "q23",
    "q31",
    "q32"
  ),
  quasse = list(
    "lambda",
    "mu",
    "char"
  ),
  geosse = list(
    "sA",
    "sB",
    "sAB",
    "xA",
    "xB",
    "dA",
    "dB"
  ),
  bisseness = list(
    "lambda0",
    "lambda1",
    "mu0",
    "mu1",
    "q01",
    "q10",
    "p0c",
    "p0a",
    "p1c",
    "p1a" 
  ),
  classe = list(
    "lambda111",
    "lambda112",
    "lambda122",
    "lambda211",
    "lambda212",
    "lambda222",
    "mu1",
    "mu2",
    "q12",
    "q21"
  )
)

distribution.parameters <- function(distributionString, distributionInputID) {
  distributionParameters <- list()
  for(i in 1:length(distributions[[distributionString]])) {
    distributionParameters[[i]] <- paste0(names(distributions[[distributionString]])[[i]], "=", input[[paste0(distributionInputID, input[[distributionInputID]], distributions[[input[[distributionInputID]]]][[i]])]])
  }
  return(paste0(distributionParameters, collapse = ","))
}

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
      style = "background-color: #ffffff; overflow-y:scroll; max-height:90vh",
      wellPanel(
        # Row for newick text/file input 
        fluidRow(
          h3(strong(em("Newick Input"))),
          textInput(inputId = "newickString", label = "Enter a Newick String"), 
          fileInput(inputId = "newickFile", label = "Choose a Newick File"),
          actionButton(inputId = "processString", label = "Process String"),
          actionButton(inputId = "processFile", label = "Process File")
        )
      ),
      wellPanel(
        # Row for configuration creation
        fluidRow(
          h3(strong(em("SMC Settings Initialization"))),
          numericInput(inputId = "particleNumber", label = "Number of Particles", value = 1000),
          numericInput(inputId = "sampleNumber", label = "Number of Samples", value = 5),
          numericInput(inputId = "ESSTolerance", label = "Effective Sample Size (ESS) Tolerance", value = 1.5),
          numericInput(inputId = "finalEpsilon", label = "Final Epsilon", value = 0.01),
          numericInput(inputId = "finalAcceptanceRate", label = "Final Acceptance Rate", value = 0.015),
          numericInput(inputId = "quality", label = "Quality", value = 0.95),
          numericInput(inputId = "stepTolerance", label = "Step Tolerance", value = 1e-5),
          actionButton(inputId = "initializeSMCSettings", label = "Initialize SMC Settings")
        )
      ),
      wellPanel(
        # Row for model selection and initialization
        fluidRow(
          h3(strong(em("Model Selection and Initialization"))),
          selectInput("generalModel", "General Model", names(models)),
          selectInput("specificModel", "Specific Model", models[[1]]),
          tabsetPanel(
            # Tab for priors
            tabPanel(
              title = "Priors",
              uiOutput("priorsTabs")
            ),
            # Tab for proposals
            tabPanel(
              title = "Proposals",
              uiOutput("proposalsTabs")
            )
          ),
          actionButton(inputId = "initializePriors&Proposals", label = "Initialize Priors & Proposals")
        )
      ),
      wellPanel(
        # Row for running simulation
        fluidRow(
          h3(strong(em("Run Kaphi"))),
          actionButton(inputId = "runKaphi", label = "Run Kaphi")
        )
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
  
  config <- list(
    params=NA,
    priors=list(),
    prior.densities=list(),
    constraints=NULL,
    proposals=list(),
    proposal.densities=list(),
    model=NA,
    
    # SMC settings
    nparticle=1000,
    nsample=5,
    ess.tolerance=1.5,
    final.epsilon=0.01,
    final.accept.rate=0.015,
    quality=0.95,
    step.tolerance=1e-5,
    
    # Distance settings: kernel, sackin, tree.width, etc
    dist='kernel.dist(x,y,decay.factor=0.2,rbf.variance=100.0,sst.control=1.0)',
    
    # Cached kernel settings, left alone if not specified in user-provided yaml/distance string
    decay.factor=0.2,
    rbf.variance=100.0,
    sst.control=1.0,
    norm.mode='NONE'
  )
  
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
  
  # Initializing SMC settings
  observeEvent(
    input$initializeSMCSettings,
    {
      config$nparticle <- input$particleNumber
      config$nsample <- input$sampleNumber
      config$ess.tolerance <- input$ESSTolerance
      config$final.epsilon <- input$finalEpsilon
      config$final.accept.rate <- input$finalAcceptanceRate
      config$quality <- input$quality
      config$step.tolerance <- input$stepTolerance
    }
  )
  
  # Handling general and specific model selection
  observe({
    updateSelectInput(session, "specificModel", choices = models[[input$generalModel]])
  })
  
  # Displaying priors for a specific model in tabs
  output$priorsTabs <- renderUI({
    nTabs = length(parameters[[input$specificModel]])
    tabs = lapply(seq_len(nTabs), function(i) {
      distribution = paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]], "Distribution")
      tabPanel(
        paste0(parameters[[input$specificModel]][[i]]),
        uiOutput(paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]])),
        uiOutput(paste0(distribution, "Parameters"))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Creating a distribution  drop down menu input for each specific prior
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      output[[paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]])]] <- renderUI({
        distribution = paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]], "Distribution")
        selectInput(inputId = distribution, label = "Distribution",  choices = names(distributions))
      })
    }),
    priority = 100
  )
  
  # Creating a series of numeric inputs for each prior's distribution parameters
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      distribution = paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]], "Distribution")
      output[[paste0(distribution, "Parameters")]] <- renderUI({
        nNumericInputs = length(distributions[[input[[distribution]]]])
        numericInputs = lapply(seq_len(nNumericInputs), function(i) {
          numericInput(
            inputId = paste0(distribution, input[[distribution]], distributions[[input[[distribution]]]][[i]]),
            label = paste0(names(distributions[[input[[distribution]]]])[[i]]),
            value = distributions[[input[[distribution]]]][[i]][[3]],
            max = distributions[[input[[distribution]]]][[i]][[2]],
            min = distributions[[input[[distribution]]]][[i]][[1]]
          )
        })
        do.call(wellPanel, numericInputs)
      })
    }),
    priority = 99
  )
    
  # Displaying proposals for a specific model in tabs
  output$proposalsTabs <- renderUI({
    nTabs = length(parameters[[input$specificModel]])
    tabs = lapply(seq_len(nTabs), function(i) {
      distribution = paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]], "Distribution")
      tabPanel(
        paste0(parameters[[input$specificModel]][[i]]),
        uiOutput(paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]])),
        uiOutput(paste0(distribution, "Parameters"))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Creating a distribution  drop down menu input for each specific proposal
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      output[[paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]])]] <- renderUI({
        distribution = paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]], "Distribution")
        selectInput(inputId = distribution, label = "Distribution",  choices = names(distributions))
      })
    }),
    priority = 100
  )
  
  
  # Creating a series of numeric inputs for each proposal's distribution parameters
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      distribution = paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]], "Distribution")
      output[[paste0(distribution, "Parameters")]] <- renderUI({
        nNumericInputs = length(distributions[[input[[distribution]]]])
        numericInputs = lapply(seq_len(nNumericInputs), function(i) {
          numericInput(
            inputId = paste0(distribution, input[[distribution]], distributions[[input[[distribution]]]][[i]]),
            label = paste0(names(distributions[[input[[distribution]]]])[[i]]),
            value = distributions[[input[[distribution]]]][[i]][[3]],
            max = distributions[[input[[distribution]]]][[i]][[2]],
            min = distributions[[input[[distribution]]]][[i]][[1]]
          )
        })
        do.call(wellPanel, numericInputs)
      })
    }),
    priority = 99
  )
  
  # Initializing Priors & Proposals
  observeEvent(
    input$initializePriors&Proposals,
    {
      for(i in 1:length(parameters[[input$specificModel]])) {
        parameter <- toString(parameters[[input$specificModel]][[i]])
        priorDistribution <- paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]], "Distribution")
        proposalDistribution <- paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]], "Distribution")
        configTest$params[[i]] <- parameter
        configTest$priors[[parameter]] = paste0("r", input[[priorDistribution]], "(n=1,", distribution.parameters(input[[priorDistribution]], priorDistribution), ")")
        configTest$prior.densities[[parameter]] = paste0("d", input[[priorDistribution]], "(arg.prior,", distribution.parameters(input[[priorDistribution]], priorDistribution), ")")
        configTest$proposals[[parameter]] = paste0("r", input[[proposalDistribution]], "(n=1,", distribution.parameters(input[[proposalDistribution]], proposalDistribution), ")")
        configTest$proposal.densities[[parameter]] = paste0("d", input[[proposalDistribution]], "(arg.delta,", distribution.parameters(input[[proposalDistribution]], proposalDistribution), ")")
      }
    }
  )
  
}

shinyApp(ui = ui, server = server)
