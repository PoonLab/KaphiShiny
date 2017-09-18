## TODO: assign a smaller variable dist <- [[input[[distribution]]]] to clean up code

library(shiny)
library(shinyLP)
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

ui <- fluidPage(
    
  navbarPage(
    
    # Page title
    windowTitle = "Kaphi - Kernel-embedded ABC-SMC for phylodynamic inference",
      
    # Navbar title
    title = strong("Kaphi"),
      
    tabPanel(
      title = "Home",
      jumbotron(header = "Kaphi", content = "Kernel-embedded ABC-SMC for phylodynamic inference", button = FALSE)
    ),
    
    tabPanel(
      title = "Tree Input",
      sidebarLayout(
        
        sidebarPanel(
          # Newick text/file input 
          fluidRow(
            h3(strong(em("Newick Input"))),
            textInput(inputId = "newickString", label = "Enter a Newick String"), 
            fileInput(inputId = "newickFile", label = "Choose a Newick File"),
            actionButton(inputId = "processString", label = "Process String"),
            actionButton(inputId = "processFile", label = "Process File")
          )
        ),
        
        mainPanel(
          # Tree visualization
          selectInput(
            inputId = "treeFormat", 
            label = "Tree Format", 
            choices = c(Rectangular = "rectangular", Circular = "circular", Radial = "radial", Hierarchical = "hierarchical")
          ),
          fluidRow(
            column(width = 6, sliderInput("width", "Panel Width (px)", min = 1, max = 10000, value = 1000)),
            column(width = 6, sliderInput("height", "Panel Height (px)", min = 1, max = 10000, value = 1000))
          ),
          uiOutput("treeVisualization")
        )
        
      )
    ),
    
    tabPanel(
      title = "SMC Settings"
    ),
    
    tabPanel(
      title = "Priors Settings"
    ),
    
    tabPanel(
      title = "Run & Review"
    ),
    
    uiOutput(outputId = "resultsPage")
    
  )
  
)  

server <- function(input, output, session) {
  
}

shinyApp(ui = ui, server = server)