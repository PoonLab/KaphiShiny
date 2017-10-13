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
          tabsetPanel(
            # Tree plot
            tabPanel(
              title = "Tree Plot",
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
        )
      )
    ),
    tabPanel(
      title = "SMC Settings",
      wellPanel(
        # SMC settings
        fluidRow(
          h3(strong(em("SMC Settings"))),
          numericInput(inputId = "particleNumber", label = "Number of Particles", value = 100),
          numericInput(inputId = "sampleNumber", label = "Number of Samples", value = 5),
          numericInput(inputId = "ESSTolerance", label = "Effective Sample Size (ESS) Tolerance", value = 50.0),
          numericInput(inputId = "finalEpsilon", label = "Final Epsilon", value = 0.05),
          numericInput(inputId = "finalAcceptanceRate", label = "Final Acceptance Rate", value = 0.05),
          numericInput(inputId = "quality", label = "Quality", value = 0.95),
          numericInput(inputId = "stepTolerance", label = "Step Tolerance", value = 1e-4)
        )
      )
    ),
    tabPanel(
      title = "Priors Settings",
      sidebarLayout(
        sidebarPanel(
          # Model selection and config settings
          fluidRow(
            h3(strong(em("Model Selection and Config Settings"))),
            selectInput("generalModel", "General Model", names(models)),
            selectInput("specificModel", "Specific Model", models[[1]]),
            tabsetPanel(
              # Priors
              tabPanel(
                title = "Priors",
                uiOutput("priorsTabs")
              ),
              # Proposals
              tabPanel(
                title = "Proposals",
                uiOutput("proposalsTabs")
              )
            ),
            actionButton(inputId = "runKaphi", label = "Run Kaphi")
          )
        ),
        mainPanel(
          tabsetPanel(
            # Prior distributions
            tabPanel(
              title = "Priors Distributions",
              uiOutput(outputId = "priorsDistributionsPlots")
            )
          )
        )
      )
    ),
    tabPanel(
      title = "Results",
      wellPanel(
        tabsetPanel(
          # Plots of tsv file
          tabPanel(
            title = "Results",
            tabsetPanel(
              tabPanel(
                title = "Means Trajectories",
                uiOutput("meansTrajectories")
              ),
              tabPanel(
                title = "Posteriors Approximations",
                uiOutput("posteriorsApproximations")
              )
            )
          )
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
    dist="1*Kaphi::kernel.dist(x, y, decay.factor=0.2, rbf.variance=100, sst.control=1, norm.mode=NONE)",
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
            inputId = paste0(distribution, input[[distribution]], i),
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
            inputId = paste0(distribution, input[[distribution]], i),
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
  # Function for creating string expressions of distribution parameters that correspond to config formatting
  distribution.parameters <- function(distributionString, distributionInputID) {
    distributionParameters <- list()
    for(i in seq_len(length(distributions[[distributionString]]))) {
      distributionParameters[[i]] <- paste0(names(distributions[[distributionString]])[[i]], "=", input[[paste0(distributionInputID, input[[distributionInputID]], i)]])
    }
    return(paste0(distributionParameters, collapse = ","))
  }
  # Initializing config, plotting priors distributions, and running Kaphi
  uniqueTraceFileName <- Sys.time()
  trace <- reactiveValues()
  observeEvent(
    input$runKaphi,
    {
      # Setting config class
      class(config) <- "smc.config"
      # Populating config with SMC settings
      config$nparticle <- input$particleNumber
      config$nsample <- input$sampleNumber
      config$ess.tolerance <- input$ESSTolerance
      config$final.epsilon <- input$finalEpsilon
      config$final.accept.rate <- input$finalAcceptanceRate
      config$quality <- input$quality
      config$step.tolerance <- input$stepTolerance
      # Populating config with priors and proposals
      for(i in seq_len(length(parameters[[input$specificModel]]))) {
        parameter <- toString(parameters[[input$specificModel]][[i]])
        priorDistribution <- paste0(input$specificModel, "Prior", parameters[[input$specificModel]][[i]], "Distribution")
        proposalDistribution <- paste0(input$specificModel, "Proposal", parameters[[input$specificModel]][[i]], "Distribution")
        config$params[[i]] <- parameter
        config$priors[[parameter]] <- paste0("r", input[[priorDistribution]], "(n=1,", distribution.parameters(input[[priorDistribution]], priorDistribution), ")")
        config$prior.densities[[parameter]] <- paste0("d", input[[priorDistribution]], "(arg.prior,", distribution.parameters(input[[priorDistribution]], priorDistribution), ")")
        config$proposals[[parameter]] <- paste0("r", input[[proposalDistribution]], "(n=1,", distribution.parameters(input[[proposalDistribution]], proposalDistribution), ")")
        config$proposal.densities[[parameter]] <- paste0("d", input[[proposalDistribution]], "(arg.delta,", distribution.parameters(input[[proposalDistribution]], proposalDistribution), ")")
      }
      # Setting config model
      config <- set.model(config, input$specificModel)
      # Plotting prior distributions (heavily inspired by plot.smc.config)
      y <- rbind(sapply(1:1000, function(x) sample.priors(config)))
      if (nrow(y) == 1){
        rownames(y)[1] <- names(config$priors)
      }
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
      # Loading tree input
      if (is.null(newickInput$data)) return()
      obs.tree <- newickInput$data
      obs.tree <- parse.input.tree(obs.tree, config)
      # Initializing workspace
      ws <- init.workspace(obs.tree, config)
      # Running ABC-SMC 
      res <- run.smc(ws, trace.file = sprintf("tmp/%s.tsv", uniqueTraceFileName), model=input$specificModel)
      # Examining the content of the trace file
      trace <- read.table(sprintf("tmp/%s.tsv", uniqueTraceFileName), header=T, sep='\t')
      # Rendering means trajectories in separate tabs
      output$meansTrajectories <- renderUI({
        nTabs = length(parameters[[input$specificModel]])
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(parameters[[input$specificModel]][[i]]),
            plotOutput(outputId = paste0("meanTrajectoryOf", parameters[[input$specificModel]][[i]]))
          )
        })
        do.call(tabsetPanel, tabs)
      })
      observe(
        lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
          output[[paste0("meanTrajectoryOf", parameters[[input$specificModel]][[i]])]] <- renderPlot(
            plot(
              sapply(split(trace[[parameters[[input$specificModel]][[i]]]]*trace$weight, trace$n), sum),
              ylim=c(0, 2),
              type='o',
              xlab='Iteration',
              ylab=paste0('Mean', parameters[[input$specificModel]][[i]]),
              cex.lab=1,
              main=paste0('Trajectory of Mean ',  names(parameters[[input$specificModel]][[i]]), ' (',  names(input$specificModel), ' Model, ', input$particleNumber, ' Particles)')
            )
          )
        })
      )
      # Rendering posteriors approximations in separate tabs
      output$posteriorsApproximations <- renderUI({
        nTabs = length(parameters[[input$specificModel]])
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(parameters[[input$specificModel]][[i]]),
            plotOutput(outputId = paste0("posteriorApproximationsOf", parameters[[input$specificModel]][[i]]))
          )
        })
        do.call(tabsetPanel, tabs)
      })
      observe(
        lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
          pal = rainbow(n=6, start=0, end=0.5, v=1, s=1)
          output[[paste0("posteriorApproximationsOf", parameters[[input$specificModel]][[i]])]] <- renderPlot(
            plot(density
                 (trace[[parameters[[input$specificModel]][[i]]]][trace$n==1], weights=trace$weight[trace$n==1]),
                 col=pal[1],
                 lwd=2,
                 main=paste0(names(parameters[[input$specificModel]]), " ", names(parameters[[input$specificModel]][[i]])),
                 xlab=paste0(names(parameters[[input$specificModel]]), ' rate parameter (', names(parameters[[input$specificModel]][[i]]), ')'),
                 cex.lab=1.2
            )
            # ,
            # lapply(seq_len(5), function(i) {
            #   temp <- trace[trace$n==i*10,]
            #   lines(density(temp[[parameters[[input$specificModel]][[i]]]], weights=temp$weight), col=pal[i+1], lwd=1.5)
            # }),
            # lines(density(trace[[parameters[[input$specificModel]][[i]]]][trace$n==max(trace$n)], weights=trace$weight[trace$n==max(trace$n)]), col='black', lwd=2)
          )
        })
      )
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
