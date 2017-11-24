# Import needed libraries
library(shiny)
library(shinyLP)
library(Kaphi)
library(phylocanvas)

# Import needed files
source("dataStructs.R")

# UI section of the app
ui <- fluidPage(
  
  # Multi-page layout
  navbarPage(
    
    # Page title
    windowTitle = "Kaphi - Kernel-embedded ABC-SMC for phylodynamic inference",
    
    # Navbar title
    title = strong("Kaphi"),
    
    # Splash screen page
    tabPanel(
      title = "Home",
      jumbotron(header = "Kaphi", content = "Kernel-embedded ABC-SMC for phylodynamic inference", button = FALSE)
    ),
    
    # Tree input page
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
    
    # SMC settings page
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
    
    # Priors settings page
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
    
    # Results page
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

# Server section of the app
server <- function(input, output, session) {
  
  # Initialize newickInput to store the nwk string or file inputted by the user
  newickInput <- reactiveValues(data = NULL)
  
  # Initialize config
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
      req(inFile$datapath)
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
    modelParameters = parameters[[input$specificModel]]
    nTabs = length(modelParameters)
    tabs = lapply(seq_len(nTabs), function(i) {
      distribution = paste0(input$specificModel, "Prior", modelParameters[[i]], "Distribution")
      tabPanel(
        paste0(modelParameters[[i]]),
        uiOutput(paste0(input$specificModel, "Prior", modelParameters[[i]])),
        uiOutput(paste0(distribution, "Parameters"))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Creating a distribution  drop down menu input for each specific prior
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      modelParameters = parameters[[input$specificModel]]
      output[[paste0(input$specificModel, "Prior", modelParameters[[i]])]] <- renderUI({
        distribution = paste0(input$specificModel, "Prior", modelParameters[[i]], "Distribution")
        selectInput(inputId = distribution, label = "Distribution",  choices = names(distributions))
      })
    }),
    priority = 100
  )
  
  # Creating a series of numeric inputs for each prior's distribution parameters
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      modelParameters = parameters[[input$specificModel]]
      distributionID = paste0(input$specificModel, "Prior", modelParameters[[i]], "Distribution")
      chosenDistribution = input[[distributionID]]
      output[[paste0(distributionID, "Parameters")]] <- renderUI({
        nNumericInputs = length(distributions[[chosenDistribution]])
        numericInputs = lapply(seq_len(nNumericInputs), function(i) {
          numericInput(
            inputId = paste0(distributionID, chosenDistribution, i),
            label = paste0(names(distributions[[chosenDistribution]])[[i]]),
            value = distributions[[chosenDistribution]][[i]][[3]],
            max = distributions[[chosenDistribution]][[i]][[2]],
            min = distributions[[chosenDistribution]][[i]][[1]]
          )
        })
        do.call(wellPanel, numericInputs)
      })
    }),
    priority = 99
  )
  
  # Displaying proposals for a specific model in tabs
  output$proposalsTabs <- renderUI({
    modelParameters = parameters[[input$specificModel]]
    nTabs = length(modelParameters)
    tabs = lapply(seq_len(nTabs), function(i) {
      distribution = paste0(input$specificModel, "Proposal", modelParameters[[i]], "Distribution")
      tabPanel(
        paste0(modelParameters[[i]]),
        uiOutput(paste0(input$specificModel, "Proposal", modelParameters[[i]])),
        uiOutput(paste0(distribution, "Parameters"))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Creating a distribution  drop down menu input for each specific proposal
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      modelParameters = parameters[[input$specificModel]]
      output[[paste0(input$specificModel, "Proposal", modelParameters[[i]])]] <- renderUI({
        distribution = paste0(input$specificModel, "Proposal", modelParameters[[i]], "Distribution")
        selectInput(inputId = distribution, label = "Distribution",  choices = names(distributions))
      })
    }),
    priority = 100
  )
  
  # Creating a series of numeric inputs for each proposal's distribution parameters
  observe(
    lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
      modelParameters = parameters[[input$specificModel]]
      distributionID = paste0(input$specificModel, "Proposal", modelParameters[[i]], "Distribution")
      chosenDistribution = input[[distributionID]]
      output[[paste0(distributionID, "Parameters")]] <- renderUI({
        nNumericInputs = length(distributions[[chosenDistribution]])
        numericInputs = lapply(seq_len(nNumericInputs), function(i) {
          numericInput(
            inputId = paste0(distributionID, chosenDistribution, i),
            label = paste0(names(distributions[[chosenDistribution]])[[i]]),
            value = distributions[[chosenDistribution]][[i]][[3]],
            max = distributions[[chosenDistribution]][[i]][[2]],
            min = distributions[[chosenDistribution]][[i]][[1]]
          )
        })
        do.call(wellPanel, numericInputs)
      })
    }),
    priority = 99
  )
  
  # Function for creating string expressions of distribution parameters that correspond to config formatting
  distribution.parameters <- function(distributionString, distributionID) {
    distributionParameters <- list()
    for(i in seq_len(length(distributions[[distributionString]]))) {
      distributionParameters[[i]] <- paste0(names(distributions[[distributionString]])[[i]], "=", input[[paste0(distributionID, input[[distributionID]], i)]])
    }
    return(paste0(distributionParameters, collapse = ","))
  }
  
  # Running Kaphi
  observeEvent(
    input$runKaphi,
    {
      
      # Initializing variables needed when running Kaphi
      uniqueTraceFileName <- Sys.time()
      trace <- reactiveValues()
      
      # Setting config class
      class(config) <- "smc.config"
      
      # Requiring SMC settings inputs for populating the config
      req(input$particleNumber)
      req(input$sampleNumber)
      req(input$ESSTolerance)
      req(input$finalEpsilon)
      req(input$finalAcceptanceRate)
      req(input$quality)
      req(input$stepTolerance)
      
      # Populating config with SMC settings
      config$nparticle <- input$particleNumber
      config$nsample <- input$sampleNumber
      config$ess.tolerance <- input$ESSTolerance
      config$final.epsilon <- input$finalEpsilon
      config$final.accept.rate <- input$finalAcceptanceRate
      config$quality <- input$quality
      config$step.tolerance <- input$stepTolerance
      
      # Populating config with priors and proposals
      modelParameters = parameters[[input$specificModel]]
      for(i in seq_len(length(modelParameters))) {
        parameter <- toString(modelParameters[[i]])
        priorDistribution <- paste0(input$specificModel, "Prior", modelParameters[[i]], "Distribution")
        proposalDistribution <- paste0(input$specificModel, "Proposal", modelParameters[[i]], "Distribution")
        # Requiring data needed to populate the config
        req(parameter)
        req(input[[priorDistribution]])
        req(input[[proposalDistribution]])
        req(distribution.parameters(input[[priorDistribution]], priorDistribution))
        req(distribution.parameters(input[[proposalDistribution]], proposalDistribution))
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
      res <- run.smc(ws, trace.file = sprintf("tmp/%s.tsv", uniqueTraceFileName), model=input$specificModel, nthreads = 10)
      
      # Examining the content of the trace file
      trace <- read.table(sprintf("tmp/%s.tsv", uniqueTraceFileName), header=T, sep='\t')
      
      # Rendering means trajectories tabs
      output$meansTrajectories <- renderUI({
        nTabs = length(modelParameters)
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(modelParameters[[i]]),
            plotOutput(outputId = paste0("meanTrajectoryOf", modelParameters[[i]]))
          )
        })
        do.call(tabsetPanel, tabs)
      })
      
      # Rendering means trajectories plots in separate tabs
      observe(
        lapply(seq_len(length(modelParameters)), function(i) {
          output[[paste0("meanTrajectoryOf", modelParameters[[i]])]] <- renderPlot(
            plot(
              sapply(split(trace[[modelParameters[[i]]]]*trace$weight, trace$n), sum),
              type = 'o',
              xlab='Iteration',
              ylab=paste0('Mean ', modelParameters[[i]]),
              cex.lab=1,
              main=paste0('Trajectory of Mean ',  modelParameters[[i]], ' (',  input$specificModel, ' Model, ', input$particleNumber, ' Particles)')
            )
          )
        })
      )
      
      # Rendering posteriors approximations tabs
      output$posteriorsApproximations <- renderUI({
        nTabs = length(modelParameters)
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(modelParameters[[i]]),
            plotOutput(outputId = paste0("posteriorApproximationsOf", modelParameters[[i]]))
          )
        })
        do.call(tabsetPanel, tabs)
      })
      
      # Rendering posteriors approximations plots in separate tabs
      observe(
        lapply(seq_len(length(modelParameters)), function(i) {
          nIterations = length(unique(trace$n)) %/% 10
          nColours = nIterations + 1
          pal = rainbow(n=nColours, start=0, end=0.5, v=1, s=1)
          output[[paste0("posteriorApproximationsOf", modelParameters[[i]])]] <- renderPlot({
            plot.new()
            plot.window(xlim=c(0, 3), ylim=c(0, 20))
            axis(1)
            axis(2)
            title(main=paste0(input$specificModel, " ", modelParameters[[i]]))
            title(xlab=paste0(input$specificModel, ' rate parameter (', modelParameters[[i]], ')'))
            title(ylab="Density")
            box()
            lines(density(trace[[modelParameters[[i]]]][trace$n==1], weights=trace$weight[trace$n==1]))
            for (j in 1:nIterations) {
              temp <- trace[trace$n==j*10,]
              lines(density(temp[[modelParameters[[i]]]], weights=temp$weight), col=pal[j+1], lwd=1.5)
            }
            lines(density(trace[[modelParameters[[i]]]][trace$n==max(trace$n)], weights=trace$weight[trace$n==max(trace$n)]), col='black', lwd=2)
            # Show the prior distribution
            x <- sort(replicate(1000, eval(parse(text=config$priors[[modelParameters[[i]]]]))))
            y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[[modelParameters[[i]]]]))}
            lines(x, y(x), lty=5)
          })
        })
      )
      
    }
  )
  
}

shinyApp(ui = ui, server = server)
