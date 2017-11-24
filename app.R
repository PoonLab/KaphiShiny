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
  # Setting config class
  class(config) <- "smc.config"
  
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
  
  ## SHINY function derived from run.smc in Kaphi
  run.smc.shiny <- function(ws, trace.file='', regex=NA, seed=NA, nthreads=1, verbose=FALSE, model='', modelsList = list(), parametersList = list(), distributionsList = list(), ...) {
    # @param ws: workspace
    # @param obs.tree: object of class 'phylo'
    # @param trace.file: (optional) path to a file to write outputs
    # @param seed: (optional) integer to set random seed
    # @param nthreads: (optional) for running on multiple cores
    # @param ...: additional arguments to pass to config@generator via
    #   simulate.trees()
    
    config <- ws$config
    
    # clear file and write header row
    write.table(t(c(
      'n', 'part.num', 'weight', config$params, paste0('dist.', 1:config$nsample)
    )), file=trace.file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    
    ## SHINY data frame to be populated                                                                                                
    shiny.df <- data.frame(n=numeric(), 
                           part.num=numeric(), 
                           weight=numeric(), 
                           sapply(config$params, function(x) {x=numeric()}), 
                           sapply(sapply(1:config$nsample, function(y) {paste0('dist.', y)}), function(z) {z=numeric()}))
    ind <- 1    # keeps track of row in shiny,df to add to
    
    # space for returned values
    result <- list(niter=0, theta=list(), weights=list(), accept.rate={}, epsilons={})
    
    # draw particles from prior distribution, assign weights and simulate data
    ptm <- proc.time()  # start timer
    cat ("Initializing SMC-ABC run with", config$nparticle, "particles\n")
    ws <- initialize.smc(ws, model, seed=seed, ...)
    
    niter <- 0
    ws$epsilon <- .Machine$double.xmax
    
    # report stopping conditions
    if (verbose) {
      cat ("ws$epsilon: ", ws$epsilon, "\n");
      cat ("config$final.epsilon: ", config$final.epsilon, "\n");
    }
    
    while (ws$epsilon != config$final.epsilon) {
      niter <- niter + 1
      
      # update epsilon
      ws <- .next.epsilon(ws)
      
      # provide some feedback
      lap <- proc.time() - ptm
      cat ("Step ", niter, " epsilon:", ws$epsilon, " ESS:", .ess(ws$weights),
           "accept:", result$accept.rate[length(result$accept.rate)],
           "elapsed:", round(lap[['elapsed']],1), "s\n")
      
      # resample particles according to their weights
      if (.ess(ws$weights) < config$ess.tolerance) {
        ws <- .resample.particles(ws)
      }
      
      # perturb particles
      ws$accept <- vector()    # vector to keep track of which particles were accepted through parallelization in .perturb.particles
      ws$alive <- 0
      ws <- .perturb.particles(ws, model, n.threads=nthreads)  # Metropolis-Hastings sampling
      
      # record everything
      result$theta[[niter]] <- ws$particles
      result$weights[[niter]] <- ws$weights
      result$epsilons <- c(result$epsilons, ws$epsilon)
      result$accept.rate <- c(result$accept.rate, ws$accepted / ws$alive)      # changed ws$accept to ws$accepted; didn't want dual behaviour of ws$accept switching back and forth between vector and int
      
      # write output to file if specified
      for (i in 1:config$nparticle) {
        write.table(
          x=t(c(niter, i, round(ws$weights[i],10), round(ws$particles[i,],5), round(ws$dists[,i], 5))),
          file=trace.file,
          append=TRUE,
          sep="\t",
          row.names=FALSE,
          col.names=FALSE
        )
        
        ## SHINY data frame being populated
        shiny.df[ind,] <- t(c(niter, i, round(ws$weights[i],10), round(ws$particles[i,],5), round(ws$dists[,i], 5)))
        ind <- ind + 1
      }
      
      # report stopping conditions
      if (verbose) {
        cat("run.smc niter: ", niter, "\n")
        cat ("ws$epsilon: ", ws$epsilon, "\n");
        cat ("config$final.epsilon: ", config$final.epsilon, "\n");
        cat ("result$accept.rate: ", result$accept.rate, "\n");
        cat ("config$final.accept.rate: ", config$final.accept.rate, "\n");
      }
      
      
      ## SHINY function for param trajectories and updated distributions --> update delay of ten iterations
      userParams = parametersList[[model]]
      if (niter %% 10 == 0) {
        observe(
        lapply(seq_len(length(userParams)), function(i) {
          
          # param trajectory
          output[[paste0("meanTrajectoryOf", userParams[[i]])]] <- renderPlot(
            plot(
              sapply(split(shiny.df[[userParams[[i]]]]*shiny.df$weight, shiny.df$n), sum),
              type = 'o',
              xlab='Iteration',
              ylab=paste0('Mean ', userParams[[i]]),
              cex.lab=1,
              main=paste0('Trajectory of Mean ',  userParams[[i]], ' (',  input$specificModel, ' Model, ', input$particleNumber, ' Particles)')
            )
          )
          
          # use denistiies to visualize posterior approximations
          nIterations = length(unique(shiny.df$n)) %/% 10
          nColours = nIterations + 1
          pal = rainbow(n=nColours, start=0, end=0.5, v=1, s=1)
          output[[paste0("posteriorApproximationsOf", userParams[[i]])]] <- renderPlot({
            plot(density
                 (shiny.df[[userParams[[i]]]][shiny.df$n==1], 
                   weights=shiny.df$weight[shiny.df$n==1]), 
                 col=pal[1], 
                 lwd=2, 
                 main=paste0('SIR ', config$priors[[userParams[[i]]]]), 
                 xlab=paste0('SIR rate parameter (', userParams[[i]], ')',
                             '\nMean: ',
                             mean(shiny.df[[userParams[[i]]]][shiny.df$n==max(shiny.df$n)]), 
                             '    Median: ', 
                             median(shiny.df[[userParams[[i]]]][shiny.df$n==max(shiny.df$n)]),
                             '\n95% CI (',
                             quantile(shiny.df[[userParams[[i]]]][shiny.df$n==max(shiny.df$n)], c(0.025, 0.975))[1],
                             ' , ',
                             quantile(shiny.df[[userParams[[i]]]][shiny.df$n==max(shiny.df$n)], c(0.025, 0.975))[2],
                             ')'), 
                 cex.lab=0.8
            )
            
            for (j in 1:nIterations) {
              temp <- shiny.df[shiny.df$n==j*10,]
              lines(density(temp[[userParams[[i]]]], weights=temp$weight), col=pal[j+1], lwd=1.5)
            }
            lines(density
                  (shiny.df[[userParams[[i]]]][shiny.df$n==max(shiny.df$n)], 
                    weights=shiny.df$weight[shiny.df$n==max(shiny.df$n)]), 
                  col='black', 
                  lwd=2)
            # Show the prior distribution
            x <- sort(replicate(1000, eval(parse(text=config$priors[[userParams[[i]]]]))))
            y <- function(x) {arg.prior <- x; eval(parse(text=config$prior.densities[[userParams[[i]]]]))}
            lines(x, y(x), lty=5)
          })
        })
        )                                                                                
      }
      
      
      # if acceptance rate is low enough, we're done
      if (result$accept.rate[niter] <= config$final.accept.rate) {
        ws$epsilon <- config$final.epsilon
        break  # FIXME: this should be redundant given loop condition above
      }
    }
    
    # finally sample from the estimated posterior distribution
    ws <- .resample.particles(ws)
    result$theta[[niter]] <- ws$particles
    result$weights[[niter]] <- ws$weights
    result$niter <- niter
    
    for (i in 1:config$nparticle) {
      write.table(
        x=t(c((niter + 1), i, round(ws$weights[i],10), round(ws$particles[i,],5), round(ws$dists[,i], 5))),
        file=trace.file,
        append=TRUE,
        sep="\t",
        row.names=FALSE,
        col.names=FALSE
      )
    }  
    # pack ws and result into one list to be returned
    ret <- list(workspace=ws, result=result)
    return (ret)
  }    # end of modified function
  
  # Initializing variables needed when running Kaphi
  uniqueTraceFileName <- Sys.time()
  trace <- reactiveValues()
  
  # Running Kaphi
  observeEvent(
    input$runKaphi,
    {
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
        config$params[[i]] <- parameter
        config$priors[[parameter]] <- paste0("r", input[[priorDistribution]], "(n=1,", distribution.parameters(input[[priorDistribution]], priorDistribution), ")")
        config$prior.densities[[parameter]] <- paste0("d", input[[priorDistribution]], "(arg.prior,", distribution.parameters(input[[priorDistribution]], priorDistribution), ")")
        config$proposals[[parameter]] <- paste0("r", input[[proposalDistribution]], "(n=1,", distribution.parameters(input[[proposalDistribution]], proposalDistribution), ")")
        config$proposal.densities[[parameter]] <- paste0("d", input[[proposalDistribution]], "(arg.delta,", distribution.parameters(input[[proposalDistribution]], proposalDistribution), ")")
      }
      
      # Setting config model
      config <- set.model(config, input$specificModel)
      
      # Plotting prior distributions (derived from plot.smc.config)
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
      
      res <- run.smc.shiny(ws, trace.file = sprintf("tmp/%s.tsv", uniqueTraceFileName), model=input$specificModel, modelsList = models, parametersList = parameters, distributionsList = distributions)
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
