# Import needed libraries
library(shiny)
library(shinyLP)
library(Kaphi)
library(phylocanvas)

# Import needed files
source("dataStructs.R")

# UI section of the app
ui <- fluidPage(
  
  # tags$style(type="text/css",
  #            ".shiny-output-error { visibility: hidden; }",
  #            ".shiny-output-error:before { visibility: hidden; }"
  # ),
  
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
            actionButton(inputId = "initializeWS", label = "Initialize Workspace"),
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
    modelParams = parameters[[input$specificModel]]
    nTabs = length(modelParams)
    tabs = lapply(seq_len(nTabs), function(i) {
      priorStr <- paste0(input$specificModel, "Prior", modelParams[[i]])
      distr = paste0(priorStr, "Distribution")
      tabPanel(
        paste0(modelParams[[i]]),
        uiOutput(priorStr),
        uiOutput(paste0(distr, "Parameters"))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Creating a distribution  drop down menu input for each specific prior
  observe(
    {
      modelParams = parameters[[input$specificModel]]
      lapply(seq_len(length(modelParams)), function(i) {
        priorStr <- paste0(input$specificModel, "Prior", modelParams[[i]])
        output[[priorStr]] <- renderUI({
          distr = paste0(priorStr, "Distribution")
          selectInput(inputId = distr, label = "Distribution",  choices = names(distributions))
        })
      })
    }
  )
  
  # Creating a series of numeric inputs for each prior's distribution parameters
  observe(
    {
      modelParams = parameters[[input$specificModel]]
      lapply(seq_len(length(modelParams)), function(i) {
        priorStr <- paste0(input$specificModel, "Prior", modelParams[[i]])
        distributionID = paste0(priorStr, "Distribution")
        chosenDistr = input[[distributionID]]
        output[[paste0(distributionID, "Parameters")]] <- renderUI({
          input.vals <- distributions[[chosenDistr]]
          nNumericInputs = length(input.vals)
          numericInputs = lapply(seq_len(nNumericInputs), function(i) {
            numericInput(
              inputId = paste0(distributionID, chosenDistr, i),
              label = paste0(names(input.vals)[[i]]),
              value = input.vals[[i]][[3]],
              max = input.vals[[i]][[2]],
              min = input.vals[[i]][[1]]
            )
          })
          do.call(wellPanel, numericInputs)
        })
      })
    }
  )
  
  
  
  
  
  # Displaying proposals for a specific model in tabs
  output$proposalsTabs <- renderUI({
    modelParams = parameters[[input$specificModel]]
    nTabs = length(modelParams)
    tabs = lapply(seq_len(nTabs), function(i) {
      propStr <- paste0(input$specificModel, "Proposal", modelParams[[i]])
      distr = paste0(propStr, "Distribution")
      tabPanel(
        paste0(modelParams[[i]]),
        uiOutput(propStr),
        uiOutput(paste0(distr, "Parameters"))
      )
    })
    do.call(tabsetPanel, tabs)
  })
  
  # Creating a distribution  drop down menu input for each specific proposal
  observe(
    {
      modelParams <- parameters[[input$specificModel]]
      lapply(seq_len(length(modelParams)), function(i) {
        propStr <- paste0(input$specificModel, "Proposal", modelParams[[i]])
        output[[propStr]] <- renderUI({
          distr = paste0(propStr, "Distribution")
          selectInput(inputId = distr, label = "Distribution",  choices = names(distributions))
        })
      })
    }
  )
  
  # Creating a series of numeric inputs for each proposal's distribution parameters
  observe(
    {
      modelParams = parameters[[input$specificModel]]
      lapply(seq_len(length(parameters[[input$specificModel]])), function(i) {
        propStr <- paste0(input$specificModel, "Proposal", modelParams[[i]])
        distributionID = paste0(propStr, "Distribution")
        chosenDistr = input[[distributionID]]
        output[[paste0(distributionID, "Parameters")]] <- renderUI({
          input.vals <- distributions[[chosenDistr]]
          nNumericInputs = length(input.vals)
          numericInputs = lapply(seq_len(nNumericInputs), function(i) {
            numericInput(
              inputId = paste0(distributionID, chosenDistr, i),
              label = paste0(names(input.vals)[[i]]),
              value = input.vals[[i]][[3]],
              max = input.vals[[i]][[2]],
              min = input.vals[[i]][[1]]
            )
          })
          do.call(wellPanel, numericInputs)
        })
      })
    }
  )
  
  
  
  
  
  # Function for creating string expressions of distribution parameters that correspond to config formatting
  create.str <- function(distributionString, distributionID) {
    distributionParameters <- list()
    for(i in seq_len(length(distributions[[distributionString]]))) {
      distributionParameters[[i]] <- paste0(names(distributions[[distributionString]])[[i]], "=", input[[paste0(distributionID, input[[distributionID]], i)]])
    }
    return(paste0(distributionParameters, collapse = ","))
  }
  


  # Running Kaphi
  observeEvent(
    input$initializeWS,
    {
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
      
      # Setting config model
      config <- set.model(config, input$specificModel)
      
      # Populating config with priors and proposals
      modelParams = parameters[[input$specificModel]]
      for(i in seq_len(length(modelParams))) {
        indiv <- toString(modelParams[[i]])
        priorDistr <- paste0(input$specificModel, "Prior", modelParams[[i]], "Distribution")
        propDistr <- paste0(input$specificModel, "Proposal", modelParams[[i]], "Distribution")
        
        # Requiring data needed to populate the config
        req(indiv)
        req(input[[priorDistr]])
        req(input[[propDistr]])
        req(create.str(input[[priorDistr]], priorDistr))
        req(create.str(input[[propDistr]], propDistr))
        
        config$params[[i]] <- indiv
        config$priors[[indiv]] <- paste0("r", input[[priorDistr]], "(n=1,", create.str(input[[priorDistr]], priorDistr), ")")
        config$prior.densities[[indiv]] <- paste0("d", input[[priorDistr]], "(arg.prior,", create.str(input[[priorDistr]], priorDistr), ")")
        config$proposals[[indiv]] <- paste0("r", input[[propDistr]], "(n=1,", create.str(input[[propDistr]], propDistr), ")")
        config$proposal.densities[[indiv]] <- paste0("d", input[[propDistr]], "(arg.delta,", create.str(input[[propDistr]], propDistr), ")")
      }
      
      
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
      
      # Rendering means trajectories tabs
      output$meansTrajectories <- renderUI({
        nTabs = length(modelParams)
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(modelParams[[i]]),
            plotOutput(outputId = paste0("meanTrajectoryOf", modelParams[[i]]))
          )
        })
        do.call(tabsetPanel, tabs)
      })
      
      # Rendering posteriors approximations tabs
      output$posteriorsApproximations <- renderUI({
        nTabs = length(modelParams)
        tabs = lapply(seq_len(nTabs), function(i) {
          tabPanel(
            paste0(modelParams[[i]]),
            plotOutput(outputId = paste0("posteriorApproximationsOf", modelParams[[i]]))
          )
        })
        do.call(tabsetPanel, tabs)
      })

    }
  )
  
  
  
  
  
  observeEvent(
    input$runKaphi,
    {
      
      # Initializing variables needed when running Kaphi
      uniqueTraceFileName <- Sys.time()
      trace <- reactiveValues()
      
      # this is a chunk of duplicated code frin initializeWS actionButton evaluation
      # reason for this is that we want the user to visualize the priors as often as possible rather than waiting for them to run Kaphi first
      
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
      
      # Setting config model
      config <- set.model(config, input$specificModel)
      
      # Populating config with priors and proposals
      modelParams = parameters[[input$specificModel]]
      for(i in seq_len(length(modelParams))) {
        indiv <- toString(modelParams[[i]])
        priorDistr <- paste0(input$specificModel, "Prior", modelParams[[i]], "Distribution")
        propDistr <- paste0(input$specificModel, "Proposal", modelParams[[i]], "Distribution")
        
        # Requiring data needed to populate the config
        req(indiv)
        req(input[[priorDistr]])
        req(input[[propDistr]])
        req(create.str(input[[priorDistr]], priorDistr))
        req(create.str(input[[propDistr]], propDistr))
        
        config$params[[i]] <- indiv
        config$priors[[indiv]] <- paste0("r", input[[priorDistr]], "(n=1,", create.str(input[[priorDistr]], priorDistr), ")")
        config$prior.densities[[indiv]] <- paste0("d", input[[priorDistr]], "(arg.prior,", create.str(input[[priorDistr]], priorDistr), ")")
        config$proposals[[indiv]] <- paste0("r", input[[propDistr]], "(n=1,", create.str(input[[propDistr]], propDistr), ")")
        config$proposal.densities[[indiv]] <- paste0("d", input[[propDistr]], "(arg.delta,", create.str(input[[propDistr]], propDistr), ")")
      }
      
      
      # Loading tree input
      if (is.null(newickInput$data)) return()
      obs.tree <- newickInput$data
      obs.tree <- parse.input.tree(obs.tree, config)
      
      # Initializing workspace
      ws <- init.workspace(obs.tree, config)
      trace.file <- sprintf("tmp/%s.tsv", uniqueTraceFileName)
      model <- input$specificModel
      nthreads <- 1                                                      # hard code a max?
      verbose <- FALSE                                                   # user can probably modify this later
      
      config <- ws$config
      
      # clear file and write header row
      write.table(t(c(
        'n', 'part.num', 'weight', config$params, paste0('dist.', 1:config$nsample)
      )), file=trace.file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
      
      
      # draw particles from prior distribution, assign weights and simulate data
      ptm <- proc.time()  # start timer
      cat ("Initializing SMC-ABC run with", config$nparticle, "particles\n")
      ws <- initialize.smc(ws, input$specificModel)
      
      #result$niter <- 0
      ws$epsilon <- .Machine$double.xmax
      
      
      ## create reactiveValues objects where we can track the elements in results, shiny.df, and ws
      result <- reactiveValues(niter=0, theta=list(), weights=list(), accept.rate={}, epsilons={},
                               shiny.df = data.frame(n=numeric(), 
                                                     part.num=numeric(), 
                                                     weight=numeric(), 
                                                     sapply(config$params, function(x) {x=numeric()}), 
                                                     sapply(sapply(1:config$nsample, function(y) {paste0('dist.', y)}), function(z) {z=numeric()})),
                               ind=1, 
                               ws=ws)   #tracking ws and updating is very important, otherwise will 'restart' the simulation every 10 iterations
      
      # this section of code will continue to repeat until the stopping condition is met
      observe({
        isolate({
          # this is where we do the expensive computing
          for (iteration in 1:5) {             # chunk of x iterations delay for plot update
            
            result$niter <- result$niter + 1
            
            # update epsilon
            result$ws <- .next.epsilon(result$ws)
            
            # provide some feedback
            lap <- proc.time() - ptm
            cat ("Step ", result$niter, " epsilon:", result$ws$epsilon, " ESS:", .ess(result$ws$weights),
                 "accept:", result$accept.rate[length(result$accept.rate)],
                 "elapsed:", round(lap[['elapsed']],1), "s\n")
            
            # resample particles according to their weights
            if (.ess(result$ws$weights) < config$ess.tolerance) {
              result$ws <- .resample.particles(result$ws)
            }
            
            # perturb particles
            result$ws$accept <- vector()    # vector to keep track of which particles were accepted through parallelization in .perturb.particles
            result$ws$alive <- 0
            result$ws <- .perturb.particles(result$ws, model, nthreads=nthreads)  # Metropolis-Hastings sampling
            
            # record everything
            result$theta[[result$niter]] <- result$ws$particles
            result$weights[[result$niter]] <- result$ws$weights
            result$epsilons <- c(result$epsilons, result$ws$epsilon)
            result$accept.rate <- c(result$accept.rate, result$ws$accepted / result$ws$alive)      # changed result$ws$accept to result$ws$accepted; didn't want dual behaviour of result$ws$accept switching back and forth between vector and int
            
            # write output to file if specified
            for (i in 1:config$nparticle) {
              write.table(
                x=t(c(result$niter, i, round(result$ws$weights[i],10), round(result$ws$particles[i,],5), round(result$ws$dists[,i], 5))),
                file=trace.file,
                append=TRUE,
                sep="\t",
                row.names=FALSE,
                col.names=FALSE
              )
              
              ## SHINY data frame being populated
              result$shiny.df[result$ind,] <- t(c(result$niter, i, round(result$ws$weights[i],10), round(result$ws$particles[i,],5), round(result$ws$dists[,i], 5)))
              result$ind <- result$ind + 1
            }
            
            # report stopping conditions
            if (verbose) {
              cat("run.smc result$niter: ", result$niter, "\n")
              cat ("result$ws$epsilon: ", result$ws$epsilon, "\n");
              cat ("config$final.epsilon: ", config$final.epsilon, "\n");
              cat ("result$accept.rate: ", result$accept.rate, "\n");
              cat ("config$final.accept.rate: ", config$final.accept.rate, "\n");

            }
            
            ## SHINY function for param trajectories and updated distributions --> update delay of ten iterations
            userParams = parameters[[model]]
            if (iteration %% 5 == 0) {
              observe(
                lapply(seq_len(length(userParams)), function(i) {
                  
                  # param trajectory
                  output[[paste0("meanTrajectoryOf", userParams[[i]])]] <- renderPlot(
                    plot(
                      sapply(split(result$shiny.df[[userParams[[i]]]]*result$shiny.df$weight, result$shiny.df$n), sum),
                      type = 'o',
                      xlab='Iteration',
                      ylab=paste0('Mean ', userParams[[i]]),
                      cex.lab=1,
                      main=paste0('Trajectory of Mean ',  userParams[[i]], ' (',  input$specificModel, ' Model, ', input$particleNumber, ' Particles)')
                    )
                  )
                  
                  # use denistiies to visualize posterior approximations
                  nIterations = length(unique(result$shiny.df$n)) %/% 10
                  nColours = nIterations + 1
                  pal = rainbow(n=nColours, start=0, end=0.5, v=1, s=1)
                  output[[paste0("posteriorApproximationsOf", userParams[[i]])]] <- renderPlot({
                    plot(density
                         (result$shiny.df[[userParams[[i]]]][result$shiny.df$n==1], 
                           weights=result$shiny.df$weight[result$shiny.df$n==1]), 
                         col=pal[1], 
                         lwd=2, 
                         main=paste0(model, ' ', config$priors[[userParams[[i]]]]), 
                         xlab=paste0(model, ' rate parameter (', userParams[[i]], ')',
                                     '\nMean: ',
                                     mean(result$shiny.df[[userParams[[i]]]][result$shiny.df$n==max(result$shiny.df$n)]), 
                                     '    Median: ', 
                                     median(result$shiny.df[[userParams[[i]]]][result$shiny.df$n==max(result$shiny.df$n)]),
                                     '\n95% CI (',
                                     quantile(result$shiny.df[[userParams[[i]]]][result$shiny.df$n==max(result$shiny.df$n)], c(0.025, 0.975))[1],
                                     ' , ',
                                     quantile(result$shiny.df[[userParams[[i]]]][result$shiny.df$n==max(result$shiny.df$n)], c(0.025, 0.975))[2],
                                     ')'), 
                         cex.lab=0.8
                    )
                    
                    for (j in 1:nIterations) {
                      temp <- result$shiny.df[result$shiny.df$n==j*10,]
                      lines(density(temp[[userParams[[i]]]], weights=temp$weight), col=pal[j+1], lwd=1.5)
                    }
                    lines(density
                          (result$shiny.df[[userParams[[i]]]][result$shiny.df$n==max(result$shiny.df$n)], 
                            weights=result$shiny.df$weight[result$shiny.df$n==max(result$shiny.df$n)]), 
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
            if (result$accept.rate[result$niter] <= config$final.accept.rate) {
              result$ws$epsilon <- config$final.epsilon
              break  # FIXME: this should be redundant given loop condition above
            }
            
          }
          
          
        })
        
        ## if we're not done yet, then schedule this block to execute again ASAP
        # note that we can be interrupted by other reactive updates to, for instance, update a text output
        if (isolate(result$ws$epsilon != config$final.epsilon)) {   # stopping condition for Kaphi::run.smc
          invalidateLater(0, session)
        }
        
        
    })
  })
  
  
  
  
}

shinyApp(ui = ui, server = server)
