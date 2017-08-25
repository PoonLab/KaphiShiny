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

configTest <- list(
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
  
  # distance settings
  # kernel, sackin, tree.width, etc
  dist="kernel.dist(x, y, decay.factor=0.2, rbf.variance=100.0, sst.control=1.0)",
  
  # cached kernel settings, left alone if not specified in user-provided yaml/distance string
  decay.factor=0.2,
  rbf.variance=100.0,
  sst.control=1.0,
  norm.mode='NONE'
)

specificModel = "yule"

for(i in seq_len(length(parameters[[specificModel]]))) {
  configTest$params[[i]] <- parameters[[specificModel]][[i]]
  configTest$priors$parameters[[specificModel]][[i]] = paste0('r', "(n=1,", ")")
  configTest$prior.densities$parameters[[specificModel]][[i]] = paste0('d', "(arg.prior,", ")")
  configTest$proposals$parameters[[specificModel]][[i]] = paste0('r', "Distribution", "(n=1,", ")")
  configTest$proposal.densities$parameters[[specificModel]][[i]] = paste0('d', "Distribution", "(arg.delta,", ")")
}
