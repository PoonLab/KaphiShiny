# This file contains the list data structures used to dynamically generate UI elements in the KaphiShiny app

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
    "t.end",
    "N",
    "beta",
    "gamma",
    "mu", 
    "alpha"
  ),
  sir.nondynamic = list(
    "t.end",
    "N",
    "beta",
    "gamma",
    "mu", 
    "alpha"
  ),
  seir= list(
    "t.end",
    "N",
    "beta",
    "gamma",
    "mu", 
    "alpha"
  ),
  sis = list(
    "t.end",
    "N",
    "beta",
    "gamma",
    "mu", 
    "alpha"
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
