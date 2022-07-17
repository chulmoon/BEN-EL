library(MASS)
library(MCMCpack)
library(SuppDists)
library(lars)
library(elasticnet)
library(gsl)
library(gmm)
library(GIGrvg)
library(snow)
library(maxLik)
library(invgamma)

n.sim=100
nvar = 30

# EN
ENSIM = function(n.sim, nvar, datnum) {
  start_time = Sys.time()
  
  # load datasets and functions
  load(paste0("./simulation/data/sim_2_",datnum,".Rdata"))
  source("./functions.R")
  beta.en = matrix(NA, n.sim, nvar)
  lambda.en = matrix(NA, n.sim, 2)
  resid.en = rep(NA, n.sim)
  
  for (sim in 1:n.sim) {
    lapply(seq_along(simdata[[sim]]), function(i) assign(names(simdata[[sim]])[i], simdata[[sim]][[i]], envir = .GlobalEnv))
    
    result.en = ENLV(x.train, y.train)
    
    beta.en[sim,] = result.en$beta
    lambda.en[sim,] = c(result.en$s, result.en$lambda2)
    resid.en[sim] = (mean(((y.test-mean(y.test))  - scale(x.test) %*% result.en$beta)^2))
  }
  
  save(beta.en, lambda.en, resid.en,
       file=paste("./simulation/sim_2_", datnum, "_EN.Rdata", sep="") )
}

## normal
### n=100
ENSIM(n.sim, nvar, "normal_n100")
### n=200
ENSIM(n.sim, nvar, "normal_n200")

## t
### n=100
ENSIM(n.sim, nvar, "t_n100")
### n=200
ENSIM(n.sim, nvar, "t_n200")

## skew t
### n=100
ENSIM(n.sim, nvar, "st_n100")
### n=200
ENSIM(n.sim, nvar, "st_n200")
