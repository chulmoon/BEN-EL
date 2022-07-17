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

# BENEL

BENELSIM = function(sim,datnum) {
  start_time <- Sys.time()
  
  # dataset
  load(paste("./simulation/data/sim_2_",datnum,".Rdata", sep=""))
  # functions
  source("./functions.R")
  
  # load datasets to global environment
  lapply(seq_along(simdata[[sim]]), function(i) assign(names(simdata[[sim]])[i], simdata[[sim]][[i]], envir = .GlobalEnv))
  
  # parameter settings
  epsilontemp = 0.5 # starting value of epsilon
  epsilonstep = 0.5 
  scounter = 0 # counter for smaller
  conditer = 1 # counter for iterations
  lambdas = c(1,1) # initial lambda
  
  while(conditer <= 10){
    result.benel = BENELP(scale(rbind(x.train)),c(y.train) - mean(c(y.train)), nsim=2000,
                          epsilon=epsilontemp,L=10,
                          L1=lambdas[1],L2=lambdas[2],a=10,b=10,itermax=5,tol=0.1)
    
    if (result.benel$acceptance.rate <= 0.7 &  result.benel$acceptance.rate >= 0.6){
      break
    } else if (result.benel$acceptance.rate < 0.6){ # need to set smaller step size
      scounter = scounter + 1
      epsilonstep = epsilonstep/2
      epsilontemp = epsilontemp-epsilonstep # decrease epsilon by epsilonstep
    } else if (result.benel$acceptance.rate > 0.7) { # need to set larger step size
      if (scounter == 0) { # constantly getting larger
        epsilonstep = epsilonstep
      } else { # there was a step that made epsilon smaller
        epsilonstep = epsilonstep/2
      }
      epsilontemp = epsilontemp+epsilonstep # increase epsilon by epsilonstep
    }
    
    conditer = conditer+1
    lambdas = result.benel$lambda[dim(result.benel$lambda)[1],]
  }
  return(result=list(acceptance=result.benel$acceptance.rate,epsilon=epsilontemp,lambdas=lambdas))
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")

## normal
### n=100
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal_n100")
save(result.benel, file = paste("./simulation/sim_2_normal_n100_BENEL_parameter.Rdata", sep=""))
### n=200
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal_n200")
save(result.benel, file = paste("./simulation/sim_2_normal_n200_BENEL_parameter.Rdata", sep=""))

## t
### n=100
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "t_n100")
save(result.benel, file = paste("./simulation/sim_2_t_n100_BENEL_parameter.Rdata", sep=""))
### n=200
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal2_n200")
save(result.benel, file = paste("./simulation/sim_2_t_n200_BENEL_parameter.Rdata", sep=""))

## skew t
### n=100
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "st_n100")
save(result.benel, file = paste("./simulation/sim_2_st_n100_BENEL_parameter.Rdata", sep=""))
### n=200
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "st_n200")
save(result.benel, file = paste("./simulation/sim_2_st_n200_BENEL_parameter.Rdata", sep=""))

stopCluster(cl)