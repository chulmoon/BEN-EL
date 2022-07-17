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


# r1: 0.25 to 10 by 0.25
r1=seq(0.25,10,by=0.25)
# d1: 0.25 to 10 by 0.25
d1=seq(0.25,10,by=0.25)
hyp=expand.grid(r1,d1)
n.sim=nrow(hyp)

BENELSIM = function(combnum,datnum) {
  # l1 hyperprior: gamma
  r1=seq(0.25,10,by=0.25)
  # d1 0.25 to 10 by 0.25
  d1=seq(0.25,10,by=0.25)
  hyp=expand.grid(r1,d1)
  
  # l2 hyperprior
  chi2=1 # fixed
  nu2=1
  psi2=1
  
  # dataset
  load(paste("./simulation/data/sim_1_",datnum,".Rdata", sep=""))
  # functions
  source("./functions.R")
  
  # load datasets to global environment
  lapply(seq_along(simdata[[1]]), function(i) assign(names(simdata[[1]])[i], simdata[[1]][[i]], envir = .GlobalEnv))
  
  # parameter settings
  epsilontemp = 0.25 # starting value of epsilon
  epsilonstep = 0.25 
  scounter = 0 # counter for smaller
  conditer = 1 # counter for iterations
  lambdas = c(1,1) # initial lambda
  
  while(conditer <= 10){
    result.benel = BENELF_raw(x=scale(rbind(x.train)),y=c(y.train) - mean(c(y.train)), 
                                               nsim=2000,epsilon=epsilontemp,L=10,
                                               L1=lambdas[1],L2=lambdas[2],
                                               a=10,b=10,r1=hyp[combnum,1],d1=hyp[combnum,2],
                                               nu2=1,chi2=1,psi2=1)
    
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
    #lambdas = result.benel$lambda[dim(result.benel$lambda)[1],]
  }
  
  return(result=list(acceptance=result.benel$acceptance.rate,epsilon=epsilontemp))
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")

result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal_n50")
save(result.benel, file = "./simulation/sim_3_BENEL_l1_parameter.Rdata")

stopCluster(cl)
