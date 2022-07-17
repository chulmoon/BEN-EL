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

  # dataset
  load(paste("./simulation/data/sim_2_",datnum,".Rdata", sep=""))
  # parameters
  load(paste("./simulation/sim_2_",datnum,"_BENEL_parameter.Rdata", sep=""))
  # functions
  source("./functions.R")

  # load datasets to global environment
  lapply(seq_along(simdata[[sim]]), function(i) assign(names(simdata[[sim]])[i], simdata[[sim]][[i]], envir = .GlobalEnv))
  
  # parameter settings
  epsilon = result.benel[[sim]]$epsilon # starting value of epsilon
  lambdas = result.benel[[sim]]$lambdas # initial lambda
  
  result.benel = BENELP_chain(x=scale(rbind(x.train)),y=c(y.train) - mean(c(y.train)),nchain=4,nsim=2000,
                              epsilon=epsilon,L=10,
                              L1=lambdas[1],L2=lambdas[2],a=10,b=10,itermax=5,tol=0.05)
  return(result=list(beta.benel=result.benel$theta,lambda.benel=result.benel$lambda,
                     y.test=y.test,y.train=y.train,
                     x.train=x.train,x.test=x.train,time=time,
                     acceptance.rate = result.benel$acceptance.rate, rhat=result.benel$rhat))
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")

## normal
### n=100
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal_n100")
save(result.benel, file = paste("./simulation/sim_2_normal_n100_BENEL.Rdata", sep=""))
### n=200
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal_n200")
save(result.benel, file = paste("./simulation/sim_2_normal_n200_BENEL.Rdata", sep=""))

## normal2
### n=100
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "t_n100")
save(result.benel, file = paste("./simulation/sim_2_t_n100_BENEL.Rdata", sep=""))
### n=200
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "t_n200")
save(result.benel, file = paste("./simulation/sim_2_t_n200_BENEL.Rdata", sep=""))

## skew t
### n=100
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "st_n100")
save(result.benel, file = paste("./simulation/sim_2_st_n100_BENEL.Rdata", sep=""))
### n=200
result.benel=parLapply(cl, 1:n.sim, BENELSIM, "st_n200")
save(result.benel, file = paste("./simulation/sim_2_st_n200_BENEL.Rdata", sep=""))

stopCluster(cl)