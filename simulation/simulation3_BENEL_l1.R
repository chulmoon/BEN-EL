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
  start_time <- Sys.time()
  
  # l1 hyperprior: gamma
  r1=seq(0.25,10,by=0.25)
  # d1 0.25 to 10 by 0.25
  d1=seq(0.25,10,by=0.25)
  hyp=expand.grid(r1,d1)
  
  # l2 hyperprior
  chi2=1 # fixed
  nu2=1
  psi2=1
  
  hyp=expand.grid(r1,d1)
  
  # data
  load(paste("./simulation/data/sim_1_",datnum,".Rdata", sep=""))
  # parameters
  load("./simulation/sim_3_BENEL_l1_parameter.Rdata")
  # functions
  source("./functions.R")
  
  # load datasets to global environment
  lapply(seq_along(simdata[[1]]), function(i) assign(names(simdata[[1]])[i], simdata[[1]][[i]], envir = .GlobalEnv))
  
  result.benel = BENELF_chain(x=scale(rbind(x.train)),y=c(y.train) - mean(c(y.train)), 
                        nsim=2000,nchain=4,epsilon=result.benel[[combnum]]$epsilon,L=10,
                        L1=1,L2=1,
                        a=10,b=10,r1=hyp[combnum,1],d1=hyp[combnum,2],
                        nu2=1,chi2=1,psi2=1)
  
  end_time <- Sys.time()
  time = end_time - start_time
  
  return(result=list(beta.benel=apply(result.benel$theta,2,median),lambda.benel=apply(result.benel$penalty,2,median),
                     y.test=y.test,y.train=y.train,
                     x.train=x.train,x.test=x.test,time=time,
                     acceptance.rate = result.benel$acceptance.rate, rhat=result.benel$rhat))
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")

result.benel=parLapply(cl, 1:n.sim, BENELSIM, "normal_n50")
save(result.benel, file = "./simulation/sim_3_BENEL_l1.Rdata")

stopCluster(cl)
