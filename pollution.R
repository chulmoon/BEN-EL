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
library(lasso2)
library(caret)
library(ggcorrplot)

# air pollution data
load("./pollution.Rdata")

# correlation plot
corr = round(cor(pollution[,1:15]), 1)
ggcorrplot(corr, hc.order = TRUE, type = "upper", lab = TRUE)

# load functions
source("./functions.R")

n.sim=100
nvar = 15

###################################################################
# BENEL
###################################################################

BENEL_parameter = function(nseed) {
  # load functions
  source("./functions.R")
  # load datasets to global environment
  load("./pollution.Rdata")
  
  set.seed(nseed)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  # parameter settings
  epsilontemp = 0.3 # starting value of epsilon
  epsilonstep = 0.3 
  scounter = 0 # counter for smaller
  conditer = 1 # counter for iterations
  lambdas = c(0.1,0.1) # initial lambda
  
  while(conditer <= 10){
    result.benel = BENELP(scale(rbind(x.train)),c(y.train) - mean(c(y.train)), nsim=2000,
                          epsilon=epsilontemp,L=10,
                          L1=lambdas[1],L2=lambdas[2],a=5,b=5,itermax=5,tol=0.1)
    
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

BENEL_pollution = function(nseed) {
  # load functions
  source("./functions.R")
  # load datasets 
  load("./pollution.Rdata")
  # load parameters
  load("./pollution_parameter.Rdata")
  
  set.seed(nseed)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  # BENEL
  # parameter settings
  epsilon = result.benel.parameter[[nseed]]$epsilon # starting value of epsilon
  lambdas = result.benel.parameter[[nseed]]$lambdas # initial lambda
  
  result.benel = BENELP_chain(x=scale(rbind(x.train)),y=c(y.train) - mean(c(y.train)),
                                        nchain=4,nsim=2000,
                                        epsilon=epsilon,L=10,
                                        L1=lambdas[1],L2=lambdas[2],a=5,b=5,itermax=10,tol=0.05)
  return(result.benel)
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")
## parameters
result.benel.parameter=parLapply(cl, 1:n.sim, BENEL_parameter)
save(result.benel.parameter, file = paste("./pollution_parameter.Rdata", sep=""))
## BENEL
result.benel=parLapply(cl, 1:n.sim, BENEL_pollution)
stopCluster(cl)

###################################################################
# BL
###################################################################
BL = function(nseed) {
  
  source("./functions.R")
  # load datasets to global environment
  load("./pollution.Rdata")
  
  set.seed(nseed)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  # BL
  lambda.ini = GSIN(scale(rbind(x.train)),c(y.train) - mean(c(y.train)))
  result.bl = GSL_BL(scale(rbind(x.train)),c(y.train) - mean(c(y.train)),
                     n.burnin=2000,n.sampler=12000,lambda1=lambda.ini$lambda1,itermax=10, tol=0.01)
  # BEN
  result.ben = GSL(scale(rbind(x.train)),c(y.train) - mean(c(y.train)),n.burnin=2000,n.sampler=12000,
                   lambda1=lambda.ini$lambda1,lambda2=1, itermax=20, tol = 0.05)
  
  return(result.bl)
}
# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")
result.bl=parLapply(cl, 1:n.sim, BL)
stopCluster(cl)

###################################################################
# BEN
###################################################################
BEN = function(nseed) {

  source("./functions.R")
  # load datasets to global environment
  load("./pollution.Rdata")
  
  set.seed(nseed)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  # BEN
  lambda.ini = GSIN(scale(rbind(x.train)),c(y.train) - mean(c(y.train)))  # BEN
  result.ben = GSL(scale(rbind(x.train)),c(y.train) - mean(c(y.train)),n.burnin=2000,n.sampler=12000,
                   lambda1=lambda.ini$lambda1,lambda2=1, itermax=20, tol = 0.05)
  
  return(result.ben)
}
# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")
result.bl=parLapply(cl, 1:n.sim, BEN)
stopCluster(cl)

###################################################################
# Results of Bayesian methods
###################################################################
p=15
method =1
d=0.5

resid.bl = rep(NA,nfold)	# to record the median squared errors
bsm.bl = matrix(NA,nfold,p)		# short for beta.select.m.blel
beta.bl.cor = matrix(NA,nfold,p)

resid.ben = rep(NA,nfold)	# to record the median squared errors
bsm.ben = matrix(NA,nfold,p)		# short for beta.select.m.blel
beta.ben.cor = matrix(NA,nfold,p)

resid.benel = rep(NA,nfold)	# to record the median squared errors
bsm.benel = matrix(NA,nfold,p)		# short for beta.select.m.blel
beta.benel.cor = matrix(NA,nfold,p)

for (nseed in 1:100){
  # create data
  set.seed(nseed)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  n.train=nrow(x.train)
  n.test=nrow(x.test)
  
  x.test.bar = scale(x.test)
  y.test.bar = y.test - mean(y.test)
  
  ########################################
  # BL
  ########################################
  beta.bl = result.benel[[nseed]]$result.bl$beta
  
  # beta
  beta.temp = beta.bl
  
  intervals = NULL
  interval.prob = NULL
  for (ii in 1:p) {
    if (method == 1) intervals = rbind(intervals, quantile(beta.temp[,ii], probs = c(d/2, 1-d/2)))
    if (method == 2) intervals = rbind(intervals, c(-1,1))
    if (method == 3) intervals = rbind(intervals, c(-1,1) * sqrt(var(beta.temp[,ii])))
  }
  for (ii in 1:p) {
    temp = ecdf(beta.temp[,ii])
    interval.prob = c(interval.prob, temp(intervals[ii,2]) - temp(intervals[ii,1]))
  }
  if (method == 1) bsm.bl[nseed,] = (intervals[,1] < 0 & intervals[,2] > 0)
  if (method == 3) bsm.bl[nseed,] = (interval.prob > d)
  
  beta.upt = beta.temp %*% (diag(1-bsm.bl[nseed,]))
  beta.bl.cor[nseed,] = apply(beta.upt,2,median)
  
  resid.bl[nseed] = mean((y.test.bar - x.test.bar %*% beta.bl.cor[nseed,])^2)
  
  ########################################
  # BEN
  ########################################
  beta.ben = result.benel[[nseed]]$result.ben$beta
  
  # beta
  beta.temp = beta.ben
  
  intervals = NULL
  interval.prob = NULL
  for (ii in 1:p) {
    if (method == 1) intervals = rbind(intervals, quantile(beta.temp[,ii], probs = c(d/2, 1-d/2)))
    if (method == 2) intervals = rbind(intervals, c(-1,1))
    if (method == 3) intervals = rbind(intervals, c(-1,1) * sqrt(var(beta.temp[,ii])))
  }
  for (ii in 1:p) {
    temp = ecdf(beta.temp[,ii])
    interval.prob = c(interval.prob, temp(intervals[ii,2]) - temp(intervals[ii,1]))
  }
  if (method == 1) bsm.ben[nseed,] = (intervals[,1] < 0 & intervals[,2] > 0)
  if (method == 3) bsm.ben[nseed,] = (interval.prob > d)
  
  beta.upt = beta.temp %*% (diag(1-bsm.ben[nseed,]))
  beta.ben.cor[nseed,] = apply(beta.upt,2,median)
  
  resid.ben[nseed] = mean((y.test.bar - x.test.bar %*% beta.ben.cor[nseed,])^2)
  
  ########################################
  # BENEL
  ########################################
  
  beta.benel = result.benel[[ifold]]$result.benel$theta
  
  n.train=nrow(x.train)
  n.test=nrow(x.test)
  
  x.test.bar = scale(x.test)
  y.test.bar = y.test - mean(y.test)
  
  # beta
  beta.temp = beta.benel
  
  intervals = NULL
  interval.prob = NULL
  for (ii in 1:p) {
    if (method == 1) intervals = rbind(intervals, quantile(beta.temp[,ii], probs = c(d/2, 1-d/2)))
    if (method == 2) intervals = rbind(intervals, c(-1,1))
    if (method == 3) intervals = rbind(intervals, c(-1,1) * sqrt(var(beta.temp[,ii])))
  }
  for (ii in 1:p) {
    temp = ecdf(beta.temp[,ii])
    interval.prob = c(interval.prob, temp(intervals[ii,2]) - temp(intervals[ii,1]))
  }
  if (method == 1) bsm.benel[ifold,] = (intervals[,1] < 0 & intervals[,2] > 0)
  if (method == 3) bsm.benel[ifold,] = (interval.prob > d)
  
  beta.upt = beta.temp %*% (diag(1-bsm.benel[ifold,]))
  beta.benel.cor[ifold,] = apply(beta.upt,2,median)
  
  resid.benel[ifold] = mean((y.test.bar - x.test.bar %*% beta.benel.cor[ifold,])^2)
}

###################################################################
# LADL
###################################################################

beta.ladlasso = matrix(NA, n.sim, nvar)
resid.ladlasso = rep(NA, n.sim)

for (ii in 1:100){
  set.seed(ii)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  result.ladlasso = slim(scale(rbind(x.train)), c(y.train) - mean(c(y.train)),method="lq",q=1,verbose = F)
  beta.ladlasso[ii,] = result.ladlasso$beta0[[length(result.ladlasso$beta0)]]
  resid.ladlasso[ii] = (mean(((y.test-mean(y.test))  - scale(x.test) %*% beta.ladlasso[ii,])^2))
  
}

###################################################################
# EN
###################################################################
beta.en = matrix(NA, n.sim, nvar)
resid.en = rep(NA, n.sim)

for (ii in 1:100){
  set.seed(ii)
  
  # create data
  x = as.matrix(pollution[,-16])
  y = pollution[,16]
  trainind = sample(1:60,30)
  x.train = x[trainind,]
  x.test = x[-trainind,]
  y.train = y[trainind]
  y.test = y[-trainind]
  
  x.test.bar = scale(x.test)
  y.test.bar = y.test - mean(y.test)
  
  result.en = ENLV(x.train, y.train)
  
  beta.en[ii,] = result.en$beta
  resid.en[ii] = mean((y.test.bar - x.test.bar %*% beta.en[ii,])^2)
}

###################################################################
# Summary
###################################################################
resdat = data.frame(resid.benel,resid.ben,resid.bl,resid.ladlasso)
colnames(resdat)=c("BEN-EL","BEN","BL","LADL")

par(mar = c(2, 4, 1, 1))
boxplot(resdat,ylab="Prediction Error")
