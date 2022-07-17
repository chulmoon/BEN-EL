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

# BEN

BENSIM = function(sim,datnum) {
	start_time = Sys.time()
	
	# load datasets and functions
	load(paste0("./simulation/data/sim_1_",datnum,".Rdata"))
	source("./functions.R")
	
	lapply(seq_along(simdata[[sim]]), function(i) assign(names(simdata[[sim]])[i], simdata[[sim]][[i]], envir = .GlobalEnv))
	
	# use BEN
	lambda.ini = GSIN(scale(rbind(x.train)),c(y.train) - mean(c(y.train)))
	result.ben = GSL(scale(rbind(x.train)),c(y.train) - mean(c(y.train)),n.burnin=2000,n.sampler=12000,
	                 lambda1=lambda.ini$lambda1,lambda2=1, itermax=20, tol = 0.05)
	
	beta.ben = result.ben$beta
	lambda.ben = result.ben$lambda
	y.ben = rbind(y.train,y.test)
	x.ben = rbind(x.train,x.test)
	
	end_time = Sys.time()
	time = end_time - start_time
	
	return(result=list(beta.ben=beta.ben,lambda.ben=lambda.ben,
										 y.test=y.test,y.train=y.train,
										 x.train=x.train,x.test=x.test,time=time))
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")

## normal
### n=50
result.ben=parLapply(cl, 1:n.sim, BENSIM, "normal_n50")
save(result.ben, file = paste("./simulation/sim_1_normal_n50_BEN.Rdata", sep=""))
### n=100
result.ben=parLapply(cl, 1:n.sim, BENSIM, "normal_n100")
save(result.ben, file = paste("./simulation/sim_1_normal_n100_BEN.Rdata", sep=""))
### n=200
result.ben=parLapply(cl, 1:n.sim, BENSIM, "normal_n200")
save(result.ben, file = paste("./simulation/sim_1_normal_n200_BEN.Rdata", sep=""))

## normal2
### n=50
result.ben=parLapply(cl, 1:n.sim, BENSIM, "normal2_n50")
save(result.ben, file = paste("./simulation/sim_1_normal2_n50_BEN.Rdata", sep=""))
### n=100
result.ben=parLapply(cl, 1:n.sim, BENSIM, "normal2_n100")
save(result.ben, file = paste("./simulation/sim_1_normal2_n100_BEN.Rdata", sep=""))
### n=200
result.ben=parLapply(cl, 1:n.sim, BENSIM, "normal2_n200")
save(result.ben, file = paste("./simulation/sim_1_normal2_n200_BEN.Rdata", sep=""))

## skew t
### n=50
result.ben=parLapply(cl, 1:n.sim, BENSIM, "st_n50")
save(result.ben, file = paste("./simulation/sim_1_st_n50_BEN.Rdata", sep=""))
### n=100
result.ben=parLapply(cl, 1:n.sim, BENSIM, "st_n100")
save(result.ben, file = paste("./simulation/sim_1_st_n100_BEN.Rdata", sep=""))
### n=200
result.ben=parLapply(cl, 1:n.sim, BENSIM, "st_n200")
save(result.ben, file = paste("./simulation/sim_1_st_n200_BEN.Rdata", sep=""))

stopCluster(cl)
