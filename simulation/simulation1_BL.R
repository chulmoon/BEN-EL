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

# BL

BLSIM = function(sim,datnum) {
	start_time = Sys.time()
	
	# load datasets and functions
	load(paste0("./simulation/data/sim_1_",datnum,".Rdata"))
	source("./functions.R")
	
	lapply(seq_along(simdata[[sim]]), function(i) assign(names(simdata[[sim]])[i], simdata[[sim]][[i]], envir = .GlobalEnv))
	
	# use BL
	lambda.ini = GSIN(scale(rbind(x.train)),c(y.train) - mean(c(y.train)))
	result.bl = GSL_BL(scale(rbind(x.train)),c(y.train) - mean(c(y.train)),n.burnin=2000,n.sampler=12000,lambda1=lambda.ini$lambda1,itermax=13, tol=0.01)
	
	beta.bl = result.bl$beta
	lambda.bl = result.bl$lambda
	
	end_time <- Sys.time()
	time = end_time - start_time
	
	return(result=list(beta.bl=beta.bl,lambda.bl=lambda.bl,
	                   y.test=y.test,y.train=y.train,
	                   x.train=x.train,x.test=x.train,time=time))
}

# run in parallel
corenum = detectCores()
cl = makeCluster(corenum, type="SOCK")

## normal
### n=50
result.bl=parLapply(cl, 1:n.sim, BLSIM, "normal_n50")
save(result.bl, file = paste("./simulation/sim_1_normal_n50_BL.Rdata", sep=""))
### n=100
result.bl=parLapply(cl, 1:n.sim, BLSIM, "normal_n100")
save(result.bl, file = paste("./simulation/sim_1_normal_n100_BL.Rdata", sep=""))
### n=200
result.bl=parLapply(cl, 1:n.sim, BLSIM, "normal_n200")
save(result.bl, file = paste("./simulation/sim_1_normal_n200_BL.Rdata", sep=""))

## normal2
### n=50
result.bl=parLapply(cl, 1:n.sim, BLSIM, "normal2_n50")
save(result.bl, file = paste("./simulation/sim_1_normal2_n50_BL.Rdata", sep=""))
### n=100
result.bl=parLapply(cl, 1:n.sim, BLSIM, "normal2_n100")
save(result.bl, file = paste("./simulation/sim_1_normal2_n100_BL.Rdata", sep=""))
### n=200
result.bl=parLapply(cl, 1:n.sim, BLSIM, "normal2_n200")
save(result.bl, file = paste("./simulation/sim_1_normal2_n200_BL.Rdata", sep=""))

## skew t
### n=50
result.bl=parLapply(cl, 1:n.sim, BLSIM, "st_n50")
save(result.bl, file = paste("./simulation/sim_1_st_n50_BL.Rdata", sep=""))
### n=100
result.bl=parLapply(cl, 1:n.sim, BLSIM, "st_n100")
save(result.bl, file = paste("./simulation/sim_1_st_n100_BL.Rdata", sep=""))
### n=200
result.bl=parLapply(cl, 1:n.sim, BLSIM, "st_n200")
save(result.bl, file = paste("./simulation/sim_1_st_n200_BL.Rdata", sep=""))

stopCluster(cl)
