library(flare)

n.sim=100
nvar = 8

datnum=c("normal_n50","normal_n100","normal_n200",
         "normal2_n50","normal2_n100","normal2_n200",
         "st_n50","st_n100","st_n200")

for (sim.i in 1:9) {	
  load(paste("./simulation/data/sim_1_",datnum[sim.i],".Rdata", sep=""))
  
  beta.ladlasso = matrix(NA, n.sim, nvar)
  resid.ladlasso = rep(NA, n.sim)
  
  for (sim in 1:n.sim) {
    lapply(seq_along(simdata[[sim]]), function(i) assign(names(simdata[[sim]])[i], simdata[[sim]][[i]], envir = .GlobalEnv))
    
    result.ladlasso = slim(scale(rbind(x.train)), c(y.train) - mean(c(y.train)),method="lq",q=1,verbose = F)
    
    beta.ladlasso[sim,] = result.ladlasso$beta0[[length(result.ladlasso$beta0)]]
    resid.ladlasso[sim] = (mean(((y.test-mean(y.test))  - scale(x.test) %*% beta.ladlasso[sim,])^2))
  }
  
  save(beta.ladlasso, resid.ladlasso,
       file=paste("./simulation/sim_1_", datnum[sim.i], "_LADL.Rdata", sep="") )
}