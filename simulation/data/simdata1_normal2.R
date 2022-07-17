# Generating simulation data for simulation 1
## error ~ N(3,1) or N(-3,1)

datagen = function(n.train,n.test,n.sim){
  simdata=list()
  
  for (sim in 1:n.sim) {
    p = 8					# number of parameters
    
    # set beta.true and sigma.true
    beta.true = c(3, 1.5, 0, 0, 2, 0, 0, 0)
    sigma.true = 3
    
    # generate X
    x.sigma = matrix(NA, p, p)
    rho = 0.5
    for(i in 1:p) {
      for (j in 1:p) {
        x.sigma[i,j] = rho^abs(i-j)
      }
    }
    x.train = mvrnorm(n=n.train, mu=rep(0,p), Sigma=x.sigma, empirical=T)		# This x is standardized.
    x.test = mvrnorm(n=n.test, mu=rep(0,p), Sigma=x.sigma, empirical=T)		# This x is standardized.
    x.train      = scale(x.train)
    x.test       = scale(x.test)
    
    # generate y
    y.train      = x.train %*% beta.true + rnorm(n.train,3,1)*(2*(runif(n.train)>0.5)-1)
    y.test       = x.test %*% beta.true + rnorm(n.train,3,1)*(2*(runif(n.train)>0.5)-1)
    
    y.train      = y.train - mean(y.train)
    y.test       = y.test - mean(y.test)
    
    simdata[[sim]]=list(x.train=x.train, x.test=x.test, y.train=y.train, y.test=y.test)
  }
  return(simdata)
}

# parameters
n.test = 400
n.sim = 100

## n=50
n.train  = 50
simdata=datagen(n.train,n.test,n.sim)
save(simdata, file = "./simulation/data/sim_1_normal2_n50.Rdata")
## n=100
n.train  = 100
simdata=datagen(n.train,n.test,n.sim)
save(simdata, file = "./simulation/data/sim_1_normal2_n100.Rdata")
## n=200
n.train  = 200
simdata=datagen(n.train,n.test,n.sim)
save(simdata, file = "./simulation/data/sim_1_normal2_n200.Rdata")
