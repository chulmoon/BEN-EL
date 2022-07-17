# Generating simulation data for simulation 2
## error ~ N(0,15^2)

datagen = function(n.train,n.test,n.sim){
  
  simdata=list()
  
  for (sim in 1:n.sim) {
    beta.true = c(rep(3, 15), rep(0, 15))
    sigma.true = 15
    
    # generate X	
    z1 = rnorm(1)
    z2 = rnorm(1)
    z3 = rnorm(1)
    x.train = x.test = NULL
    for (i in 1:n.train) {
      x.train      = rbind(x.train,
                           c(rnorm(5,mean=z1,sd=0.1),
                             rnorm(5,mean=z2,sd=0.1),
                             rnorm(5,mean=z3,sd=0.1),
                             rnorm(15)))}
    for (i in 1:n.test) {
      x.test       = rbind(x.test,
                           c(rnorm(5,mean=z1,sd=0.1),
                             rnorm(5,mean=z2,sd=0.1),
                             rnorm(5,mean=z3,sd=0.1),
                             rnorm(15)))}
    x.train      = scale(x.train)
    x.test       = scale(x.test)
    
    # generate y
    y.train      = x.train %*% beta.true + sigma.true * rnorm(n.train)
    y.test       = x.test %*% beta.true + sigma.true * rnorm(n.test)
    
    y.train      = y.train - mean(y.train)
    y.test       = y.test - mean(y.test)
    
    simdata[[sim]]=list(z=c(z1,z2,z3), x.train=x.train, x.test=x.test, y.train=y.train, y.test=y.test)
  }
  return(simdata)
}

# parameters
n.test = 400
n.sim = 100

## n=100
n.train  = 100
simdata=datagen(n.train,n.test,n.sim)
save(simdata, file = "./simulation/data/sim_2_normal_n100.Rdata")
## n=200
n.train  = 200
simdata=datagen(n.train,n.test,n.sim)
save(simdata, file = "./simulation/data/sim_2_normal_n200.Rdata")

