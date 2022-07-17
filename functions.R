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

#############################################################
# HMC sampler
#############################################################
HMC = function (log_U, grad_U, epsilon, L, q0){
  qq             <-  q0
  pp             <-  rnorm(length(qq),0,1)  
  current_p      <-  pp
  pp             <-  pp - epsilon * grad_U(qq) / 2
  for (i in 1:L){
    qq           <-  qq + epsilon * pp
    if (i!=L) pp <-  pp - epsilon * grad_U(qq)}
  pp             <-  pp - epsilon*grad_U(qq) / 2
  pp             <-  -pp
  current_U      <-  log_U(q0)
  current_K      <-  sum(current_p^2) / 2
  proposed_U     <-  log_U(qq)
  proposed_K     <-  sum(pp^2) / 2
  r              <-  exp(current_U-proposed_U+current_K-proposed_K)
  if(is.nan(r))  r <- 0
  if (runif(1) < r)  { 
    th_new=qq
    accepted=TRUE
  } else { 
    th_new=q0
    accepted=FALSE
  }
  return(list(th=th_new,accepted=accepted))
}


#############################################################
# To compute the Lagrange multiplier:                       #
# this is a modified version of owen's lagrange optimizer   #
# for mean. We added some twicks so it will work for linear #
# models.                                                   #
#############################################################

lag = function (x,y, theta, lam, maxit = 25, gradtol = 1e-07, svdtol = 1e-09){
  x <- as.matrix(x);y <- as.vector(y);n <- nrow(x);  p <- ncol(x)
  llog  = function(z){
    eps = 1/n ;ans = z;lo= (z<eps)
    ans[ lo  ] = log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
    ans[ !lo ] = log( z[!lo] )
    ans}
  llogp  = function(z){
    eps = 1/n; ans = z; lo = (z<eps)
    ans[ lo  ] = 2.0/eps - z[lo]/eps^2
    ans[ !lo ] = 1/z[!lo]
    ans}
  llogpp = function(z){
    eps = 1/n;   ans = z;  lo = (z<eps)
    ans[ lo  ] = -1.0/eps^2
    ans[ !lo ] = -1.0/z[!lo]^2
    ans}
  logelr = function(z,theta,lam){
    if(n <= p)                  stop("Need more observations than variables in logelr.")
    theta = as.vector(theta)
    if(length(theta)!= p)       stop("Length of parameters doesn't match number of variables in logelr.")
    arg = 1 + z %*% lam
    return(-sum( llog(arg)))}
  theta <- as.vector(theta)
  if (length(theta) != p)     stop("theta must have same dimension as observation vectors.")
  if (n <= p)                 stop("Need more observations than length(theta) in el.test().")
  z     <-  x*drop(y-x%*%theta)
  TINY  <- sqrt(.Machine$double.xmin)
  scale <- mean(abs(z)) + TINY
  z     <- z/scale
  if (!missing(lam)) {
    lam <- as.vector(lam)
    lam <- lam * scale
    if (logelr(z, rep(0, p), lam) > 0)  lam     <- rep(0, p)
  }
  if (missing(lam))                   lam     <- rep(0, p)
  if (svdtol < TINY)                  svdtol  <- TINY
  if (gradtol < TINY)                 gradtol <- TINY
  nwts        <- c(3^-c(0:3), rep(0, 12))
  gwts        <- 2^(-c(0:(length(nwts) - 1)))
  gwts        <- (gwts^2 - nwts^2)^0.5
  gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
  nits        <- 0
  gsize       <- gradtol + 1
  while (nits < maxit && gsize > gradtol) {
    arg          <- 1 + z %*% lam
    wts1         <- as.vector(llogp(arg))
    wts2         <- as.vector(-llogpp(arg))^0.5
    grad         <- as.matrix(-z * wts1)
    grad         <- as.vector(rowsum(grad, rep(1, nrow(grad))))
    gsize        <- mean(abs(grad))
    hess         <- z * wts2
    svdh         <- svd(hess)
    if (min(svdh$d) < max(svdh$d) * svdtol)     svdh$d <- svdh$d + max(svdh$d) * svdtol
    nstep        <- svdh$v %*% (t(svdh$u)/svdh$d)
    nstep        <- as.vector(nstep %*% matrix(wts1/wts2, n, 1))
    gstep        <- -grad
    if (sum(nstep^2) < sum(gstep^2))            gstep <- gstep * (sum(nstep^2)^0.5/sum(gstep^2)^0.5)
    ologelr      <- -sum(llog(arg))
    ninner       <- 0
    for (i in 1:length(nwts)) {
      nlogelr <- logelr(z, rep(0, p), lam + nwts[i] * nstep + gwts[i] * gstep)
      if (nlogelr < ologelr) {
        lam      <- lam + nwts[i] * nstep + gwts[i] * gstep
        ninner   <- i
        break}}
    nits  <- nits + 1
    if (ninner == 0)  nits <- maxit}
  list(`-2LLR` = -2 * nlogelr,
       Pval    = 1 - pchisq(-2 * nlogelr, df = p),
       lambda  = lam/scale,
       grad    = grad * scale,
       hess    = t(hess) %*%hess * scale^2,
       weights = wts1/n)}

#############################################################
# BENEL
#############################################################
BENEL <- function(x,y,nsim=2000,nwarm=floor(nsim/2),
                  seed=sample.int(.Machine$integer.max, 1),epsilon,L,L1,L2,a,b){
  set.seed(seed)
  x                 <- as.matrix(x)
  y                 <- as.vector(y)
  n                 <- nrow(x)
  p                 <- ncol(x)
  theta.array       <- array(NA, dim = c(nsim,p))
  tau_1.array       <- array(NA, dim = c(nsim,p))
  sigma.array       <- array(NA, dim = nsim)
  acceptance        <- rep(NA, nsim - 1)
  
  ind.tau           <- rep(0,nsim)
  ind.sigma         <- rep(0,nsim)
  
  psi.tau           <- rep(0,nsim)
  chi.tau           <- array(NA, dim = c(nsim,p))
  taumean           <- array(NA, dim = c(nsim,p))
  
  if(L1 == 0)    L1 = L1 + .Machine$double.eps
  if(L2 == 0)    L2 = L2 + .Machine$double.eps
  sigma.array[1]    <- enet(x,y,lambda=L2/(L2+L1),max.steps=50)$sigma2   # lambda = L2/(L1+L2) 
  theta.array[1,]   <- coef(lm(y~. -1, data = cbind.data.frame(y = y, x = x))) 
  psi.tau0        <-  L1^2/(4*L2*sigma.array[1])
  chi.tau0        <-  (L2*theta.array[1,]^2/sigma.array[1])
  for(j in seq(p))  tau_1.array[1,j] = GIGrvg::rgig(n = 1,lambda = 0.5,chi = chi.tau0[j],psi = psi.tau0)
  #####
  for(k in 2:nsim){
    if(k%% 100 ==0) print(paste0("BENEL iteration:"," ", k))
    invD = diag((tau_1.array[k-1,]+1)/tau_1.array[k-1,])
    
    log_U <- function(th){
      lam <-  lag(x = x,y = y, theta = th, maxit=25, gradtol=1e-7, svdtol = 1e-9 )$lambda
      arg1 <- 0
      for(i in 1:n) arg1 <- arg1 + drop(log(1 + t(lam)%*%x[i,]%*%(y[i]-t(x[i,])%*%th)))
      return(arg1 + 0.5*L2*drop(t(th)%*%invD%*%th)/sigma.array[k-1])}
    
    
    grad_U <- function(th){
      lam  <- lag(x = x,y = y, theta = th, maxit=25, gradtol=1e-7, svdtol = 1e-9 )$lambda  
      arg1 <- 0
      for(i in 1:n) arg1 <- arg1 - t(lam)%*%x[i,]%*%t(x[i,])/drop(1 + t(lam)%*%x[i,]%*%(y[i]-t(x[i,])%*%th))
      return(t(arg1+L2*t(th)%*%invD/sigma.array[k-1]))}
    
    next.sample = HMC(log_U = log_U, grad_U = grad_U, epsilon = epsilon, L = L, q0 = theta.array[k-1,])
    
    theta.array[k,] = next.sample$th
    
    acceptance[k-1] <- next.sample$accepted
    
    # Sampling tau
    psi.tau_1        <-  L1^2/(4*L2*sigma.array[k-1])
    chi.tau_1        <-  (L2*theta.array[k,]^2/sigma.array[k-1])
    
    for(j in seq(p)) {
      tau_1.array[k,j] = GIGrvg::rgig(n = 1,lambda = 0.5,chi = chi.tau_1[j],psi = psi.tau_1)
    }
    
    # Sampling Sigma(^2)
    shape.sigma    <- a +p
    rate.sigma     <- b +0.5*sum(L2*theta.array[k,]^2*(tau_1.array[k,]+1)/tau_1.array[k,]+
                                   0.25*L1^2*(tau_1.array[k,]+1)/L2)
    
    sigma.array[k] <- invgamma::rinvgamma(1,shape = shape.sigma, rate= rate.sigma)
    
  }
  
  acceptance.rate <- sum(acceptance[-(1:nwarm)]) / (nsim - nwarm - 1)
  
  list(theta = theta.array[-(1:nwarm),], tau_1 = tau_1.array[-(1:nwarm),], 
       sigma = sigma.array[-(1:nwarm)],acceptance.rate=acceptance.rate)
}


BENEL_chain <- function(x,y,nsim=2000,nwarm=floor(nsim/2),nchain=4,
                        seed=sample.int(.Machine$integer.max, 1),epsilon,L=10,L1,L2,a=10,b=10){
  result.chain = list()
  thetalist = list()
  theta = tau_1 = sigma = acceptance = NULL
  for (ii in 1:nchain){
    result.chain[[ii]]=BENEL(x,y,nsim=nsim,seed=seed+ii,epsilon=epsilon,L=L,L1=L1,L2=L2,a=a,b=b)
    thetalist[[ii]]=result.chain[[ii]]$theta
    theta = rbind(theta,result.chain[[ii]]$theta)
    tau_1 = rbind(tau_1,result.chain[[ii]]$tau_1)
    sigma = c(sigma,result.chain[[ii]]$sigma)
    acceptance = c(acceptance,result.chain[[ii]]$acceptance.rate)
  }
  rhat = splitrhat(ndraw=( (nsim-nwarm)*nchain),theta=thetalist)
  return(list(theta = theta, tau_1 = tau_1, 
              sigma = sigma, acceptance.rate=acceptance, rhat=rhat))
}

# examine convergence
colVar <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=var, na.rm=na.rm)
## split-R-hat
splitrhat = function(ndraw=8000,theta){
  thetabar=t(sapply(theta,colMeans)) # theta-bar
  thetavar=t(sapply(theta,colVar)) # variance of theta, s^2_m
  B=ndraw*colVar(thetabar)
  W=colMeans(thetavar)
  rhat=sqrt( ( (ndraw-1)/(ndraw)*W + B/ndraw )/W)
  return(rhat)
}


#############################################################
# BENEL HMC chains
#############################################################
BENELP_chain <- function(x,y,nsim,nwarm=floor(nsim/2),nchain=4,L1,L2,epsilon,L,a,b,itermax=10,tol=0.1){
  n = dim(x)[1]
  p = dim(x)[2]
  flag = 1
  counter = 1
  lambda.record = c(L1, L2)
  f.record = NULL

  while(flag) {
    print(counter)
    print(paste("BENEL f.record=",f.record))
    print(lambda.record)
    result = BENEL_chain(x=x,y=y,nsim=nsim,nchain=nchain,L1=L1,L2=L2,epsilon=epsilon,L=L,a=a,b=b)
    print(paste("acceptance rate=",result$acceptance.rate))
    # E step
    # conditional likelihood function 
    # a: c(L1, L2)
    f = function(a) {
      a1 = a[1]
      a2 = a[2]

      temp1 = p*log(a1)

      temp2 = 0
      for (j in 1:p) {
        temp2 = temp2 + mean( (result$tau_1[,j]+1) / (result$tau_1[,j]) * result$theta[,j]^2 / result$sigma)
        #temp2 = temp2 + median( (result$tau_1[,j]+1) / (result$tau_1[,j]) * result$theta[,j]^2 / result$sigma)
      }
      temp2 = -a2/2*temp2
      
      temp3 = 0
      for (j in 1:p) {
        temp3 = temp3 + mean( (result$tau_1[,j]+1) / result$sigma)
        #temp3 = temp3 + median( (result$tau_1[,j]+1) / result$sigma)
      }
      temp3 = -a1^2/(8*a2) * temp3
      # return the value of function Q
      return(- temp1 - temp2 - temp3)
    }
    
    # gradient functions for M-step
    ## a is same as above
    gf = function(a) {
      a1 = a[1]
      a2 = a[2]
      # gradient of lambda1
      temp1.1 = p/a1
      
      temp1.2 = 0
      for (j in 1:p) {
        temp1.2 = temp1.2 + mean( (result$tau_1[,j]+1) / result$sigma)
      }
      result.1 = temp1.1 - a1/(4*a2) * temp1.2
      
      # gradient of lambda2
      temp2 = 0
      for (j in 1:p) {
        temp2 = temp2 + mean( (result$tau_1[,j]+1) / (result$tau_1[,j]) * result$theta[,j]^2 / result$sigma)
      }
      temp2 = -0.5 * temp2
      
      temp3 = 0
      for (j in 1:p) {
        temp3 = temp3 + mean( (result$tau_1[,j]+1) / result$sigma )
      }
      temp3 = a1^2/(8*a2^2) * temp3
      result.2 = temp2 + temp3
      return(c(-result.1, -result.2))
    }
    
    lambda.current = c(L1,L2)
    # M step
    lambda.update = constrOptim(lambda.current, f=f, grad=gf, ui=diag(c(1,1)), ci=.Machine$double.eps*c(1,1))$par
    flag = (counter < itermax) & 
      ( sum((lambda.current - lambda.update)^2) >= tol )
    convergence = ( sum((lambda.current - lambda.update)^2) < tol )
    lambda.record = rbind(lambda.record, lambda.update)
    f.current = f(lambda.current)
    f.record = c(f.record, f.current)
    L1 = lambda.update[1]
    L2 = lambda.update[2]
    counter = counter + 1
  }
  out = list(theta = result$theta, sigma = result$sigma, 
             tau_1 = result$tau_1, lambda = lambda.record, 
             acceptance.rate = result$acceptance.rate, rhat=result$rhat,
             convergence = convergence, f=f.record)
  return(out)
  
}


#############################################################
# BENEL Full Bayes (FB)
#############################################################
BENELF_base <- function(x,y,nsim,nwarm=floor(nsim/2),
                       epsilon,L=10,L1,L2,a=10,b=10,
                       r1,d1,nu2,psi2,chi2){
  
  x                       <- as.matrix(x)
  y                       <- as.vector(y)
  n                       <- nrow(x)
  p                       <- ncol(x)
  theta.array             <- array(NA,dim=c(nsim,p))
  tau_1.array             <- array(NA,dim=c(nsim,p)) 
  sigma.array             <- array(NA,dim=nsim)
  acceptance              <- rep(NA, nsim - 1)
  
  if(L1==0)                  L1 = L1 + .Machine$double.eps
  if(L2==0)                  L2 = L2 + .Machine$double.eps
  
  sigma.array[1]    <- elasticnet::enet(x,y,lambda=L2/(L2+L1),max.steps=50)$sigma2   # lambda = L2/(L1+L2) 
  theta.array[1,]   <- coef(lm(y~. -1, data = cbind.data.frame(y = y, x = x))) 
  psi.tau0        <-  L1^2/(4*L2*sigma.array[1])
  chi.tau0        <-  (L2*theta.array[1,]^2/sigma.array[1])
  for(j in seq(p))  tau_1.array[1,j] = GIGrvg::rgig(n = 1,lambda = 0.5,chi = chi.tau0[j],psi = psi.tau0)
  penalty         <- array(NA,dim=c(nsim,2))
  
  penalty[1,]=c(L1,L2)
  
  for (k in 2:nsim){
    if(k%% 100 ==0) print(paste0("BENEL iteration:"," ", k))
    
    # Sample theta's
    invD = diag((tau_1.array[k-1,]+1)/tau_1.array[k-1,])
    
    log_U <- function(th){
      lam <-  lag(x = x,y = y, theta = th, maxit=55, gradtol=1e-7, svdtol = 1e-9 )$lambda
      arg1 <- 0
      for(i in 1:n) arg1 <- arg1 + drop(log(1 + t(lam)%*%x[i,]%*%(y[i]-t(x[i,])%*%th)))
      return(arg1 + 0.5*L2*drop(t(th)%*%invD%*%th)/sigma.array[k-1])}
    
    grad_U <- function(th){
      lam  <- lag(x = x,y = y, theta = th, maxit=55, gradtol=1e-7, svdtol = 1e-9 )$lambda  
      arg1 <- 0
      for(i in 1:n) arg1 <- arg1 - t(lam)%*%x[i,]%*%t(x[i,])/drop(1 + t(lam)%*%x[i,]%*%(y[i]-t(x[i,])%*%th))
      return(t(arg1+L2*t(th)%*%invD/sigma.array[k-1]))}
    
    next.sample = HMC(log_U = log_U, grad_U = grad_U, epsilon = epsilon, L = L, q0 = theta.array[k-1,])
    
    theta.array[k,] = next.sample$th
    acceptance[k-1] <- next.sample$accepted

    # Sampling tau
    psi.tau_1        <-  L1^2/(4*L2*sigma.array[k-1])
    chi.tau_1        <-  (L2*theta.array[k,]^2/sigma.array[k-1])
    
    for(j in seq(p)) {
      tau_1.array[k,j] = GIGrvg::rgig(n = 1,lambda = 0.5,chi = chi.tau_1[j],psi = psi.tau_1)
    }
    
    # Sampling Sigma(^2)
    shape.sigma    <- a +p
    rate.sigma     <- b +0.5*sum(L2*theta.array[k,]^2*(tau_1.array[k,]+1)/tau_1.array[k,]+
                                   0.25*L1^2*(tau_1.array[k,]+1)/L2)
    
    sigma.array[k] <- invgamma::rinvgamma(1,shape = shape.sigma, rate= rate.sigma)
    
    # Sample L2
    ## gamma prior
    ##chi.l2 = sum(tau_1.array[k,]+1)*L1^2/(4*sigma.array[k])
    ##psi.l2 = sum(theta.array[k,]^2*(tau_1.array[k,]+1)/tau_1.array[k,])/sigma.array[k] + 2*d2
    ##L2=GIGrvg::rgig(n = 1,lambda = r2,chi = chi.l2,psi = psi.l2)
    
    ## GIG prior
    chi.l2 = sum(tau_1.array[k,]+1)*L1^2/(4*sigma.array[k]) + chi2
    psi.l2 = sum(theta.array[k,]^2*(tau_1.array[k,]+1)/tau_1.array[k,])/sigma.array[k] + psi2
    L2=GIGrvg::rgig(n = 1,lambda = nu2,chi = chi.l2,psi = psi.l2)
    
    # Sample L1
    L1=sqrt(rgamma(1, p/2+r1 , sum(tau_1.array[k,]+1)/(8*L2*sigma.array[k])+d1 ) )
    
    penalty[k,]=c(L1,L2)
    
  }      
  acceptance.rate <- sum(acceptance[-(1:nwarm)]) / (nsim - nwarm - 1)
  return(result=list(theta = theta.array[-(1:nwarm),], tau_1 = tau_1.array[-(1:nwarm),], 
                     sigma = sigma.array[-(1:nwarm)], penalty=penalty[-(1:nwarm),],
                     acceptance.rate = acceptance.rate))
}

#############################################################
# BENEL Full Bayes (FB) chains
#############################################################
BENELF_chain <- function(x,y,nsim,nwarm=floor(nsim/2),
                         epsilon,nchain=4,L=10,
                         L1,L2,a=10,b=10,r1,d1,nu2,psi2,chi2){
  result.chain = list()
  thetalist = list()
  theta = tau_1 = sigma = acceptance = penalty = NULL
  for (ii in 1:nchain){
    result.chain[[ii]]=BENELF_base(x=x,y=y,nsim,nwarm=nwarm,
                                  epsilon=epsilon,L=L,L1=L1,L2=L2,a=10,b=10,
                                  r1=r1,d1=d1,nu2=nu2,psi2=psi2,chi2=chi2)
    thetalist[[ii]]=result.chain[[ii]]$theta
    theta = rbind(theta,result.chain[[ii]]$theta)
    tau_1 = rbind(tau_1,result.chain[[ii]]$tau_1)
    penalty = rbind(penalty,result.chain[[ii]]$penalty)
    sigma = c(sigma,result.chain[[ii]]$sigma)
    acceptance = c(acceptance,result.chain[[ii]]$acceptance.rate)
  }
  rhat = splitrhat(ndraw=( (nsim-nwarm)*nchain),theta=thetalist)
  return(list(theta = theta, tau_1 = tau_1, 
              sigma = sigma, acceptance.rate=acceptance, 
              penalty=penalty, rhat=rhat))
}

#############################################################
# BEN
#############################################################

GSIN = function(x, y) {
  
  #	This function picks the initial value of lambda1 and lambda2 
  # to be used later in the EM algorithm in GSL.
  
  x = scale(x)
  y = y - mean(y)
  n = dim(x)[1]
  p = dim(x)[2]
  
  ls = lm(y ~ 0 + x)
  sigma2hat = var(ls$resid)
  betahat = na.omit(ls$coef)
  
  result = list(lambda1 = p*sqrt(sigma2hat)/sum(abs(betahat)), lambda2 = 0.0001)
  return(result)
}


phi = function(t) {
  
  #	The function to be used in GSL
  
  return(t^(-1/2)*exp(-t))
}

GS = function(x,y,lambda1,lambda2,n.burnin = 1000,n.sampler=11000) {
  
  #	This is the Gibbs sampler for BEN. 
  # In this function, we are not updating lambda1 and lambda2. 
  
  #	n.burnin is the length of the burn-in period
  #	n.sampler is the length of the Markov chain
  #	n.lambda is the number of iterations we want for updating lambda
  
  x = scale(x)
  y = y - mean(y)
  n = dim(x)[1]
  p = dim(x)[2]
  
  #	The parameters with ".c" are the temporary ones that we use for updating.
  #	The parameters with ".p" are the recorded ones.
  
  #---	Initialization
  beta.c = rep(1, p)
  sigma2.c = 10000
  tau.c = rep(1.5, p)
  
  #---	Iteration
  beta.p = sigma2.p = tau.p = NULL
  xtx = t(x) %*% x
  ytx = t(y) %*% x
  for (iter in 1:n.sampler) {
    if (iter/1000 == as.integer(iter/1000)) {print(iter)}
    
    a.temp = lambda2 * tau.c / (tau.c - 1) + diag(xtx)
    for (j in 1:p) {
      bj.temp = sum((xtx)[-j,j] * beta.c[-j]) - 2 * ytx[j]
      beta.c[j] = rnorm(1, mean=-bj.temp/(2*a.temp[j]), sd=sqrt(sigma2.c/a.temp[j]))
    }
    
    temp.nu = sqrt(lambda1) / (2*lambda2*abs(beta.c))
    for (j in 1:p) {
      flag = 1
      while (flag) {
        temp = SuppDists::rinvGauss(n=1, nu=temp.nu[j], lambda=lambda1/(4*lambda2*sigma2.c))
        flag = (temp <= 0)
      }
      tau.c[j] = 1 + 1 / temp
    }
    
    a.temp = n/2 + p
    b.temp = 1/2*sum( (y - x%*%beta.c)^2 ) + 1/2*lambda2*sum((tau.c / (tau.c - 1)) * beta.c^2) + 1/8*lambda1^2/lambda2*sum(tau.c)
    #--- The MH algorithm
    #	Z.temp = rinvgamma(n=1, shape = a.temp, scale = b.temp)
    #	prob.temp = min(1, (gamma_inc(1/2, lambda1^2/(8*sigma2.c*lambda2)) / gamma_inc(1/2, lambda1^2/(8*Z.temp*lambda2)))^p)
    #	u.temp = runif(1)
    #	if (u.temp <= prob.temp) {sigma2.c = Z.temp} 
    
    #--- The AR algorithm
    flag.temp = 1
    while(flag.temp) {
      z.temp = MCMCpack::rinvgamma(n=1, shape = a.temp, scale = b.temp)
      u.temp = runif(1)
      if (log(u.temp) <= p*log(gamma(0.5))-p*log(gsl::gamma_inc(1/2, lambda1^2/(8*z.temp*lambda2)))) {
        sigma2.c = z.temp
        flag.temp = 0
      }
    }
    
    beta.p = rbind(beta.p, beta.c)
    tau.p = rbind(tau.p, tau.c)
    sigma2.p = c(sigma2.p, sigma2.c)
  }
  
  result = list(beta = beta.p[-(1:n.burnin),], sigma2 = sigma2.p[-(1:n.burnin)], tau = tau.p[-(1:n.burnin),])
  return(result)
}

GSL = function(x,y,lambda1,lambda2,n.burnin=1000,n.sampler=11000,itermax=10,tol=0.1) {
  
  #	This function is based on the function GS. 
  # Here we can update the penalty parameter lambda1 and lambda2. 
  # This is done by EM. 
  # The result contains all lambda's we have tried and 
  # the Markov Chain using the last lambda.
  
  #	itermax: The maximum number of EM iterations
  #	tol: The tolerance criterion for convergence
  
  n = dim(x)[1]
  p = dim(x)[2]
  flag = 1
  counter = 1
  lambda.record = c(lambda1, lambda2)
  f.record = NULL
  
  while(flag) {
    print(counter)
    print(paste("BEN f.record=",f.record))
    print(lambda.record)
    result = GS(x=x,y=y,lambda1=lambda1,lambda2=lambda2,n.burnin=n.burnin,n.sampler=n.sampler)
    
    # E step
    # conditional likelihood function 
    ## a = [lambda1, lambda2]
    f = function(a) {
      a1 = a[1]
      a2 = a[2]
      # page 168 function Q, the first two elements
      temp1 = p*log(a1) - 
        p*mean(gsl::gamma_inc(1/2, a1^2/(8*result$sigma2*a2)))
      temp2 = 0
      # page 168 function Q, the third element
      for (j in 1:p) {
        temp2 = temp2 + mean(result$tau[,j] / (result$tau[,j]-1) * result$beta[,j]^2 / result$sigma2)
      }
      temp2 = -a2/2 * temp2
      # page 168 function Q, the last element
      temp3 = 0
      for (j in 1:p) {
        temp3 = temp3 + mean(result$tau[,j] / result$sigma2)
      }
      temp3 = -a1^2/(8*a2) * temp3
      # return the value of function Q
      return(- temp1 - temp2 - temp3)
    }
    
    # gradient functions for M-step
    ## a is same as above
    df = function(a) {
      a1 = a[1]
      a2 = a[2]
      # page 168, gradient of lambda1
      temp1 = p/a1 + p*a1/(4*a2)*mean(gsl::gamma_inc(1/2, a1^2/(8*result$sigma2*a2))^(-1) * phi(a1^2/(8*result$sigma2*a2)) / result$sigma2)
      temp2 = 0
      for (j in 1:p) {
        temp2 = temp2 + mean(result$tau[,j] / result$sigma2)
      }
      temp2 = -a1/(4*a2) * temp2
      result.1 = temp1 + temp2
      
      # page 168, gradient of lambda2
      temp1 = -p*a1^2/(8*a2^2)*mean(gsl::gamma_inc(1/2, a1^2/(8*result$sigma2*a2))^(-1) * phi(a1^2/(8*result$sigma2*a2)) / result$sigma2)
      temp2 = 0
      for (j in 1:p) {
        temp2 = temp2 + mean(result$tau[,j] / (result$tau[,j]-1) * result$beta[,j]^2 / result$sigma2)
      }
      temp2 = -1/2 * temp2
      temp3 = 0
      for (j in 1:p) {
        temp3 = temp3 + mean(result$tau[,j] / result$sigma2)
      }
      temp3 = a1^2/(8*a2^2) * temp3
      result.2 = temp1 + temp2 + temp3
      return(c(-result.1, -result.2))
    }
    
    lambda.current = c(lambda1, lambda2)
    # M step
    lambda.update = constrOptim(lambda.current, f, grad=df, ui=diag(c(1,1)), ci=.Machine$double.eps*c(1,1))$par
    # lambda.update = constrOptim(lambda.current, f, grad=NULL,ui=diag(c(1,1)), ci=.Machine$double.eps*c(1,1))$par
    flag = (counter < itermax) & 
      ( sum((lambda.current - lambda.update)^2) >= tol )
    convergence = ( sum((lambda.current - lambda.update)^2) < tol )
    lambda.record = rbind(lambda.record, lambda.update)
    f.current = f(lambda.current)
    f.record = c(f.record, f.current)
    lambda1 = lambda.update[1]
    lambda2 = lambda.update[2]
    counter = counter + 1
  }
  out = list(beta = result$beta, sigma2 = result$sigma2, 
             tau = result$tau, lambda = lambda.record, 
             convergence = convergence, f=f.record)
  return(out)
}

GSIN = function(x, y) {
  
  #	This function picks the initial value of lambda1 and lambda2 
  # to be used later in the EM algorithm in GSL.
  
  x = scale(x)
  y = y - mean(y)
  n = dim(x)[1]
  p = dim(x)[2]
  
  ls = lm(y ~ 0 + x)
  sigma2hat = var(ls$resid)
  betahat = na.omit(ls$coef)
  
  result = list(lambda1 = p*sqrt(sigma2hat)/sum(abs(betahat)), lambda2 = 0.0001)
  return(result)
}

#############################################################
# Elastic net
#############################################################

ENLV = function(x.train, y.train,lambda2 = c(0,0.01,0.1,1,10,100,1000), mode = "fraction", max.steps = 200) {
  
  #	The elastic net algorithm with L1 penalty lambda1 chosen by the validation data set.
  #	The loss function is MSE. The default abscissa is on the fractions seq(0,1,length=100). This 
  #	is easy to control. If we want to use lambda1, it is not straightforward how to
  #	decide the abscissa. The L2 penalty term lambda2 is chosen as the best one in the pool.
  
  n.train = dim(x.train)[1]
  p = dim(x.train)[2]
  s = err = NULL	# to record, for each lambda2, the L1 norm fraction chosen by the validation set 
  # and the corresponding prediction error.
  if (mode == "fraction") 
    abscissa = seq(0,1,length=100)
  if (mode == "step")
    abscissa = 1:max.steps
  for (index2 in 1:length(lambda2)) {
    temp = elasticnet::enet(x.train, y.train, lambda=lambda2[index2])
    pred.temp = predict(temp, x.train, type="fit", s=abscissa,mode=mode)$fit
    mse = apply((as.numeric(y.train) - pred.temp)^2, 2, mean)
    s = c(s, abscissa[which(mse == min(mse))])
    err = c(err, min(mse))
  }
  lambda2.chosen = lambda2[which(err == min(err))[1]]
  s.chosen = s[which(err == min(err))[1]]
  
  temp = elasticnet::enet(scale(rbind(x.train)),c(y.train) - mean(c(y.train)),lambda=lambda2.chosen, max.steps = max.steps)
  if (lambda2.chosen == 0) temp = lars::lars(scale(rbind(x.train)),c(y.train) - mean(c(y.train)), max.steps = max.steps)
  result.en = predict(temp, s = s.chosen, type = "coefficients", mode = mode)
  result = list(s = s.chosen, lambda2 = lambda2.chosen, beta = result.en$coef, err = min(err))
  return(result)
}

#############################################################
# LARS
#############################################################

LARSL = function(x,y,K = 10) {
  abscissa = seq(0,1,length=100)
  temp = lars::cv.lars(x,y,fraction=abscissa,K=K, plot.it=F)
  s = abscissa[which(temp$cv == min(temp$cv))]
  err = min(temp$cv)
  
  temp = lars(x,y)
  result.lars = predict(temp, s = s, type = "coefficients", mode = "fraction")
  result = list(s = s, beta = result.lars$coef, err = err)
  return(result)
}

LARSLV = function(x.train, y.train, x.validation, y.validation) {
  abscissa = seq(0,1,length=100)
  temp = lars::lars(x.train,y.train)
  pred.temp = predict(temp, x.validation, type="fit", s=abscissa,mode="fraction")$fit
  mse = apply((as.numeric(y.validation) - pred.temp)^2, 2, mean)
  s = abscissa[which(mse == min(mse))]
  err = min(mse)
  
  temp = lars::lars(scale(rbind(x.train, x.validation)),c(y.train, y.validation) - mean(c(y.train,y.validation)))
  result.lars = predict(temp, s = s, type = "coefficients", mode = "fraction")
  result = list(s = s, beta = result.lars$coef, err = err)
  return(result)
}

#############################################################
# BL
#############################################################

GS_BL = function(x,y,lambda1,n.burnin = 1000,n.sampler=11000) {
  
  #	This is the Gibbs sampler for BL. In this function, we are not updating lambda1.
  
  #	n.burnin is the length of the burn-in period
  #	n.sampler is the length of the Markov chain
  #	n.lambda is the number of iterations we want for updating lambda
  
  x = scale(x)
  y = y - mean(y)
  n = dim(x)[1]
  p = dim(x)[2]
  
  #	The parameters with ".c" are the temporary ones that we use for updating.
  #	The parameters with ".p" are the recorded ones.
  
  #---	Initialization
  beta.c = rep(1, p)
  sigma2.c = 1
  tau2.c = rep(1.5, p)
  iD = diag(1/tau2.c)
  
  #---	Iteration
  beta.p = sigma2.p = tau2.p = NULL
  xtx = t(x) %*% x
  xty = t(x) %*% y
  
  for (iter in 1:n.sampler) {
    if (iter/1000 == as.integer(iter/1000)) {print(iter)}
    
    A = xtx + iD
    iA = solve(A)
    beta.c = MASS::mvrnorm(n=1, mu = iA%*%xty, Sigma=sigma2.c*iA)
    
    temp.nu = lambda1*sqrt(sigma2.c) / abs(beta.c)
    for (j in 1:p) {
      flag = 1
      while (flag) {
        temp = SuppDists::rinvGauss(n=1, nu=temp.nu[j], lambda=lambda1^2)
        flag = (temp <= 0)
      }
      tau2.c[j] = 1 / temp
    }
    iD = diag(1/tau2.c)
    
    temp =  1/2 * ( sum(y^2) - 2*sum(xty*beta.c) + sum(beta.c*(A%*%beta.c)) )
    sigma2.c = MCMCpack::rinvgamma(n=1, shape = (n+p)/2, scale = temp)
    
    beta.p = rbind(beta.p, beta.c)
    tau2.p = rbind(tau2.p, tau2.c)
    sigma2.p = c(sigma2.p, sigma2.c)
  }
  
  result = list(beta = beta.p[-(1:n.burnin),], sigma2 = sigma2.p[-(1:n.burnin)], tau2 = tau2.p[-(1:n.burnin),])
  return(result)
}

GSL_BL = function(x,y,lambda1,n.burnin=1000,n.sampler=11000,itermax=10,tol=0.05) {
  
  n = dim(x)[1]
  p = dim(x)[2]
  flag = 1
  counter = 1
  lambda.record = lambda1
  
  while(flag) {
    print(counter)
    print(lambda.record)
    result = GS_BL(x,y,lambda1,n.burnin,n.sampler)
    
    lambda.current = lambda1
    lambda.update = sqrt(2*p / (sum(apply(result$tau2, 2, mean))))
    
    flag = (counter < itermax) & ( sum((lambda.current - lambda.update)^2) >= tol )
    convergence = ( sum((lambda.current - lambda.update)^2) < tol )
    lambda.record = c(lambda.record, lambda.update)
    lambda1 = lambda.update
    counter = counter + 1
  }
  out = list(beta = result$beta, sigma2 = result$sigma2, tau2 = result$tau2, lambda = lambda.record, convergence = convergence)
  return(out)
}

