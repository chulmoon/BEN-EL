# setup
n.sim = 100
p=8
method = 3
d=0.5

datnum=c("normal_n50","normal_n100","normal_n200",
         "normal2_n50","normal2_n100","normal2_n200",
         "st_n50","st_n100","st_n200")

#####################################################################
# BEN
#####################################################################

medse.ben = matrix(NA,length(datnum),2) # median of residual and sd

for (sim in 1:length(datnum)) {
	load(paste0("./simulation/sim_1_",datnum[sim],"_BEN.Rdata"))
	load(paste0("./simulation/data/sim_1_",datnum[sim],".Rdata"))
	
	resid.ben = rep(NA,n.sim)	# to record the median squared errors
	bsm.ben = matrix(NA,n.sim,p)		# short for beta.select.m.benel
	beta.ben.cor = matrix(NA,n.sim,p)
	
	for (i in 1:n.sim) {
		x.train = result.ben[[i]]$x.train
		x.test = simdata[[i]]$x.test
		y.train = result.ben[[i]]$y.train
		y.test = result.ben[[i]]$y.test
		
		beta.ben = result.ben[[i]]$beta.ben

		n.train=nrow(x.train)
		n.test=nrow(x.test)
		
		x.test.bar = scale(x.test)
		y.test.bar = y.test - mean(y.test)
		
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
		if (method == 1) bsm.ben = rbind(bsm.ben, (intervals[,1] < 0 & intervals[,2] > 0))
		if (method == 3) bsm.ben[i,] = (interval.prob > d)
		
		beta.upt = beta.temp %*% (diag(1-bsm.ben[i,]))
		beta.ben.cor[i,] = apply(beta.upt,2,median)
		
		resid.ben[i] = mean((y.test.bar - x.test.bar %*% beta.ben.cor[i,])^2)
	}
	
	median.bsample = NULL
	for (b in 1:1000) {
		bsample = sample(resid.ben, replace = T)
		median.bsample = c(median.bsample, median(bsample))
	}
	sd.ben = sqrt(var(median.bsample))
	
	medse.ben[sim,]=c(median(resid.ben),sd.ben)
}
medse.ben

#####################################################################
# BL
#####################################################################
medse.bl = matrix(NA,length(datnum),2) # median of residual and sd

for (sim in 1:length(datnum)) {
  load(paste0("./simulation/sim_1_",datnum[sim],"_BL.Rdata"))
  load(paste0("./simulation/data/sim_1_",datnum[sim],".Rdata"))
  
	resid.bl = rep(NA,n.sim)	# to record the median squared errors
	bsm.bl = matrix(NA,n.sim,p)		# short for beta.select.m.blel
	beta.bl.cor = matrix(NA,n.sim,p)
	
	for (i in 1:n.sim) {
		x.train = result.bl[[i]]$x.train
		x.test = simdata[[i]]$x.test
		y.train = result.bl[[i]]$y.train
		y.test = result.bl[[i]]$y.test
		
		beta.bl = result.bl[[i]]$beta.bl
		
		n.train=nrow(x.train)
		n.test=nrow(x.test)
		
		x.test.bar = scale(x.test)
		y.test.bar = y.test - mean(y.test)
		
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
		if (method == 1) bsm.bl = rbind(bsm.bl, (intervals[,1] < 0 & intervals[,2] > 0))
		if (method == 3) bsm.bl[i,] = (interval.prob > d)
		
		beta.upt = beta.temp %*% (diag(1-bsm.bl[i,]))
		beta.bl.cor[i,] = apply(beta.upt,2,median)
		
		resid.bl[i] = mean((y.test.bar - x.test.bar %*% beta.bl.cor[i,])^2)
	}
	
	median.bsample = NULL
	for (b in 1:1000) {
		bsample = sample(resid.bl, replace = T)
		median.bsample = c(median.bsample, median(bsample))
	}
	sd.bl = sqrt(var(median.bsample))
	
	medse.bl[sim,]=c(median(resid.bl),sd.bl)
}
medse.bl


#####################################################################
# BENEL
#####################################################################
medse.benel = matrix(NA,length(datnum),2) # median of residual and sd

for (sim in 1:length(datnum)) {
  load(paste0("./simulation/sim_1_",datnum[sim],"_BENEL.Rdata"))
  conv=NULL
  for (ii in 1:100){
    conv=conv+sum(result.benel[[ii]]$rhat>1.01)
  }
  print(conv)
}

for (sim in 1:length(datnum)) {
  load(paste0("./simulation/sim_1_",datnum[sim],"_BENEL.Rdata"))
  load(paste0("./simulation/data/sim_1_",datnum[sim],".Rdata"))
  
	resid.benel = rep(NA,n.sim)	# to record the median squared errors
	bsm.benel = matrix(NA,n.sim,p)		# short for beta.select.m.benelel
	beta.benel.cor = matrix(NA,n.sim,p)
	
	for (i in 1:n.sim) {
		x.train = result.benel[[i]]$x.train
		x.test = simdata[[i]]$x.test
		y.train = result.benel[[i]]$y.train
		y.test = result.benel[[i]]$y.test
		
		beta.benel = result.benel[[i]]$beta.benel
		
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
		if (method == 1) bsm.benel = rbind(bsm.benel, (intervals[,1] < 0 & intervals[,2] > 0))
		if (method == 3) bsm.benel[i,] = (interval.prob > d)
		
		beta.upt = beta.temp %*% (diag(1-bsm.benel[i,]))
		beta.benel.cor[i,] = apply(beta.upt,2,median)
		
		resid.benel[i] = mean((y.test.bar - x.test.bar %*% beta.benel.cor[i,])^2)
	}
	
	median.bsample = NULL
	for (b in 1:1000) {
		bsample = sample(resid.benel, replace = T)
		median.bsample = c(median.bsample, median(bsample))
	}
	sd.benel = sqrt(var(median.bsample))
	
	medse.benel[sim,]=c(median(resid.benel),sd.benel)
	
	print(apply(bsm.benel,2,sum))
}
medse.benel

#####################################################################
# EN
#####################################################################
medse.en = matrix(NA,length(datnum),2) # median of residual and sd

for (sim in 1:length(datnum)) {
  load(paste0("./simulation/sim_1_",datnum[sim],"_EN.Rdata"))
  median.bsample = NULL
	for (b in 1:1000) {
		bsample = sample(resid.en, replace = T)
		median.bsample = c(median.bsample, median(bsample))
	}
	sd.en = sqrt(var(median.bsample))
	
	medse.en[sim,]=c(median(resid.en),sd.en)
}
medse.en

bsm.en = (beta.en!=0)

#####################################################################
# LADL
#####################################################################
medse.ladl = matrix(NA,length(datnum),2) # median of residual and sd

for (sim in 1:length(datnum)) {
  load(paste0("./simulation/sim_1_",datnum[sim],"_LADL.Rdata"))
  median.bsample = NULL
  for (b in 1:1000) {
    bsample = sample(resid.ladlasso, replace = T)
    median.bsample = c(median.bsample, median(bsample))
  }
  sd.ladlasso = sqrt(var(median.bsample))
  
  medse.ladlasso[sim,]=c(median(resid.ladlasso),sd.ladlasso)
  bsm.ladlasso = (beta.ladlasso!=0)
}
medse.ladlasso

#############################################
# number of nonzero coefficients
#############################################
apply(bsm.ben,2,sum)
apply(bsm.benel,2,sum)
apply(bsm.bl,2,sum)
100-apply(bsm.en,2,sum)
100-apply(bsm.ladlasso,2,sum)