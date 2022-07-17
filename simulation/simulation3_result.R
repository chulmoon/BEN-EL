library(tidyverse)
library(viridis)

# l1
load("./simulation/sim_3_BENEL_l1.Rdata")

# r1: 0.25 to 10 by 0.25
r1=seq(0.25,10,by=0.25)
# d1: 0.25 to 10 by 0.25
d1=seq(0.25,10,by=0.25)
hyp=expand.grid(r1,d1)
n.sim=nrow(hyp)

conv = NULL
theta = matrix(NA,n.sim,8)
for (ii in 1:n.sim) {
	theta[ii,]=result.benel[[ii]]$beta.benel
	conv=conv+sum(result.benel[[ii]]$rhat>1.01)
}
print(conv)

theta.l1 = data.frame(cbind(hyp,theta))

ggplot(theta.l1)+
	geom_raster(aes(Var1,Var2,fill=X1)) +
	labs(x=expression(r[1]),y=expression(delta[1]),fill=expression(theta[1])) +
	theme_minimal()+
	scale_fill_viridis()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank())

ggplot(theta.l1,aes(Var1,Var2,fill=X4))+
  geom_raster() +
  labs(x=expression(r[1]),y=expression(delta[1]),fill=expression(theta[4])) +
  theme_minimal()+
  scale_fill_viridis()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


# l2
load("./simulation/sim_3_BENEL_l2.Rdata")

# nu: -5 to 5 by 0.25
nu2=seq(-5,5,by=0.25)
# psi: 0.25 to 10 by 0.25
psi2=seq(0.25,10,by=0.25)

hyp=expand.grid(nu2,psi2)
n.sim=nrow(hyp)

conv = NULL
theta = matrix(NA,n.sim,8)
for (ii in 1:n.sim) {
  theta[ii,]=result.benel[[ii]]$beta.benel
  conv=conv+sum(result.benel[[ii]]$rhat>1.01)
}
print(conv)

theta.l2 = data.frame(cbind(hyp,theta))

ggplot(theta.l2,aes(Var1,Var2,fill=X1))+
	geom_raster() +
	labs(x=expression(nu[2]),y=expression(psi[2]),fill=expression(theta[1])) +
	theme_minimal()+
	scale_fill_viridis()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank())

ggplot(theta.l2,aes(Var1,Var2,fill=X4))+
	geom_raster() +
	labs(x=expression(nu[2]),y=expression(psi[2]),fill=expression(theta[4])) +
	theme_minimal()+
	scale_fill_viridis()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank())
