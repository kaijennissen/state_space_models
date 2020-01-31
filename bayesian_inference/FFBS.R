# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
# p. 166
library(dlm)

a1 <- 2
b1 <- 0.0001
a2 <- 2
b2 <- 0.0001

psi1 <- 1
psi2 <- 2
mod_level <- dlmModPoly(order = 1, dV = 1/psi1, dW = 1/psi2)
mc <- 30000

psi1_save <- numeric(mc)
psi2_save <- numeric(mc)
theta_save <- matrix(0, nrow = n+1, ncol=mc)

sh1 <- a1+n/2
sh2 <- a2+n/2

set.seed(10)
now <- Sys.time()
for (it in 1:mc){
  ## draw the states: FFBS
  filt <- dlmFilter(Nile, mod_level)
  theta <- dlmBSample(filt)
  
  browser()
  ## draw observation precision psi1
  rate <- b1+crossprod(Nile-theta[-1])/2
  psi1 <- rgamma(1, shape = sh1, rate = rate)
  
  ## draw system precision psi2
  rate <- b2+crossprod(theta[-1]-theta[-n])/2
  psi2 <- rgamma(1, shape = sh2, rate = rate)
  
  ## update and save
  V(mod_level) <- 1/psi1
  W(mod_level) <- 1/psi2
  
  theta_save[,it] <- c(theta)
  psi1_save[it] <- psi1
  psi2_save[it] <- psi2
  
}

Sys.time()-now



plot.ts(cumsum(psi1_save)/1:mc)
plot.ts(cumsum(psi2_save)/1:mc)

mean(psi1_save[-(1:5000)])
mean(psi2_save[-(1:5000)])

theta_hat <- rowMeans(theta_save[-1, -(1:5000)])
plot.ts(cbind(c(Nile), theta_hat), plot.type = "single")















