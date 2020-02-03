# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
# p. 166
library(dlm)

# Example 1
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

# Example 2
y <- c(AirPassengers)
n <- length(y)
ay <- 2
by <- 0.0001
a_psi1 <- 2
b_psi1 <- 0.0001
a_psi2 <- 2
b_psi2 <- 0.0001
a_psi3 <- 2
b_psi3 <- 0.0001

psiy <- 1
psi1 <- 1
psi2 <- 1
psi3 <- 1
V <- 1/psiy
W <- diag(rep(0, 13))
diag(W)[1:3] <- c(1/psi1, 1/psi2, 1/psi3)
mod_level <- dlmModPoly(order = 2)+dlmModSeas(frequency = 12)
V(mod_level) <- V
W(mod_level) <- W

mc <- 30

psiy_save <- numeric(mc)
psi1_save <- numeric(mc)
psi2_save <- numeric(mc)
psi3_save <- numeric(mc)
theta_save <- matrix(0, nrow = n+1, ncol=mc)

shy <- ay+n/2
sh1 <- a_psi1+n/2
sh2 <- a_psi2+n/2
sh3 <- a_psi3+n/2

set.seed(10)
now <- Sys.time()
for (it in 1:mc){
  
  FF <- FF(mod_level)
  GG <- GG(mod_level)
  ## draw the states: FFBS
  filt <- dlmFilter(y, mod_level)
  theta <- dlmBSample(filt)
  
  ## draw observation precision psiy
  rate <- by + crossprod(y-theta[-1,]%*%t(FF))/2
  psiy <- rgamma(1, shape = sh1, rate = rate)
  
  SS_theta <- colSums(theta[-1,]-theta[-n,]%*%t(GG))
  ## draw system precision psi2
  rate <- b_psi1 + SS_theta[1]/2
  psi1 <- rgamma(1, shape = sh2, rate = rate)
  
  ## draw system precision psi2
  rate <- b_psi2 + SS_theta[2]/2
  psi2 <- rgamma(1, shape = sh2, rate = rate)
  
  ## draw system precision psi2
  rate <- b_psi3 + SS_theta[3]/2
  psi3 <- rgamma(1, shape = sh2, rate = rate)
  
  ## update and save
  V(mod_level) <- 1/psiy
  diag(W(mod_level))[1:3] <- c(1/psi1, 1/psi2, 1/psi3)
  theta_save[,it] <- c(theta)
  psiy_save[it] <- psi1
  psi1_save[it] <- psi2
  psi2_save[it] <- psi2
  psi3_save[it] <- psi2
  
}

Sys.time()-now



plot.ts(cumsum(psi1_save)/1:mc)
plot.ts(cumsum(psi2_save)/1:mc)

mean(psi1_save[-(1:5000)])
mean(psi2_save[-(1:5000)])

theta_hat <- rowMeans(theta_save[-1, -(1:5000)])
plot.ts(cbind(c(Nile), theta_hat), plot.type = "single")













