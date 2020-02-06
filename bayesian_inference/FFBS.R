# Petris, G. & Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
# p. 166
library(dlm)

# Example 1
# a1 <- 2
# b1 <- 0.0001
# a2 <- 2
# b2 <- 0.0001
# 
# psi1 <- 1
# psi2 <- 2
# mod_level <- dlmModPoly(order = 1, dV = 1/psi1, dW = 1/psi2)
# mc <- 30000
# 
# psi1_save <- numeric(mc)
# psi2_save <- numeric(mc)
# theta_save <- matrix(0, nrow = n+1, ncol=mc)
# 
# sh1 <- a1+n/2
# sh2 <- a2+n/2
# 
# set.seed(10)
# now <- Sys.time()
# for (it in 1:mc){
#   ## draw the states: FFBS
#   filt <- dlmFilter(Nile, mod_level)
#   theta <- dlmBSample(filt)
#   
#   browser()
#   ## draw observation precision psi1
#   rate <- b1+crossprod(Nile-theta[-1])/2
#   psi1 <- rgamma(1, shape = sh1, rate = rate)
#   
#   ## draw system precision psi2
#   rate <- b2+crossprod(theta[-1]-theta[-n])/2
#   psi2 <- rgamma(1, shape = sh2, rate = rate)
#   
#   ## update and save
#   V(mod_level) <- 1/psi1
#   W(mod_level) <- 1/psi2
#   
#   theta_save[,it] <- c(theta)
#   psi1_save[it] <- psi1
#   psi2_save[it] <- psi2
#   
# }
# 
# Sys.time()-now
# 
# 
# 
# plot.ts(cumsum(psi1_save)/1:mc)
# plot.ts(cumsum(psi2_save)/1:mc)
# 
# mean(psi1_save[-(1:5000)])
# mean(psi2_save[-(1:5000)])
# 
# theta_hat <- rowMeans(theta_save[-1, -(1:5000)])
# plot.ts(cbind(c(Nile), theta_hat), plot.type = "single")

# Example 2
y <- log(c(AirPassengers))
plot.ts(y)
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

mc <- 30000

psiy_save <- numeric(mc)
psi1_save <- numeric(mc)
psi2_save <- numeric(mc)
psi3_save <- numeric(mc)
theta_save <- array(0, dim=c(n+1, nrow(W), mc))

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
  
  SS_theta <- diag(crossprod(theta[-1,]-theta[-n,]%*%t(GG)))
  ## draw system precision psi2
  rate <- b_psi1 + SS_theta[1]/2
  psi1 <- rgamma(1, shape = sh1, rate = rate)
  
  ## draw system precision psi2
  rate <- b_psi2 + SS_theta[2]/2
  psi2 <- rgamma(1, shape = sh2, rate = rate)
  
  ## draw system precision psi2
  rate <- b_psi3 + SS_theta[3]/2
  psi3 <- rgamma(1, shape = sh3, rate = rate)
  
  ## update and save
  V(mod_level) <- 1/psiy
  diag(W(mod_level))[1:3] <- c(1/psi1, 1/psi2, 1/psi3)
  theta_save[,,it] <- theta
  psiy_save[it] <- psiy
  psi1_save[it] <- psi1
  psi2_save[it] <- psi2
  psi3_save[it] <- psi3
  
}
Sys.time()-now


plot.ts(cumsum(psiy_save)/1:mc)
plot.ts(cumsum(psi1_save)/1:mc)
plot.ts(cumsum(psi2_save)/1:mc)
plot.ts(cumsum(psi3_save)/1:mc)

psiy_hat <- mean(psiy_save[-(1:20000)])
psi1_hat <- mean(psi1_save[-(1:20000)])
psi2_hat <- mean(psi2_save[-(1:20000)])
psi3_hat <- mean(psi3_save[-(1:20000)])

theta_hat <- apply(theta_save[,, -(1:20000)], MARGIN=c(1,2), FUN=mean)

plot.ts(cbind(y, (theta_hat%*%t(FF))[-1,]), plot.type = "single", col=c("black", "red"))
plot.ts(cbind(y, theta_hat[-1,1]), plot.type = "single", col=c("black", "red"))
plot.ts(theta_hat[-1,3])

# Diagnostics
C <- array(0, c(13, 13, 145))
for (i in 1:145){
  C[,,i] <- filt$U.C[[i]] %*% diag(filt$D.C[i,])^2 %*% t(filt$U.C[[i]])
}

thetaF <- apply(theta_save[,,], MARGIN=3, FUN = function(x){x%*%t(FF)})
SS_y<- crossprod(thetaF)
plot.ts(diag(SS_y)[-(1:100)])
median(diag(SS_y)[-(1:100)])

thetaG <- array(0, c(144, 13, mc))
for (i in 1:mc){
  thetaG[,,i] <- theta_save[-1,,i] - theta_save[-145,,i]%*%t(GG)
}

SS_theta1 <- array(0, mc)
for (i in 1:mc){
  SS_theta1[i] <- crossprod(thetaG[,1,i])
}
plot.ts(SS_theta1[-(1:100)])
median(SS_theta1[-(1:100)])

SS_theta2 <- array(0, mc)
for (i in 1:mc){
  SS_theta2[i] <- crossprod(thetaG[,2,i])
}
plot.ts(SS_theta2[-(1:100)])

SS_theta3 <- array(0, mc)
for (i in 1:mc){
  SS_theta3[i] <- crossprod(thetaG[,3,i])
}
plot.ts(SS_theta3[-(1:100)])








# Results from SVD-Kalman-Algorithm
mod <- dlmModPoly(order = 2)+dlmModSeas(frequency = 12)
V(mod) <- 1.0/22814.31
diag(W(mod))[1:3] <- c(1.0/2472.64,1.0/55947.63,1.0/4309.207)

filt <- dlmFilter(y, mod)
samp <- dlmBSample(filt)


plot.ts(cbind(y, filt$m[-1,1]), plot.type = "single",col=c("black", "red"))
filt$m%*%t(FF(mod))
filt$m%*%c(1, rep(0, 12))

samp <- dlmBSample(filt)
plot.ts(cbind(y, samp[-1,1]), plot.type = "single",col=c("black", "red"))
plot.ts(cbind(y, samp[-1,]%*%t(FF(mod))), plot.type = "single",col=c("black", "red"))


tmp1 <- svd(X)
X
tmp1$u%*%diag(tmp1$d)%*%t(tmp1$v)
all.equal(tmp1$u%*%diag(tmp1$d)%*%t(tmp1$v), X)

tmp2 <- La.svd(X)
tmp2$u%*%diag(tmp2$d)%*%tmp2$v
all.equal(tmp2$u%*%diag(tmp2$d)%*%tmp2$vt, X)



V <- 1.0/22814.31
W <- matrix(0, 13, 13)
diag(W)[1:3] <- c(1.0/2472.64,1.0/55947.63,1.0/4309.207)

tmp <- La.svd(V,nu=0)
Uv <- t(tmp$vt);
Dv <- sqrt(tmp$d)
Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
sqrtVinv <- Dv.inv * tmp$vt

svdW <- La.svd(W,nu=0)
sqrtW <- sqrt(svdW$d) * svdW$vt 
