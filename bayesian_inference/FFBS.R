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

simulate_y <- function(n, psiy, psi1, psi2, psi3){
  
  # model
  FF <- matrix(c(1, 0, 1,rep(0,10)), nrow=1)
  GG <- matrix(0, 13, 13)
  GG[1,1] <- 1
  GG[2,2] <- 1
  GG[1,2] <- 1
  GG[3, 3:13] <- -1
  GG[4:13, 3:12] <- diag(10)

  V <- matrix(1/psiy)
  W <- diag(c(1/psi1, 1/psi2, 1/psi3, rep(0,10)))
  
  seas <-  3*sin(seq(-1,1, length.out=12))[-12]
  theta0 <-  matrix(c(4,0.1, seas))
  
  # simulate
  y <- vector("numeric", 144)
  
  theta <-  theta0
  tmp <- La.svd(W)
  sqrtW <- diag(sqrt(tmp$d)) %*% tmp$vt
  tmp <- La.svd(V)
  sqrtV <- sqrt(tmp$d) %*% tmp$vt
  
  for (i in 1:n){
    theta <- GG%*%theta+sqrtW%*%rnorm(dim(sqrtW)[1])
    y[i] <- FF%*%theta + sqrtV%*%rnorm(dim(sqrtV)[1])
  }
  
  return(y)
  
}


# Example 2
#y <- log(c(AirPassengers))
#plot.ts(y)
ffbs <- function(y, mc) {
  n <- length(y)
  ay <- 2
  by <- 0.0001
  a_psi1 <- 2
  b_psi1 <- 0.0001
  a_psi2 <- 2
  b_psi2 <- 0.0001
  a_psi3 <- 2
  b_psi3 <- 0.0001
  
  set.seed(123)
  psiy <- 1
  psi1 <- 1
  psi2 <- 1
  psi3 <- 1
  V <- 1 / psiy
  W <- diag(rep(0, 13))
  diag(W)[1:3] <- c(1 / psi1, 1 / psi2, 1 / psi3)
  mod_level <- dlmModPoly(order = 2) + dlmModSeas(frequency = 12)
  V(mod_level) <- V
  W(mod_level) <- W
  
  
  
  psiy_save <- numeric(mc)
  psi1_save <- numeric(mc)
  psi2_save <- numeric(mc)
  psi3_save <- numeric(mc)
  theta_save <- array(0, dim = c(n + 1, nrow(W), mc))
  
  shy <- ay + n / 2
  sh1 <- a_psi1 + n / 2
  sh2 <- a_psi2 + n / 2
  sh3 <- a_psi3 + n / 2
  
  set.seed(10)
  
  for (it in 1:mc) {
    FF <- FF(mod_level)
    GG <- GG(mod_level)
    ## draw the states: FFBS
    filt <- dlmFilter(y, mod_level)
    theta <- dlmBSample(filt)
    
    ## draw observation precision psiy
    rate <- by + crossprod(y - theta[-1, ] %*% t(FF)) / 2
    psiy <- rgamma(1, shape = sh1, rate = rate)
    
    SS_theta <- diag(crossprod(theta[-1, ] - theta[-n, ] %*% t(GG)))
    ## draw system precision psi2
    rate <- b_psi1 + SS_theta[1] / 2
    psi1 <- rgamma(1, shape = sh1, rate = rate)
    
    ## draw system precision psi2
    rate <- b_psi2 + SS_theta[2] / 2
    psi2 <- rgamma(1, shape = sh2, rate = rate)
    
    ## draw system precision psi2
    rate <- b_psi3 + SS_theta[3] / 2
    psi3 <- rgamma(1, shape = sh3, rate = rate)
    
    ## update and save
    V(mod_level) <- 1 / psiy
    diag(W(mod_level))[1:3] <- c(1 / psi1, 1 / psi2, 1 / psi3)
    theta_save[, , it] <- theta
    psiy_save[it] <- psiy
    psi1_save[it] <- psi1
    psi2_save[it] <- psi2
    psi3_save[it] <- psi3
    
    
  }
  return(list(psiy_save, psi1_save, psi2_save, psi3_save))
}




y <- simulate_y(n=250, psiy=200, psi1=10, psi2=1e5, psi3=1e4)
plot.ts(y)

now <- Sys.time()
resu <- ffbs(y, 30000)
Sys.time()-now


plot.ts(cumsum(psiy_save)/1:mc)
plot.ts(cumsum(psi1_save)/1:mc)
plot.ts(cumsum(psi2_save)/1:mc)
plot.ts(cumsum(psi3_save)/1:mc)

# psiy_hat = 21881.95; 22814.31; 23091.88
# psi1_hat = 2669.815; 2472.64; 2468.197
# psi2_hat = 56132.51; 55947.63; 56553.41
# psi3_hat = 4565.963; 4309.207; 4571.657
psiy_hat <- mean(psiy_save[-(1:10000)])
psi1_hat <- mean(psi1_save[-(1:10000)])
psi2_hat <- mean(psi2_save[-(1:10000)])
psi3_hat <- mean(psi3_save[-(1:10000)])

theta_hat <- apply(theta_save[,, -(1:20000)], MARGIN=c(1,2), FUN=mean)
theta_max <- apply(theta_save[,, -(1:20000)], MARGIN=c(1,2), FUN=max)
theta_min <- apply(theta_save[,, -(1:20000)], MARGIN=c(1,2), FUN=min)

plot.ts(cbind(y, (theta_hat%*%t(FF))[-1,]), plot.type = "single", col=c("black", "red"))
plot.ts(cbind(y, theta_hat[-1,1], theta_max[-1,1], theta_min[-1,1] ), plot.type = "single", col=c("black", "blue", "red", "red"))
plot.ts(theta_hat[-1,3])

# Diagnostics
C <- array(0, c(13, 13, 145))
for (i in 1:145){
  C[,,i] <- filt$U.C[[i]] %*% diag(filt$D.C[i,])^2 %*% t(filt$U.C[[i]])
}

thetaF <- apply(theta_save[-1,,], MARGIN=3, FUN = function(x){y-x%*%t(FF)})
SS_y<- crossprod(thetaF)
plot.ts(diag(SS_y)[1:100])
plot.ts(diag(SS_y)[1001:30000])
mean(diag(SS_y)[5001:30000])

thetaG <- array(0, c(144, 13, mc))
for (i in 1:mc){
  thetaG[,,i] <- theta_save[-1,,i] - theta_save[-145,,i]%*%t(GG)
}

SS_theta1 <- array(0, mc)
for (i in 1:mc){
  SS_theta1[i] <- crossprod(thetaG[,1,i])
}
plot.ts(SS_theta1[-(1:100)])
plot.ts(SS_theta1[1:100])
plot.ts(SS_theta1[1001:10000])
mean(SS_theta1[1001:30000])


SS_theta2 <- array(0, mc)
for (i in 1:mc){
  SS_theta2[i] <- crossprod(thetaG[,2,i])
}
plot.ts(SS_theta2[-(1:100)])
plot.ts(SS_theta2[1:100])
plot.ts(SS_theta2[1001:10000])
mean(SS_theta2[1001:30000])

SS_theta3 <- array(0, mc)
for (i in 1:mc){
  SS_theta3[i] <- crossprod(thetaG[,3,i])
}
plot.ts(SS_theta3[-(1:100)])
plot.ts(SS_theta3[1:100])
plot.ts(SS_theta3[1001:10000])
mean(SS_theta3[1001:30000])




# Results from SVD-Kalman-Algorithm
mod <- dlmModPoly(order = 2)+dlmModSeas(frequency = 12)
V(mod) <- 1.0
diag(W(mod))[1:3] <- c(1.0, 1.0, 1.0)
#V(mod) <- 1.0/22814.31
#diag(W(mod))[1:3] <- c(1.0/2472.64,1.0/55947.63,1.0/4309.207)

filt <- dlmFilter(y, mod)
samp <- dlmBSample(filt)

eps <- .Machine$double.eps^.4
i=27
filt$m[i,]
filt$a[i,]
filt$U.C[[i]] %*% diag(filt$D.C[i,])
filt$D.C[i,]

tmp <- La.svd(filt$mod$W,nu=0)
Dw <- sqrt(tmp$d)
Dw <- pmax(Dw, eps)
Dw.inv <- 1/Dw
sqrtWinv <- Dw.inv * tmp$vt # t()%*%() = W^(-1)
sqrtWinv[abs(sqrtWinv)==Inf] <- 0

D.inv <- 1/filt$D.C[i,]
D.inv[abs(D.inv)==Inf] <- 0
La.svd(rbind(sqrtWinv %*% filt$mod$GG %*% filt$U.C[[i]],
             diag(x=D.inv,nrow=length(D.inv))), nu=0)

tG.Winv <- t(filt$mod$GG) %*% crossprod(sqrtWinv)
D.inv <- 1/filt$D.C[i,]; D.inv[abs(D.inv)==Inf] <- 0
tmp <- La.svd(rbind(sqrtWinv %*% filt$mod$GG %*% filt$U.C[[i]],
                    diag(x=D.inv,nrow=length(D.inv))), nu=0)
U.H <- filt$U.C[[i]] %*% t(tmp$vt)
D.H <- 1/tmp$d;
D.H[abs(D.H)==Inf] <- 0
h <- filt$m[i,] + crossprod(D.H*t(U.H)) %*%
  tG.Winv %*% (t(theta[i+1,,drop=F])-filt$a[i,])
theta[i,] <- h + U.H %*% matrix(D.H*rnorm(p))
U.H %*% diag(D.H)

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
