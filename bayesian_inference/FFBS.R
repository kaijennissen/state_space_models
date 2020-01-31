library(dlm)


y <- AirPassengers
y_log <- log(y)


build_dlm <- function(param) {
  dlm1 <- dlmModPoly(order=2) + dlmModSeas(frequency = 12)
  V(dlm1) <- exp(param[1])
  diag(W(dlm1))[c(1:3)] <- exp(param[2:4])
  return(dlm1)
}


fit <- dlmMLE(y = y_log, parm = rnorm(4), build = build_dlm)
fit$convergence
exp(fit$par)

dlm1 <- build_dlm(fit$par)
dlm1_filtered <- dlmFilter(y = y_log, mod = dlm1)
plot.ts(cbind(dlm1_filtered$m[-(1:21),1], y_log[-(1:20)]), plot.type = "single", col= c("red", "blue"))
plot.ts(cbind(dlm1_filtered$m[-(1:21),1]+dlm1_filtered$m[-(1:21),3], y_log[-(1:20)]), plot.type = "single", col= c("red", "blue"))


# p 166
write.csv(Nile, "Nile.csv")

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















