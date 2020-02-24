#--------------------------------------------------------------------------
# Petris, G. & Petrone, S. & Campagnoli, P.,                             
# Dynamic Linear Models with R,                                          
# Springer (2009)                                                        
#--------------------------------------------------------------------------
library(dlm)

# Example 1 ---------------------------------------------------------------
# Forward Filtering Backward Sampling using the dlm package
y <- Nile

ffbs_nile <- function(y, mc) {
  N <- length(y)
  a1 <- 2
  b1 <- 0.0001
  a2 <- 2
  b2 <- 0.0001
  
  psi1 <- 1
  psi2 <- 2
  mod_level <- dlmModPoly(order = 1,
                          dV = 1 / psi1,
                          dW = 1 / psi2)
  
  psi1_save <- numeric(mc)
  psi2_save <- numeric(mc)
  theta_save <- matrix(0, nrow = N + 1, ncol = mc)
  
  sh1 <- a1 + N / 2
  sh2 <- a2 + N / 2
  
  set.seed(10)
  
  for (it in 1:mc) {
    # draw the states: FFBS
    filt <- dlmFilter(y, mod_level)
    theta <- dlmBSample(filt)
    
    # draw observation precision psi1
    rate <- b1 + crossprod(y - theta[-1]) / 2
    psi1 <- rgamma(1, shape = sh1, rate = rate)
    
    # draw system precision psi2
    rate <- b2 + crossprod(theta[-1] - theta[-N]) / 2
    psi2 <- rgamma(1, shape = sh2, rate = rate)
    
    # update and save
    V(mod_level) <- 1 / psi1
    W(mod_level) <- 1 / psi2
    
    theta_save[, it] <- c(theta)
    psi1_save[it] <- psi1
    psi2_save[it] <- psi2
  }
  return(list(theta_save, psi1_save, psi2_save))
}

# MCMC estimation
# approx. 6 sec. per 1000 draws
mc <- 5e4
now <- Sys.time()
resu <- ffbs_nile(y, mc)
Sys.time() - now

theta <- resu[[1]]
psi1 <- resu[[2]]
psi2 <- resu[[3]]

# Diagnostic
nburn <- 2.5e4
plot.ts(cumsum(psi1) / 1:mc)
plot.ts(cumsum(psi2) / 1:mc)

mean(psi1[-(1:nburn)])
mean(psi2[-(1:nburn)])

theta_hat <- rowMeans(theta[-1, -(1:nburn)])
plot.ts(cbind(c(y), theta_hat), plot.type = "single")


# Example 2 ---------------------------------------------------------------
# Forward Filtering Backward Sampling using the dlm package

y <- log(c(AirPassengers))

ffbs_airpassengers <- function(y, mc) {
  N <- length(y)
  ay <- 1e-2
  by <- 1e-3
  a_psi1 <- 2.5e-1
  b_psi1 <- 5e-3
  a_psi2 <- 2.5e7
  b_psi2 <- 500
  a_psi3 <- 2.5e7
  b_psi3 <- 50
  
  set.seed(123)
  psiy <- 1
  psi1 <- 50
  psi2 <- 5e5
  psi3 <- 5e5
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
  theta_save <- array(0, dim = c(N + 1, nrow(W), mc))
  
  shy <- ay + N / 2
  sh1 <- a_psi1 + N / 2
  sh2 <- a_psi2 + N / 2
  sh3 <- a_psi3 + N / 2
  
  set.seed(10)
  
  for (it in 1:mc) {
    FF <- FF(mod_level)
    GG <- GG(mod_level)
    
    # draw the states: FFBS
    filt <- dlmFilter(y, mod_level)
    theta <- dlmBSample(filt)
    
    # draw observation precision psiy
    rate <- by + crossprod(y - theta[-1,] %*% t(FF)) / 2
    psiy <- rgamma(1, shape = shy, rate = rate)
    
    SS_theta <- diag(crossprod(theta[-1,] - theta[-N,] %*% t(GG)))
    # draw system precision psi2
    rate <- b_psi1 + SS_theta[1] / 2
    psi1 <- rgamma(1, shape = sh1, rate = rate)
    
    # draw system precision psi2
    rate <- b_psi2 + SS_theta[2] / 2
    psi2 <- rgamma(1, shape = sh2, rate = rate)
    
    # draw system precision psi2
    rate <- b_psi3 + SS_theta[3] / 2
    psi3 <- rgamma(1, shape = sh3, rate = rate)
    
    # update and save
    V(mod_level) <- 1 / psiy
    diag(W(mod_level))[1:3] <- c(1 / psi1, 1 / psi2, 1 / psi3)
    theta_save[, , it] <- theta
    psiy_save[it] <- psiy
    psi1_save[it] <- psi1
    psi2_save[it] <- psi2
    psi3_save[it] <- psi3
    
    if (10 * it %% mc == 0) {
      comp_perc <- 100 * it / mc
      print(paste("completion: ", comp_perc))
    }
  }
  return(list(psiy_save, psi1_save, psi2_save, psi3_save, theta_save))
}

# MCMC estimation
# approx. 34 sec per 1000 draws
now <- Sys.time()
resu <- ffbs_airpassengers(y, 5e4)
Sys.time() - now

psiy <- resu[[1]]
psi1 <- resu[[2]]
psi2 <- resu[[3]]
psi3 <- resu[[4]]
theta <- resu[[5]]

# Diagnostics
nburn <- 2.5e4
plot.ts(cumsum(psiy) / 1:mc)
plot.ts(cumsum(psi1) / 1:mc)
plot.ts(cumsum(psi2) / 1:mc)
plot.ts(cumsum(psi3) / 1:mc)
mean(psiy[-(1:nburn)])
mean(psi1[-(1:nburn)])
mean(psi2[-(1:nburn)])
mean(psi3[-(1:nburn)])

FF <- matrix(c(1, 0, 1, rep(0, 10)), nrow = 1)
theta_hat <- apply(theta[, , -(1:nburn)], c(1, 2), mean)
plot.ts(cbind(c(y), theta_hat[-1, 1]),
        plot.type = "single",
        col = c("black", "red"))

# Example 3 ---------------------------------------------------------------
# Forward Filtering Backward Sampling using the dlm package
# p. 182

y <- log(UKgas)
set.seed(4521)
now <- Sys.time()
MCMC <- 10500
mod <- dlmModPoly(2) + dlmModSeas(4)
gibbsOut <- dlmGibbsDIGt(
  y = y,
  mod = mod,
  A_y = 10000,
  B_y = 10000,
  n.sample = MCMC,
  save.states = TRUE,
  thin = 2
)
Sys.time()-now

burn <- 1:500
nuRange <- c(1:10, seq(20, 100, by = 10))


omega_y <- ts(colMeans(gibbsOut$omega_y[-burn, ]),
              start = start(y),
              freq = 4)

omega_theta <- ts(apply(gibbsOut$omega_theta[,,-burn],
                        MARGIN = c(1:2),
                        FUN = mean),
              start = start(y),
              freq = 4)

layout(matrix(c(1, 2, 3, 4), 4, 1, TRUE))
par(mar = c(5.1, 4.1, 2.1, 2.1))

plot(omega_y,
     type = "p",
     ylim = c(0, 1.2),
     pch = 16,
     xlab = "",
     ylab = expression(omega[list(y, t)])
     )
abline(h=1, lty = "dashed")
for (i in 1:3){
  plot(omega_theta[, i],
     type = "p",
     ylim = c(0, 1.2),
     pch = 16,
     xlab = "",
     ylab = bquote(omega[list(theta, t * .(i))])
     )
     abline(h=1, lty = "dashed")

  
}

par(mfrow=c(1,1))

thetaMean <- ts(apply(gibbsOut$theta, MARGIN = c(1, 2), FUN = mean),
                start = start(y),
                freq = frequency(y)
                )
LprobLim <- ts(apply(gibbsOut$theta, MARGIN = c(1, 2), FUN = quantile, probs = 0.025),
               start = start(y),
               freq = frequency(y)
)
UprobLim <- ts(apply(gibbsOut$theta, MARGIN = c(1, 2), FUN = quantile, probs = 0.975),
               start = start(y),
               freq = frequency(y)
)
par(mfrow = c(2, 1), mar = c(5.1, 4.1, 2.1, 2.1))

plot(thetaMean[, 1], xlab = "", ylab = "Trend")
lines(LprobLim[, 1], lty = 2)
lines(UprobLim[, 1], lty = 2)

plot(thetaMean[, 3], xlab = "", ylab = "Seasonal", type = "o")
lines(LprobLim[, 3], lty = 2)
lines(UprobLim[, 3], lty = 2)

