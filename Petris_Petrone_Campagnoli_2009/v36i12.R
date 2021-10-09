###################################################
### Preliminaries
###################################################
library("dlm")


###################################################
### How to define a local level model
### (two equivalent ways)
###################################################
myModel <- dlm(FF = 1, V = 3.1, GG = 1, W = 1.2, m0 = 0, C0 = 100)
myModel <- dlmModPoly(order = 1, dV = 3.1, dW = 1.2, C0 = 100)


###################################################
### A local trend plus quarterly seasonal component
###################################################
myModel <- dlmModPoly(2) + dlmModSeas(4)


###################################################
### Same as above, for a bivariate series
###################################################
bivarMod <- myModel %+% myModel


###################################################
### Maximum Likelihood Estimation
###################################################
buildTemp <- function(x) {
     L <- matrix(0, 2, 2)
     L[upper.tri(L, TRUE)] <- x[1 : 3]
     diag(L) <- exp(diag(L))
     modTemp <- dlm(FF = matrix(1, 2, 1),
                    V = crossprod(L),
                    GG = 1,
                    W = exp(x[4]),
                    m0 = 0,
                    C0 = 1e7)
     return(modTemp)
}
y1 <- scan("HL.dat") # land
#or:# y1 <- scan("http://www.stat.pitt.edu/stoffer/tsa2/data/HL.dat") # land
y2 <- scan("Folland.dat") # sea
#or:# y2 <- scan("http://www.stat.pitt.edu/stoffer/tsa2/data/Folland.dat") # sea

y <- ts(cbind(y1, y2), start = 1880)
fitTemp <- dlmMLE(y, parm = rep(0, 4), build = buildTemp,
                  hessian = TRUE,
                  control = list(maxit = 500))

modTemp <- buildTemp(fitTemp$par) # fitted model

V(modTemp) # MLE of 'V'
drop(W(modTemp)) # MLE of 'W'


###################################################
### MLE: obtaining standard errors
###################################################
buildNile <- function(x) dlmModPoly(1, dV = x[1], dW = x[2])
fitNile <- dlmMLE(Nile, parm = rep(100, 2), build = buildNile,
                  lower = rep(1e-8, 2), hessian = TRUE)
fitNile$par

aVar <- solve(fitNile$hessian)
aVar
sqrt(diag(aVar)) # estimated standard errors


###################################################
### Forward Filtering Backward Sampling
###################################################
modNile <- dlmModPoly(1, dV = 15100, dW = 1470)
nileFilt <- dlmFilter(Nile, modNile)
eta <- replicate(1000, max(diff(dlmBSample(nileFilt))))


###################################################
### Gibbs sampler for d-Inverse Gamma model
###################################################
set.seed(999)
mcmc <- 1000; burn <- 500
outGibbsIRW <- dlmGibbsDIG(y = log(UKgas),
                           mod = dlmModPoly(2) + dlmModSeas(4),
                           shape.y = 1e-3, rate.y = 1e-3,
                           shape.theta = 1e-3, rate.theta = 1e-3,
                           n.sample = mcmc + burn, ind = c(2, 3))

### Trace plot of sampler's output
par(mar = c(3.1,4.1,1.1,2.1), cex = 0.7)
matplot(with(outGibbsIRW, sqrt(cbind(dV[-(1:burn)], dW[-(1:burn), ]))),
        ylab = "", type = "l")
legend("topleft", lty = "solid", col = 1 : 3, bty = "n",
       legend = c(expression(sigma[y]), expression(sigma[beta]),
       expression(sigma[s])))

### Bayes estimates of V and W, with Monte Carlo standard errors
mcmcMean(with(outGibbsIRW, sqrt(cbind(V = dV[-(1 : burn)],
                                      dW[-(1 : burn), ]))))


###################################################
### Gibbs sampler for outliers and structural breaks
### (d-Inverse Gamma model with t-distributed innovations)
### Warning: fairly slow!!
###################################################
source("dlmGibbsDIGt.R")
y <- ts(read.table("IPCONGD.txt",
                   skip = 11, header = TRUE)[, 2],
        frequency = 12, start = c(1939, 1))
yM <- log(y)
plot(yM, main = "Data")
mc <- 2000
burn <- 1000
gibbsOut <- dlmGibbsDIGt(yM, mod = dlmModPoly(2),
                         A_y = 5e4, B_y = 5e4, ind = 1 : 2,
                         save.states = TRUE,
                         n.sample = mc + burn, thin = 4)

### Output analysis
burn <- 1 : burn

### Posterior means of the omega's
omega_y <- ts(colMeans(gibbsOut$omega_y[-burn, ]),
              start = start(y), freq = 12)
omega_theta <- ts(apply(gibbsOut$omega_theta[, , -burn], 1 : 2,
                        mean), end = end(y), freq = 12)
## Posterior means of the states
thetaMean <- ts(apply(gibbsOut$theta[, 1 : 2, -burn], 1 : 2, mean),
                end = end(y),
                freq = frequency(y))
## Posterior means of residuals
resid <- yM - dropFirst(thetaMean[, 1])
## Plots
pdf("Rplot-lambdaOmegaT.pdf",
   width = 6 * 0.618, height = 5, pointsize = 8)
layout(matrix(1 : 6), heights = rep(c(5, 3), 3))
par(oma = c(2, 0, 2, 0) + 0.1, mar = c(0, 4, 0, 4), xpd = NA)
plot(thetaMean[, 1], axes = FALSE, xlab = "", ylab = "Trend")
axis(2); axis(3)
usr <- par("usr")
segments(x0 = c(usr[1], usr[2]),
         y0 = rep(usr[3], 2),
         x1 = c(usr[1], usr[2]),
         y1 = rep(usr[4], 2))
abline(h = usr[4], xpd = FALSE)
plot(omega_theta[, 1], type = "p", ylim = c(0, 1.2), pch = 20,
     cex = 0.5, xlab = "", ylab = "", axes = FALSE)
abline(h = 1, lty = "dashed", col = "grey", xpd = FALSE)
axis(4); mtext(expression(omega[list(theta * 1, t)]), side = 4, line = 3)
usr <- par("usr")
segments(x0 = c(usr[1], usr[2]),
         y0 = rep(usr[3], 2),
         x1 = c(usr[1], usr[2]),
         y1 = rep(usr[4], 2))
abline(h = usr[3], col = "grey", xpd = FALSE)
plot(thetaMean[, 2], axes = FALSE, xlab = "", ylab = "Slope")
axis(2)
abline(h = 0, lty = "dashed", col = "grey", xpd = FALSE)
usr <- par("usr")
segments(x0 = c(usr[1], usr[2]),
         y0 = rep(usr[3], 2),
         x1 = c(usr[1], usr[2]),
         y1 = rep(usr[4], 2))
plot(omega_theta[, 2], type = "p", ylim = c(0, 1.2), pch = 20,
     cex = 0.5, xlab = "", ylab = "", axes = FALSE)
abline(h = 1, lty = "dashed", col = "grey", xpd = FALSE)
axis(4); mtext(expression(omega[list(theta * 2, t)]), side = 4, line = 3)
usr <- par("usr")
segments(x0 = c(usr[1], usr[2]),
         y0 = rep(usr[3], 2),
         x1 = c(usr[1], usr[2]),
         y1 = rep(usr[4], 2))
abline(h = usr[3], col = "grey", xpd = FALSE)
plot(resid, axes = FALSE, xlab = "", ylab = "Residuals")
abline(h = 0, lty = "dashed", col = "grey", xpd = FALSE)
axis(2)
usr <- par("usr")
segments(x0 = c(usr[1], usr[2]),
         y0 = rep(usr[3], 2),
         x1 = c(usr[1], usr[2]),
         y1 = rep(usr[4], 2))
plot(omega_y, type = "p", ylim = c(0, 1.2), pch = 20,
     cex = 0.5, xlab = "", ylab = "", axes = FALSE)
abline(h = 1, lty = "dashed", col = "grey", xpd = FALSE)
axis (1); axis(4); mtext(expression(omega[list(y, t)]), side = 4, line = 3)
box(bty = "u")
dev.off()

