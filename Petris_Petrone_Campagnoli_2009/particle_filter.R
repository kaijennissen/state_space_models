#--------------------------------------------------------------------------
# Petris, G. & Petrone, S. & Campagnoli, P.,                             
# Dynamic Linear Models with R,                                          
# Springer (2009)                                                        
#--------------------------------------------------------------------------
library(dlm)

log_pdf_mvnorm <- function(x, mu, Sigma=NULL, P = NULL, L = NULL) {
  
  eps <- .Machine$double.eps^.4
  if (is.null(P)) {
      tmp <- La.svd(Sigma)
      D <- sqrt(tmp$d)
      D <- pmax(D, eps)
      D.inv <- 1 / D
      sqrtP <- D.inv * tmp$vt
      sqrtP[abs(sqrtP) == Inf] <- 0
  }
  
  d <- dim(sqrtP)[1]
  K <- -(d/2)*log((2*pi))
  p <- -0.5*log(det(Sigma))-0.5*crossprod(sqrtP %*% (x - mu))
  return(p+K)
}

# Particle Filter ---------------------------------------------------------
y <- log(AirPassengers)
plot.ts(y)
N <- 1000
wt <- matrix(NA_real_, N)
mod <- dlmModPoly(2) + dlmModSeas(12)
q <- dim(mod$W)[1]
thetat <- array(NA_real_, c(q, N, 145))
m0 <- matrix(0, nrow=q)
C0 <- 1e3*diag(q)
L0 <- chol(C0)

# 0 
for (i in 1:N){
  thetat[, i, 1] <- m0 + L0 %*% rnorm(n = q)
}
wt[] <- 1/N

# TODO use svd-algorithm to compute covariance
XX <-  mod$W %*% t(mod$FF) %*% solve(mod$FF %*% mod$W %*% t(mod$FF) + mod$V)
Sigma_hat <-  mod$W - XX %*% mod$FF %*% mod$W
#L <- chol(Sigma_hat)
eps <- .Machine$double.eps^.4
tmp <- La.svd(Sigma_hat)
D <- sqrt(tmp$d)
#D <- pmax(D, eps)
sqrtSigmahat <- D * tmp$vt
L <- sqrtSigmahat 
# 1
#for (t in 2:TT){

  # 1.1
  for (it in 1:N){
    # a) transition density
    # b) optimal importance kernel
      # draw from Normal
    mu_hat <- mod$G %*% thetat[, it, t-1]  + XX %*% (y[t-1]-mod$FF %*% mod$GG %*% thetat[, it, t-1])
    thetat[, it, t] <- mu_hat + L %*% rnorm(q)
    x <- dnorm(x = y[t], mean = c(mod$FF %*% thetat[, it, t]), sd = c(mod$V), log=TRUE)
    + log_pdf_mvnorm(x = thetat[, it, t], mu = mod$GG %*% thetat[, it, t-1], Sigma = mod$W)
    - log_pdf_mvnorm(x = thetat[, it, t], mu = mu_hat, Sigma=tcrossprod(L))
    
    wt[i] <- wt[i]*exp(x)

  # 1.2
  wt <- wt/sum(wt)
  
    # 1.3
  N_eff <- 1/(crossprod(wt[,i]))

  # 1.4
  if (N_eff < N_0){
    sample(x = N, N, replace = TRUE, prob = wt[ut, ])
  }
  
  # 1.5
  }
#}
