#--------------------------------------------------------------------------
# Petris, G. & Petrone, S. & Campagnoli, P.,
# Dynamic Linear Models with R,
# Springer (2009)
#--------------------------------------------------------------------------
library(dlm)

log_pdf_mvnorm <- function(x,
                           mu,
                           Sigma = NULL,
                           P = NULL,
                           L = NULL) {
  eps <- .Machine$double.eps ^ .4
  if (is.null(P)) {
    tmp <- La.svd(Sigma)
    D <- sqrt(tmp$d)
    D <- pmax(D, eps)
    D.inv <- 1 / D
    sqrtP <- D.inv * tmp$vt
    sqrtP[abs(sqrtP) == Inf] <- 0
  }
  
  d <- dim(sqrtP)[1]
  K <- -(d / 2) * log((2 * pi)) - 0.5 * log(det(Sigma))
  p <- -0.5 * crossprod(sqrtP %*% (x - mu))
  return(p)
}

vec_log_pdf_mvnorm <- function(x,
                               mu,
                               Sigma) {
  eps <- .Machine$double.eps ^ .4
  tmp <- La.svd(Sigma)
  D <- sqrt(tmp$d)
  D <- pmax(D, eps)
  D.inv <- 1 / D
  sqrtP <- D.inv * tmp$vt
  sqrtP[abs(sqrtP) == Inf] <- 0
  
  
  
  n <- dim(mu)[2]
  log_pdf <- vector("numeric", n)
  d <- dim(sqrtP)[1]
  K <- -(d / 2) * log((2 * pi)) - 0.5 * log(det(Sigma))
  
  if (dim(x)[2] > 1) {
    if (dim(mu)[2] > 1) {
      p <- -0.5 * diag(crossprod(sqrtP %*% (x - mu)))
    } else {
      p <- -0.5 * diag(crossprod(sqrtP %*% (x - kronecker(rep(1, n), mu))))
    }
  } else {
    if (dim(mu)[2] > 1) {
      p <- -0.5 * diag(crossprod(sqrtP %*% (kronecker(rep(1, n), x) - mu)))
    } else{
      p <- -0.5 * diag(crossprod(sqrtP %*% (x - mu)))
    }
  }
  
  return(p)
}


# Particle Filter ---------------------------------------------------------
particle_filter <- function(y, N, mod) {
  FF <- FF(mod)
  GG <- GG(mod)
  W <- W(mod)
  V <- V(mod)
  
  TT <- length(y)
  q <- dim(W)[1]
  theta_pf <- array(NA_real_, c(q, N, TT + 1))
  wt <- matrix(NA_real_, N)
  m0 <- m0(mod) 
  C0 <- C0(mod) 
  L0 <- chol(C0)

  # 0
  theta_pf[, , 1] <- m0 + L0 %*% matrix(rnorm(q*N), nrow = q) # initial draws
  wt[] <- 1 / N
  
  N_0 <-  0.5 * N
  
  # TODO use svd-algorithm to compute covariance
  Sigma <- FF %*% W %*% t(FF) + V
  Sigma_inv <- solve(Sigma, diag(dim(Sigma)[1]))
  importanceSD <- W - W %*% t(FF) %*% Sigma_inv %*% FF %*% W
  
  eps <- .Machine$double.eps ^ .4
  tmp <- La.svd(importanceSD)
  D <- sqrt(tmp$d)
  
  sqrtSigmahat <- tmp$u %*% diag(D)
  L <- sqrtSigmahat
  C <- chol(importanceSD)
  
  # 1
  for (t in 2:(TT)) {
    # 1.1
    # a) transition density
    # b) optimal importance kernel
    # draw from Normal
    mu_hat <-
      (GG %*% theta_pf[, , t - 1]  + W %*% t(FF) %*%
         Sigma_inv %*% (y[t - 1] - FF %*% GG %*% theta_pf[, , t - 1]))
    theta_pf[, , t] <-
      mu_hat + C %*% matrix(rnorm(N * q), nrow = 2)
    log_target <-
      (
        dnorm(
          x = y[t - 1],
          mean = c(FF %*% theta_pf[, , t]),
          sd = c(V),
          log = TRUE
        ) +
          vec_log_pdf_mvnorm(
            x = theta_pf[, , t],
            mu = GG %*% theta_pf[, , t - 1],
            Sigma = W
          ) -
          vec_log_pdf_mvnorm(
            x = theta_pf[, , t],
            mu = mu_hat,
            Sigma = tcrossprod(L)
          )
      )
    
    wt <- wt * exp(log_target)
    
    # 1.2
    wt <- wt / sum(wt)
    
    # 1.3
    N_eff <- 1 / crossprod(wt)
    
    # 1.4
    if (N_eff < N_0  & N_eff > 5 * N_0) {
      idx <- sample(x = N,
                    N,
                    replace = TRUE,
                    prob = wt[, 1])
      theta_pf[, , t] <- theta_pf[, idx, t]
      wt[] <- 1 / N
    }
  }
  return(theta_pf)
}

# Simulate ----------------------------------------------------------------
n <- 200
mod <-
  dlmModPoly(
    2,
    dV = 2,
    dW = c(0.1, 0.5),
    m0 = c(100, 0.05),
    C0 = diag(c(3, 3))
  )
GG(mod)[2, 2] <- 0.5 # AR-Trend
simData <- dlmForecast(mod, nAhead = n, sampleNew = 1)
y <- simData$newObs[[1]]
plot.ts(y)


# Particle Filter
mod_pf <-
  dlmModPoly(
    2,
    dV = 1,
    dW = c(1, 1),
    m0 = c(100, 0.05),
    C0 = diag(c(100, 100))
  )
theta_pf <- particle_filter(y, 2000, mod_pf)
theta_pf_hat <- t(apply(theta_pf, MARGIN = c(1, 3), FUN = mean))

# Kalman Filter
y_filt <- dlmFilter(y, mod)
theta_kf <- y_filt$m

# Plot
plot.ts(
  cbind(c(y), # observed
        theta_pf_hat[-1, 1], # particle filter
        #theta_gibbs_hat[1, -1],
        theta_kf[-1, 1]),
  plot.type = "single",
  col = c("red", "blue", "darkgreen")
)
