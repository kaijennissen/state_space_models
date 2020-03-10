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
  thetaPF <- array(NA_real_, c(q, TT + 1, N))
  wt <- matrix(NA_real_, TT+1, N)
  m0 <- m0(mod) 
  C0 <- C0(mod) 
  L0 <- chol(C0)

  # 0
  thetaPF[, 1, ] <- m0 + L0 %*% matrix(rnorm(q*N), nrow = q) # initial draws
  wt[1, ] <- 1 / N
  
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
      (GG %*% thetaPF[, t - 1, ]  + W %*% t(FF) %*%
         Sigma_inv %*% (y[t - 1] - FF %*% GG %*% thetaPF[ , t - 1, ]))
    thetaPF[, t, ] <-
      mu_hat + C %*% matrix(rnorm(N * q), nrow = 2)
    log_target <-
      (
        dnorm(
          x = y[t - 1],
          mean = c(FF %*% thetaPF[, t, ]),
          sd = c(V),
          log = TRUE
        ) +
          vec_log_pdf_mvnorm(
            x = thetaPF[, t, ],
            mu = GG %*% thetaPF[ , t - 1, ],
            Sigma = W
          ) -
          vec_log_pdf_mvnorm(
            x = thetaPF[, t, ],
            mu = mu_hat,
            Sigma = tcrossprod(L)
          )
      )
    
    wt[t, ] <- wt[t-1, ] * exp(log_target)
    
    # 1.2
    wt[t, ] <- wt[t, ]  / sum(wt[t, ])
    
    # 1.3
    N_eff <- 1 / crossprod(wt[t, ])
    
    # 1.4
    #browser()
    if (N_eff < N_0 ) {
      idx <- sample(x = N,
                    N,
                    replace = TRUE,
                    prob = wt[t, ])
      thetaPF[, t, ] <- thetaPF[, t, idx]
      wt[t, ] <- 1 / N
    }
  }
  return(list(theta=thetaPF, wt=wt))
}

# Simulate ----------------------------------------------------------------
n <- 100
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
# mod_pf <-
#   dlmModPoly(
#     2,
#     dV = 1,
#     dW = c(1, 1),
#     m0 = c(100, 0.05),
#     C0 = diag(c(100, 100))
#   )
resuPF <- particle_filter(y, 1000, mod)
thetaPF <- resuPF$theta[, -1, ]
wt <- resuPF$wt[-1,  ]
thetaHatPF <- t(sapply(1:TT, FUN = function(i){thetaPF[, i, ]%*%wt[i, ]} ))

# Kalman Filter
y_filt <- dlmFilter(y, mod)
thetaHatKF <- y_filt$m[-1, ]

# Plot
plot.ts(
  cbind(c(y), # observed
        thetaHatPF[, 1], # particle filter
        #theta_gibbs_hat[1, -1],
        thetaHatKF[, 1]),
  plot.type = "single",
  col = c("red", "blue", "darkgreen")
)



# Petris ------------------------------------------------------------------
mod <- dlmModPoly(1, dV=2, dW=1, m0=10, C0=9)
n <- 100
set.seed(23)
simData <- dlmForecast(mod, nAhead = n, sampleNew = 1)
y <- simData$newObs[[1]]
plot.ts(y)

# particle filter
N <- 1000
N_0 <- N/2
pfOut <- matrix(NA_real_, n+1, N)
wt <- matrix(NA_real_, n+1, N)
importanceSD <- sqrt(drop(W(mod) - W(mod)^2 / (W(mod)+V(mod))))
predSd <- sqrt(drop(W(mod) + V(mod)))

pfOut[1, ] <- rnorm(N, mean=m0(mod), sd=sqrt(C0(mod)))
wt[1, ] <- rep(1/N, N)
for (it in 2:(n+1)){
  means <- pfOut[it-1, ] + drop(W(mod)) * (y[it-1] -pfOut[it-1, ]) / drop(W(mod) + V(mod))
  pfOut[it, ] <- rnorm(N, mean = means, sd = importanceSD)
  wt[it, ] <- dnorm(y[it-1], mean = pfOut[it-1, ], sd =  predSd) * wt[it-1, ]
  
  wt[it, ] <- wt[it, ] / sum(wt[it, ])
  N.eff <- 1/crossprod(wt[it, ])
if (N.eff < N_0){
  
  index <- sample(N, N, replace=TRUE, prob = wt[it, ])
  pfOut[it, ] <- pfOut[it, index]
  wt[it, ] <- 1/N
  }
}



modFilt <- dlmFilter(y, mod)
thetaHatKF <- modFilt$m[-1]
sdKF <- with(modFilt, sqrt(unlist(dlmSvd2var(U.C, D.C))))[-1]
pfOut <- pfOut[-1, ]
wt <- wt[-1, ]
thetaHatPF <- sapply(1:n, function(i){
  weighted.mean(pfOut[i, ], wt[i, ])
})
sdPF <- sapply(1:n, function(i){
  sqrt(weighted.mean((pfOut[i, ] - thetaHatPF[i])^2, wt[i, ]))
})

plot.ts(cbind(thetaHatKF, thetaHatPF), plot.type = "s", lty = c("dotted", "longdash"), 
        xlab="", ylab=expression(m[t]))
legend("topleft", c("Kalman", "Particle"), lty=c("dotted", "longdash"), bty = "n")

plot.ts(cbind(sdKF, sdPF), plot.type = "s", lty = c("dotted", "longdash"), 
        xlab="", ylab=expression(sqrt(C[t])))
legend("topleft", c("Kalman", "Particle"), lty=c("dotted", "longdash"), bty = "n")
