# Petris, G., Petrone, S. & Campagnoli, P. (2009).
# Dynamic Linear Models with R.
# p. 166
library(dlm)

# Simulate from a given state space model
simulate_y <- function(n, mod){

  # model
  FF <- matrix(c(1, 0, 1,rep(0,2)), nrow=1)
  GG <- matrix(0, 5, 5)
  GG[1,1] <- 1
  GG[2,2] <- 1
  GG[1,2] <- 1
  GG[3, 3:5] <- -1
  GG[4:5, 3:4] <- diag(2)

  #V <- matrix(1/psiy)
  #W <- diag(c(1/psi1, 1/psi2, 1/psi3, rep(0,2)))
  
  
  V <- mod$V
  W <- mod$W
  
  FF <- mod$FF
  GG <- mod$GG
  
  #seas <-  4*sin(seq(-1,1, length.out=4))[-4]
  L <- chol(mod$C0)
  theta0 <-  mod$m0 + L%*% matrix(rnorm(dim(L)[1]))
  
  # simulate
  y <- vector("numeric", n)
  
  theta <-  theta0
  tmp <- La.svd(W)
  sqrtW <- tmp$u %*% diag(sqrt(tmp$d)) 
  tmp <- La.svd(V)
  sqrtV <- tmp$u %*% sqrt(tmp$d)
  
  for (i in 1:n){
    theta <- GG %*% theta+sqrtW %*% rnorm(dim(sqrtW)[1])
    y[i] <- FF %*% theta + sqrtV %*% rnorm(dim(sqrtV)[1])
  }
  
  return(y)
  
}

# calculate gamma hyperparameter
gamma_prior <- function(a, b){
  alpha <- a^2/b
  beta <- a/b
  
  return(list(alpha, beta))
}
