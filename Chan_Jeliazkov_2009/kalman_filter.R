#------------------------------------------------------------------------------
# R Code for a fast Kalman Filter based on 
# Chan, J.C.C. and Jeliazkov, I. (2009). Efficient Simulation and
# Integrated Likelihood Estimation in State Space Models,
# International Journal of Mathematical Modelling and Numerical
# Optimisation, 1, 101-120.
#------------------------------------------------------------------------------

library(Matrix)
library(SuppDists)
library(RcppZiggurat)
library(Rcpp)


# TVP - VAR #---------------------------------------------------------------
# VAR:
# y[t] = mu[t] + Gamma[t] * y[t-1] + epsilon[t],   epsilon[t] ~ N(0, Omega_11)
# 
# State Space Form
# y[t] = X[t] * beta[t] + epsilon[t]   epsilon[t] ~ N(0, Omega_11)
# beta[t]=beta[t-1] + nu[t]     nu[t] ~ N(0, Omega_22)
#
# beta[t] = vec(mu[t], Gamma[t]')'
# 

rm(list = ls())


set.seed(123)

nsim <- 100
nburn <- 0
total_runs <- nsim + nburn

data <- read.csv("./papers/Chan_Jeliazkov_2009/USdata.csv", header = FALSE)
colnames(data) <- c("gpd_growth", "unemp", "inter", "inf")
y <- ts(data = data,
        start = c(1948, 1),
        freq = 4)

#y <- c(log(y)) * 100
y <- as.matrix(y)
colnames(y) <- NULL
y0 <- y[1:3, ]
#y <- y[-c(1:3), ]

Y <- matrix(c(t(y[-c(1:3),])))

tt <- nrow(y)
nn <- ncol(y)
TT <- nrow(Y)/nn
qq <- nn*(nn+1)
TTqq <- TT*qq
# priors #---------------------------------------------------------------------

# Omega_11
nu1 <- nn + 3 # nn?
S1 <- diag(nn)

# Omega_22
DD <- 5 * diag(qq)
nu2 <- rep(6, qq)
S2 <- rep(0.01, qq)


# initial values #-------------------------------------------------------------

# Omega_11
Omega11 <- Matrix(cov(y))
Omega11_inv <- Matrix(solve(Omega11, diag(ncol(Omega11))))

# H
diags_H <- list(rep(1, TT*qq), rep(-1, (TT - 1)*qq))
H <- bandSparse(n = TT*qq,
               k = c(0, -(qq+1)),
               diag = diags_H,
               symm = FALSE
               )

# S
DD_inv <- solve(DD)
Omega22 <- 0.01*diag(qq) # initial values for variance of betas?
Omega22_inv <- solve(Omega22)
S <- bdiag(replicate((TT-1), Omega22, simplify = F))
S <- bdiag(DD, S)  

S_inv <- bdiag(replicate((TT-1), Omega22_inv, simplify = F))
S_inv <-  bdiag(DD_inv, S_inv)  

# bigG = X
#G <- bdiag(replicate(qq, matrix(c(1, y[3:(TT-1),]), ncol=nn+1), simplify = F))

G_list <- vector("list", 245)
for (i in 3:247){
    G_list[[i-2]] <- kronecker(diag(4), matrix(c(1, y[i,]), ncol=nn+1))
}
G <- bdiag(G_list)

new_nu2 <- (nu2+(TT-1))/2

# store
store_beta <- array(NA, dim=c(nsim, TT*qq))
store_Omega11 <- array(NA, dim=c(nn, nn, nsim))
store_Omega22 <- array(NA, dim=c(nsim, qq))
now <- Sys.time()
# Gibbs Sampler -----------------------------------------------------------
for (ii in 1:total_runs) {

# beta --------------------------------------------------------------------
    S_inv <- kronecker(Diagonal(TT), Omega22_inv)
    S_inv[1:qq, 1:qq] <- DD_inv 
    K <- Matrix::forceSymmetric(crossprod(H, S_inv) %*% H)
    G_Omega11_inv <- crossprod(G, kronecker(diag(TT), Omega11_inv))
    G_Omega11_inv_G <-  Matrix::forceSymmetric(G_Omega11_inv %*% G)
    P <- K + G_Omega11_inv_G
    L <- chol(P)
    beta_hat <- Matrix(backsolve(L, forwardsolve(L, G_Omega11_inv  %*% Y, upper.tri = FALSE),
                          upper.tri = TRUE, transpose = TRUE))
    beta <- beta_hat + solve(L, zrnorm(n = length(beta_hat)))
    
# Omega_11 ----------------------------------------------------------------
    e1 <- matrix(Y-G%*%beta, nrow=nn)
    new_S1 <- S1 + tcrossprod(e1) 
    Omega11_inv <- rWishart(n = 1, df = nu1 + TT, Sigma =  solve(new_S1, diag(nn)))[,,1]
    Omega11 <- solve(Omega11, diag(nn))

# Omega_22 |---------------------------------------------------------------
    e2 <- matrix(H %*% beta, ncol=qq)
    new_S2 <- (S2 + colSums(e2[-1,]**2))/2
    diag(Omega22) <- 1/rgamma(n=20,
                 shape = new_nu2,
                 scale = 1/new_S2
                     )

# store |------------------------------------------------------------------
    if (ii > nburn) {
        run <- ii - nburn
        store_beta[run, ] <- beta@x
        store_Omega11[ , , run] <- matrix(Omega11@x, nrow=nn)
        store_Omega22[run, ] <- diag(Omega22)
    }
}
print(Sys.time() - now)


# Results #--------------------------------------------------------------------
beta_hat <- colMeans(store_beta)
Omega11_hat <- apply(store_Omega11, MARGIN = c(1,2), FUN = mean)
Omega22_hat <- colMeans(store_Omega22)

y_hat <- matrix(beta_hat, nrow=qq)
XX <- array(NA, c(4,5,245))
for (i in 1:245){
    XX[,,i] <- matrix(y_hat[,i], nrow=nn)
}


dimnames(XX) <- list(NULL, c("mu", "gdp", "une", "int", "inf"), NULL)
gdp <- t(XX[1,,])
une <- t(XX[2,,])
int <- t(XX[3,,])
inf <- t(XX[4,,])
#reshape2::melt(gdp)

ggplot2::ggplot(reshape2::melt(inf),
    ggplot2::aes(x = Var1, y = value, col=Var2)) +
    ggplot2::geom_line()+
    ggplot2::ylim(c(-0.4,0.4))






