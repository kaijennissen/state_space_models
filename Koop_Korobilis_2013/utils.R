.transform <- function(ydata, tcode, yearlab) {
  d <- dim(ydata)
  Yraw <-  matrix(0, d[1], d[2])
  for (i in 1:ncol(ydata)) {
    Yraw[, i] = .transx(ydata[, i], tcode[i])
  }
  return(Yraw)
}


.transx <- function(x, tcode) {
  #    Transform x
  #    Return Series with same dimension and corresponding dates
  #    Missing values where not calculated
  #    -- Tcodes:
  #             1 Level
  #             2 First Difference
  #             3 Second Difference
  #             4 Log-Level
  #             5 Log-First-Difference
  #             6 Log-Second-Difference
  #             #7 Detrend Log Using 1-sided HP detrending for Monthly data
  #             #8 Detrend Log Using 1-sided HP detrending for Quarterly data
  #            #16 Log-Second-Difference
  #            #17 (1-L)(1-L^12)
  
  small <- 1e-4
  relvarm <- .00000075
  relvarq <- .000625    #HP parameter
  #.00000075 for monthly data
  #.000625 for quarterly data, see Harvey/Jeager (1993), page 234 @
  n <-  length(x)
  y <- rep(0, n)      #storage space for y
  
  if (tcode == 1) {
    # Level
    y <- x
  } else if (tcode == 2) {
    # First Difference
    y[2:n] <- x[2:n] - x[1:(n - 1)]
  } else if (tcode == 3) {
    # Second Difference
    y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
  } else if (tcode == 4) {
    # Log-Level
    if (min(x) < small) {
      y <- NA
    }
    x <- log(x)
    y <- x
    
  } else if (tcode == 5) {
    # Log-First-Difference
    if (min(x) < small) {
      y <- NA
    }
    x <- log(x)
    y[2:n] <- x[2:n] - x[1:(n - 1)]
  } else {
    #tcode == 6
    # Log-Second-Difference
    if (min(x) < small) {
      y <- NA
    }
    x <- log(x)
    y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
  }
  return(y)
}


.standardize <- function(x) {
  # Function to make your data have mean 0 and variance 1.
  # Data in x are Txp, i.e. T time series observations times p variables
  y <-
    (x - kronecker(matrix(1, nrow = nrow(x)), matrix(colMeans(x), nrow =
                                                       1))) / kronecker(matrix(1, nrow = nrow(x)), matrix(apply(x, MARGIN =
                                                                                                                  2, FUN = sd), nrow =
                                                                                                            1))
  return(y)
}


.standardize1 <- function(x, T_thres) {
  # standardize using mean and variance from training sample
  # Data in x are Txp, i.e. T time series observations times p variables
  y <-
    (x - kronecker(matrix(1, nrow = nrow(x)), matrix(colMeans(x[1:T_thres,]), nrow =
                                                       1))) / kronecker(matrix(1, nrow = nrow(x)), matrix(apply(x[1:T_thres,], MARGIN =
                                                                                                                  2, FUN = sd), nrow =
                                                                                                            1))
  return(y)
}


.create_RHS <- function(YY, M, p, t) {
  #YY <- ylag
  # M <- M[[ss]]
  p <- p
  t <- t
  # K is the number of elements in the state vector
  K <-  M + p * (M ^ 2)
  # Create x_t matrix.
  # first find the zeros in matrix x_t
  x_t <-  .zeros((t - p) * M, K)
  for (i in 1:(t - p)) {
    ztemp = .eye(M)
    for (j in 1:p) {
      xtemp <-  YY[i, ((j - 1) * M + 1):(j * M), drop = FALSE]
      xtemp <-  kronecker(.eye(M), xtemp)
      ztemp <-  cbind(ztemp, xtemp)
    }
    x_t[((i - 1) * M + 1):(i * M),] = ztemp
  }
  return(list(x_t, K))
}


.zeros <- function(n, k) {
  if (!is.null(n) & !is.null(k)) {
    return(matrix(0, nrow = n, ncol = k))
  } else if (any(!is.null(n), !is.null(k))) {
    return(rep(0, max(n, k)))
  }
}

.ones <- function(n, k) {
  if (!is.null(n) & !is.null(k)) {
    return(matrix(1, nrow = n, ncol = k))
  } else if (any(!is.null(n), !is.null(k))) {
    return(rep(1, max(n, k)))
  }
}


.eye <- function(d) {
  return(diag(d))
}


.mlag2 <- function(X, p) {
  size = dim(X)
  Traw = size[1]
  N = size[2]
  Xlag = .zeros(n = Traw, k = N * p)
  for (ii in 1:p) {
    Xlag[((p + 1):Traw), (N * (ii - 1) + 1):(N * ii)] <-
      X[(p + 1 - ii):(Traw - ii), 1:N]
  }
  return(Xlag)
}



.Minn_prior_KOOP <- function(alpha_bar, gamma, M, p, K) {
  # This is the version of the Minnesota prior with no dependence on the
  # standard deviations of the univariate regressions. This prior allows
  # online estimation and forecasting of the large TVP-VAR.
  
  # 1. Minnesota Mean on VAR regression coefficients
  A_prior <- cbind(.zeros(M, 1), 0 * .eye(M), .zeros(M, M * (p - 1)))
  a_prior <- c(A_prior)
  
  
  # 2. Minnesota Variance on VAR regression coefficients
  
  # Create an array of dimensions K x M, which will contain the K diagonal
  # elements of the covariance matrix, in each of the M equations.
  V_i <- .zeros(K / M, M)
  
  for (i in 1:M) {
    # for each i-th equation
    for (j in 1:K / M) {
      # for each j-th RHS variable
      if (j == 1) {
        # if there is constant, use this code
        V_i[j, i] <-
          alpha_bar # variance on intercept
      } else{
        V_i[j, i] <-
          gamma / ((ceiling((j - 1) / M)) ^ 2) # variance on own lags
        # Note: the "ceil((j-1)/M^2)" command finds the associated lag number for each parameter
      }
    }
  }
  
  # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'
  V_i_T <-  t(V_i)
  V_prior <-
    diag(c(V_i_T))  # this is the prior variance of the vector alpha
  
  return(list(a_prior, V_prior))
}

# 
# .Minn_prior_LITT <- function(Y, Ylag, alpha_bar, gamma, M, p, K, t) {
#   # 1. Minnesota Mean on VAR regression coefficients
#   
#   A_prior <-
#     A_prior <- cbind(.zeros(M, 1), 0 * .eye(M), .zeros(M, M * (p - 1)))
#   #A_prior = [mean(Y(1:40,:)); 0.9*eye(M); zeros((p-1)*M,M)]';
#   #A_prior = [zeros(1,M) ; 0.9*repmat(eye(M),p,1)]';
#   a_prior <-  c(A_prior)
#   
#   # 2. Minnesota Variance on VAR regression coefficients
#   # Now get residual variances of univariate p-lag autoregressions. Here
#   # we just run the AR(p) model on each equation, ignoring the constant
#   # and exogenous variables (if they have been specified for the original
#   # VAR model)
#   p_MIN <-  p
#   sigma_sq <-  .zeros(M, 1) # vector to store residual variances
#   for (i in 1:M) {
#     Ylag_i <- Ylag
#     X_i <-
#       rbind(.ones(t - p_MIN + p, 1), Ylag_i[, seq(i, M * p_MIN, by = M)])#
#     Y_i <- Y[, i]
#     # OLS estimates of i-th equation
#     alpha_i <-  solve(crossprod(X_i)) %*% crossprod(X_i, Y_i)
#     sigma_sq[i, 1] <-
#       (1. / (t - p_MIN + 1)) %*% t(Y_i - X_i * alpha_i) %*% (Y_i - X_i%*%alpha_i)
#   }
#   
#   # Now define prior hyperparameters.
#   # Create an array of dimensions K x M, which will contain the K diagonal
#   # elements of the covariance matrix, in each of the M equations.
#   V_i <-  zeros(K / M, M)
#   
#   # index in each equation which are the own lags
#   ind <-  zeros(M, p)
#   
#   for (i in 1:M) {
#     ind[i,] <-  seq(i + 1, K / M, M)
#   }
#   for (i in 1:M) {
#     # for each i-th equation
#     for (j in 1:K / M) {
#       # for each j-th RHS variable
#       if (j == 1) {
#         # if there is constant, use this code
#         V_i(j, i) <-
#           alpha_bar * sigma_sq(i, 1)
#         # variance on intercept
#       }elseif (which(j == ind[i, ]) > 0){
#         V_i[j, i] <-  gamma / ((ceiling(j - 1) / M) ^ 2)
#         # variance on own lags
#         # Note: the "ceil((j-1)/M)" command finds the associated lag
#         # number for each parameter
#       } else{
#         for (kj in 1:M) {
#           if (which(j == ind[kj, ]) > 0) {
#             ll = kj
#             
#           }
#         }    # variance on other lags
#         V_i[j, i] <-  (gamma * sigma_sq(i, 1)) / (((ceiling(j - 1) / M) ^ 2) * sigma_sq(ll, 1))
#         
#       }
#     }
#   }
#   
#   # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'
#   V_i_T <-  t(V_i)
#   V_prior <- diag(c(V_i_T))
#   # this is the prior variance of the vector alpha
#   
#   Sigma_0 <- diag(sigma_sq)
#   
#   return(list(a_prior, V_prior, Sigma_0))
# }

.mvnpdf <- function(x, mu = NULL, sigma = NULL) {
  #browser()
  if (!is.matrix(x)) {
    x <- matrix(x)
  }
  
  if (is.matrix(x) & ncol(x) > 1) {
    x <- t(x)
  }
  
  
  if (is.null(mu)) {
    mu <- matrix(0, nrow = nrow(x))
  } else if (!is.matrix(mu)) {
    mu <- matrix(mu)
  }
  
  if (is.null(sigma)) {
    sigma <- diag(nrow(x))
  } else if (!is.matrix(sigma)) {
    sigma <- matrix(sigma)
  }
  
  K <- solve(sigma, diag(dim(sigma)[1]))
  
  resu <-
    sqrt(2 * pi) ** -1 * det(sigma) ** (-0.5) * exp(-0.5 * t(x - mu) %*% K %*% (x -
                                                                                  mu))
  return (resu)
}



.mvnpdfs <- function(X, Mu = NULL, Sigma = NULL) {
  pp <- nrow(X)
  qq <- ncol(X)
  
  if (is.null(Mu)) {
    Mu <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  }
  
  if (is.null(Sigma)) {
    Sigma <- array(rep(eye(ncol(X)), nrow(X)), dim = c(qq, qq, pp))
  }
  
  if (nrow(X) == 1) {
    .mvnpdf(X, Mu, Sigma)
  } else {
    # pp <- matrix(apply(X, MARGIN=1, FUN=mvnpdf))
    pp <- vector(mode = "numeric", length = nrow(X))
    
    for (i in 1:nrow(X)) {
      pp[i] <-
        .mvnpdf(t(X[i, , drop = FALSE]), t(Mu[i, , drop = FALSE]), Sigma[, , i,  drop = TRUE])
      
      
    }
    
  }
}
