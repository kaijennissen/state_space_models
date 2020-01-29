# TVP_VAR_DPS_DMA_sim.m - Forecasting with Large TVP-VAR using forgetting factors
# MULTIPLE MODEL CASE / DYNAMIC PRIOR SELECTION (DPS) AND DYNAMIC MODEL
# AVERAGING (DMA)
#------------------------------------------------------------------------------
# The model is:
#
#	 y[t] = theta[t] x[t] + e[t]
#	 theta[t] = theta[t-1] + u[t]
#
# where x[t] = I x (y[t-1],...,y[t-p]) (Kronecker product), and e[t]~N(0,V[t])
# and u[t]~N(0,Q[t]).
#
# Additionally:
#
#  V[t] = kappa V[t-1] + (1-kappa) e[t-1]e[t-1]'
#  Q[t] = (1 - 1/lambda) S[t-1|t-1]
#
# This code estimates lambda and allows it to be time-varying. The specification is:
#
#  lambda[t] = lambda[min] + (1-lambda[min]) LL^(e[t]e[t]')
#
#------------------------------------------------------------------------------
#  - This code allows to calculate ONLY iterated forecasts
#  - This code does predictive simulation by fixing theta[T+1],
#  theta[T+2],... from Normals with mean theta[T] (no drifting parameters
#  out-of-sample).
#  - This code does "online" forecasting, i.e. the Minnesota prior should not be
#  dependent on the data, so that the Kalman filter runs once for 1:T.
#------------------------------------------------------------------------------
library(Matrix)
rm(list = ls())

# --------------------| HELPER FUNCTIONS |-------------------------------------

# zeros |----------------------------------------------------------------------
ones <- function(n = NULL, k = NULL) {
    if (!is.null(n) & !is.null(k)) {
        return(matrix(1, nrow = n, ncol = k))
    } else if (any(!is.null(n), !is.null(k))) {
        return(rep(1, max(n, k)))
    }
}

zeros <- function(n, k) {
    if (!is.null(n) & !is.null(k)) {
        return(matrix(0, nrow = n, ncol = k))
    } else if (any(!is.null(n), !is.null(k))) {
        return(rep(0, max(n, k)))
    }
}

eye <- function(d) {
    return(diag(d))
}


# mlag2 |----------------------------------------------------------------------
mlag2 <- function(X, p) {
    size = dim(X)
    Traw = size[1]
    N = size[2]
    Xlag = zeros(n = Traw, k = N * p)
    for (ii in 1:p) {
        Xlag[((p + 1):Traw), (N * (ii - 1) + 1):(N * ii)] <-
            X[(p + 1 - ii):(Traw - ii), 1:N]
    }
    return(Xlag)
}

# create_RHS |-----------------------------------------------------------------
create_RHS <- function(YY, M, p, t) {
    #YY <- ylag
    M <- M[[ss]]
    p <- p
    t <- t
    # K is the number of elements in the state vector
    K <-  M + p * (M ^ 2)
    # Create x_t matrix.
    # first find the zeros in matrix x_t
    x_t <-  zeros((t - p) * M, K)
    for (i in 1:(t - p)) {
        ztemp = eye(M)
        for (j in 1:p) {
            xtemp <-  YY[i, ((j - 1) * M + 1):(j * M), drop = FALSE]
            xtemp <-  kronecker(eye(M), xtemp)
            ztemp <-  cbind(ztemp, xtemp)
        }
        x_t[((i - 1) * M + 1):(i * M),] = ztemp
    }
    return(list(x_t, K))
}

# mvnnpdf |-----------------------------------------------------------------
mvnpdf <- function(x, mu = NULL, sigma = NULL) {
    
    #browser()
    if (!is.matrix(x)){
        x <- matrix(x)
    }

    if (is.matrix(x) & ncol(x) >1 ){
        x <- t(x)
    }
    
    
    if (is.null(mu)) {
        mu <- matrix(0, nrow = nrow(x))
    } else if (!is.matrix(mu)) {
        mu<- matrix(mu)
    }
    
    if (is.null(sigma)) {
        sigma <- diag(nrow(x))
    } else if (!is.matrix(sigma)) {
        sigma <- matrix(sigma)
    }
    
    K <- solve(sigma, diag(dim(sigma)[1]))
    
    resu <-
        sqrt(2 * pi) ** -1 * det(sigma) ** (-0.5) * exp(-0.5 * t(x - mu) %*% K %*% (x -mu))
    return (resu)
}



mvnpdfs <- function(X, Mu=NULL, Sigma=NULL) {
   
    pp <- nrow(X)
    qq <- ncol(X)
    
     if (is.null(Mu)) {
        Mu <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    }
    
    if (is.null(Sigma)) {
        Sigma <- array(rep(eye(ncol(X)), nrow(X)), dim=c(qq, qq, pp))
    }
    
    if (nrow(X) == 1) {
        mvnpdf(X, Mu, Sigma)
    } else {
       # pp <- matrix(apply(X, MARGIN=1, FUN=mvnpdf))
       pp <- vector(mode = "numeric", length = nrow(X))
        
       for (i in 1:nrow(X)){
            pp[i] <-  mvnpdf(t(X[i, , drop = FALSE]), t(Mu[i, , drop = FALSE]), Sigma[, , i,  drop = TRUE])
            
            
        }
       
       }
}
# --------------------| END HELPER FUNCTIONS |---------------------------------

# Add path of data and functions
data_path <- "data"
func_path <- "functions"

#-------------------------------PRELIMINARIES----------------------------------
# Choose grids for major tuning parameters
gamma  <-  c(1e-10, 1e-5, 0.001, 0.01, 0.05, 0.1, 1, 5)
lambda <- 0.99
kappa  <-  0.94

eta <- 0.99   # Forgetting factor for DPS (dynamic prior selection) and DMA

# Please choose:
p <-  2            # p is number of lags in the VAR part
N <-  19            # Number of cross-sections (countries)

prior <- 1        # 1: Use Koop-type Minnesota prior
# 2: Use Litterman-type Minnesota prior

# Variables not to include in DMA
n_DMA <-  3      # Variables always included
varsN <-  0 + 1  # Additional variables in model averaging

# Forecasting
nfore <-  12       # Forecast horizon (note: forecasts are iterated)
t0 <- "2005M12"    # Set last observation of initial estimation period
nsim <- 1000      # Number of times to simulate from the predictive density

# Choose which results to print
# NOTE: CHOOSE ONLY 0/1 (FOR NO/YES) VALUES!
print_fore <- 1           # summary of forecasting results
print_coefficients <- 0   # plot volatilities and lambda_t (but not theta_t which is huge)
print_pred <- 0           # plot predictive likelihoods over time
print_Min <- 0            # print the Minnesota prior over time

#----------------------------------LOAD DATA----------------------------------------
# Observations of 
# INF = Inflation,
# IP = Industrial Production
# UN = Unemplyoemnt
# COMP = Real Effective Exchange Rate 
# FS1Y = Financial Situation
# GES1Y = General Economic Situation
# IE1Y = Inflation Expectation
# for 19 countries 
colnames(.data) == "INF"
.data = readr::read_csv2(
    "./papers/Koop_Korobolis_2019/data/data.csv",
    col_names = TRUE,
    col_types = paste0("c", paste0(rep("n", 133), collapse = ""), collapse =
                           "")
)
Y <- as.matrix(.data[-1, -1][]) #as(as.matrix(.data[-1,-1][]), "dgeMatrix")
yearlab <- .data[, 1, drop = TRUE][-c(1:(p + 1))]

#plot.ts(Y[,grep("INF", colnames(Y))], plot.type = "single")
ggplot2::ggplot(reshape2::melt(Y[,grep("INF", colnames(Y))]), ggplot2::aes(x=Var1, y=value, col=Var2))+
    ggplot2::geom_line(alpha=0.5)

# initial observations
T_thres <- which(yearlab == t0)  # Convert t0 to numeric value
country_index = 0

# for (country_index in 1:N) {
country_index = country_index + 1
Y1 <-vector("list", varsN) # TODO: decide whether to use arrays, lists or matrices
M <- vector("list", varsN)

for (ss in 1:varsN) {
    G <- ss + n_DMA    # Number of macro fundamentals + the exchange rate
    SS1 <- 1:G
    SS <- vector(mode = "numeric")
    GG <- SS1 + 7 * (country_index - 1)
    SS <- c(SS, GG)
    Y1[[ss]] <- Y[, SS] #Y1[ss, 1] <- Y[, SS]
    M[[ss]] <- length(SS) # M is the dimensionality of Y
}
t <- nrow(Y1[[1]])

nos <- varsN

# Inflation is the variable of interest for forecasting
nfocus = 1

# ===================================| VAR EQUATION |======================
# Generate lagged Y matrix. This will be part of the X matrix
x_t <-  vector("list", nos)
x_f <-  vector("list", nos)
y_t <-  vector("list", nos)
K <-  zeros(n = nos, k = 1)

for (ss in 1:nos) {
    ylag <-  mlag2(Y1[[ss]], p)
    ylag <-  ylag[-(1:p),]
    resu_list <-  create_RHS(ylag, M[[ss]], p, t)
    temp <-  resu_list[[1]]
    kk <- resu_list[[2]]
    x_t[[ss]] <-  cbind(ones(dim(ylag)[1], 1), ylag)
    K[[ss]] <-  kk
    x_f[[ss]] <-  ylag
    y_t[[ss]] <-  Y1[[ss]][-(1:p),]
}

yearlab <-  yearlab[-c(1:p)]
# Time series observations
t <-  dim(y_t[[1]])[1] ##ok<*NASGU>

#----------------------------| PRELIMINARIES |-----------------------------

anumber <-  t - T_thres + 1
y_fore <-  vector("list", length = nos)
y_fore[[ss]] <- array(NA, dim = c(nfore, 4, nsim))

if (country_index == 1) {
    LOG_PL_VAR <-  zeros(anumber, nfore)
    MSFE_VAR <-  array(0, c(anumber, N * nfocus, nfore))
    MAFE_VAR <-  array(0, c(anumber, N * nfocus, nfore))
    logpl_VAR <-  array(0, c(anumber, N * nfocus, nfore))
}

offset <-  1e-9  # just a constant for numerical stability

#-----------------------------| END OF PRELIMINARIES |---------------------

#=======================| BEGIN KALMAN FILTER ESTIMATION |=================

#for (irep in T_thres:t) {
irep = T_thres # initial observations
# if (irep %/% ceiling(t / 40) == 0) {
#     disp t([num2str(100 * (irep / t))]) # completed'
#     toc
# }

beta_OLS  <- solve(crossprod(x_t[[1]][1:irep,]),  crossprod(x_t[[1]][1:irep,],  y_t[[1]][1:irep,]))
sigma_OLS <- crossprod(y_t[[1]][1:irep,] - x_t[[1]][1:irep,] %*% beta_OLS) / (irep - M[[1]])

#if (irep >= T_thres) {
# Start predictive simulation
chol_S <-  chol(sigma_OLS)
Yraw_f <- vector("list", length = 1)


# simulate forecasts 
for (sim in 1:nsim) {
    Y_hat <-  0
    # Now create forecast for h=1
    X_FORE <- c(1, y_t[[ss]][irep,], x_f[[ss]][irep, 1:M[[ss]] * (p - 1)])
    Y_hat <- X_FORE %*% beta_OLS + rnorm(M[[ss]]) %*% chol_S
    y_fore[[ss]][1, , sim] <-  Y_hat
    
    # Now do forecasts for h>1
    
    for (ii in 1:(nfore - 1)) {
        if (ii <= p) {
            # if h<=p (number of lags)
            
            if (ii == 1) {
                X_new_temp <- c(1, Y_hat, X_FORE[2:(M[[ss]] * (p - ii) + 1)])
            } else{
                X_new_temp <- c(1, Y_hat)
            }
            
            Y_temp <- X_new_temp %*% beta_OLS + rnorm(M[[ss]]) %*% chol_S
            Y_hat <- c(Y_temp, Y_hat)
        } else {
            # if h>p (number of lags)
            X_new_temp <-  c(1, Y_hat[1:(M[[ss]] * p)])
            Y_temp <- X_new_temp %*% beta_OLS + rnorm(M[[ss]]) %*% chol_S
            Y_hat <-  c(Y_temp, Y_hat)
        }
        
        # This cell array saves the draws from all forecast
        # horizons, from all model sizes.
        y_fore[[ss]][ii + 1, , sim] <-  Y_temp
        
    }
} 
# TODO Continue
#-----------------------------| 09.01.2020 |---------------------
#apply(y_fore[[ss]], MARGIN = c(1, 2), FUN=mean)


# Find "observed" out-of-sample data for MSFE and MAFE calculations
if (irep <= t - nfore) {
    Yraw_f[[1]] = y_t[[1]][(irep + 1):(irep + nfore),] #Pseudo out-of-sample observations
} else {
    Yraw_f[[1]] = rbind(y_t[[1]][(irep + 1):t,],  matrix(NA, nrow = nfore - (t - irep), ncol = M[[1]]))
}

# Now we have the predictions for each model & the associated model
# probabilities
y_t_VAR <- vector(mode="numeric") #array(NA, dim = c(nfore, 100, 100))
ii <- 0
# TODO Calc of out-of-sample metrics
#for (ii in 1:nfore) {
ii <- ii+1
    focus_vars = 1
    
    y_t_VAR <- append(y_t_VAR, mean(y_fore[[ss]][ii, focus_vars, ], na.rm =
                                                 TRUE)
    )
    # y_t_VAR[ii, , irep - T_thres + 1] = mean(y_fore[[ss]][ii, focus_vars, ], na.rm =
    #                                              TRUE)
    variance_VAR = var(c(y_fore[[ss]][ii, focus_vars, ,drop=TRUE]))
    #browser()
    # TODO fix mvnpdf to cover 3D-case
    LOG_PL_VAR[irep - T_thres + 1, ii] <-  log(mvnpdfs(t(Yraw_f[[1]][ii, focus_vars]), t(y_t_VAR[ii]), variance_VAR) + offset)
    #LOG_PL_VAR[irep - T_thres + 1, ii] <-  log(mvnpdfs(t(Yraw_f[[1]][ii, focus_vars]), t(y_t_VAR[ii, , irep -T_thres + 1]), variance_VAR) + offset)
    
    MAFE_VAR[irep - T_thres + 1, (country_index - 1) * nfocus +
                 1:country_index * nfocus, ii] <-  abs(Yraw_f[[1]][ii, focus_vars] - squeeze(y_t_VAR[ii, , irep -
                                                                                                       T_thres + 1]))
    MSFE_VAR[irep - T_thres + 1, (country_index - 1) * nfocus +
                 1:country_index * nfocus, ii] <-  (Yraw_f[[1]][ii, focus_vars] - squeeze(y_t_VAR[ii, , irep -
                                                                                                    T_thres + 1])) ^ 2
    
    MAFE_RW[irep - T_thres + 1, (country_index - 1) * nfocus +
                1:country_index * nfocus, ii] <-  abs(t(y_t[[1]][irep, focus_vars]) - Yraw_f[[1]][ii, focus_vars])
    MSFE_RW[irep - T_thres + 1, (country_index - 1) * nfocus +
                1:country_index * nfocus, ii] <-  t(y_t[[1]][irep, focus_vars]) - Yraw_f[[1]][ii, focus_vars] ^
        2
    
    j_in = 0
    
    for (j in 1:focus_vars) {
        j_in <-  j_in + 1
        logpl_VAR[irep - T_thres + 1,
                  (country_index - 1) * nfocus + 1:country_index * nfocus,
                  ii] <-  log(mvnpdfs(t(Yraw_f[[ss]][ii, j]), t(y_t_VAR[ii, j_in, irep -
                                                                          T_thres + 1]), variance_VAR[j_in, j_in]) + offset)
    }
}
}

}
}

# END OF KALMAN FILTER ESTIMATION |============================================

# PRINT RESULTS |==============================================================

# format short g
# summary of forecasting results
# if (print_fore == 1) {
#     disp('MSFE for the key variables of interest')
#     msfe_focus = mean(MSFE_VAR(1:end - 1, :, 1))
#     for (iii = 2:nfore) {
#         msfe_focus = [msfe_focus mean(MSFE_VAR(1:end - iii, :, iii))]
#     }
#
#     disp(['    Horizon    GDP       INFL      INTR']) ##ok<*NBRAK>
#     disp([(1:nfore)' msfe_focus])
#
#     disp('MAFE for the key variables of interest')
#     mafe_focus = mean(MAFE_VAR(1:end-1,:,1))
#
#     for (iii in 2:nfore) {
#         mafe_focus=[mafe_focus mean(MAFE_VAR(1:end-iii,:,iii))]
#     }
#
#     disp(['    Horizon    GDP       INFL      INTR'])
#     disp([(1:nfore)' mafe_focus])
#
#     disp('sum of log pred likes for the key variables of interest individually')
#     lpls_focus = sum(logpl_VAR(1:end - 1, :, 1))
#
#     for (iii in 2:nfore) {
#         lpls_focus = [lpls_focus sum(logpl_VAR(1:end - iii, :, iii))]
#     }
#
#     disp(['    Horizon    GDP      INFL      INTR'])
#     disp([(1:nfore)' lpls_focus])
#     disp('sum of log pred likes for the key variables of interest jointly')
#     lpl_focus = sum(LOG_PL_VAR(1:end-1,1))
#
#     for (iii in 2:nfore){
#         lpl_focus=[lpl_focus sum(LOG_PL_VAR(1:end-iii,iii))]
#     }
#
#     disp(['    Horizon   TOTAL'])
#     disp([(1:nfore)' lpl_focus])
#     disp('                      ')
#     disp('                      ')
#
# }


# plot volatilities and lambda_t (but not theta_t which is huge)
# if (print_coefficients == 1) {
#     prob_pl = omega_update{
#         1, 1
#     }(index_DMA(:, 1))
#     for (ss in 2:nos) {
#         prob_pl = [prob_pl omega_update{
#             ss, 1
#         }(index_DMA(:, ss))] ##ok<*AGROW>
#     }
#
#     final_prob = 0 * prob_pl
#
#     for (sss in 1:nos) {
#         final_prob(:, sss) = prob_pl(:, sss). / sum(prob_pl, 2)
#     }
#     figure
#     plot(yearlab, final_prob)
#     warning off
#     legend('small VAR', 'medium VAR', 'large VAR')
#     title('Time-varying probabilities of small/medium/large VARs')
# }

# plot predictive likelihoods over time
# if (print_pred == 1) {
#     figure
#     subplot(2, 2, 1)
#     plot(yearlab(T_thres + 1:end), cumsum(LOG_PL_DMS(1:end - 1, 1)))
#     title('Cumulative sum of total predictive likelihood h=1')
#     subplot(2, 2, 2)
#     plot(yearlab(T_thres + 2:end), cumsum(LOG_PL_DMS(1:end - 2, 2)))
#     title('Cumulative sum of total predictive likelihood h=2')
#     subplot(2, 2, 3)
#     plot(yearlab(T_thres + 3:end), cumsum(LOG_PL_DMS(1:end - 3, 3)))
#     title('Cumulative sum of total predictive likelihood h=3')
#     subplot(2, 2, 4)
#     plot(yearlab(T_thres + 4:end), cumsum(LOG_PL_DMS(1:end - 4, 4)))
#     title('Cumulative sum of total predictive likelihood h=4')
# }
# print the Minnesota prior over time
# if (print_Min == 1) {
#     print(c(
#         "=========MINNESOTA PRIOR==========",
#         "\n",
#         "Best gammas for each time period"
#     ))
#     vars_Min = t(T_thres:t)
#     for (ss in 1:nos) {
#         vars_Min = cbind(vars_Min  Minn_gamms[ss, 1])
#     }
#     lll = c("period", "smallVAR", "mediumVAR", "largeVAR")
#     print(lll(1:nos + 1))
#     print(vars_Min)
# }

#save(sprintf('#s.mat','competing_forecasts'),'-mat')