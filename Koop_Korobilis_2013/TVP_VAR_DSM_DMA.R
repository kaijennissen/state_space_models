library(rlist)
library(Matrix)
source("./Koop_Korobilis_2013/utils.R")

# Preliminaries -----------------------------------------------------------
forgetting <- 1     # 1: use constant factor;
                    # 2: use variable factor
lambda <- 0.99
kappa <- 0.96       # Decay factor for measurement error variance
eta <- 0.99         # Forgetting factor for DPS (dynamic prior selection) and DMA

# Please choose:
p <- 4             # p is number of lags in the VAR part
nos <- 3           # number of subsets to consider (default is 3, i.e. 3, 7, and 25 variable VARs)

# if nos=1 you might want a single model. Which one is this?
single <- 1         # 1: 3 variable VAR  - 
                    # 2: 7 variable VAR
                    # 3: 25 variable VAR

prior <- 1          # 1: Use Koop-type Minnesota prior
                    # 2: Use Litterman-type Minnesota prior

# Forecasting
first_sample_ends <- 1974.75 # The end date of the first sample in the 
# recursive exercise (default value: 1969:Q4)
nfore <- 8                   # Select forecast horizon
nsim <- 50             # Number of times to simulate from the predictive density

# Choose which results to print
# NOTE: CHOOSE ONLY 0/1 (FOR NO/YES) VALUES!
print_fore <- 1           # summary of forecasting results
print_coefficients <- 1   # plot volatilities and lambda_t (but not theta_t which is huge)
print_pred <- 1           # plot predictive likelihoods over time
print_Min <- 1            # print the Minnesota prior over time


# Load Data ---------------------------------------------------------------

ydata <- read.delim("./Koop_Korobilis_2013/data/ydata.dat", header = FALSE, sep="\t") 
tcode <- read.delim("./Koop_Korobilis_2013/data/tcode.dat", header = FALSE, sep="\t") 
tcode <- tcode[,1,  drop=TRUE]
vars <- list.load("./Koop_Korobilis_2013/vars.rdata")
ynames <- list.load("./Koop_Korobilis_2013/ynames.rdata")

# Create dates variable
start_date <- 1959.00 #1959.Q1
end_date <- 2010.25   #2010.Q2
yearlab <- seq(1959, 2010.25, 0.25)
T_thres <- which(yearlab == first_sample_ends) # find tau_0 (first sample)

# Transform data to stationarity
# Y: standard transformations (for iterated forecasts, and RHS of direct forecasts)
Y <- .transform(ydata, tcode, yearlab)

# Select a subset of the data to be used for the VAR
if (nos > 3){
  print('DMA over too many models, memory concerns...')
}

Y1 <- vector("list", nos)
Ytemp <-  .standardize1(Y,T_thres)
M <-  rep(0, nos)
for (ss in 1:nos) {
  if (nos != 1) {
    single <-  ss
  }
select_subset <-  vars[[single]] # vars stores which variables should be included in which model, 
Y1[[ss]] <-  Ytemp[, select_subset] # select only vars included in the model
M[ss] <-  length(select_subset) # M is the dimensionality of Y
}

t <-nrow(Y1[[1]]) # number of observations?

# The first nfocus variables are the variables of interest for forecasting
nfocus <- 3;

# VAR ---------------------------------------------------------------------
# Generate lagged Y matrix. This will be part of the X matrix
x_t <- vector("list", nos)
x_f <- vector("list", nos)
y_t <- vector("list", nos)
K <- .zeros(nos, 1)

for (ss in 1:nos){
  ylag <- .mlag2(Y1[[ss]], p)
  ylag <- ylag[(p+1):nrow(ylag), ]
  tmp <- .create_RHS(ylag, M[ss], p, t) # returns matrix
  x_t[[ss]] <- tmp[[1]]
  K[ss, 1] <-  tmp[[2]]
  x_f[[ss]] <- ylag
  y_t[[ss]] <- Y1[[ss]][(p+1):nrow(Y1[[ss]]), ]
}

yearlab <- yearlab[(p+1):length(yearlab)]

# Time series observations
t <- nrow(y_t[[1]]) 


# Preliminaries -----------------------------------------------------------
# Set the alpha_bar and the set of gamma values
alpha_bar <-  10
gamma <-  c(1e-10, 1e-5, 0.001, 0.005, 0.01, 0.05, 0.1)
nom <-  max(length(gamma))  # This variable defines the number of DPS models - DYNAMIC PRIOR SHRINKAGE

# prior means and variances (_prmean / _prvar)
theta_0_prmean <-  vector("list", nos)  # TODO create list matrices prior to assignment
theta_0_prvar <-  vector("list", nos)
Sigma_0 <-  vector("list", nos)

for (ss in 1:nos) {
  theta_0_prmean[[ss]] <- matrix(NA, nrow = M[ss] * (1 + p * M[ss]), ncol = nom)
  theta_0_prvar[[ss]] <- array(NA, c(M[ss] * (1 + p * M[ss]), M[ss] * (1 + p * M[ss]), nom))
  Sigma_0[[ss]] <- matrix(NA, nrow =  M[ss], ncol =  M[ss])
  if (prior == 1) {# 1) "No dependence" prior
    for (i in 1:nom) {
      tmp <- .Minn_prior_KOOP(alpha_bar, gamma[i], M[ss], p, K[ss])
      prior_mean <- tmp[[1]]
      prior_var <- tmp[[2]]
      theta_0_prmean[[ss]][, i] <- prior_mean # TODO create list matrices prior to assignment
      theta_0_prvar[[ss]][, , i] <- prior_var     # TODO create list matrices prior to assignment
    }
    Sigma_0[[ss]] <-  cov(y_t[[ss]][1:T_thres,])
  } else if (prior == 2) { # 2) Full Minnesota prior
    print("Minnesota prior is under development")
    #   for (i in 1:nom){
    #     tmp <-  .Minn_prior_LITT(
    #       y_t[[ss]][1:T_thres, ],
    #       x_f[[ss]][1:T_thres, ],
    #       alpha_bar,
    #       gamma(i),
    #       M[ss],
    #       p,
    #       K(ss),
    #       T_thres)
    #     prior_mean <- tmp[[1]]
    #     prior_var <- tmp[[2]]
    #     sigma_var <- tmp[[3]]
    #     theta_0_prmean[[ss]][, i] <-  prior_mean
    #     theta_0_prvar[[ss]][, , i] <-  prior_var
    # }
  # Sigma_0[[ss]] <-  sigma_var # Initialize the measurement covariance matrix (Important!)
  }
}


# Define forgetting factor lambda:
lambda_t <-   vector("list", nos)
for (ss in 1:nos){
  if (forgetting == 1){# CASE 1: Choose the forgetting factor   
    inv_lambda = 1/lambda
    lambda_t[[ss]] <-  lambda*.ones(t, nom)
  }else if (forgetting == 2){# CASE 2: Use a variable (estimated) forgetting factor
    lambda_min <-  0.97
    inv_lambda <-  1/0.99
    alpha <-  1
    LL <-  1.1
    lambda_t[[ss]] <- .zeros(t, nom) # TODO
  }else{
    print('Wrong specification of forgetting procedure')
  }
}

# Initialize matrices
theta_pred <- vector("list", nos)  
theta_update <- vector("list", nos)
R_t <- vector("list", nos)
S_t <- vector("list", nos)
y_t_pred <- vector("list", nos)
e_t <-vector("list", nos)
A_t <- vector("list", nos)
V_t <- vector("list", nos)
y_fore <- vector("list", nos)
omega_update <- vector("list", nos)
omega_predict <- vector("list", nos)
ksi_update <- array(0, c(t,nos))
ksi_predict <- array(0, c(t,nos))
w_t <- vector("list", nos)
w2_t <- array(0, c(t,nos))
f_l <- array(0, c(nom,1))
max_prob <- rep(0, nos)
k_max<- rep(0, nos)
max_prob_DMS <- array(0, c(t,1))
index_best <- array(0, c(t,1))
index_DMA <- array(0, c(t,nos))
sum_prob_omega <- vector("list", nos)
  
anumber <- t-T_thres+1
y_t_DMA <- array(0, c(nfore,nfocus,anumber))
y_t_DMS <- array(0, c(nfore,nfocus,anumber))
LOG_PL_DMA <- array(0, c(anumber,nfore))
MSFE_DMA <- array(0, c(anumber,nfocus,nfore))
MAFE_DMA <- array(0, c(anumber,nfocus,nfore))
LOG_PL_DMS <- array(0, c(anumber,nfore))
MSFE_DMS <- array(0, c(anumber,nfocus,nfore))
MAFE_DMS <- array(0, c(anumber,nfocus,nfore))
logpl_DMA <- array(0, c(anumber,nfocus,nfore))
logpl_DMS <- array(0, c(anumber,nfocus,nfore))
offset <- 1e-9  # just a constant for numerical stability

# TODO create
index_DMS <- rep(0, t)
Yraw_f <-  vector("list", nos)


for(ss in 1:nos){
  theta_pred[[ss]] <- array(NA, c( M[ss] * (1 + p * M[ss]), t, nom))
  theta_update[[ss]] <- array(NA, c( M[ss] * (1 + p * M[ss]), t, nom))
  R_t[[ss]] <- array(NA, c( M[ss] * (1 + p * M[ss]),  M[ss] * (1 + p * M[ss]), nom))
  omega_predict[[ss]] <- matrix(NA, nrow=t, ncol=nom)
  omega_update[[ss]] <- matrix(NA, nrow=t, ncol=nom)
  y_t_pred[[ss]] <- array(NA, c( M[ss], t, nom))
  e_t[[ss]] <-  array(NA, c( M[ss], t, nom))
  #lambda_t[[ss]] <- matrix(NA, nrow=t, ncol=nom)
  V_t[[ss]] <- array(NA, c(M[ss], M[ss], t, nom))
  S_t[[ss]] <- array(NA, c(M[ss] * (1 + p * M[ss]), M[ss] * (1 + p * M[ss]), nom))
  w_t[[ss]] <- array(NA, c(1, t, nom))
  Yraw_f[[ss]] <- matrix(NA, nrow=nfore, ncol=M[ss])
  y_fore[[ss]] <- array(NA, c(nfore, M[ss], nom, nsim))
}




now <- Sys.time()
# Kalman Filter -----------------------------------------------------------
for (irep in 1:t){# Loop over time axis
    print(paste("irep:", irep, "of",  t, "===================================="))

    for(ss in 1:nos){# LOOP FOR 1 TO NOS VAR MODELS OF DIFFERENT DIMENSIONS - DIFFERENT SIZE OF VARS
        #
      print(paste("ss:", ss))
        # Find sum of probabilities for DPS
        if (irep>1){
            sum_prob_omega[ss] <- sum((omega_update[[ss]][irep-1, ])^eta) # this is the sum of the nom model probabilities (all in the power of the forgetting factor 'eta')
        }
        for (k in 1:nom){ # LOOP FOR 1 TO NOM VAR MODELS WITH DIFFERENT DEGREE OF SHRINKAGE - DIFFERENT DEGREE OF SHRINKAGE
            #print(paste("k:", k))  
            # Predict   
            if (irep == 1) { # start recursion
              theta_pred[[ss]][, irep, k] <-  theta_0_prmean[[ss]][, k]
              R_t[[ss]][, , k] <-  theta_0_prvar[[ss]][ , , k]
              omega_predict[[ss]][irep, k] <-  1/nom
            } else  {             
              theta_pred[[ss]][, irep, k] <-  theta_update[[ss]][, irep-1, k]
              R_t[[ss]][ , , k] <-  (1/lambda_t[[ss]][irep-1, k]) * S_t[[ss]][, , k]
              omega_predict[[ss]][irep, k] <- ((omega_update[[ss]][irep-1, k])^eta + offset)/(sum_prob_omega[[ss]] + offset)
            }
      
            xx <-  x_t[[ss]][((irep-1)*M[ss]+1):(irep*M[ss]), , drop=FALSE]
            y_t_pred[[ss]][, irep, k] <-  xx%*%theta_pred[[ss]][, irep, k]  # this is one step ahead prediction
      
            # Prediction error
            e_t[[ss]][, irep, k] <-  y_t[[ss]][irep, , drop=TRUE] - y_t_pred[[ss]][, irep, k, drop=TRUE]  # this is one step ahead prediction error
          
            # Update forgetting factor
            if (forgetting == 2){ # update factors in CASE 2
                  lambda_t[[ss]][irep, k] <- lambda_min + (1-lambda_min)*(LL^(-round(alpha*crossprod(e_t[[ss]][1:nfocus, irep, k]))))
              }

            # first update V[t], see the part below equation (10)
            A_t <- tcrossprod(e_t[[ss]][,irep,k]) 
            if (irep==1){
              V_t[[ss]][, , irep, k] <-  kappa*Sigma_0[[ss]]
            }else{
              V_t[[ss]][, , irep, k] <-  kappa*V_t[[ss]][ , , irep-1, k, drop=TRUE] + (1-kappa)*A_t
            }
            #          if all(eig(squeeze(V_t(:,:,irep,k))) < 0)
            #              V_t(:,:,irep,k[ <-V_t(:,:,irep-1,k);       
            #          }
            # 
            # % update theta[t] and S[t] - what is S == Sigma?
      
            Rx <-  tcrossprod(R_t[[ss]][, , k], xx)
            KV <-  V_t[[ss]][, , irep, k] + xx%*%Rx
            KG <- t(solve(t(KV), t(Rx))) #  all.equal(t(solve(t(KV), t(Rx))), Rx%*%solve(KV))
            theta_update[[ss]][, irep, k] <-  theta_pred[[ss]][, irep, k] + (KG%*%e_t[[ss]][, irep, k])
            S_t[[ss]][ , ,k] <-  R_t[[ss]][, , k] - KG%*%(xx%*%R_t[[ss]][, , k])
    
            # Find predictive likelihood based on Kalman filter and update DPS weights     
            if (irep == 1){
              variance <-  Sigma_0[[ss]][1:nfocus,1:nfocus] + xx[1:nfocus, ]%*%Rx[,1:nfocus]
            }else{
              variance <-  V_t[[ss]][1:nfocus,1:nfocus,irep,k] + xx[1:nfocus, ]%*%Rx[,1:nfocus]
            }

            if (any(eigen(variance)$values <= 0)) { # replace negative eigenvalues for invertability
              variance <-  abs(diag(diag(variance)))       
            }
      
            ymean <-  y_t_pred[[ss]][1:nfocus, irep, k]
            ytemp <-  t(y_t[[ss]][irep, 1:nfocus])
            f_l[k] <-  .mvnpdfs(ytemp, ymean, variance)  # likelihood calculation for DMA / DMS
            w_t[[ss]][, irep, k] <- omega_predict[[ss]][irep, k]*f_l[k]
      

            # Individual model forecasts for DPS
            if (irep >= T_thres){
                #if (ss==2){
                  #browser()
                #}
                chol_var_beta <-  t(chol(S_t[[ss]][ , ,k]))
                chol_var_y <-  chol(V_t[[ss]][, , irep, k])
             
                # Start predictive simulation
                for (sim in 1:nsim){
                    #print(paste("sim: ", sim))
                    Y_hat <-  0
                    # Draw from the posterior of the coefficients, theta
                    beta_sim  <-  theta_update[[ss]][, irep, k] + chol_var_beta%*%rnorm(K[ss])
                    bbtemp <-  beta_sim[(M[ss]+1):K[ss]]  # get the draw of B(t) at time i=1,...,T  (exclude intercept)
                    splace <- 0
                    biga <- matrix(0, nrow=M[ss], ncol=p*M[ss])
      
                    for (ii in 1:p){
                        for (iii in 1:M[ss]){
                            biga[iii, ((ii-1)*M[ss]+1):(ii*M[ss])] <- bbtemp[(splace+1):(splace+M[ss])]
                            splace <-  splace + M[ss]
                        }
                    }
      
                  beta_fore <-  rbind(t(beta_sim[1:M[ss]]), t(biga))
                  # Now create forecast for h=1
                  X_FORE <-  c(1, y_t[[ss]][irep, ], x_f[[ss]][irep, 1:(M[ss]*(p-1))])
                  Y_hat <-  X_FORE%*%beta_fore + rnorm(M[ss])%*%chol_var_y
                  y_fore[[ss]][1, , k, sim] <-  Y_hat
      
                  # Now do forecasts for h>1
                  for (ii in 1:(nfore-1)){
                        if (ii < p) {  # if h<=p (number of lags)
                            X_new_temp <-  c(1, Y_hat, X_FORE[2:(M[ss]*(p-ii)+1)])
                            Y_temp <-  X_new_temp%*%beta_fore + rnorm(M[ss])%*%chol_var_y
                            Y_hat <-  c(Y_temp, Y_hat)
                        }else if (ii==p){
                          X_new_temp <-  c(1, Y_hat)
                          Y_temp <-  X_new_temp%*%beta_fore + rnorm(M[ss])%*%chol_var_y
                          Y_hat <-  c(Y_temp, Y_hat)
                        }else{  #% if h>p (number of lags)
                            X_new_temp <-  c(1, Y_hat[1:(M[ss]*p)])
                            Y_temp <-  X_new_temp%*%beta_fore + rnorm(M[ss])%*%chol_var_y
                            Y_hat <-  c(Y_temp, Y_hat)
                        }
                      # This cell array saves the draws from all forecast
                      # horizons, from all model sizes.
                      y_fore[[ss]][ii+1, , k, sim] <-  Y_temp
                  }
                } # End predictive simulation
              
            } # End cycling through nom models with different shrinkage factors
        }
      # First calculate the denominator of Equation (19) (the sum of the w's)
      sum_w_t <-  0   
      for (k_2 in 1:nom){
        sum_w_t <-  sum_w_t + w_t[[ss]][, irep, k_2]
      }

      # Then calculate the DPS probabilities  
      for (k_3 in 1:nom){
        omega_update[[ss]][irep, k_3] <- (w_t[[ss]][, irep, k_3] + offset) / (sum_w_t + offset)  # this is Equation (19)
      }
      max_prob[ss] <- max(omega_update[[ss]][irep,])
      k_max[ss] <- which(omega_update[[ss]][irep,] ==  max_prob[ss])
      index_DMA[irep, ss] <- k_max[ss]

      # Use predictive likelihood of best (DPS) model, and fight the weight for DMA     
      w2_t[irep, ss] <- omega_predict[[ss]][irep, k_max[ss]]*f_l[k_max[ss],1]
  
      # Find "observed" out-of-sample data for MSFE and MAFE calculations
      if (irep <= t-nfore){
        Yraw_f[[ss]] <- y_t[[ss]][(irep+1):(irep+nfore), ] #Pseudo out-of-sample observations                   
      }else if (irep>=t){
        Yraw_f[[ss]] <-  matrix(NA, nrow=nfore-(t-irep), ncol= M[ss])
        }else{
        # TODO fix out of bound error
        Yraw_f[[ss]] <- rbind(y_t[[ss]][(irep+1):t, ] , matrix(NA, nrow=nfore-(t-irep), ncol= M[ss]) )
      }
    }
  
      # First calculate the denominator of Equation (19) (the sum of the w's)
      sum_w2_t <-  0
      for (k_2 in 1:nos){
          sum_w2_t <-  sum_w2_t + w2_t[irep, k_2]
      }
  
        # Then calculate the DPS probabilities
      for (k_3 in 1:nos){
          ksi_update[irep, k_3] <-  (w2_t[irep, k_3] + offset) / (sum_w2_t + offset)  # this is Equation (19)
      }

      # Find best model for DMS
      max_prob_DMS[irep, ] <- max(ksi_update[irep, ])
      ss_max <- which(ksi_update[irep, ] == max_prob_DMS[irep,])
      index_DMS[irep] <- k_max[ss_max]
      
      # Now we cycled over NOM and NOS models, do DMA-over-DPS
      # if (irep >= T_thres){
      #     # Now we have the predictions for each model & the associated model probabilities
      #     for (ii in 1:nfore){
      #         weight_pred <- 0*y_fore[[ss]]ii, 1:nfocus, k_max[[ss]], ]
      #         # DPS-DMA prediction
      #         for (ss in 1:nos){
      #             temp_predict = y_fore[[ss]](ii,1:nfocus,k_max{ss},:).*ksi_update(irep,ss);
      #             weight_pred = weight_pred + temp_predict;
      #             Minn_gamms[[ss]](irep-T_thres+1,1[ <- gamma(1,k_max[[ss]]);
      #         }
      #         
      #         y_t_DMA[ii,:,irep-T_thres+1] <-  mean(weight_pred,4);
      #         variance_DMA <- cov(squeeze(weight_pred)#')
      # 
      #         y_t_DMS[ii, , irep-T_thres+1] = mean(y_fore[[ss_max]][ii, 1:nfocus, k_max[[ss_max]],], 4)
      #         variance_DMS <- cov(squeeze(y_fore{ss_max,1}(ii,1:nfocus,k_max{ss_max},:))#')
      #         
      #         
      #         LOG_PL_DMS[irep-T_thres+1,ii] <- log(mvnpdfs(Yraw_f{ss}(ii,1:nfocus)#',y_t_DMS(ii,:,irep-T_thres+1)',variance_DMS) + offset);
      #         MAFE_DMS[irep-T_thres+1,:,ii] <- abs(Yraw_f{ss}(ii,1:nfocus) - squeeze(y_t_DMS(ii,:,irep-T_thres+1)));
      #         MSFE_DMS[irep-T_thres+1,:,ii] <- (Yraw_f{ss}(ii,1:nfocus) - squeeze(y_t_DMS(ii,:,irep-T_thres+1))).^2;
      #         
      #         
      #         LOG_PL_DMA[irep-T_thres+1,ii] <-log(.mvnpdfs(Yraw_f{ss}(ii,1:nfocus)#',y_t_DMA(ii,:,irep-T_thres+1)',variance_DMA) + offset);
      #         MAFE_DMA[irep-T_thres+1,:,ii] <-abs(Yraw_f{ss}(ii,1:nfocus) - squeeze(y_t_DMA(ii,:,irep-T_thres+1)));
      #         MSFE_DMA[irep-T_thres+1,:,ii] <-(Yraw_f{ss}(ii,1:nfocus) - squeeze(y_t_DMA(ii,:,irep-T_thres+1))).^2;
      #         for (j in 1:nfocus){
      #             logpl_DMA[irep-T_thres+1,j,ii] <- log(.mvnpdfs(Yraw_f{ss}(ii,j)#',y_t_DMA(ii,j,irep-T_thres+1)',variance_DMA(j,j)) + offset);
      #             logpl_DMS[irep-T_thres+1,j,ii] <- log(.mvnpdfs(Yraw_f{ss}(ii,j)#',y_t_DMS(ii,j,irep-T_thres+1)',variance_DMS(j,j)) + offset);
      #         }      
      #     }
      # 
      
}



print(Sys.time()-now)

# Plots -------------------------------------------------------------------

prob_pl <- matrix(NA, nrow=t, ncol=nos)
for (ss in 1:nos){
  for (i in 1:t){
    prob_pl[i,ss] <-  omega_update[[ss]][i,index_DMA[i,ss]]
  }
}
final_prop <- prob_pl / rowSums(prob_pl)


plot.ts(final_prop,
        plot.type = "single",
        col=c("red", "blue", "green"))
legend("topleft",
       legend=c("small VAR", "medium VAR", "large VAR"),
       col=c("red", "blue", "green"),
       lty=1, 
       cex=0.6)




# % summary of forecasting results
# if print_fore == 1;
# disp('==============| DMA RESULTS |==================')
# disp('MSFE for the key variables of interest')
# msfe_focus = mean(MSFE_DMA(1:end-1,1:nfocus,1));
# for iii = 2:nfore
#     msfe_focus = [msfe_focus; mean(MSFE_DMA(1:end-iii,1:nfocus,iii))];
# }
# disp(['    Horizon    GDP       INFL      INTR']) %#ok<*NBRAK>
# disp([(1:nfore)' sqrt(msfe_focus)])
# 
# disp('MAFE for the key variables of interest')
# mafe_focus = mean(MAFE_DMA(1:end-1,1:nfocus,1));
# for iii = 2:nfore
# mafe_focus=[mafe_focus; mean(MAFE_DMA(1:end-iii,1:nfocus,iii))];
# end
# disp(['    Horizon    GDP       INFL      INTR'])
# disp([(1:nfore)' mafe_focus])
# 
# disp('sum of log pred likes for the key variables of interest individually')
# lpls_focus = sum(logpl_DMA(1:end-1,:,1));
# for iii = 2:nfore
#     lpls_focus = [lpls_focus; sum(logpl_DMA(1:end-iii,:,iii))];
# end
# disp(['    Horizon    GDP      INFL      INTR'])
# disp([(1:nfore)' lpls_focus])
# disp('sum of log pred likes for the key variables of interest jointly')
# lpl_focus = sum(LOG_PL_DMA(1:end-1,1));
# for iii=2:nfore
# lpl_focus=[lpl_focus; sum(LOG_PL_DMA(1:end-iii,iii))];
# end
# disp(['    Horizon   TOTAL'])
# disp([(1:nfore)' lpl_focus])
# disp('                      ')
# disp('                      ')
# 
# disp('==============| DMS RESULTS |==================')
# disp('MSFE for the key variables of interest')
# msfe_focus = mean(MSFE_DMS(1:end-1,1:nfocus,1));
# for iii = 2:nfore
#     msfe_focus = [msfe_focus; mean(MSFE_DMS(1:end-iii,1:nfocus,iii))];
# end
# disp(['    Horizon    GDP       INFL      INTR'])
# disp([(1:nfore)' msfe_focus])
# 
# disp('MAFE for the key variables of interest')
# mafe_focus = mean(MAFE_DMS(1:end-1,1:nfocus,1));
# for iii = 2:nfore
# mafe_focus=[mafe_focus; mean(MAFE_DMS(1:end-iii,1:nfocus,iii))];
# end
# disp(['    Horizon    GDP       INFL      INTR'])
# disp([(1:nfore)' mafe_focus])
# 
# disp('sum of log pred likes for the key variables of interest individually')
# lpls_focus = sum(logpl_DMS(1:end-1,:,1));
# for iii = 2:nfore
#     lpls_focus = [lpls_focus; sum(logpl_DMS(1:end-iii,:,iii))];
# end
# disp(['    Horizon    GDP      INFL      INTR'])
# disp([(1:nfore)' lpls_focus])
# disp('sum of log pred likes for the key variables of interest jointly')
# lpl_focus = sum(LOG_PL_DMS(1:end-1,1));
# for iii=2:nfore
# lpl_focus=[lpl_focus; sum(LOG_PL_DMS(1:end-iii,iii))];
# end
# disp(['    Horizon   TOTAL'])
# disp([(1:nfore)' lpl_focus])
# end
# 
# 
# % plot volatilities and lambda_t (but not theta_t which is huge)
# if print_coefficients == 1;
#     prob_pl = omega_update{1,1}(index_DMA(:,1));
#     for ss=2:nos
#         prob_pl = [ prob_pl omega_update[[ss]](index_DMA(:,ss))]; %#ok<*AGROW>
#     end
#     final_prob = 0*prob_pl;
#     for sss = 1:nos
#         final_prob(:,sss[ <-prob_pl(:,sss)./sum(prob_pl,2);
#     end
#     figure
#     plot(yearlab,final_prob)
#     warning off;
#     legend('small VAR','medium VAR','large VAR')
#     title('Time-varying probabilities of small/medium/large VARs')
#     if forgetting == 2;
#         lam_pl = lambda_t{1,1}(index_DMA(:,1));
#         for ss=2:nos
#             lam_pl = [lam_pl lambda_t[[ss]](index_DMA(:,ss))];
#         end
#         figure
#         plot(yearlab,lam_pl)
#         legend('small VAR','medium VAR','large VAR')
#         title('This is the value attained by the time varying forgetting factor \lambda')
#     end
# end
# % plot predictive likelihoods over time
# if print_pred == 1;
#     figure
#     subplot(2,2,1)
#     plot(yearlab(T_thres+1:end),cumsum(LOG_PL_DMS(1:end-1,1)))
#     title('Cumulative sum of total predictive likelihood h=1')
#     subplot(2,2,2)
#     plot(yearlab(T_thres+2:end),cumsum(LOG_PL_DMS(1:end-2,2)))
#     title('Cumulative sum of total predictive likelihood h=2')
#     subplot(2,2,3)
#     plot(yearlab(T_thres+3:end),cumsum(LOG_PL_DMS(1:end-3,3)))       
#     title('Cumulative sum of total predictive likelihood h=3')
#     subplot(2,2,4)
#     plot(yearlab(T_thres+4:end),cumsum(LOG_PL_DMS(1:end-4,4)))    
#     title('Cumulative sum of total predictive likelihood h=4')
# end
# % print the Minnesota prior over time
# if print_Min == 1;
#     disp('=========MINNESOTA PRIOR==========')
#     disp('Best gammas for each time period')
#     vars_Min = (T_thres:t)';
#       for ss=1:nos
#       vars_Min = [vars_Min  Minn_gamms[[ss]]];
#       end
#       lll={'period', 'smallVAR', 'mediumVAR', 'largeVAR'};
#       disp(lll(1:nos+1));
#       disp(vars_Min);
#       end
#       
#       
#       save(sprintf('%s_%g_%g_%g_%g_%g.mat','TVP_VAR_DPS_DMA_sim',nos,single,forgetting,kappa,lambda),'-mat');