## LOAD LIBRARIES =============================================================
library("tidyquant")
library("tidyverse")
library("tsibble")
library("tibbletime")
library("timetk")
library("dlm")
library("forecast")
library("xts")
library("tictoc")
# library("foreach")
# library("doMC")
# library("doParallel")


## LOAD DATA ==================================================================
symbols <- c("XOM", "CVX", "SLB", "COP", "OXY", "VLO")
prices <- getSymbols(symbols,
  src = "yahoo",
  from = "1990-07-31",
  to = "2018-12-31",
  auto.assign = TRUE, warnings = FALSE
) %>%
  map(~ Ad(get(.))) %>%
  reduce(merge) %>%
  `colnames<-`(symbols)

w <- c(0.3, 0.3, 0.10, 0.10, 0.10, 0.10)

asset_returns <- prices %>%
  to.monthly(indexAt = "lastof", OHLC = FALSE) %>%
  tk_tbl(preserve_index = TRUE, rename_index = "date") %>%
  gather(asset, returns, -date) %>%
  group_by(asset) %>%
  mutate(returns = (log(returns) - log(lag(returns)))) %>%
  na.omit()

portfolio_returns <- asset_returns %>%
  tq_portfolio(
    assets_col = asset,
    returns_col = returns,
    weights = w,
    col_rename = "returns",
    rebalance_on = "months"
  )



temp <- tempfile()

base <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/"

factor <- "Global_3_Factors"

format <- "_CSV.zip"

full_url <- str_c(base, factor, format, sep = "")

download.file(
  full_url,
  temp,
  quiet = TRUE
)


Global_3_Factors <- read_csv(unz(temp, "Global_3_Factors.csv"), skip = 6) %>%
  rename(DATE = X1) %>%
  mutate_at(vars(-DATE), as.numeric) %>%
  filter(nchar(date) == 6) %>%
  mutate(DATE = rollback(ymd(parse_date_time(DATE, "%Y%m") + months(1)))) %>%
  filter(DATE >= first(portfolio_returns$DATE) &
    DATE <= last(portfolio_returns$DATE))

ff_portfolio_returns <- portfolio_returns %>%
  left_join(Global_3_Factors, by = "DATE") %>%
  rename(MKT_RF = `Mkt-RF`, R = returns) %>%
  mutate(
    MKT_RF = MKT_RF / 100,
    SMB = SMB / 100,
    HML = HML / 100,
    RF = RF / 100,
    R_EX = round(R - RF, 4)
  )


# define a rolling ff model with tibbletime


## ROLLING BETAS ==============================================================

roll_fama_x <- function(tbl_df, window = 12) {
  rolling_lm <- rollify(.f = function(R_EX, MKT_RF, SMB, HML) {
    lm(R_EX ~ MKT_RF + SMB + HML)
  }, window = window, unlist = FALSE)

  rolling_ff_betas <- ff_portfolio_returns %>%
    mutate(rolling_ff = rolling_lm(R_EX, MKT_RF, SMB, HML)) %>%
    slice(-1:-c(window - 1)) %>%
    mutate(tidied = map(rolling_ff, tidy, conf.int = F)) %>%
    unnest(tidied) %>%
    select(date, term, estimate) %>%
    filter(term != "(Intercept)") %>%
    rename(beta = estimate, factor = term) # %>%
  # spread(factor, beta)

  return(rolling_ff_betas)
}

args <- list(window = c(12, 24, 48), tbl_df = ff_portfolio_returns)
ff_roll_betas <- args %>%
  pmap(roll_fama_x) %>%
  reduce(left_join, by = c("date", "factor")) %>%
  rename(
    beta_12 = beta.x,
    beta_24 = beta.y,
    beta_48 = beta
  )

ff_roll_betas %>%
  gather("window", "beta", -c(date, factor)) %>%
  separate(window, c("a", "window")) %>%
  select(-a) %>%
  ggplot(aes(x = date, y = beta, col = factor)) +
  geom_line() +
  facet_wrap(~window) +
  labs(title = "Rolling FF Factor Betas", x = "rolling betas") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 30)
  )

ff_roll_betas %>%
  gather("window", "beta", -c(date, factor)) %>%
  separate(window, c("a", "window")) %>%
  select(-a) %>%
  ggplot(aes(x = date, y = beta, col = window)) +
  geom_line() +
  facet_wrap(~factor) +
  labs(title = "Rolling FF Factor Betas", x = "rolling betas") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 30)
  )

## DLM BETAS ==================================================================

y <- ff_portfolio_returns %>% pull(R_EX)
X <- ff_portfolio_returns %>%
  select(MKT_RF, SMB, HML) %>%
  as.matrix()

build_reg_seas <- function(param) {
  dlm1 <- dlmModReg(X = X, addInt = FALSE) + dlmModSeas(frequency = 12)
  V(dlm1) <- exp(param[1])
  diag(W(dlm1))[c(1:4)] <- exp(param[2:5])
  return(dlm1)
}

build_reg <- function(param) {
  dlm1 <- dlmModReg(X = X, addInt = FALSE)
  V(dlm1) <- exp(param[1])
  diag(W(dlm1)) <- exp(param[2:4])
  return(dlm1)
}

fit <- dlmMLE(y = y, parm = rnorm(4), build = build_reg)
fit$convergence
fit$par

mod1 <- build_reg(fit$par)
W(mod1)
V(mod1)

yFilter <- dlmFilter(y, mod1)
ySmooth <- dlmFilter(y, mod1)
res <- residuals(yFilter, sd = FALSE)

checkresiduals(res)
shapiro.test(res)
sapply(4:42, function(x) Box.test(x = res, lag = x, type = "Ljung", fitdf = 4)$p.value)

res <- y - c(mod1$FF %*% mod1$GG %*% t(yFilter$m))[-1]
plot.ts(res)

colnames(yFilter$m)[c(1:3)] <- colnames(X)
dlm_betas <- as_tibble(yFilter$m[-1, ]) %>% mutate(date = ff_portfolio_returns$date)

dlm_betas %>%
  gather("factor", "beta", -date) %>%
  left_join(ff_roll_betas, by = c("date", "factor")) %>%
  rename(dlm = beta) %>%
  gather("method", "beta", -c(date, factor)) %>%
  separate(method, c("del", "method"), fill = "left") %>%
  select(-del) %>%
  ggplot(aes(x = date, y = beta, col = method)) +
  geom_line() +
  facet_wrap(~factor)

plot.ts(yFilter$m[-c(1:50), ], plot.type = "s", col = c(1:3))

ggplot(aes(x = DATE, y = beta, col = factor)) +
  geom_line()

## MULTIVARIATE DLM ===========================================================

select_series <- c(9, 17, 21)
ids <- M4_test[select_series] %>% map(pluck("st")) %>% flatten_dbl()
mv_ts <- M4_test[select_series] %>% 
  map(pluck("x_ts")) %>% 
  reduce(ts.union) 
colnames(mv_ts) <- paste0("id_",ids)
y <- mv_ts
ccf(y[,1], y[,2])

y_sub <- subset(y, end = 60)
ndim <- dim(y_sub)[2]

param_length <- 
  
  build_SUTSE <- function(param) {
    m <- 3
    n_poly <- 2
    
    dlm1 <- dlmModPoly(order = n_poly) + dlmModSeas(frequency = 12)
    # dlm1%+%dlm1
    dlm1$FF <- dlm1$FF %x% diag(m)
    dlm1$GG <- dlm1$GG %x% diag(m)
    dlm1$m0 <- rep(0, 2*m)
    dlm1$C0 <- dlm1$C0 %x% diag(m) 
    diag(dlm1$C0) <- 1e7
    
    param_seq <- seq(0, length(param))%/%m
    ## V
    dlm1$V <- dlm1$V %x% diag(m) 
    U0 <- matrix(0, nrow = 3, ncol = 3)
    U0[upper.tri(U0)] <- param[param_seq == 0] # m elements; m+1:2m
    diag(U0) <- exp(0.5 * param[param_seq == 7])# m elments 2*m+1:3*m
    dlm1$V  <- crossprod(U0)
    # diag(dlm1$V)[1:m] <- exp(param[seq(0, length(param))%/%m==0]) # m elemnts; 1:m
    
    ## W
    ## cov theta
    U1 <- matrix(0, nrow = 3, ncol = 3)
    U1[upper.tri(U1)] <- param[param_seq == 1] # m elements; m+1:2m
    diag(U1) <- exp(0.5 * param[param_seq == 2])# m elments 2*m+1:3*m
    U1 <- crossprod(U1)
    
    ## cov beta
    U2 <- matrix(0, nrow = 3, ncol = 3)
    U2[upper.tri(U2)] <- param[param_seq == 3] # m elements; 3m+1:4m
    diag(U2) <- exp(0.5 * param[param_seq == 4]) # m elments 4m+1:5*m
    U2 <- crossprod(U2)
    
    ## cov seasonality
    U3 <- matrix(0, nrow = 3, ncol = 3)
    U3[upper.tri(U3)] <- param[param_seq == 5] # m elements; 5m+1:6m
    diag(U3) <- exp(0.5 * param[param_seq == 6]) # m elements; 6m+1:7m
    U3 <- crossprod(U3)
    
    dlm1$W <- dlm1$W %x% diag(m)  
    dlm1$W[1:(3*m), 1:(3*m)] <- bdiag(U1, U2, U3)
    
    
    return(dlm1)
  }
build_dlm <- function(param) {
  m <- 3
  n_poly <- 2
  
  dlm1 <- dlmModPoly(order = n_poly) + dlmModSeas(frequency = 12)
  # dlm1%+%dlm1
  dlm1$FF <- dlm1$FF %x% diag(m)
  dlm1$GG <- dlm1$GG %x% diag(m)
  dlm1$m0 <- rep(0, 2*m)
  dlm1$C0 <- dlm1$C0 %x% diag(m) 
  diag(dlm1$C0) <- 1e7
  
  param_seq <- seq(0, length(param))%/%m
  div_seq <- 0
  ## V
  dlm1$V <- dlm1$V %x% diag(m) 
  # U0 <- matrix(0, nrow = 3, ncol = 3)
  # U0[upper.tri(U0)] <- param[param_seq == div_seq]
  # div_seq <- div_seq+1
  # diag(U0) <- exp(0.5 * param[param_seq == div_seq])
  # dlm1$V  <- crossprod(U0)
  # diag(dlm1$V)[1:m] <- exp(param[seq(0, length(param))%/%m==0]) # m elemnts; 1:m
  
  ## W
  ## cov theta
  U1 <- matrix(0, nrow = 3, ncol = 3)
  # U1[upper.tri(U1)] <- param[param_seq == 1] # m elements; m+1:2m
  # diag(U1) <- exp(0.5 * param[param_seq == 2])# m elments 2*m+1:3*m
  # U1 <- crossprod(U1)
  
  ## cov beta
  U2 <- matrix(0, nrow = 3, ncol = 3)
  div_seq <- div_seq+1
  U2[upper.tri(U2)] <- param[param_seq == div_seq]
  div_seq <- div_seq+1
  diag(U2) <- exp(0.5 * param[param_seq == div_seq])
  U2 <- crossprod(U2)
  
  ## cov seasonality
  U3 <- matrix(0, nrow = 3, ncol = 3)
  div_seq <- div_seq+1
  U3[upper.tri(U3)] <- param[param_seq == div_seq]
  div_seq <- div_seq+1
  diag(U3) <- exp(0.5 * param[param_seq == div_seq])
  U3 <- crossprod(U3)
  dlm1$W <- dlm1$W %x% diag(m)  
  dlm1$W[1:(3*m), 1:(3*m)] <- bdiag(U1, U2, U3)
  
  return(dlm1)
}

par_dlm1 <- dlmMLE(y = y_sub, parm = rnorm(15), build = build_dlm, method = "SANN")
if (par_dlm1$convergence == 0) "MLE converged" else ("did not converge")
par_dlm1$par
dlm1 <- build_dlm(par_dlm1$par)
W(dlm1)[1:9, 1:9]
cov2cor(V(dlm1))
cov2cor(W(dlm1)[4: 6, 4 : 6])
y_Filter <- dlmFilter(y_sub, mod = dlm1)
y_Smooth <- dlmSmooth(y_sub, mod = dlm1)

## forecast 
fcast <- dlmForecast(mod = y_Filter, nAhead = 12)
plot_cov <- plot_dlm(y_Filter, fcast, y_true = y, nseries = 3)
plot_cov

## BAYESIAN DLM ===============================================================

model <- Arima(rnorm(100), order = c(0, 0, 0), seasonal = list(c(1, 1, 0), period = 12))

foo <- simulate(model, nsim = 1000)
fit <- Arima(foo, order = c(1,1,1), seasonal=c(1,1,1))

## Gamma-Prior & Gibbs-Sampler
y <- prices[,1]

set.seed(5678)
MCMC <- 12000
gibbsOut <- dlmGibbsDIG(y, mod = dlmModPoly(order = 2),
                        a.y = 1, b.y = 1000, a.theta = 10, b.theta = 1000,
                        n.sample = MCMC, thin = 1, save.states = TRUE)
# run again with save.states = TRUE
# variance required for smoothing/filtering?
# how to handle changing variance?
# how to handle multiple series with covariance?



## prior hyperparameters
delta0 <- delta2 <- 3
delta1 <- 100
V0 <- (delta0 - 2) * diag(c(10^2, 500^2))
Wmu0 <- (delta1 - 2) * diag(0.01^2, 2)
Wbeta0 <- (delta2 - 2) * diag(c(5^2, 100^2))

## Gibbs sampling
MC <- 3000
TT <- nrow(y)
gibbsTheta <- array(0, dim = c(TT + 1, 4, MC - 1))
gibbsV <- array(0, dim = c(2, 2, MC))
gibbsWmu <- array(0, dim = c(2, 2, MC))
gibbsWbeta <- array(0, dim = c(2, 2, MC))
mod <- dlm(
  FF = matrix(c(1, 0), nrow = 1) %x% diag(2),
  V = diag(2),
  GG = matrix(c(1, 0, 1, 1), 2, 2) %x% diag(2),
  W = bdiag(diag(2), diag(2)),
  m0 = rep(0, 4),
  C0 = diag(x = 1e7, nrow = 4)
)

# starting values
mod$V <- gibbsV[, , 1] <- V0 / (delta0 - 2)
gibbsWmu[, , 1] <- Wmu0 / (delta1 - 2)
gibbsWbeta[, , 1] <- Wbeta0 / (delta2 - 2)
mod$W <- bdiag(gibbsWmu[, , 1], gibbsWbeta[, , 1])

tic()
set.seed(3420)
for (it in 1:(MC - 1)) {
  
  # generate states - FFBS
  modFilt <- dlmFilter(y, mod, simplify = TRUE)
  gibbsTheta[, , it] <- theta <- dlmBSample(modFilt)
  
  # update V
  S <- crossprod(y - theta[-1, ] %*% t(mod$FF)) + V0
  gibbsV[, , it + 1] <- solve(rwishart(df = delta0 + 1 + TT, p = 2, Sigma = solve(S)))
  mod$V <- gibbsV[, , it + 1]
  
  # update Wmu and Wbeta
  theta.center <- theta[-1, ] - (theta[-(TT + 1), ] %*% t(mod$GG))
  SS1 <- crossprod(theta.center)[1:2, 1:2] + Wmu0
  SS2 <- crossprod(theta.center)[3:4, 3:4] + Wbeta0
  gibbsWmu[, , it + 1] <- solve(rwishart(df = delta1 + 1 + TT, Sigma = solve(SS1)))
  gibbsWbeta[, , it + 1] <- solve(rwishart(df = delta2 + 1 + TT, Sigma = solve(SS2)))
  mod$W <- bdiag(gibbsWmu[, , it + 1], gibbsWbeta[, , it + 1])
}
toc()

## MCMC diagnostics
burn <- 1:2000
# traces
ts.plot(sqrt(gibbsV[1, 1, -burn]), xlab = "iteration", ylab = expression(sigma[1]))
ts.plot(sqrt(gibbsV[2, 2, -burn]), xlab = "iteration", ylab = expression(sigma[2]))
ts.plot(gibbsV[2, 1, -burn], xlab = "iteration", ylab = expression(sigma[1][2]))
# V
par(mar = c(2, 4, 1, 1) + 0.1, cex = 1)
par(mfrow = c(3, 2))
plot(ergMean(sqrt(gibbsV[1, 1, -burn])),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[1]), xlab = "MCMC iteration"
)
acf(sqrt(gibbsV[1, 1, -burn]), main = "")
plot(ergMean(sqrt(gibbsV[2, 2, -burn])),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[2]), xlab = "MCMC iteration"
)
acf(sqrt(gibbsV[2, 2, -burn]), main = "")
plot(ergMean(gibbsV[2, 1, -burn]),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[1][2]), xlab = "MCMC iteration"
)
acf(gibbsV[2, 1, -burn], main = "")
# Wmu
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.8)
par(mfrow = c(3, 2))
plot(ergMean(sqrt(gibbsWmu[1, 1, -burn])),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[mu][1]), xlab = "MCMC iteration"
)
acf(sqrt(gibbsWmu[1, 1, -burn]), main = "")
plot(ergMean(sqrt(gibbsWmu[2, 2, -burn])),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[mu][2]), xlab = "MCMC iteration"
)
acf(sqrt(gibbsWmu[2, 2, -burn]), main = "")
plot(ergMean(gibbsWmu[2, 1, -burn]),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[mu][1][2]), xlab = "MCMC iteration"
)
acf(gibbsWmu[2, 1, -burn], main = "")
# Wbeta
par(mar = c(2, 4, 1, 1) + 0.1, cex = 1)
par(mfrow = c(3, 2))
plot(ergMean(sqrt(gibbsWbeta[1, 1, -burn])),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[beta][1]), xlab = "MCMC iteration"
)
acf(sqrt(gibbsWbeta[1, 1, -burn]), main = "")
plot(ergMean(sqrt(gibbsWbeta[2, 2, -burn])),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[beta][2]), xlab = "MCMC iteration"
)
acf(sqrt(gibbsWbeta[2, 2, -burn]), main = "")
plot(ergMean(gibbsWbeta[2, 1, -burn]),
  type = "l", main = "", cex.lab = 1.5,
  ylab = expression(sigma[beta][1][2]), xlab = "MCMC iteration"
)
acf(gibbsWbeta[2, 1, -burn], main = "")

## Synthesis of MCMC output
estV <- cbind(
  mcmcMean(gibbsV[1, 1, -burn]), mcmcMean(gibbsV[2, 2, -burn]),
  mcmcMean(gibbsV[2, 1, -burn])
)
estV
estWmu <- cbind(
  mcmcMean(gibbsWmu[1, 1, -burn]),
  mcmcMean(gibbsWmu[2, 2, -burn]),
  mcmcMean(gibbsWmu[2, 1, -burn])
)
round(estWmu, 4)
estWbeta <- cbind(
  mcmcMean(gibbsWbeta[1, 1, -burn]),
  mcmcMean(gibbsWbeta[2, 2, -burn]),
  mcmcMean(gibbsWbeta[2, 1, -burn])
)
estWbeta

## Plot: data and estimated level
meanGibbsTheta <- apply(gibbsTheta[, , -burn], c(1, 2), mean)
qinf <- function(x) {
  quantile(x, prob = .025)
}
qinfGibbsTheta <- apply(gibbsTheta[, , -burn], c(1, 2), qinf)
qsup <- function(x) {
  quantile(x, prob = .975)
}
qsupGibbsTheta <- apply(gibbsTheta[, , -burn], c(1, 2), qsup)
par(mar = c(2, 4, 1, 1) + 0.1, cex = 1)
par(mfrow = c(2, 1))
require(zoo)
plot(as.zoo(y[, 1]),
  main = "", xlab = "year", cex.lab = 0.7,
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  col = "darkgray", ylab = "Investments - Denmark", type = "o", pch = 20
)
lines(as.zoo(ts(meanGibbsTheta[-1, 1], freq = frequency(y), start = start(y))),
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  type = "o", pch = 20, lty = 3, cex = .8
)
lines(as.zoo(ts(qinfGibbsTheta[-1, 1], freq = frequency(y), start = start(y))),
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  type = "l", lty = 2
)
lines(as.zoo(ts(qsupGibbsTheta[-1, 1], freq = frequency(y), start = start(y))),
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  type = "l", lty = 2
)
plot(as.zoo(y[, 2]),
  main = "", xlab = "year", cex.lab = 0.7,
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  col = "darkgray", ylab = "Investments - Spain", type = "o", pch = 20
)
lines(as.zoo(ts(meanGibbsTheta[-1, 2], freq = frequency(y), start = start(y))),
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  type = "o", pch = 20, lty = 3, cex = .8
)
lines(as.zoo(ts(qinfGibbsTheta[-1, 2], freq = frequency(y), start = start(y))),
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  type = "l", lty = 2
)
lines(as.zoo(ts(qsupGibbsTheta[-1, 2], freq = frequency(y), start = start(y))),
  oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1),
  type = "l", lty = 2
)


## sequentiel - ca. 60 sec
# tic()
# n <- 1e3
# x <- foreach(i=1:n, mu = c(1:n), sd  = c(1:n), .combine = "+") %do% {
#   rnorm(1e7, mean = mu, sd = sd)
# }
# toc()


## parallel - sec
# tic()
# registerDoParallel(4)
# n <- 1e3
# x <- foreach(i=1:n, mu = c(1:n), sd  = c(1:n), .combine = "+") %dopar% {
#   rnorm(1e7, mean = mu, sd = sd)
# }
# toc()

registerDoSEQ()

## Trash

# diesel <- read.zoo("diesel.csv", header = TRUE)
# dim(diesel)
# miss_obs <- 1
# set_count <- 0
# while (miss_obs == 1) {
#   ts <- as.xts(diesel[, sample(c(1:588), 2)])
#   set_count <- set_count + 1
#   if (sum(is.na(ts)) == 0) {
#     miss_obs <- 0
#   }
#   if (set_count >= 100) {
#     break
#   }
# }
# 
# tbl_diesel <- tibble(
#   PRICE = coredata(ts),
#   INDEX = as.POSIXct(.index(ts), origin = "1970-01-01")
# )
# 
# 
# tbl_diesel <- tbl_diesel %>%
#   mutate(DATE = lubridate::date(INDEX)) %>%
#   select(DATE = DATE, PRICE = PRICE) %>%
#   as_tsibble(index = DATE) %>%
#   mutate(
#     YEAR = year(DATE), MONTH = month(DATE, label = T),
#     WDAY = wday(DATE, label = T)
#   )
# 
# price <- tbl_diesel %>% pull(PRICE) * 0.1
