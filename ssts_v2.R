## load libraries and data ====================================================
library("tidyverse")
# library("tidybayes")
library("lubridate")
#library("anytime")
library("tsibble")
library("ggplot2")
library("stats")
# library("MARSS")
library("forecast")
# library("datasets")
# library("KFAS")
library("dlm")
library("xts")


diesel <- read.zoo("./diesel.csv", header = TRUE)
dim(diesel)
miss_obs <- 1
set_count <- 0
while (miss_obs == 1) {
  ts <- as.xts(diesel[, sample(c(1:588), 10)])
  set_count <- set_count + 1
  if (sum(is.na(ts)) == 0) {
    miss_obs <- 0
  }
  if (set_count >= 1000) {
    print(set_count)
    break
  }
}

# days <- ymd(zoo::index(ts))
tbl_diesel <-
  as_tibble(coredata(ts)) %>%
  bind_cols(DATE = ymd(zoo::index(ts))) %>%
  select(DATE, everything()) %>%
  rename_if(is.numeric, funs(str_c("STATION_", 1:10)))

tsb_d <- as_tsibble(tbl_diesel, index = DATE)

tsb_d <-
  tsb_d %>%
  gather(-DATE, key = "STATION", value = "PRICE", factor_key = T) %>%
  mutate(
    YEAR = year(DATE), MONTH = month(DATE, label = T),
    WDAY = wday(DATE, label = T)
  ) %>%
  separate(STATION, c("DROP", "STATION")) %>%
  select(-DROP) %>%
  mutate_at(vars("STATION"), as.factor)

tsb_d %>%
  ggplot(aes(x = DATE, y = PRICE, col = STATION)) +
  geom_line() +
  facet_wrap(~YEAR, scales = "free_x")

price <- tbl_diesel %>%
  transmute(ST1 = STATION_1 * 0.1, ST2 = STATION_2 * 0.1) %>%
  as.matrix()
price1 <- price[, 1]
price1 <- price[, 2]
# mod_ex25 <- dlmModPoly(order = 1, dV = 0.25, dW = 25, m0 = 17, C0=1)
# y <- c(13, 16.6, 16.3, 16.1, 17.1, 16.9, 16.8, 17.4, 17.1, 17)
# (dlmFilter(y = y, mod = mod_ex25)$m)


## Constant Level =============================================================
## W = 0

build <- function(param) {
  dlmModPoly(order = 1, dV = exp(param[1]), dW = 0)
}
dlm1 <- dlmMLE(y = price, parm = rep(0, 2), build = build)
fit$convergence
unlist(build(fit$par)[c("V", "W")])
mod1 <- dlmModPoly(order = 1, dV = 45.89112, dW = 0)
yFilter <- dlmFilter(y = price, mod = mod1)
ySmooth <- dlmSmooth(y = price, mod = mod1)
plot.ts(price, type = "o", col = "darkgrey")
lines(yFilter$m[-1], col = "red")
lines(ySmooth$s[-1], col = "blue")

cbind(price, yFilter$m[-1], ySmooth$s[-1]) %>%
  as_tibble() %>%
  rename(filter = V2, smooth = V3) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line()


## 1. Local Level =============================================================
## W - 1x1

dlm1 <- dlmModPoly(order = 1)
buildFun <- function(x) {
  V(dlm1) <- exp(x[1])
  W(dlm1) <- exp(x[2])
  return(dlm1)
}

fit <- dlmMLE(y = price, parm = rep(0, 2), build = buildFun)
fit$convergence
unlist(buildFun(fit$par)[c("V", "W")])

dlm1 <- buildFun(fit$par)
yFilter <- dlmFilter(y = price, mod = dlm1)
ySmooth <- dlmSmooth(y = price, mod = dlm1)

cbind(price, yFilter$m[-1], ySmooth$s[-1]) %>%
  as_tibble() %>%
  rename(filter = V2, smooth = V3) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line()

res <- residuals(yFilter, sd = FALSE)
sapply(1:42, function(x) Box.test(res, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res)


## 2. Local Linear Trend ======================================================
## W - 2x2

dlm2 <- dlmModPoly(order = 2)
buildFun2 <- function(x) {
  V(dlm2) <- exp(x[1])
  diag(W(dlm2))[1:2] <- exp(x[2:3])
  return(dlm2)
}

fit2 <- dlmMLE(y = price, parm = rep(0, 3), build = buildFun)
fit2$convergence
unlist(buildFun2(fit2$par)[c("V", "W")])

dlm2 <- buildFun2(fit2$par)
yFilter2 <- dlmFilter(y = price, mod = dlm2)
ySmooth2 <- dlmSmooth(y = price, mod = dlm2)

cbind(price, ySmooth$s[-1, c(1, 2)]) -> x
colnames(x) <- c("price", "level", "trend")
plot.ts(x, type = "o") # %>% as_tibble %>%

res2 <- residuals(yFilter2, sd = FALSE)
plot.ts(res2)
abline(a = 0, b = 0, col = "red")
checkresiduals(res2)
sapply(1:42, function(x) Box.test(res2, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res2)

## 3.1 Local Level + Daily Seasonal ===========================================
## Fixed Seasonal Component
## W - 7x7

dlm31 <- dlmModPoly(order = 1) + dlmModSeas(frequency = 7)

buildFun31 <- function(x) {
  V(dlm31) <- exp(x[1])
  W(dlm31)[2, 2] <- 0
  diag(W(dlm31))[1] <- exp(x[2])
  return(dlm31)
}

fit31 <- dlmMLE(y = price, parm = rep(0, 2), build = buildFun31)
fit31$convergence
unlist(buildFun31(fit31$par)[c("V", "W")])
dlm31 <- buildFun31(fit31$par)
yFilter31 <- dlmFilter(price, mod = dlm31)
ySmooth31 <- dlmSmooth(price, mod = dlm31)

cbind(price, ySmooth31$s[-1, c(1, 2)]) %>%
  as_tibble() %>%
  rename(a_price = price, b_level = V2, c_seasonal = V3) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line() +
  facet_grid(series ~ ., scales = "free_y")

res31 <- residuals(yFilter31, sd = F)
plot.ts(res31)
abline(a = 0, b = 0, col = "red")
checkresiduals(res31)
sapply(1:42, function(x) Box.test(res31, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res31)



## 3.1 Local Level + Daily Seasonal ===========================================
## Random Seasonal Component
## W - 7x7

dlm32 <- dlmModPoly(order = 1) + dlmModSeas(frequency = 7)

buildFun32 <- function(x) {
  V(dlm32) <- exp(x[1])
  diag(W(dlm32))[1:2] <- exp(x[2:3])
  return(dlm32)
}

fit32 <- dlmMLE(y = price, parm = rep(0, 3), build = buildFun32)
fit32$convergence
unlist(buildFun32(fit32$par)[c("V", "W")])
dlm32 <- buildFun32(fit32$par)
yFilter32 <- dlmFilter(price, mod = dlm32)
ySmooth32 <- dlmSmooth(price, mod = dlm32)

cbind(price, ySmooth32$s[-1, c(1, 2)]) %>%
  as_tibble() %>%
  rename(a_price = price, b_level = V2, c_seasonal = V3) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line() +
  facet_grid(series ~ ., scales = "free_y")

res32 <- residuals(yFilter32, sd = F)
plot.ts(res32)
abline(a = 0, b = 0, col = "red")
checkresiduals(res32)
sapply(1:42, function(x) Box.test(res32, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res32)

## 4.1 Local Linear Trend + Daily Seasonal ====================================
## Fixed Seasonal Component
## W - 8x8

dlm41 <- dlmModPoly(order = 2) + dlmModSeas(frequency = 7)

buildFun41 <- function(x) {
  V(dlm41) <- exp(x[1])
  diag(W(dlm41))[2:3] <- c(exp(x[2]), 0)
  return(dlm41)
}

fit41 <- dlmMLE(y = price, parm = rep(0, 2), build = buildFun41)
fit41$convergence
unlist(buildFun41(fit41$par)[c("V", "W")])
dlm41 <- buildFun41(fit41$par)
yFilter41 <- dlmFilter(price, mod = dlm41)
ySmooth41 <- dlmSmooth(price, mod = dlm41)

cbind(price, ySmooth41$s[-1, c(1, 2, 3)]) %>%
  as_tibble() %>%
  rename(a_price = price, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line() +
  facet_grid(series ~ ., scales = "free_y")

res41 <- residuals(yFilter41, sd = F)
plot.ts(res41)
abline(a = 0, b = 0, col = "red")
checkresiduals(res41)
sapply(1:42, function(x) Box.test(res41, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res41)


## 4.2 Local Linear Trend + Daily Seasonal ====================================
## Random Seasonal Component
## W - 8x8

dlm42 <- dlmModPoly(order = 2) + dlmModSeas(frequency = 7)

buildFun42 <- function(x) {
  V(dlm42) <- exp(x[1])
  diag(W(dlm42))[2:3] <- exp(x[2:3])
  return(dlm42)
}

fit42 <- dlmMLE(y = price, parm = rep(0, 3), build = buildFun42)
fit42$convergence
unlist(buildFun42(fit42$par)[c("V", "W")])
dlm42 <- buildFun42(fit42$par)
yFilter42 <- dlmFilter(price1, mod = dlm42)
ySmooth42 <- dlmSmooth(price1, mod = dlm42)

cbind(price1, ySmooth42$s[-1, c(1, 2, 3)]) %>%
  as_tibble() %>%
  rename(a_price = price1, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line() +
  facet_grid(series ~ ., scales = "free_y")

res42 <- residuals(yFilter42, sd = F)
plot.ts(res42)
abline(a = 0, b = 0, col = "red")
checkresiduals(res42)
sapply(1:42, function(x) Box.test(res42, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res42)


## 5. Local Linear Trend + Trigonometric ======================================
## W - 8x8

dlm5 <- dlmModPoly(order = 1) + dlmModTrig(s = 7, q = 3)
W(dlm5)
buildFun5 <- function(x) {
  V(dlm5) <- exp(x[1])
  diag(W(dlm5))[1] <- exp(x[2])
  return(dlm5)
}

fit5 <- dlmMLE(y = price, parm = rep(0, 2), build = buildFun5)
fit5$convergence
unlist(buildFun5(fit5$par)[c("V", "W")])
dlm5 <- buildFun5(fit5$par)
yFilter5 <- dlmFilter(price, mod = dlm5)
ySmooth5 <- dlmSmooth(price, mod = dlm5)

cbind(price, ySmooth5$s[-1, c(1, 2, 3)]) %>%
  as_tibble() %>%
  rename(a_price = price, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>%
  filter(time < 100) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line() +
  facet_grid(series ~ ., scales = "free_y")

res5 <- residuals(yFilter5, sd = F)
plot.ts(res5)
abline(a = 0, b = 0, col = "red")
checkresiduals(res5)
sapply(1:42, function(x) Box.test(res5, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res5)

## 6.1 Multivariate SSM ========================================================
## multivariate local level -- seemingly unrelated time series?
## W - 8x8

buildFun61 <- function(x) {
  dlm61 <- dlmModPoly(order = 2) # + dlmModSeas(frequency = 7)

  ## MV specification
  dlm61$GG <- dlm61$GG %x% diag(2)
  dlm61$m0 <- rep(0, 4)
  dlm61$C0 <- diag(4) * 1e7

  ## V
  Vsd <- exp(x[1:2]) # only positive values
  Vcorr <- tanh(x[3]) # values between -1 and 1
  V <- Vsd %o% Vsd
  V[1, 2] <- V[2, 1] <- V[2, 1] * Vcorr
  dlm61$V <- V

  ## W1
  W1_sd <- exp(x[4:5])
  W1_corr <- tanh(x[6])
  W1 <- W1_sd %o% W1_sd
  W1[1, 2] <- W1[2, 1] <- W1[1, 2] * W1_corr

  ## W2
  W2_sd <- exp(x[7:8])
  W2_corr <- tanh(x[9])
  W2 <- W2_sd %o% W2_sd
  W2[1, 2] <- W2[2, 1] <- W2[1, 2] * W2_corr

  dlm61$W <- bdiag(W1, W2)

  return(dlm61)
}

fit61 <- dlmMLE(
  y = price, parm = rnorm(9, 0, 1),
  build = buildFun61
)
dlm61 <- buildFun61(fit61$par)
# unlist(buildFun61(fit61$par)[c("V", "W")])
yFilter61 <- dlmFilter(price, mod = dlm61)
ySmooth61 <- dlmSmooth(price, mod = dlm61)
sdev <- residuals(yFilter61)$sd
yFilter61$m

mu <- yFilter
plot.ts(mu, plot.type = "s", col = c("red", "blue"))
plot.ts(beta)
plot.ts(beta[-1, ])

## 6.2 Multivariate SSM ========================================================
## Local Linear Trend + Fixed Seasonality
## W - 8x8

buildFun62 <- function(x) {
  dlm <- dlmModPoly(order = 2) # + dlmModSeas(frequency = 7)
  
  ## MV specification
  dlm$GG <- dlm$GG %x% diag(2)
  dlm$m0 <- rep(0, 4)
  dlm$C0 <- diag(4) * 1e7
  
  ## V
  Vsd <- exp(x[1:2]) # only positive values
  Vcorr <- tanh(x[3]) # values between -1 and 1
  V <- Vsd %o% Vsd
  V[1, 2] <- V[2, 1] <- V[2, 1] * Vcorr
  dlm62$V <- V
  
  ## W1
  W1_sd <- rep(0, 2) #exp(x[4:5])
  W1_corr <- 0#tanh(x[6])
  W1 <- W1_sd %o% W1_sd
  W1[1, 2] <- W1[2, 1] <- W1[1, 2] * W1_corr
  
  ## W2
  W2_sd <- exp(x[4:5])
  W2_corr <- tanh(x[6])
  W2 <- W2_sd %o% W2_sd
  W2[1, 2] <- W2[2, 1] <- W2[1, 2] * W2_corr
  
  dlm$W <- bdiag(W1, W2)
  
  return(dlm)
}

fit62 <- dlmMLE(
  y = price, parm = rnorm(6, 0, 1),
  build = buildFun62
)
dlm62 <- buildFun62(fit62$par)
# unlist(buildFun61(fit61$par)[c("V", "W")])
yFilter62 <- dlmFilter(price1, mod = dlm62)
ySmooth62 <- dlmSmooth(price1, mod = dlm62)

mu <- yFilter62$f
plot.ts(mu[-1,], plot.type = "s", col = c("red", "blue"))
plot.ts(beta)
plot.ts(beta[-1, ])

cbind(price, ySmooth5$s[-1, c(1, 2, 3)]) %>%
  as_tibble() %>%
  rename(a_price = price, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>%
  filter(time < 100) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series)) +
  geom_line() +
  facet_grid(series ~ ., scales = "free_y")

res62 <- residuals(yFilter62, sd = F)
plot.ts(res62)
abline(a = 0, b = 0, col = "red")
checkresiduals(res62)
sapply(1:42, function(x) Box.test(res62, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res62)



## 8. Univariate Bayesian Variance Estimation - Conjugate Prio ===========================

data(NelPlo)
### multivariate local level -- seemingly unrelated time series
buildSu <- function(x) {
  Vsd <- exp(x[1:2])
  Vcorr <- tanh(x[3])
  V <- Vsd %o% Vsd
  V[1, 2] <- V[2, 1] <- V[1, 2] * Vcorr
  Wsd <- exp(x[4:5])
  Wcorr <- tanh(x[6])
  W <- Wsd %o% Wsd
  W[1, 2] <- W[2, 1] <- W[1, 2] * Wcorr
  return(list(
    m0 = rep(0, 2),
    C0 = 1e7 * diag(2),
    FF = diag(2),
    GG = diag(2),
    V = V,
    W = W
  ))
}

suMLE <- dlmMLE(NelPlo, rep(0, 6), buildSu)
suMLE
buildSu(suMLE$par)[c("V", "W")]
StructTS(NelPlo[, 1], type = "level") ## compare with W[1,1] and V[1,1]
StructTS(NelPlo[, 2], type = "level") ## compare with W[2,2] and V[2,2]


## 8. Bayesian Variance Estimation - Sequential MCMC ==========================

## prior
## likelihood
## posterior
## new prior

## 9.1 Fama & French: Linear Regression =======================================
## 


## 9.2 Fama & French: Rolling Regression =======================================
## 


## 9.3 Fama & French: Dynamic Linear Model ====================================
## 

