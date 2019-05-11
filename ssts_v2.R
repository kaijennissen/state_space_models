#####
library("tidyverse")
library("lubridate")
library("tsibble")
library("ggplot2")
library("stats")
library("MARSS")
library("forecast")
library("datasets")
library("KFAS")
library("dlm")
library("xts")
## gasoline - dlm ####
diesel <- read.zoo("state_space_models/diesel.csv", header = TRUE)
dim(diesel)
miss_obs <- 1
set_count <- 0
while (miss_obs==1){
  ts <- as.xts(diesel[,sample(c(1:588),1)])
  set_count <- set_count + 1
if(sum(is.na(ts)) == 0){
  miss_obs <- 0
}  
if (set_count>=100){break}
}

tbl_diesel <- tibble(PRICE = as.vector(coredata(ts)), INDEX = as.POSIXct(.index(ts), origin = "1970-01-01"))
tbl_diesel

tbl_diesel %>% mutate(DATE = lubridate::date(INDEX)) %>%
  select(DATE = DATE, PRICE = PRICE) -> tbl_diesel
as_tsibble(tbl_diesel, index = DATE) -> tsb_d
# tsb_d %>% select(DATE, PRICE) -> tsb_d 
tsb_d %>% mutate(YEAR = year(DATE), MONTH = month(DATE, label = T),
                 WDAY = wday(DATE, label = T)) -> tsb_d

price <- tsb_d %>% pull(PRICE) *0.1

mod_ex25 <- dlmModPoly(order = 1, dV = 0.25, dW = 25, m0 = 17, C0=1)
y <- c(13, 16.6, 16.3, 16.1, 17.1, 16.9, 16.8, 17.4, 17.1, 17)
(dlmFilter(y = y, mod = mod_ex25)$m)


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Constant Level
## W = 0
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
build <- function(param){
  dlmModPoly(order = 1, dV = exp(param[1]), dW = 0)
}
dlm1 <- dlmMLE(y = price, parm = rep(0,2), build = build)
fit$convergence
unlist(build(fit$par)[c("V", "W")])
mod1 <- dlmModPoly(order = 1, dV = 45.89112, dW = 0)
yFilter <- dlmFilter(y = price, mod = mod1)
ySmooth <- dlmSmooth(y = price, mod = mod1)
plot.ts(price, type = "o", col = "darkgrey")
lines(yFilter$m[-1], col = "red")
lines(ySmooth$s[-1], col = "blue")

cbind(price, yFilter$m[-1], ySmooth$s[-1]) %>% as_tibble %>% 
  rename(filter = V2, smooth = V3) %>% mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 1. Local Level
## W - 1x1
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm1 <- dlmModPoly(order = 1)
buildFun = function(x) {
  V(dlm1) = exp(x[1])
  W(dlm1) = exp(x[2])
  return(dlm1)
}

fit <- dlmMLE(y = price, parm = rep(0,2), build = buildFun)
fit$convergence
unlist(buildFun(fit$par)[c("V", "W")])

dlm1 <- buildFun(fit$par)
yFilter <- dlmFilter(y = price, mod = dlm1)
ySmooth <- dlmSmooth(y = price, mod = dlm1)

cbind(price, yFilter$m[-1], ySmooth$s[-1]) %>% as_tibble %>% 
  rename(filter = V2, smooth = V3) %>% mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()

res <- residuals(yFilter, sd = FALSE)
sapply(1:42, function(x) Box.test(res, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 2. Local Linear Trend
## W - 2x2
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm2 <- dlmModPoly(order = 2)
buildFun2 = function(x) {
  V(dlm2) = exp(x[1])
  diag(W(dlm2))[1:2] = exp(x[2:3])
  return(dlm2)
}

fit2 <- dlmMLE(y = price, parm = rep(0,3), build = buildFun)
fit2$convergence
unlist(buildFun2(fit2$par)[c("V", "W")])

dlm2 <- buildFun2(fit2$par)
yFilter2 <- dlmFilter(y = price, mod = dlm2)
ySmooth2 <- dlmSmooth(y = price, mod = dlm2)

cbind(price, ySmooth$s[-1, c(1, 2)]) -> x
colnames(x) <- c("price", "level" , "trend")
plot.ts(x, type = "o")#%>% as_tibble %>% 

res2 <- residuals(yFilter2, sd = FALSE)
plot.ts(res2)
abline(a = 0, b = 0, col = "red")
checkresiduals(res2)
sapply(1:42, function(x) Box.test(res2, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res2)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 3.1 Local Level + Daily Seasonal 
## Fixed Seasonal Component
## W - 7x7
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm31 <- dlmModPoly(order = 1) + dlmModSeas(frequency = 7)

buildFun31 = function(x) {
  V(dlm31) = exp(x[1])
  W(dlm31)[2,2] <- 0
  diag(W(dlm31))[1] = exp(x[2])
  return(dlm31)
}

fit31 <- dlmMLE(y = price, parm = rep(0,2), build = buildFun31)
fit31$convergence
unlist(buildFun31(fit31$par)[c("V", "W")])
dlm31     <-  buildFun31(fit31$par)
yFilter31 <-  dlmFilter(price, mod = dlm31)
ySmooth31 <-  dlmSmooth(price, mod = dlm31)

cbind(price, ySmooth31$s[-1, c(1, 2)]) %>% as_tibble %>%
  rename(a_price = price, b_level = V2, c_seasonal = V3) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()+
  facet_grid(series~., scales = "free_y")

res31 <- residuals(yFilter31, sd = F)
plot.ts(res31)
abline(a = 0, b = 0, col = "red")
checkresiduals(res31)
sapply(1:42, function(x) Box.test(res31, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res31)


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 3.1 Local Level + Daily Seasonal 
## Random Seasonal Component
## W - 7x7
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm32 <- dlmModPoly(order = 1) + dlmModSeas(frequency = 7)

buildFun32 = function(x) {
  V(dlm32) = exp(x[1])
  diag(W(dlm32))[1:2] = exp(x[2:3])
  return(dlm32)
}

fit32 <- dlmMLE(y = price, parm = rep(0,3), build = buildFun32)
fit32$convergence
unlist(buildFun32(fit32$par)[c("V", "W")])
dlm32 = buildFun32(fit32$par)
yFilter32 = dlmFilter(price, mod = dlm32)
ySmooth32 = dlmSmooth(price, mod = dlm32)

cbind(price, ySmooth32$s[-1, c(1, 2)]) %>% as_tibble %>%
  rename(a_price = price, b_level = V2, c_seasonal = V3) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()+
  facet_grid(series~., scales = "free_y")

res32 <- residuals(yFilter32, sd = F)
plot.ts(res32)
abline(a = 0, b = 0, col = "red")
checkresiduals(res32)
sapply(1:42, function(x) Box.test(res32, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res32)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 4.1 Local Linear Trend + Daily Seasonal 
## Fixed Seasonal Component
## W - 8x8
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm41 <- dlmModPoly(order = 2) + dlmModSeas(frequency = 7)

buildFun41 = function(x) {
  V(dlm41) = exp(x[1])
  diag(W(dlm41))[2:3] = c(exp(x[2]),0)
  return(dlm41)
}

fit41 <- dlmMLE(y = price, parm = rep(0,2), build = buildFun41)
fit41$convergence
unlist(buildFun41(fit41$par)[c("V", "W")])
dlm41 = buildFun41(fit41$par)
yFilter41 = dlmFilter(price, mod = dlm41)
ySmooth41 = dlmSmooth(price, mod = dlm41)

cbind(price, ySmooth41$s[-1, c(1, 2, 3)]) %>% as_tibble %>%
  rename(a_price = price, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()+
  facet_grid(series~., scales = "free_y")

res41 <- residuals(yFilter41, sd = F)
plot.ts(res41)
abline(a = 0, b = 0, col = "red")
checkresiduals(res41)
sapply(1:42, function(x) Box.test(res41, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res41)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 4.2 Local Linear Trend + Daily Seasonal 
## Random Seasonal Component
## W - 8x8
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm42 <- dlmModPoly(order = 2) + dlmModSeas(frequency = 7)

buildFun42 = function(x) {
  V(dlm42) = exp(x[1])
  diag(W(dlm42))[2:3] = exp(x[2:3])
  return(dlm42)
}

fit42 <- dlmMLE(y = price, parm = rep(0,3), build = buildFun42)
fit42$convergence
unlist(buildFun42(fit42$par)[c("V", "W")])
dlm42     <-  buildFun42(fit42$par)
yFilter42 <-  dlmFilter(price, mod = dlm42)
ySmooth42 <-  dlmSmooth(price, mod = dlm42)

cbind(price, ySmooth42$s[-1, c(1, 2, 3)]) %>% as_tibble %>%
  rename(a_price = price, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()+
  facet_grid(series~., scales = "free_y")

res42 <- residuals(yFilter42, sd = F)
plot.ts(res42)
abline(a = 0, b = 0, col = "red")
checkresiduals(res42)
sapply(1:42, function(x) Box.test(res42, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res42)


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 5. Local Linear Trend + Trigonometric 
## W - 8x8
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm5 <- dlmModPoly(order = 1) + dlmModTrig(s = 7, q = 3)
W(dlm5)
buildFun5 = function(x) {
  V(dlm5) = exp(x[1])
  diag(W(dlm5))[1] = exp(x[2])
  return(dlm5)
}

fit5 <- dlmMLE(y = price, parm = rep(0,2), build = buildFun5)
fit5$convergence
unlist(buildFun5(fit5$par)[c("V", "W")])
dlm5     <-  buildFun5(fit5$par)
yFilter5 <-  dlmFilter(price, mod = dlm5)
ySmooth5 <-  dlmSmooth(price, mod = dlm5)

cbind(price, ySmooth5$s[-1, c(1, 2, 3)]) %>% as_tibble %>%
  rename(a_price = price, b_level = V2, c_trend = V3, d_seasonal = V4) %>%
  mutate(time = c(1:1552)) %>% filter(time < 100) %>%
  gather("series", "value", -time) %>%
  ggplot(aes(x = time, y = value, col = series))+
  geom_line()+
  facet_grid(series~., scales = "free_y")

res5 <- residuals(yFilter5, sd = F)
plot.ts(res5)
abline(a = 0, b = 0, col = "red")
checkresiduals(res5)
sapply(1:42, function(x) Box.test(res5, lag = x, type = "Ljung-Box")$p.value)
shapiro.test(res5)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 6. Multivariate SSM
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#



##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 7. Bayesian Variance Estimation - Conjugate Prio
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 8. Bayesian Variance Estimation - Sequential MCMC
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## prior
## likelihood
## posterior
## new prior

online_mean <- function(theta = 2, sig = 2, n = 100, m_0 = 2, C_0 = 10){
  # theta <- 2
  # sig <- 2
  # n <- 100
  y <- theta + rnorm(n = n, mean = 0, sd = sig)
  sig_sq <- sig^2
  
  # m_0 <- 100
  # C_0 <- 100
  
  filter_y <- vector(length = n, mode = "numeric")
  filter_ci <- vector(length = n, mode = "numeric")
  for (i in 1:n){
    
    if (i == 1) {
      m_prio <- m_0
      C_prio <- C_0
    } else{
      m_prio <- m_post
      C_prio <- C_post
    }
    y_n <- y[i]
    m_post <- m_prio+(C_prio)/(C_prio+sig_sq)*(y_n-m_prio)
    C_post <- (sig_sq*C_prio)/(sig_sq+C_prio)
    filter_y[i] <- m_post
    filter_ci[i] <- C_post
  }
  plot.ts(y, type = "b", col = "blue")
  lines(filter_y, type = "o", col = "red")
  lines(filter_y+1.96*filter_ci, type = "l", col = "red")
  lines(filter_y-1.96*filter_ci, type = "l", col = "red")
}

online_mean(theta = 10, sig = 2, n = 100, m_0 = -20, C_0 = 2000)

online_mean_var <- function(theta = 2, sig = 4, n = 100,
                            m_0 = 2, n_0 = 10, a_0 = 2, b_0 = 2){
  
  y <- theta + rnorm(n = n, mean = 0, sd = sig)
  sig_sq <- sig^2
  
  theta_hat <- vector(length = n, mode = "numeric")
  phi_hat <- vector(length = n, mode = "numeric")
  for (i in 1:n){
    
    if (i == 1) {
      m_prio <- m_0
      n_prio <- n_0
      a_prio <- a_0
      b_prio <- b_0
    } else{
      m_prio <- m_post
      n_prio <- n_post
      a_prio <- a_post
      b_prio <- b_post
    }
    y_n <- y[i]
    m_post <- m_prio+(y_n-m_prio)/(n_prio+1)
    n_post <- n_prio+1
    a_post <- a_prio+0.5
    b_post    <- b_prio+(n_prio*(y_n-m_prio)^2)/(2*n_prio+2)
    theta_hat[i] <- m_post
    phi_hat[i] <- (b_post)/(a_post)
    
  }
  plot.ts(y, type = "b", col = "blue")
  lines(theta_hat, type = "o", col = "red")
  #lines(mean_theta+1.96*var_theta, type = "l", col = "red")
  #lines(mean_theta-1.96*var_theta, type = "l", col = "red")
  return(list(y = y, theta = theta_hat, phi = phi_hat))
}


xx <- online_mean_var(theta = 2, sig = 5, n = 1000,
                      m_0 = 2, n_0 = 10, a_0 = 100, b_0 = 1)
plot.ts(xx$phi, type = "o", col = "red")
plot.ts(xx$theta, col = "red")
lines(cumsum(xx$y)/seq(1,1000), col = "blue")
#plot.ts(xx$y, type = "b", col = "blue")
plot.ts(xx$y, type = "b", col = "blue")

