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
diesel <- read.zoo("/Users/kj/Documents/02_Studium/04_Data_Masterthesis/avg_v6/xts_all_15_diesel.csv", header = TRUE)
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
## 5 Local Linear Trend + Trigonometric 
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




## bayesian online filtering ####
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


## MARSS - Fishery and Environement ######

mod.list <- list(B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0)

# simulate
q <- 0.1; r <- 0.1; n <- 100
y <- cumsum(rnorm(n,0,sqrt(q)))+rnorm(n,0,sqrt(r))

# fit
fit <- MARSS(y, model=mod.list)

mod.list$Q <- matrix(0.1)
fit <- MARSS(y,model=mod.list)

dat <- as.vector(Nile)

tibble(Year = 1871:1970, Flow = dat) %>% 
ggplot(aes(x=Year, y = Flow))+
  geom_line()+
  geom_smooth()

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 1. Local Level 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = B*x_t + u_t + w_t ; w_t ~ N(0,Q)
## x_t = Z*x_t-1 + a + v_t ;v_t ~ N(0,R)
##
## y_t = x_t + w_t ; w_t ~ N(0,q)  
## x_t = x_t-1 + v_t ;  v_t ~ N(0,0)
## x_0 = mu
## B=1; U=0; Q=0; Z=1; A=0, R="r"; x0=mu; 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

mod.nile.0 <- list( 
  B=matrix(1), U=matrix(0), Q=matrix(0),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )
fit.0 <- MARSS(dat, model=mod.nile.0)

c(coef(fit.0, type="vector"), LL=fit.0$logLik, AICc=fit.0$AICc)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 2. Local Linear Trend 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = B*x_t + a + v_t ; v_t ~ N(0,R)
## x_t = Z*x_t-1 + u_t + w_t ;w_t ~ N(0,Q)
##
## y_t = x_t + w_t ; w_t ~ N(0,r)  
## x_t = x_t-1 + u + v_t ;  v_t ~ N(0,0)
## x_0 = mu
## B=1; U=u; Q=0; Z=1; A=0, R="r"; x0=mu; 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
mod.nile.1 <- list( 
  B=matrix(1), U=matrix("u"), Q=matrix(0),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )

fit.1 <- MARSS(dat, model=mod.nile.1)

c(coef(fit.1, type="vector"), LL=fit.1$logLik, AICc=fit.1$AICc)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 3. Stochastic Level 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = B*x_t + a + v_t ; v_t ~ N(0,R)
## x_t = Z*x_t-1 + u_t + w_t ;w_t ~ N(0,Q)
## 
## y_t = x_t + w_t ; w_t ~ N(0,r)
## x_t = x_t-1 + v_t ;  v_t ~ N(0,q)
## x_0 = mu
## B=1; U=0; Q=q; Z=1; A=0, R="r"; x0=mu;
###-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
mod.nile.2 = list(
  B=matrix(1), U=matrix(0), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )

fit.2 <- MARSS(dat, model=mod.nile.2)

c(coef(fit.2, type="vector"), LL=fit.2$logLik, AICc=fit.2$AICc)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 4. Stochastic Level with Drift 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = B*x_t + a + v_t ; v_t ~ N(0,R)
## x_t = Z*x_t-1 + u_t + w_t ;w_t ~ N(0,Q)
##
## y_t = x_t + w_t ; w_t ~ N(0,r)  
## x_t = x_t-1 + u + v_t ;  v_t ~ N(0,q)
## x_0 = mu
## B=1; U=u; Q=q; Z=1; A=0, R=r; x0=mu; 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

mod.nile.3 = list(
  B=matrix(1), U=matrix("u"), Q=matrix("q"),
  Z=matrix(1), A=matrix(0), R=matrix("r"),
  x0=matrix("mu"), tinitx=0 )

fit.3 <- MARSS(dat, model=mod.nile.3)

c(coef(fit.3, type="vector"), LL=fit.3$logLik, AICc=fit.3$AICc)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 4. Stochastic Level with Drift in JAGS
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = B*x_t + a + v_t ; v_t ~ N(0,R)
## x_t = Z*x_t-1 + u_t + w_t ;w_t ~ N(0,Q)
##
## y_t = x_t + w_t ; w_t ~ N(0,r)  
## x_t = x_t-1 + u + v_t ;  v_t ~ N(0,q)
## x_0 = mu
## B=1; U=u; Q=q; Z=1; A=0, R=r; x0=mu; 
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(R2jags)
library(rjags)
library(coda)

y <- as.vector(Nile)

model.loc <- "ss_model.txt"
jagsscript <- cat(
model {  
# priors on parameters
mu ~ dnorm(Y1, 1/(Y1*100)); # normal mean = 0, sd = 1/sqrt(0.01)
tau.q ~ dgamma(0.001,0.001); # This is inverse gamma
sd.q <- 1/sqrt(tau.q); # sd is treated as derived parameter
tau.r ~ dgamma(0.001,0.001); # This is inverse gamma
sd.r <- 1/sqrt(tau.r); # sd is treated as derived parameter
u ~ dnorm(0, 0.01);
                  
# Because init X is specified at t=0
X0 <- mu
X[1] ~ dnorm(X0+u,tau.q);
Y[1] ~ dnorm(X[1], tau.r);
                  
for(i in 2:TT) {
predX[i] <- X[i-1]+u; 
X[i] ~ dnorm(predX[i],tau.q); # Process variation
Y[i] ~ dnorm(X[i], tau.r); # Observation variation
}
}                  
,file=model.loc)

jags.data <- list("Y"=y, "TT"=length(y), Y1=y[1])
jags.params <- c("sd.q", "sd.r", "X", "mu", "u")

mod_ss <- jags(jags.data, parameters.to.save=jags.params, 
               model.file=model.loc, n.chains = 3, 
               n.burnin=5000, n.thin=1, n.iter=10000, DIC=TRUE)
mod_ss
summary(mod_ss)

attach.jags(mod_ss)
par(mfrow=c(2,2))
hist(mu)
abline(v=coef(kem.3)$x0, col="red")
hist(u)
abline(v=coef(kem.3)$U, col="red")
hist(log(sd.q^2))
abline(v=log(coef(kem.3)$Q), col="red")
hist(log(sd.r^2))
abline(v=log(coef(kem.3)$R), col="red")

plotModelOutput <- function(jagsmodel, Y) {
  attach.jags(jagsmodel)
  x <- seq(1,length(Y))
  XPred <- cbind(apply(X,2,quantile,0.025), apply(X,2,mean), apply(X,2,quantile,0.975))
  ylims <- c(min(c(Y,XPred), na.rm=TRUE), max(c(Y,XPred), na.rm=TRUE))
  plot(Y, col="white",ylim=ylims, xlab="",ylab="State predictions")
  polygon(c(x,rev(x)), c(XPred[,1], rev(XPred[,3])), col="grey70",border=NA)
  lines(XPred[,2])
  points(Y)
}
plotModelOutput(mod_ss, y)

lines(kem.3$states[1,], col="red")
lines(1.96*kem.3$states.se[1,]+kem.3$states[1,], col="red", lty=2)
lines(-1.96*kem.3$states.se[1,]+kem.3$states[1,], col="red", lty=2)
title("State estimate and data from\nJAGS (black) versus MARSS (red)")



##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Residual Diagnostics
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## residuals from measurement - v_t - model residuals
## residuals from state- w_t - state residuals
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

## Local Level
par(mfrow=c(2,1))
resids <- residuals(fit.0)
plot(resids$model.residuals[1,], 
     ylab="model residual", xlab="", main="Local Level")
abline(h=0)
plot(resids$state.residuals[1,], 
     ylab="state residual", xlab="", main="Local Lecel")
abline(h=0)

## Local Level w Drift
par(mfrow=c(2,1))
resids <- residuals(fit.1)
plot(resids$model.residuals[1,], 
     ylab="model residual", xlab="", main="Local Level w Drift")
abline(h=0)
plot(resids$state.residuals[1,], 
     ylab="state residual", xlab="", main="Local Level w Drift")
abline(h=0)

## Stochastic Level 
par(mfrow=c(2,1))
resids <- residuals(fit.2)
plot(resids$model.residuals[1,], 
     ylab="model residual", xlab="", main="Stochastic Level")
abline(h=0)
plot(resids$state.residuals[1,], 
     ylab="state residual", xlab="", main="Stochastic Level")
abline(h=0)

## Stochastic Level w Drift
par(mfrow=c(2,1))
resids <- residuals(fit.3)
plot(resids$model.residuals[1,], 
     ylab="model residual", xlab="", main="Stochastic Level w Drift")
abline(h=0)
plot(resids$state.residuals[1,], 
     ylab="state residual", xlab="", main="Stochastic Level w Drift")
abline(h=0)

par(mfrow=c(2,2))
resids <- residuals(fit.0)
acf(resids$model.residuals[1,], main="flat level v(t)")
resids <- residuals(fit.1)
acf(resids$model.residuals[1,], main="linear trend v(t)")
resids <- residuals(fit.2)
acf(resids$model.residuals[1,], main="stoc level v(t)")
acf(resids$state.residuals[1,], main="stoc level w(t)", na.action=na.pass)

## dlm - Intro State Space Time Series ####
library(dlm)
data.1            <- log(read.table("/Users/kj/Desktop/SSTS/UKdriversKSI.txt",skip=1))
colnames(data.1)  <- "logKSI"
data.1            <- ts(data.1, start = c(1969), frequency = 12)
data.2            <- log(read.table("/Users/kj/Desktop/SSTS/NorwayFinland.txt",skip=1))
data.2            <- data.2[,2,drop=F]
colnames(data.2)  <- "logNorFatalities"
data.2            <- ts(data.2 , start = c(1970,1))

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 1. Deterministic Level
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
(fit <- lm(data.1[,1]~1 ))
res <- residuals(fit)
(coefs <- round(as.numeric(coef(fit)),8))
(error.var<- round(summary(fit)$sigma^2,8))

par(mfrow=c(1,1))
plot.ts(c(data.1),  col = "darkgrey", xlab="",ylab = "log KSI",pch=3,cex=0.5,
          cex.lab=0.8,cex.axis=0.7,xlim=c(0,200))
abline(h=coefs[1] , col = "blue", lwd  = 2, lty=2)
legend("topright",leg = c("log UK drivers KSI",
                            "  deterministic level"),cex = 0.6,
         lty = c(1, 2),col = c("darkgrey","blue"),
         pch=c(3,NA), bty = "y", horiz = T)

par(mfrow=c(1,1))
plot(ts(residuals(fit)),ylab="",xlab="",xlim = c(0,200), col = "darkgrey",lty=2)
abline(h=0)

##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 2. Stochastic Level
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = F_t*theta_t + v_t ; v_t ~ N(0,V_t)
## theta_t = G_t*theta_t-1 + w_t ; w_t ~ N(0,W_t)
## thata_0 ~ N(m_0,C_0)
##
## y_t = theta_t + v_t ; v_t ~ N(0,V)  
## theta_t = theta_t-1 + w_t ;  w_t ~ N(0,W)
## x_0 = mu
## 
## FF = 1; V = ; GG = 1; W; m0=0; C0=100; 
## Time Varying Models
## JFF, JV, JGG, JW, X
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

fn  <- function(params){
  dlmModPoly(order = 1, dV = exp(params[1]) , dW = exp(params[2]))
}

y                <- c(data.1)
fit              <- dlmMLE(y, rep(0,2), fn)
mod              <- fn(fit$par)
(obs.error.var   <- V(mod))
(state.error.var <- W(mod))
(mu.1            <- dropFirst(dlmFilter(y,mod)$m)[1])
res              <- residuals(dlmFilter(y,mod),sd=F)
filtered         <- dlmFilter(y,mod)
smoothed         <- dlmSmooth(filtered)
mu               <- dropFirst(smoothed$s)
mu.1             <- mu[1]


par(mfrow=c(1,1))
plot.ts(y,  col = "darkgrey", xlab="",ylab = "log KSI",pch=3,cex=0.5,
          cex.lab=0.8,cex.axis=0.7,xlim=c(0,200))
lines(mu , col = "blue", lwd  = 2, lty=2)
legend("topright",leg = c("log UK drivers KSI","  stochastic level"),
         cex = 0.6,lty = c(1, 2), col = c("darkgrey","blue"),
         pch=c(3,NA), bty = "y", horiz = T)


par(mfrow=c(1,1))
plot.ts(res,ylab="",xlab="", col = "darkgrey", main = "",cex.main = 0.8)
abline(h=0)

shapiro.test(res)
Box.test(res, lag = 15, type = "Ljung")
sapply(1:20,function(l){round(Box.test(res, lag=l, type = "Ljung-Box")$p.value,4)})



##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 3. Stochastic Level and Stochastic Slope
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = F_t*theta_t + v_t ; v_t ~ N(0,V_t)
## theta_t = G_t*theta_t-1 + w_t ; w_t ~ N(0,W_t)
## thata_0 ~ N(m_0,C_0)
##
## y_t = theta_t + v_t ; v_t ~ N(0,V)  
## theta_t = theta_t-1 + w_t ;  w_t ~ N(0,W)
## x_0 = mu
## 
## FF = 1; V = ; GG = 1; W; m0=0; C0=100; 
## Time Varying Models
## JFF, JV, JGG, JW, X
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

level0 <- y[1]
slope0 <- mean(diff(y))
fn <- function(params){
  dlmModPoly(dV = exp(params[1]), dW = exp(params[2:3]),
             m0 = c(level0, slope0),C0 = 2*diag(2))
}
fit                <- dlmMLE(y, rep(0,3), fn)
mod                <- fn(fit$par)
(obs.error.var     <- V(mod))
(state.error.var.1 <- W(mod)[1,1])
(state.error.var.2 <- W(mod)[2,2])
filtered           <- dlmFilter(y, mod)
smoothed           <- dlmSmooth(filtered)
sm                 <- dropFirst(smoothed$s)
mu.1               <- sm[1,1]
nu.1               <- sm[1,2]
mu                 <- c(sm[,1])
nu                 <- c(sm[,2])
res                <- c(residuals(filtered, sd = F))

par(mfrow=c(1,1))
temp <- window(cbind(data.1,mu),start = 1969, end = 1984)
plot(temp , plot.type="single" , col =c("darkgrey","blue"),lty=c(1,2), xlab="",ylab = "log KSI")
legend("topright",leg = c("log UK drivers KSI","  stochastic level and Slope"),
         cex = 0.7, lty = c(1, 2),col = c("darkgrey","blue"),
         pch=c(3,NA),bty = "y", horiz = T)

par(mfrow=c(1,1))
plot(ts(nu*10^4,start =1969,end=1984,frequency =12),xlab="",
col = "darkgrey",lty=2,ylab=expression(10^(-4)))

par(mfrow=c(1,1))
plot(ts(res,c(1969),frequency=12),ylab="",xlab="", col = "darkgrey",
       lty=2, main = 'Irregular component', cex.main = 0.7) > abline(h=0)

shapiro.test(res)
Box.test(res, lag = 15, type = "Ljung")
sapply(1:15,function(l){round(Box.test(res, lag=l, type = "Ljung-Box")$p.value,4)})


##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 3. Stochastic Level, Stochastic Slope and Regression Coefficient
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

intervention        <- rep(1,dim(data.1)[1])
intervention[1:169] <- 0
intervention        <- ts(intervention, start = c(1969),frequency=12)
prices <- read.table("/Users/kj/Desktop/SSTS/logUKpetrolprice.txt",skip=1)
prices <- ts(prices, start = c(1969),frequency=12)

fn <- function(params){
  mod <- dlmModPoly(order = 1, dV =exp(params[1]), dW = exp(params[2])) +
    dlmModSeas(frequency = 12, dV =exp(params[1]) , dW= rep(0,11))
  return(mod)
}
fit      <- dlmMLE(data.1, rep(0,2),fn)
mod      <- fn(fit$par)
filtered <- dlmFilter(data.1,mod)
smoothed <- dlmSmooth(filtered)
cov      <- dlmSvd2var(smoothed$U.S, smoothed$D.S)
lev.var  <- sapply(cov, function(x){x[1,1]})
mu       <- ((smoothed$s)[,1])[-1]
nu       <- ((smoothed$s)[,2])[-1]

res    <- residuals(smoothed,sd=F)
lev.ts <- ts(lev.var[-1],start = 1969, frequency  =12)
wid    <- qnorm(0.05, lower = FALSE) *sqrt(lev.ts)
temp   <- cbind(mu, mu + wid %o% c(-1, 1))
temp   <- ts(temp,start = 1969,frequency =12)

par(mfrow=c(1,1))
plot(lev.ts,xlab="",ylab = "level estimation variance")

par(mfrow=c(1,1))
plot(temp, plot.type = "s", type = "l",lty = c(1, 5, 5),
       ylab = "Level", xlab = "", ylim = range(data.1),col=c("blue","red","red"),lwd=2)
lines(data.1, type = "o", col = "darkgrey")
legend("topright",
         leg = c("log UK drivers KSI","  stochastic level +/- 1.64SE"),
         cex = 0.7, lty = c(1, 5),col = c("darkgrey","red"),
         bty = "y", horiz = F)

nu.var    <- sapply(cov, function(x){x[2,2]})
nu.var.ts <- ts(nu.var[-1], start = c(1969,1), frequency = 12)
wid       <- qnorm(0.05, lower = F) * sqrt(nu.var.ts)
temp      <- cbind(nu, nu + wid %o% c(-1,1))
temp      <- ts(temp,start = 1969,frequency =12)
par(mfrow=c(1,1))
plot(temp, plot.type = "s", type = "l",lty = c(1, 5, 5),
       ylab = "Level", xlab = "",col=c("blue","red","red"),lwd=1)
legend("topright",
         leg =   "deterministic level +/- 1.64SE",
         cex = 0.7, lty = c(5),col = c("red"),
         bty = "y", horiz = F)


## bsts ####
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## y_t = mu_t + gamma_t + beta'_t* x_t + e_t
## mu_t = mu_t-1 + delta_t-1 + u_t
## delta_t = delta_t-1 + v_t
## gamma_t = -gamma_t-1 -...- gamma_t-S+1 + w_t
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
(model_components <- list())
summary(model_components <- AddStudentLocalLinearTrend(model_components, y = y))
summary(model_components <- AddSeasonal(model_components, y = y, nseasons = 12))

fit <- bsts(y, model_components, niter = 2000)
summary(fit)
burnin <- 500 # Throw away first 500 
tibble(
  date = as.Date(time(data.1)),
  trend = colMeans(fit$state.contributions[-(1:burnin),"trend",]),
  seasonality = colMeans(fit$state.contributions[-(1:burnin),"seasonal.12.1",])) %>%
  gather("component", "value", trend, seasonality) %>%
  ggplot(aes(x = date, y = value)) + 
  geom_line() + theme_bw() + 
  theme(legend.title = element_blank()) + ylab("") + xlab("") +
  facet_grid(component ~ ., scales="free") + guides(colour=FALSE) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

pred <- predict(fit, horizon = 100, burn = burnin, quantiles = c(.05, .95))
plot(pred)
  

residuals(fit)
errors <- bsts.prediction.errors(fit, burn = 1000, )
PlotDynamicDistribution(errors)






## dlm - Vignette ####
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## 2. Stochastic Level
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
dlm(FF = 1, V = NA, GG = 1, W = NA, m0 = 0, C0 = 100)
dlmModPoly(order = 1, dV = NA , dW = exp(params[2]))

dlmMLE

## MLE Example 1
buildFun <- function(x) {
  dlmModPoly(order = 1, dV = exp(x[1]), dW = exp(x[2]))
    }
fit <- dlmMLE(Nile, parm = c(0,0), build = buildFun)
fit$conv
dlmNile <-  buildFun(fit$par)
V(dlmNile)
W(dlmNile)

## MLE Example 2
buildFun <- function(x) {
m <- dlmModPoly(1, dV = exp(x[1]))
m$JW <- matrix(1)
m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
j <- which(time(Nile) == 1899)
m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
return(m)
}


fit <- dlmMLE(Nile, parm = c(0,0,0), build = buildFun) 
fit$conv
dlmNileJump <-  buildFun(fit$par)
V(dlmNileJump)
dlmNileJump$X[c(1,which(time(Nile) == 1899)),1]
plot(Nile)

nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
plot(Nile, type = 'o', col = "seagreen")
lines(dropFirst(nileJumpFilt$m),
      type = 'o', pch = 20, col = "brown")

attach(nileJumpFilt)
v <- unlist(dlmSvd2var(U.C, D.C))
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")

nileJumpSmooth <- dlmSmooth(nileJumpFilt)
plot(Nile, type = 'o', col = "seagreen")
attach(nileJumpSmooth)
lines(dropFirst(s), type = 'o', pch = 20, col = "brown") > v <- unlist(dlmSvd2var(U.S, D.S))
pl <- dropFirst(s) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(s) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown")
lines(pu, lty = 2, col = "brown")


## MLE Example 3
y = log(UKgas)
par(mfrow=c(1,1))
plot(log(UKgas))

## Free-form seasonal factor DLM
##-#-#-#-#-#-#-#-#
## y_t = F*theta_t+v_t
## theta_t = G*theta_t-1+w_t
## theta_t = [mu_t, beta_t, alpha_1, alpha_2, alpha_3, alpha_4]'
## F = [1 0 1 0 0 ]
## G =  [ 1  1  0  0  0 
##        0  1  0  0  0
##        0  0 -1 -1 -1 
##        0  0  1  0  0
##        0  0  0  1  0 ]
##-#-#-#-#-#-#-#-#
dlm3 = dlmModPoly(order = 2) + dlmModSeas(4)

buildFun = function(x) {
  diag(W(dlm3))[2:3] = exp(x[1:2])
  V(dlm3) = exp(x[3])
  return(dlm3)
}

fit3 = dlmMLE(y, parm = c(0.1,0.1,0.1), build = buildFun)
dlm3 = buildFun(fit3$par)
ySmooth = dlmSmooth(y, mod = dlm3)
x = cbind(y, dropFirst(ySmooth$s[, c(1, 3)]))
colnames(x) = c("Gas", "Trend", "Seasonal")
par(mfrow=c(1,1))
plot(x, type = "o", main = "UK Gas Consumption")

# Fourier form representation of seasonality
dlm4 = dlmModPoly() + dlmModTrig(4)

buildFun = function(x) {
  diag(W(dlm4))[2:3] = exp(x[1:2])
  V(dlm4) = exp(x[3])
  return(dlm4)
}

fit4 = dlmMLE(y,parm=c(0.1,0.1,0.1),build=buildFun)

dlm4 = buildFun(fit4$par)

ySmooth = dlmSmooth(y, mod = dlm4)

x = cbind(y, dropFirst(ySmooth$s[, c(1, 3)]))

colnames(x) = c("Gas", "Trend", "Seasonal")

pdf(file="UKgas-harmonic.pdf",width=8,height=6)
par(mfrow=c(1,1))
plot(x, type = "o", main = "UK Gas Consumption")