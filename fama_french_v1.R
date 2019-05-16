## LOAD LIBRARIES =============================================================
library(tidyquant)
library(tidyverse)
library(tibbletime)
library(timetk)
library(dlm)
library(forecast)

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
X <- ff_portfolio_returns %>% select(MKT_RF, SMB, HML) %>% as.matrix()

build_reg_seas<- function(param) {
dlm1             <- dlmModReg(X = X, addInt = FALSE) + dlmModSeas(frequency = 12)
V(dlm1)          <- exp(param[1])
diag(W(dlm1))[c(1:4)] <- exp(param[2:5])
return(dlm1)
}

build_reg<- function(param) {
  dlm1             <- dlmModReg(X = X, addInt = FALSE)
  V(dlm1)          <- exp(param[1])
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
dlm_betas <- as_tibble(yFilter$m[-1,]) %>% mutate(date = ff_portfolio_returns$date)

dlm_betas %>% 
  gather("factor", "beta", -date) %>%
  left_join(ff_roll_betas, by = c("date", "factor")) %>%
  rename(dlm = beta) %>% 
  gather("method", "beta", -c(date, factor)) %>% 
  separate(method, c("del", "method"), fill = "left") %>%
  select(-del) %>%
  ggplot(aes(x = date, y = beta, col = method))+
  geom_line()+
  facet_wrap(~factor)

plot.ts(yFilter$m[-c(1:50),], plot.type = "s", col = c(1:3))

  ggplot(aes(x = DATE, y = beta, col = factor))+
  geom_line()






## BAYESIAN DLM ===============================================================