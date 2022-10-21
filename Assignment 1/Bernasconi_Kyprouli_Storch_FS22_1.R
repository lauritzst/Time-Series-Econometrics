library(tidyverse)
library(xts)
library(zoo)
library(forecast)
library(tsbox)
library(modeldb)

rm(list=ls())

setwd("~/Desktop/Time Series Econometrics/Assignment 1")

# Create directory to store the figures
dir.create("figures/", showWarnings=F)
dir.create("tables/", showWarnings=F)


# Import data -------------------------------------------------------------

# STG.txt
stg <- read.table("STG.txt", skip=19, header=T)
# ZH.txt
zh <- read.table("ZH.txt", skip=19, header=T)


# 4a) ---------------------------------------------------------------------

# Convert to time series and create the data column (index)
stg <- xts(stg[, 3:4], order.by=as.yearmon(stg$Year+(stg$Month-1)/12))
zh <- xts(zh[,3:4], order.by=as.yearmon(zh$Year+(zh$Month-1)/12))

names(stg) <- paste0(tolower(names(stg)), '_sg')
names(zh) <- paste0(tolower(names(zh)), '_zh')

data <- merge(stg, zh, join='outer')

# Create new columns
data$temperature_diff <- data$temperature_zh - data$temperature_sg
data$precipitation_diff <- data$precipitation_zh - data$precipitation_sg

# Omit missing values
# data <- na.omit(data)


# 4b) ---------------------------------------------------------------------

# Combine line plots
png("figures/line_plots.png", width=1500, height=750)
par(mfrow=c(2,2))

# Plot temperatures and precipitations for SG
plot(data$temperature_sg, main='Temperature (SG)')
plot(data$precipitation_sg, main='Precipitation (SG)')

# Plot differences
plot(data$temperature_diff, main='Temperature difference (ZH-SG)')
plot(data$precipitation_diff, main='Precipitation difference (ZH-SG)')

# Close line plot
dev.off()

# Combine ACF plots
png("figures/acf_plots.png", width=1500, height=750)
par(mfrow=c(2,2))

# Plot ACF temperatures and precipitations for SG
acf(as.vector(data$temperature_sg), lag.max=24, main='Temperature (SG)', na.action=na.pass)
acf(as.vector(data$precipitation_sg), lag.max=24, main='Precipitation (SG)', na.action=na.pass)

# Plot ACF differences
acf(as.vector(data$temperature_diff), lag.max=24, main='Temperature difference (ZH-SG)', na.action=na.pass)
acf(as.vector(data$precipitation_diff), lag.max=24, main='Precipitation difference (ZH-SG)', na.action=na.pass)

# Close ACF plot
dev.off()


# 4e) ---------------------------------------------------------------------

# Decompose into trend and seasonal component using STL function
stl_decomp <- stl(tsbox::ts_ts(data$temperature_sg), 'periodic')
remainder_stl <- stl_decomp$time.series[,'remainder']

# Plot charts produced by decomposition
png("figures/stl.png", width=1500, height=750)
plot(stl_decomp)
dev.off()


# 4f) ---------------------------------------------------------------------

# Estimate seasonal component with trigonometric functions
S_12 <- 12
t <- seq_along(data$temperature_sg)
ols_data <- tibble(temperature_sg = as.vector(data$temperature_sg),
                   cos_component_12 = cos(2*pi*t/S_12),
                   sin_component_12 = sin(2*pi*t/S_12),
                   t = t,
                   t2 = t^2)
s_ols <- lm(temperature_sg ~ cos_component_12 + sin_component_12, data=ols_data)
s_hat <- s_ols$fitted.values
data_seasonally_adj <- data$temperature_sg - s_hat

# Estimate seasonal component with season dummies
# m <- format(index(data), "%m")
# ols_data <- tibble(temperature_sg = as.vector(data$temperature_sg),
#                    m = factor(m))
# ols_data <- ols_data %>% modeldb::add_dummy_variables(m, values=levels(ols_data$m))
# ols_fit <- lm(temperature_sg ~ ., data=ols_data)
# s_hat <- ols_fit$fitted.values
# data_seasonally_adj <- data$temperature_sg - s_hat

# Alternative way to remove seasonality
# data_seasonally_adj <- xts::diff.xts(data$temperature_sg, lag=12, differences=1)

# Plot seasonally-adjusted TS
png("figures/plot_wo_seasonality.png", width=750, height=375)
plot(data_seasonally_adj, main="Seasonally adjusted temperature (SG)")
dev.off()

# Plot ACF and PACF
png("figures/acf_wo_seasonality.png", width=750, height=375)
acf(as.vector(data_seasonally_adj), lag.max=24, main='Seasonally adjusted temperature (SG)', na.action=na.pass)
dev.off()


# 4g) ---------------------------------------------------------------------

# Plot ACF
png("figures/acf_remainder.png", width=750, height=375)
acf(as.vector(remainder_stl), lag.max=24, main='Remainder STL decomposition')
dev.off()

# Plot PACF
png("figures/pacf_remainder.png", width=750, height=375)
pacf(as.vector(remainder_stl), lag.max=24, main='Remainder STL decomposition')
dev.off()


# 4h) ---------------------------------------------------------------------

# Define ARMA model selection routine
select_arma <- function(ts, p_vec, q_vec){
  p_vec <- 0:2
  q_vec <- 0:2
  p_q_combs <- expand.grid(list(p=p_vec,q=q_vec))
  aic_vec <- c()
  bic_vec <- c()
  for (i in 1:nrow(p_q_combs)){
    p <- p_q_combs[i, 'p']
    q <- p_q_combs[i, 'q']
    fit_model <- Arima(ts, order=c(p, 0, q))
    # aic_vec[i] <- fit_model$aic
    # bic_vec[i] <- fit_model$bic
    aic_vec[i] <- AIC(fit_model)
    bic_vec[i] <- BIC(fit_model)
  }
  table <- tibble(p = p_q_combs$p,
                  q = p_q_combs$q,
                  AIC = aic_vec,
                  BIC = bic_vec)
  aic_best <- table[which.min(table$AIC),]
  bic_best <- table[which.min(table$BIC),]
  aic_msg <- paste0("Best model according to AIC: ARMA(", aic_best$p, ",", aic_best$q, ")")
  bic_msg <- paste0("Best model according to BIC: ARMA(", bic_best$p, ",", bic_best$q, ")")
  print(aic_msg)
  print(bic_msg)
  return(table)
}


# 4i) ---------------------------------------------------------------------

# Apply routine to fit the data
arma_table <- select_arma(remainder_stl)
print(arma_table)
print(xtable::xtable(arma_table, type="latex"), file="tables/arma_table.tex", include.rownames=FALSE)

# Check results vs. auto.arima
arma_aic_opt <- auto.arima(remainder_stl, d=0, seasonal=F, max.p=2, max.q=2, ic='aic')
arma_bic_opt <- auto.arima(remainder_stl, d=0, seasonal=F, max.p=2, max.q=2, ic='bic')

# Fit an ARMA(4,4)
arma_4_4 <- Arima(remainder_stl, order=c(4,0,4))
print(paste0("AIC of ARMA(4,4): ", round(AIC(arma_4_4), 2)))
print(paste0("BIC of ARMA(4,4): ", round(BIC(arma_4_4), 2)))

