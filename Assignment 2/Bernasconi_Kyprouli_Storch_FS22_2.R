library(tseries)
library(readxl)
library(forecast)
library(urca)
library(stargazer)
library(texreg)
library(ggplot2)
library(lmtest)
library(stargazer)


rm(list=ls())

# Set working directory
setwd("~/Desktop/Time Series Econometrics/Assignment 2")

# Create directory to store the figures
dir.create("figures/", showWarnings=F)
dir.create("tables/", showWarnings=F)

# Set ggplot theme
# theme_set(theme_gray())


# 1) ----------------------------------------------------------------------

# Import the file
data <- read_xlsx("MonthlyChocolateSales2021.xlsx", sheet='Timeseries')

# Create ts
start <- c(2016, 1)
t <- ts(data, start=start, frequency=12)
t

# 1a) ---------------------------------------------------------------------

# Plot the ts
png("figures/plot_chocolate.png", width=750, height=375)
# plot(t)
autoplot(t) + 
  geom_line() + 
  labs(y='Average sales (in million CAD)', x='Year')
dev.off()


# 1b) ---------------------------------------------------------------------

# Seasonal plot of the data
png("figures/ggseasonplot_chocolate.png", width=750, height=375)
# seasonplot(t, s=12)
ggseasonplot(t, year.labels=T) + 
  labs(title="Seasonal plot",
       y="Average sales (in million CAD)")
dev.off()

# ACF plot (max 36 lags)
png("figures/acf_chocolate.png", width=750, height=375)
acf(as.vector(t), lag.max=36, main="ACF of average sales")
dev.off()


# 1c) ---------------------------------------------------------------------

# Fit ARIMA model - do not allow SARIMA
arima <- forecast::auto.arima(t, seasonal=F)
print(arima)
# texreg(arima)

# Diagnostics
checkresiduals(arima)
texreg(checkresiduals(arima))

# Fit ARIMA model - allowing SARIMA
sarima <- forecast::auto.arima(t)
print(sarima)
# texreg(sarima) # Latex output

# Diagnostics
checkresiduals(sarima)
texreg(checkresiduals(sarima))


# 1d) ---------------------------------------------------------------------

# Forecast 24 steps ahead (ARIMA)
fcast <- forecast::forecast(arima, h=24)
print(fcast)

# Plot predictions
png("figures/forecast.png", width=750, height=375)
autoplot(fcast, PI=T) + 
  labs(title = "Forecast based on ARIMA(1,1,1) model",
       y="Average sales (in million CAD)", 
       x='Year') + 
  autolayer(fcast$mean, series="Forecast", size=1) +
  guides(colour=guide_legend(title="Data series")) +
  scale_color_manual(values=c("darkblue"))
dev.off()

rm(fcast)

# Forecast 24 steps ahead (SARIMA)
fcast_sarima <- forecast::forecast(sarima, h=24)
print(fcast_sarima)

# Plot predictions
png("figures/forecast_sarima.png", width=750, height=375)
autoplot(fcast_sarima, PI=T) + 
  labs(title = "Forecast based on ARIMA(0,1,1)(0,1,1)[12] model",
       y="Average sales (in million CAD)", 
       x='Year') + 
  autolayer(fcast_sarima$mean, series="Forecast", size=1) +
  guides(colour=guide_legend(title="Data series")) +
  scale_color_manual(values=c("darkblue"))
dev.off()


# 2) ----------------------------------------------------------------------

# Import the file
data <- read_xlsx("CHMarriageAges.xlsx", sheet='px-x-0102020202_110',
                  skip=3, na="...", n_max=220, 
                  col_types=c("numeric", "numeric", "numeric"),
                  col_names=c("date", "male", "female"))

# Create ts
start <- data$date[1]
end <- data$date[nrow(data)]
t <- ts(data[,2:ncol(data)], start=start, end=end, frequency=1)

# Limit TS to 1855 to 2020
t <- window(t, 1855)


# 2a) ---------------------------------------------------------------------

# Plot the ts
png("figures/plot_marriage_age.png", width=750, height=375)
autoplot(t) +
  labs(x='Year', 
       y='Average age',
       title="Average age at 1st marriage (CH)")
dev.off()


# 2b) ---------------------------------------------------------------------

# Compute spread of the two time series
spread <- t[,'male'] - t[,'female']

# Plot the spread
png("figures/plot_spread.png", width=750, height=375)
autoplot(spread) + 
  labs(x="Year",
       y="Spread", 
       title="Spread in average age at 1st marriage (CH)")
dev.off()

# ACF spread
png("figures/acf_spread.png", width=750, height=375)
acf(spread, main="Spread in average age at 1st marriage (CH)")
dev.off()

# PACF spread
png("figures/pacf_spread.png", width=750, height=375)
pacf(spread, main="Spread in average age at 1st marriage (CH)")
dev.off()

# Take the first difference of spread
spread_diff <- diff(spread)

# ACF spread diff
png("figures/acf_spread_diff.png", width=750, height=375)
acf(spread_diff, main="Diff_spread in average age at 1st marriage (CH) - 1st difference")
dev.off()

# PACF spread diff
png("figures/pacf_spread_diff.png", width=750, height=375)
pacf(spread_diff, main="Diff_spread in average age at 1st marriage (CH) - 1st difference")
dev.off()

# Check with auto.arima
auto.arima(spread)
auto.arima(spread_diff)

# Latex tables
# texreg(list(auto.arima(spread),auto.arima(spread_diff)))


# 2c) ---------------------------------------------------------------------

# Keep TS only from 1975 onwards
t <- window(t, 1975)
png("figures/plot_1975.png", width=750, height=375)
plot(t)
dev.off()

# Plot the restricted data
png("figures/plot_marriage_age_from_1975.png", width=750, height=375)
acf(spread_diff, main="Sprad in average age \n at 1st marriage (CH) - 1st difference")
dev.off()

# Have a look at the kind of model suggested by auto.arima
auto.arima(t[, 'male'])  # (0,2,1)
auto.arima(t[, 'female'])  # (0,2,1)

# Plot ACFs to determine lags
# acf(t[, 'male'])
# acf(t[, 'female'])

# Conduct ADF tests by manual lag selection
# 1) start with many lags
adf_test_m <- ur.df(t[, 'male'], type="trend", lags=12)
adf_test_f <- ur.df(t[, 'female'], type="trend", lags=12)
# 2) select significant lags
summary(adf_test_m)
summary(adf_test_f)

# 3) run ADF again with fewer lags (in this case: none)
adf_test_m <- ur.df(t[, 'male'], type="trend", lags=0)
adf_test_f <- ur.df(t[, 'female'], type="trend", lags=0)
# 4) look at summary
summary(adf_test_m)
summary(adf_test_f)

# Alternative: 
# 1) conduct ADF test with automatic lag selection via BIC
# adf_test_m <- ur.df(t[, 'male'], type="trend", selectlags='BIC')
# adf_test_f <- ur.df(t[, 'female'], type="trend", selectlags='BIC')
# 2) look at summary of ADF tests
# summary(adf_test_m)
# summary(adf_test_f)

# --> not rejected in both cases: unit root with drift
#  because phi3 cannot be rejected!


# Plot summary ADF tests
png("figures/plot_adf_male.png", width=750, height=375)
plot(adf_test_m)
dev.off()
png("figures/plot_adf_female.png", width=750, height=375)
plot(adf_test_f)
dev.off()

# 2d) ---------------------------------------------------------------------

# Optional: Co-integration check

# 1) Run first step
ols_1 <- lm(male ~ female, data=t)
# Latex table
# stargazer(ols_1)

# 2) Get residuals
u_hat <- ols_1$residuals

# 3) Run ADF test on residuals
adf <- ur.df(u_hat, type='none', selectlags='BIC')
summary(adf)

png("figures/plot_adf_res.png", width=750, height=375)
plot(adf)
dev.off()

