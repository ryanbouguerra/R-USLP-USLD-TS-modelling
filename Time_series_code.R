setwd("~/Desktop/ST422 Project")
library(forecast)
library(lubridate)
library(tseries)
library(zoo)
library(astsa)

fuel <- read.csv("Fuel_price.csv")

ulsp <- msts(fuel[,2], seasonal=365.25/7,
             start = decimal_date(as.Date("2003-06-09")))
ulsd <- msts(fuel[,3], seasonal=365.25/7,
             start = decimal_date(as.Date("2003-06-09")))

## Question 1
# Logarithmic transformation on data

ulsp_log <- log(ulsp)
ulsd_log <- log(ulsd)

# Plotting time series and its associated ACF
plot(ulsp_log, main= "Log-prices of ultra-low sulphur petrol over time", xlab ="Time", ylab="Log-weekly prices" )
plot(ulsd_log, main= "Log-prices of ultra-low sulphur diesel over time", xlab ="Time", ylab="Log-weekly prices" )

acf2(ulsp_log, max.lag=100)
acf2(ulsd_log, max.lag=100)

# Taking log returns of data
return_ulsp <- diff(log(ulsp))
return_ulsd <- diff(log(ulsd))

# Plotting time series and its associated ACF

plot(return_ulsp, main= "Log-returns of ultra-low sulphur petrol over time", xlab ="Time", ylab="Log-weekly returns" )
acf2(return_ulsp) 

plot(return_ulsd, main= "Log-returns of ultra-low sulphur diesel over time", xlab ="Time", ylab="Log-weekly returns" )
acf2(return_ulsd)

# First order difference
dreturn_ulsp <- diff(return_ulsp)
dreturn_ulsd <- diff(return_ulsd)

# Plotting
plot(dreturn_ulsp, main= "First order difference of log-returns of ULSP over time", xlab ="Time", ylab="Difference of log returns" )
acf2(dreturn_ulsp) 

plot(dreturn_ulsd, main= "First order difference of log-returns of ULSD over time", xlab ="Time", ylab="Difference of log returns" )
acf2(dreturn_ulsd)



## Question 2
# Model for ULSP
# Start with auto.arima
auto_ulsp <- auto.arima(return_ulsp)
summary(auto_ulsp) # ARMA(1,1), AIC: -7101.69
tsdiag(auto_ulsp) # last point don't pass the test


# Best model with lowest AIC
mod_uslp <- Arima(return_ulsp, order=c(5,0,3)) # AIC -7106.01
summary(mod_uslp)
sarima(return_ulsp, 5,0,3)

# Try to add seasonality
mod_ulsp_season <- Arima(return_ulsp, order=c(5,0,3), seasonal = list(order=c(0,0,1), period = 52))
summary(mod_ulsp_season) # ARIMA(5,0,3)(0,0,1)[52], AIC: -7104.56
sarima(return_ulsp, 5,0,3,0,0,1,52)

## Predictions using mod_ulsp
# Log-return predictions
pred_ulsp_return <- forecast(mod_uslp, 5)
plot(pred_ulsp_return, 52, main="Forecasted log-returns for ULSP time series", xlab="Time", ylab="Log-returns")

# Log-prices predictions
pred_ulsp_price <- forecast(Arima(ulsp_log, order=c(5,1,3)), 5)
plot(pred_ulsp_price, 52, main="Forecasted log-prices for ULSP", xlab="Time", ylab="Log-prices")


# Model for ULSD

# Start with auto.arima
auto_ulsd <- auto.arima(return_ulsd) # limit the search to 10 models because of lengthy running time
summary(auto_ulsd) # suggests ARMA(1,0,1) with AIC: -7417.5
sarima(return_ulsd, 1,0,1)

# Best model with lowest AIC
mod_ulsd <- Arima(return_ulsd, order=c(3,0,2))
summary(mod_ulsd) # ARMA(3,2) with AIC: -7418.78 
sarima(return_ulsd, 3,0,2)

# Try to add seasonality
mod_ulsd_s <- Arima(dreturn_ulsd, order=c(3,0,2), seasonal = list(order=c(0,0,1), period = 52))
summary(mod_ulsd_s) # ARIMA(3,0,2)(0,0,1)[52], AIC: -7399.04
sarima(return_ulsd, 3,0,2,0,0,1,52)


## Predictions using mod_ulsd
# Log-return predictions
pred_ulsd_return <- forecast(mod_ulsd, 5)
plot(pred_ulsd_return, 52, main="Forecasted log-returns for ULSD time series", xlab="Time", ylab="Log-returns")

# Log-prices predictions
pred_ulsd_price <- forecast(Arima(ulsd_log, order=c(3,1,2)), 5)
plot(pred_ulsd_price, 52, main="Forecasted log-prices for ULSD", xlab="Time", ylab="Log-prices")



## Question 3

# Splitting the data into training and validation data
n <- length(fuel[,2])
training_length <- ceiling(n*0.60)

ulsp_training <- msts(fuel[,2][1:training_length], seasonal=365.25/7,
                      start = decimal_date(as.Date("2003-06-09")))

ulsp_validation <- msts(fuel[,2][(training_length+1):n], seasonal=365.25/7,
                        start = decimal_date(as.Date("2003-06-09")))

ulsd_training <- msts(fuel[,3][1:training_length], seasonal=365.25/7,
                      start = decimal_date(as.Date("2003-06-09")))

ulsd_validation <- msts(fuel[,3][(training_length+1):n], seasonal=365.25/7,
                        start = decimal_date(as.Date("2003-06-09")))

return_ulsp_training <- diff(log(ulsp_training))
return_ulsp_valid <- diff(log(ulsp_validation))
return_ulsd_training <- diff(log(ulsd_training))
return_ulsd_valid <- diff(log(ulsd_validation))

## Rolling window functions

pred.roll <- function(data, D, A) {
  # Function that predicts the next A data points, using a window of length D each time 
  # moving the window forward A data points. The fitted model for prediction is derived 
  # from auto.arima(). Computes the prediction error for each window and outputs the final 
  # aggregate error term, as well as, a list of fitted model.
  n <- length(data)
  N <- ceiling((n-D)/A) + 1
  fit <- list()
  pred_error <- 0
  for (i in 1:N){
    data.ts  <- msts(data[((i-1)*A+1):(min((i-1)*A+D,n))], seasonal=365.25/7 )
    fit[[i]] = auto.arima(data.ts)
    #tsdiag(fit[[i]])
    pred <- forecast(fit[[i]], A)
    if ((min((i-1)*A+D,n)) + A < n) {
      error <- sum(  (data[((min((i-1)*A+D,n))+1):((min((i-1)*A+D,n)) + A)] - pred$mean )^2 ) # predict next A data points
      pred_error <- pred_error + as.numeric(error) # aggregate error component
    }
  }
  return(list(fit, pred_error))
}


automated_pred.roll <- function(data, rep) {
  # Function that iterates the fit.roll (rep-3) times and outputs the window length 
  # with lowest prediction error. Take A = ceiling(D/3)
  min_error <- 1000 # some large integer that will never be exceeded
  D_min <- 0
  for (D in 4:rep) {
    A <- ceiling(D/3)
    temp <- pred.roll(data, D, A)
    if (temp[[2]] < min_error) {
      min_error <- temp[[2]]
      D_min <- D
    }
  }
  return(list(D_min, min_error))
}

# Optimal window length for ULSP
automated_pred.roll(return_ulsp_training, 32) #  D = 6, with error = 0.03391066
automated_pred.roll(return_ulsp_valid, 32) # D = 6, with error = 0.01729922
# Same window length

# Optimal window length for ULSD
automated_pred.roll(return_ulsd_training, 32) # D = 6, with error = 0.02670696
automated_pred.roll(return_ulsd_valid, 32) # D = 6, with error = 0.01426493
# Same window length



## Forecasting the next 5 data points
## Predictions
ulsp_wind_fcst <- forecast(auto.arima(return_ulsp[(n-4):n]), 5)
plot(ulsp_wind_fcst, main="Forecasted log-prices for ULSP using rolling window", xlab="Time", ylab="Log-prices")
ulsd_wind_fcst <- forecast(auto.arima(return_ulsd[(n-4):n]), 5)
plot(ulsd_wind_fcst, main="Forecasted log-prices for ULSD using rolling window", xlab="Time", ylab="Log-prices")

# Prediction points from question 2
pred_ulsp_return$mean
pred_ulsd_return$mean


## Question 4
# Perform decomposition of the ULSP and ULSD log-returns time series
ulsp_decomp <- mstl(return_ulsp, iterate = 100)
ulsd_decomp <- mstl(return_ulsd, iterate = 100)

# Plot each component
plot(ulsp_decomp, main="Decomposition of the log-return ULSP time series")
plot(ulsd_decomp, main="Decomposition of the log-return ULSD time series")

# We can forecast the data using the stlf() function
ulsp_decomp_fcst <- forecast.mstl(ulsp_decomp, h=5)
plot(ulsp_decomp_fcst, main="Forecasted log-returns for ULSP using multiple seasonality model", xlab="Time", ylab="Log-prices", include = 50)

ulsd_decomp_fcst <- forecast.mstl(ulsd_decomp, h=5)
plot(ulsd_decomp_fcst, main="Forecasted log-returns for ULSD using multiple seasonality model", xlab="Time", ylab="Log-prices", include = 50)

# Compare to previous prediction point
ulsp_decomp_fcst
pred_ulsp_return$mean

ulsd_decomp_fcst
pred_ulsd_return$mean

