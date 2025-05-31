# Load the required libraries
library(dplyr)
require(zoo)
require(tseries)
library(fUnitRoots)
library(astsa)
library(ellipse)

# -----------------------------------------------------------
# PART 1: Data Import and Preprocessing
# - Load monthly series data, clean it, and convert to time series object.
# - Visualize the series to check overall structure and trend.
# -----------------------------------------------------------
data <- read.csv("~/work/Séries Temp/valeurs_mensuelles.csv", sep=";", skip=3)
data_clean <- data[, -3]
colnames(data_clean) <- c("Périod", "Value")
dates <- as.yearmon(seq(from = 1990, to = 2025+2/12, by = 1/12))
serie_temp <- zoo(data_clean$Value, order.by = dates)
serie_temp <- rev(serie_temp)

head(serie_temp)
tail(serie_temp)
plot(serie_temp, main = "Time Series", xlab = "Time", ylab = "Value", col = "red")

# -----------------------------------------------------------
# Check for deterministic trend via linear regression.
# -----------------------------------------------------------
temps <- 1:length(serie_temp)
reg_tendance <- lm(coredata(serie_temp) ~ temps)
summary(reg_tendance)

# -----------------------------------------------------------
# Stationarity testing with augmented Dickey-Fuller (ADF) test.
# - Uses a custom function that iteratively chooses the number of lags.
# -----------------------------------------------------------
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK ? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05, na.rm=TRUE) == 0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("no \n")
    k <- k+1
  }
  return(adf)
}

adf <- adfTest_valid(serie_temp, 24, adftype="ct")

# Visual check of ACF and PACF of the original series
par(mfrow=c(1,2))
acf(serie_temp, 24); pacf(serie_temp, 24)
par(mfrow=c(1,1))

# -----------------------------------------------------------
# Differencing to achieve stationarity (first difference).
# -----------------------------------------------------------
serie_diff <- diff(serie_temp)
plot(serie_diff, main = "Differenced Series (d = 1)", xlab = "Time", ylab = "First Difference", col = "blue")

# Stationarity test on the differenced series
adf_diff <- adfTest_valid(serie_diff, 24, adftype="ct")
adf_diff

# -----------------------------------------------------------
# PART 2: Model Identification and Selection
# - Explore ARMA(p, q) models for p, q <= 4.
# - Check residuals for whiteness using Ljung-Box test.
# -----------------------------------------------------------
par(mfrow=c(1,2))
acf(serie_diff, 24); pacf(serie_diff, 24)
par(mfrow=c(1,1))

pmax <- 4
qmax <- 4
k_lags <- 24  

residu_OK <- matrix(NA, nrow = pmax + 1, ncol = qmax + 1)
rownames(residu_OK) <- paste0("p=", 0:pmax)
colnames(residu_OK) <- paste0("q=", 0:qmax)

for (p in 0:pmax) {
  for (q in 0:qmax) {
    modele <- try(arima(serie_diff, order = c(p, 0, q), include.mean = FALSE))
    
    if (!inherits(modele, "try-error")) {
      fitdf <- length(modele$coef)
      test_q <- Qtests(residuals(modele), k_lags, fitdf = fitdf)[, "pval"]
      residu_OK[p + 1, q + 1] <- all(test_q > 0.05, na.rm = TRUE)
    } else {
      residu_OK[p + 1, q + 1] <- NA
    }
  }
}

residu_OK

# Function to extract and summarize valid models
extract_valid_models <- function(residu_OK, series) {
  results <- list()
  for (i in 1:nrow(residu_OK)) {
    for (j in 1:ncol(residu_OK)) {
      if (isTRUE(residu_OK[i, j])) {
        p <- i - 1
        q <- j - 1
        mod <- try(arima(series, order = c(p, 0, q), include.mean = FALSE))
        if (!inherits(mod, "try-error")) {
          results[[length(results) + 1]] <- data.frame(
            p = p, q = q, AIC = mod$aic, BIC = BIC(mod)
          )
        }
      }
    }
  }
  if (length(results) == 0) {
    return(data.frame(p = integer(0), q = integer(0), AIC = numeric(0), BIC = numeric(0),
                      AIC_min = logical(0), BIC_min = logical(0)))
  }
  models_df <- do.call(rbind, results)
  models_df$AIC_min <- FALSE
  models_df$BIC_min <- FALSE
  idx_AIC_min <- which.min(models_df$AIC)
  idx_BIC_min <- which.min(models_df$BIC)
  models_df$AIC_min[idx_AIC_min] <- TRUE
  models_df$BIC_min[idx_BIC_min] <- TRUE
  return(models_df)
}

valid_models <- extract_valid_models(residu_OK, serie_diff)
valid_models

# -----------------------------------------------------------
# PART 3: Forecasting and Visualizing Results
# - Fit final ARIMA model, compute forecast error covariance.
# - Plot forecasts and 95% confidence ellipsoid.
# - Plot residual histogram.
# -----------------------------------------------------------
arima011 <- arima(serie_temp, order = c(0, 1, 1), include.mean = FALSE)
arima011

pred <- predict(arima011, n.ahead = 2)
Yhat <- pred$pred 

se <- arima011$sigma2
theta1 <- arima011$coef["ma1"]

Sigma <- se * matrix(
  c(1, theta1, theta1, 1 + theta1^2),
  nrow = 2, byrow = TRUE
)

print(Sigma)

serie_temp_ts <- ts(coredata(serie_temp), start = c(1990, 1), frequency = 12)
prev <- sarima.for(serie_temp_ts, n.ahead = 2, p = 0, d = 1, q = 1, plot = TRUE,
                   main = "2-step-ahead forecasts with 95% confidence intervals")

alpha <- 0.05
ellipse_points <- ellipse(Sigma, centre = Yhat, level = 1 - alpha)

plot(ellipse_points, type = 'l',
     xlab = expression(Y[T+1]), ylab = expression(Y[T+2]),
     main = "95% confidence ellipsoid")
points(Yhat[1], Yhat[2], col = 'blue', pch = 19)

res <- residuals(arima011)
hist(res, breaks = 20, main = "Histogram of ARIMA(0,1,1) residuals",
     xlab = "Residuals", col = "lightblue", border = "gray")
