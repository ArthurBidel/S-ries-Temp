library(dplyr)
require(zoo)
require(tseries)
library(fUnitRoots)
library(ellipse)

# Partie 1
#import et création de la série
data <- read.csv("~/work/Séries Temp/valeurs_mensuelles.csv",sep=";",skip=3)
data_clean <- data[, -3]
colnames(data_clean) <- c("Périod", "Value")
dates <- as.yearmon(seq(from = 1990, to = 2025+2/12, by = 1/12))
serie_temp <- zoo(data_clean$Value, order.by = dates)
serie_temp <- rev(serie_temp)

# visualisation série
head(serie_temp)
tail(serie_temp)
plot(serie_temp, main = "Time Series", xlab = "Time", ylab = "Value", col = "red")


# test pour valider la tendance négative 
temps <- 1:length(serie_temp)
reg_tendance <- lm(coredata(serie_temp) ~ temps)
summary(reg_tendance)
# La série présente une tendance déterministe fortement significative.
# La pente négative (
#   β<0
#   β<0) indique que la série diminue linéairement au cours du temps.
# En plus comme non démoyennisée alors il faut mettre constante + trend dans les régression 

# Test Adf 
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
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("no \n")
    k <- k+1
  }
  return(adf)
}

adf <- adfTest_valid(serie_temp,24,adftype="ct")

# lags 21 
# Tu as testé automatiquement tous les lags de 0 à 21 dans le test ADF.
# À chaque fois, tu as vérifié que les résidus du modèle ADF n’étaient pas autocorrélés, à l’aide du test de Ljung-Box (jusqu’à 24 lags).
# À chaque fois, au moins une p-value < 0.05 → ❌ résidus autocorrélés.
# À 21 lags, enfin toutes les p-values > 0.05 → ✅ résidus non autocorrélés → le test ADF devient valide.

adf # avec 21 lags
# p-value 0.6328 donc pas de stationnarité 

par(mfrow=c(1,2))
acf(serie_temp, 24);pacf(serie_temp, 24)
par(mfrow=c(1,1))

# 🔹 À gauche : ACF (Autocorrélation simple)
# Toutes les barres sont positives, décroissantes lentement, très au-dessus des bandes bleues.
# Cela suggère une forte persistance dans la série, typique d’une non-stationnarité.
# Comportement classique d’un processus AR(1) non stationnaire ou avec racine unitaire.
# 🔹 À droite : PACF (Autocorrélation partielle)
# Forte valeur au premier retard (lag = 1), puis les valeurs deviennent faibles ou insignifiantes.
# Cela suggère que si on veut modéliser cette série en AR, un AR(1) ou AR(p) avec petit p pourrait être raisonnable après stationnarisation


# différentiation première 
serie_diff <- diff(serie_temp)
plot(serie_diff, main = "Differenced Series (d = 1)", xlab = "Time", ylab = "First Difference", col = "blue")

# Test sur la série différenciée
adf_diff <- adfTest_valid(serie_diff,24,adftype="ct")
adf_diff

# la p-value est de 0.01 donc la série est stationnaire 

# Partie 2 
# Autocorrélogramme 
par(mfrow=c(1,2))
acf(serie_diff, 24);pacf(serie_diff, 24)
par(mfrow=c(1,1))


# On va tester pour p <= 4 et q <= 4 bien que quelques autres petits lags marginalement significatifs
pmax <- 4
qmax <- 4

# Résidus non auto-corrélés
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
            p = p,
            q = q,
            AIC = mod$aic,
            BIC = BIC(mod)
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





# Partie 3 
arima011 <- arima(serie_diff, order = c(0, 0, 1), include.mean = FALSE)
# estimated MA(1) coefficient is -0.257


pred <- predict(arima011, n.ahead = 2)
Xhat <- pred$pred 


se_1 <- arima011$sigma2
theta1 <- arima011$coef["ma1"]

# Construction de la matrice de variance-covariance pour horizon 2
Sigma <- se_1 * matrix(
  c(1,
    theta1,
    theta1,
    1 + theta1^2),
  nrow = 2, byrow = TRUE
)

# Affichage
print(Sigma)

# Affichage de l’ellipsoïde de confiance à 95%
alpha <- 0.05
ellipse_points <- ellipse(Sigma, centre = Xhat, level = 1 - alpha)

plot(ellipse_points, type = 'l',
     xlab = expression(X[T+1]), ylab = expression(X[T+2]),
     main = "Confidence ellipsoid 95%")

points(Xhat[1], Xhat[2], col = 'blue', pch = 19)


