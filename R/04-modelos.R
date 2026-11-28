

# ajustar_modelo(), criterios()

ajustar_modelo <- function(z, orden, lambda = 1, constante = TRUE, tendencia = FALSE, metodo = "ML") {
  y <- if (is.data.frame(z)) z$y else z
  t_vec <- if (is.data.frame(z)) z$t else seq_along(y)
  
  y_tr <- boxcox_tr(y, lambda)
  xreg <- if (tendencia) matrix(t_vec, ncol = 1, dimnames = list(NULL, "tendencia")) else NULL
  
  mod <- stats::arima(y_tr, order = orden, xreg = xreg, include.mean = constante, method = metodo)
  
  coefs <- cbind(estimacion = coef(mod), error_estandar = sqrt(diag(vcov(mod))))
  
  return(list(
    modelo = mod,
    coeficientes = coefs,
    sigma2 = mod$sigma2,
    matriz_covarianzas = vcov(mod),
    residuos = as.numeric(residuals(mod)),
    orden = orden,
    metodo = metodo
  ))
}

criterios <- function(ajuste) {
  mod <- ajuste$modelo
  n_eff <- length(ajuste$residuos) - ajuste$orden[2]
  k <- length(coef(mod)) + 1
  
  loglik <- if (ajuste$metodo == "CSS") {
    -0.5 * n_eff * (log(2 * pi) + log(ajuste$sigma2) + 1)
  } else {
    as.numeric(logLik(mod))
  }
  
  aic <- -2 * loglik + 2 * k
  aicc <- aic + (2 * k * (k + 1)) / (n_eff - k - 1)
  bic <- -2 * loglik + k * log(n_eff)
  
  return(data.frame(logLik = loglik, AIC = aic, AICc = aicc, BIC = bic))
}



# Añadido
# Traduce la especificación determinista a los argumentos de ajustar_modelo()
#   espec = "none"  : sin constante ni tendencia
#   espec = "drift" : constante (d = 0) o deriva θ0 (d = 1)
#   espec = "trend" : recta más ARMA (solo d = 0)
espec_args <- function(espec = c("none", "drift", "trend"), d) {
  espec <- match.arg(espec)
  list(
    constante = espec != "none" & d == 0,
    tendencia = (espec == "trend" & d == 0) | (espec == "drift" & d == 1)
  )
}

# Tabla de coeficientes con IC del 95 %
tabla_coef <- function(ajuste, nivel = 0.95) {
  cf <- ajuste$coeficientes
  z  <- qnorm((1 + nivel) / 2)
  tibble(
    parametro      = rownames(cf),
    estimacion     = unname(cf[, "estimacion"]),
    error_estandar = unname(cf[, "error_estandar"]),
    li             = estimacion - z * error_estandar,
    ls             = estimacion + z * error_estandar
  )
}
