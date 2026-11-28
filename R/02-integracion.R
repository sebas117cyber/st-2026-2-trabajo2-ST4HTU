

# diferencias y pruebas de Dickey y Fuller


# Diferencias de orden d de un vector
diferenciar <- function(x, d = 1) {
  if (d == 0) return(as.numeric(x))
  as.numeric(diff(x, differences = d))
}

# Valores críticos de Dickey y Fuller
# Hamilton (1994), Tabla B.6, basada en Fuller (1976).

valores_criticos_df <- function(tipo = c("none", "drift", "trend"), n) {
  tipo <- match.arg(tipo)
  tablas <- list(
    none = data.frame(
      n   = c(25, 50, 100, 250, 500, Inf),
      p01 = c(-2.66, -2.62, -2.60, -2.58, -2.58, -2.58),
      p05 = c(-1.95, -1.95, -1.95, -1.95, -1.95, -1.95),
      p10 = c(-1.60, -1.61, -1.61, -1.62, -1.62, -1.62)
    ),
    drift = data.frame(
      n   = c(25, 50, 100, 250, 500, Inf),
      p01 = c(-3.75, -3.58, -3.51, -3.46, -3.44, -3.43),
      p05 = c(-3.00, -2.93, -2.89, -2.88, -2.87, -2.86),
      p10 = c(-2.63, -2.60, -2.58, -2.57, -2.57, -2.57)
    ),
    trend = data.frame(
      n   = c(25, 50, 100, 250, 500, Inf),
      p01 = c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
      p05 = c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41),
      p10 = c(-3.24, -3.18, -3.15, -3.13, -3.13, -3.12)
    )
  )
  tb   <- tablas[[tipo]]
  fila <- which.min(abs(tb$n - n))         # fila con el tamaño de muestra más cercano
  c(n_tabla = tb$n[fila], p01 = tb$p01[fila],
    p05 = tb$p05[fila], p10 = tb$p10[fila])
}

# Prueba DF aumentada:
#   Δy_t = [c] + [β t] + γ y_{t-1} + Σ_{j=1}^p φ_j Δy_{t-j} + e_t
#   H0: γ = 0 (raíz unitaria)   H1: γ < 0 (estacionaria)
# Todos los rezagos se estiman sobre la misma muestra (t = pmax+2, ..., n),
# para que AIC/BIC sean comparables entre p.
prueba_df_aumentada <- function(x, tipo = c("trend", "drift", "none"),
                                pmax = 8, criterio = c("AIC", "BIC"),
                                p_fijo = NULL) {
  tipo     <- match.arg(tipo)
  criterio <- match.arg(criterio)
  
  y   <- as.numeric(x)
  n   <- length(y)
  dy  <- diff(y)
  idx <- (pmax + 2):n                      # índices t de la variable respuesta
  N   <- length(idx)
  Y   <- dy[idx - 1]                       # Δy_t
  
  construir <- function(p) {
    cols <- list()
    if (tipo != "none")  cols$const <- rep(1, N)
    if (tipo == "trend") cols$tend  <- idx
    cols$ylag <- y[idx - 1]
    for (j in seq_len(p)) cols[[paste0("dlag", j)]] <- dy[idx - 1 - j]
    cols$Y <- Y
    as.data.frame(cols)
  }
  
  crit <- vapply(0:pmax, function(p) {
    f <- lm(Y ~ . - 1, data = construir(p))
    if (criterio == "AIC") AIC(f) else BIC(f)
  }, numeric(1))
  
  p_opt <- if (is.null(p_fijo)) (0:pmax)[which.min(crit)] else p_fijo
  fit   <- lm(Y ~ . - 1, data = construir(p_opt))
  
  cf    <- summary(fit)$coefficients
  gamma <- cf["ylag", "Estimate"]
  ee    <- cf["ylag", "Std. Error"]
  tau   <- gamma / ee
  cv    <- valores_criticos_df(tipo, N)
  
  list(
    tipo                = tipo,
    p                   = p_opt,
    criterio            = criterio,
    N                   = N,
    gamma               = unname(gamma),
    error_estandar      = unname(ee),
    tau                 = unname(tau),
    valor_critico_05    = unname(cv["p05"]),
    tabla_cv            = cv,
    rechaza_H0          = unname(tau < cv["p05"]),
    criterio_por_rezago = data.frame(p = 0:pmax, valor = crit),
    coeficientes        = cf
  )
}

# Estrategia secuencial: con tendencia, con constante y sin nada
tabla_df_secuencial <- function(x, pmax = 8, criterio = "AIC") {
  tipos <- c("trend", "drift", "none")
  res   <- lapply(tipos, function(tp) prueba_df_aumentada(x, tp, pmax, criterio))
  data.frame(
    especificacion   = c("Con constante y tendencia", "Con constante",
                         "Sin constante ni tendencia"),
    p                = sapply(res, `[[`, "p"),
    N                = sapply(res, `[[`, "N"),
    tau              = sapply(res, `[[`, "tau"),
    valor_critico_05 = sapply(res, `[[`, "valor_critico_05"),
    rechaza_H0       = sapply(res, `[[`, "rechaza_H0"),
    row.names        = NULL
  )
}

# Verificación (bono): contraste contra urca::ur.df() con los mismos rezagos.
# Para que la muestra coincida, pmax debe ser igual a p.
verificar_df_urca <- function(x, tipo = c("trend", "drift", "none"), p) {
  tipo    <- match.arg(tipo)
  propia  <- prueba_df_aumentada(x, tipo, pmax = p, p_fijo = p)$tau
  ur      <- urca::ur.df(as.numeric(x), type = tipo, lags = p)
  externa <- as.numeric(ur@teststat)[1]          # tau1, tau2 o tau3 según el tipo
  data.frame(tau_propia = propia, tau_urca = externa,
             diferencia = abs(propia - externa))
}

# Error estándar de la media con corrección de autocorrelación
# (Newey-West con núcleo de Bartlett)
se_media_nw <- function(w, L = floor(4 * (length(w) / 100)^(2 / 9))) {
  N   <- length(w)
  u   <- w - mean(w)
  gam <- function(k) sum(u[(k + 1):N] * u[1:(N - k)]) / N
  V   <- gam(0)
  for (k in seq_len(L)) V <- V + 2 * (1 - k / (L + 1)) * gam(k)
  sqrt(V / N)
}
