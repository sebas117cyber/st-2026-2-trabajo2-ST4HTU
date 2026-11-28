


# raices_modelo(), ljung_box(), sobreajuste()


raices_modelo <- function(ajuste) {
  mod <- ajuste$modelo
  coefs <- coef(mod)
  
  ar <- coefs[grep("^ar", names(coefs))]
  ma <- coefs[grep("^ma", names(coefs))]
  
  raices_ar <- if (length(ar) > 0) polyroot(c(1, -ar)) else numeric(0)
  raices_ma <- if (length(ma) > 0) polyroot(c(1, ma)) else numeric(0)
  
  # Gráfico de inversos de las raíces en el círculo unitario
  theta <- seq(0, 2 * pi, length.out = 200)
  plot(cos(theta), sin(theta), type = "l", lty = 2, asp = 1,
       xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
       xlab = "Parte Real", ylab = "Parte Imaginaria",
       main = "Inversos de las raíces (Círculo Unitario)")
  abline(h = 0, v = 0, col = "gray")
  
  if (length(raices_ar) > 0) points(1 / raices_ar, pch = 19, col = "blue")
  if (length(raices_ma) > 0) points(1 / raices_ma, pch = 17, col = "red")
  legend("topright", legend = c("AR", "MA"), pch = c(19, 17), col = c("blue", "red"), bty = "n")
  
  return(list(
    raices_ar = raices_ar,
    modulos_ar = Mod(raices_ar),
    estacionario = all(Mod(raices_ar) > 1),
    raices_ma = raices_ma,
    modulos_ma = Mod(raices_ma),
    invertible = all(Mod(raices_ma) > 1)
  ))
}


ljung_box <- function(ajuste, rezagos = c(6, 12, 18, 24)) {
  res <- ajuste$residuos
  p <- ajuste$orden[1]
  q <- ajuste$orden[3]
  fitdf <- p + q
  
  # Filtrar rezagos válidos (m > p + q)
  rezagos <- rezagos[rezagos > fitdf]
  
  tabla <- do.call(rbind, lapply(rezagos, function(m) {
    test <- Box.test(res, lag = m, type = "Ljung-Box", fitdf = fitdf)
    data.frame(
      rezago = m,
      gl = m - fitdf,
      estadistico_Q = as.numeric(test$statistic),
      p_valor = test$p.value
    )
  }))
  
  return(tabla)
}

sobreajuste <- function(z, ajuste, lambda = 1, constante = TRUE, tendencia = FALSE) {
  ord <- ajuste$orden
  metodo <- ajuste$metodo
  
  # Modelos candidatos: aumentar p en 1 y aumentar q en 1
  ord_ar <- c(ord[1] + 1, ord[2], ord[3])
  ord_ma <- c(ord[1], ord[2], ord[3] + 1)
  
  mod_ar <- ajustar_modelo(z, orden = ord_ar, lambda = lambda, 
                           constante = constante, tendencia = tendencia, metodo = metodo)
  mod_ma <- ajustar_modelo(z, orden = ord_ma, lambda = lambda, 
                           constante = constante, tendencia = tendencia, metodo = metodo)
  
  # Función auxiliar para extraer estadístico t del parámetro añadido
  extraer_extra <- function(mod_nuevo, nombre_param) {
    coefs <- mod_nuevo$coeficientes
    est <- coefs[nombre_param, "estimacion"]
    ee  <- coefs[nombre_param, "error_estandar"]
    t_val <- est / ee
    p_val <- 2 * (1 - pnorm(abs(t_val)))
    c(estimacion = est, error_estandar = ee, t_valor = t_val, p_valor = p_val)
  }
  
  param_ar <- paste0("ar", ord_ar[1])
  param_ma <- paste0("ma", ord_ma[3])
  
  resumen_params <- rbind(
    extraer_extra(mod_ar, param_ar),
    extraer_extra(mod_ma, param_ma)
  )
  rownames(resumen_params) <- c(
    paste0("Sobreajuste_AR (", paste(ord_ar, collapse = ","), ") -> ", param_ar),
    paste0("Sobreajuste_MA (", paste(ord_ma, collapse = ","), ") -> ", param_ma)
  )
  
  comparacion_crit <- rbind(
    Base = criterios(ajuste),
    Sobreajuste_AR = criterios(mod_ar),
    Sobreajuste_MA = criterios(mod_ma)
  )
  
  return(list(
    parametros_extra = as.data.frame(resumen_params),
    criterios = comparacion_crit,
    modelo_ar = mod_ar,
    modelo_ma = mod_ma
  ))
}


# Asimetría, curtosis y prueba de Jarque-Bera (JB ~ χ²(2) bajo H0 de normalidad)
jarque_bera <- function(x) {
  n    <- length(x)
  z    <- x - mean(x)
  s2   <- mean(z^2)
  asim <- mean(z^3) / s2^(3 / 2)
  curt <- mean(z^4) / s2^2
  jb   <- n / 6 * (asim^2 + (curt - 3)^2 / 4)
  list(asimetria = asim, curtosis = curt, JB = jb, gl = 2,
       p_valor = pchisq(jb, df = 2, lower.tail = FALSE))
}
