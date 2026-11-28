

# funciones
# boxcox_tr(), boxcox_inv(), perfil_lambda()


boxcox_tr <- function(datos, lambda) {
  y <- if (is.data.frame(datos)) datos$y else datos
  if (any(y <= 0, na.rm = TRUE)) stop("Los valores deben ser positivos (> 0)")
  
  y_tr <- if (abs(lambda) < 1e-8) log(y) else (y^lambda - 1) / lambda
  
  if (is.data.frame(datos)) {
    datos$y <- y_tr
    attr(datos, "lambda") <- lambda
    return(datos)
  }
  return(y_tr)
}

boxcox_inv <- function(datos, lambda = attr(datos, "lambda")) {
  if (is.null(lambda)) stop("Debe indicar el valor de lambda")
  y <- if (is.data.frame(datos)) datos$y else datos
  
  y_inv <- if (abs(lambda) < 1e-8) exp(y) else (lambda * y + 1)^(1 / lambda)
  
  if (is.data.frame(datos)) {
    datos$y <- y_inv
    return(datos)
  }
  return(y_inv)
}


perfil_lambda <- function(datos, lambdas = seq(-2, 2, by = 0.05)) {
  y <- datos$y
  t <- datos$t
  n <- length(y)
  
  log_lik <- sapply(lambdas, function(lam) {
    z <- boxcox_tr(y, lam)
    rss <- sum(resid(lm(z ~ t))^2)
    -(n / 2) * log(rss / n) + (lam - 1) * sum(log(y))
  })
  
  optimo <- lambdas[which.max(log_lik)]
  ic <- range(lambdas[log_lik >= max(log_lik) - 0.5 * qchisq(0.95, 1)])
  
  plot(lambdas, log_lik, type = "l", xlab = "Lambda", ylab = "Log-Verosimilitud")
  abline(v = optimo, col = "red", lty = 2)
  abline(h = max(log_lik) - 0.5 * qchisq(0.95, 1), col = "blue", lty = 2)
  
  return(list(lambda_optimo = optimo, ic_95 = ic, perfil = data.frame(lambdas, log_lik)))
}



pendiente_dispersion <- function(datos, k = 5) {
  y    <- if (is.data.frame(datos)) datos$y else datos
  tam  <- length(y) %/% k                  # tamaño de cada bloque
  y    <- y[seq_len(k * tam)]              # descarta sobrantes, si los hay
  bloques <- split(y, rep(seq_len(k), each = tam))
  
  tabla <- data.frame(
    bloque   = seq_len(k),
    nivel    = unname(sapply(bloques, mean)),
    desv_dif = unname(sapply(bloques, function(b) sd(diff(b))))
  )
  
  fit  <- lm(log(desv_dif) ~ log(nivel), data = tabla)
  pend <- unname(coef(fit)[2])
  
  list(tabla = tabla, fit = fit, pendiente = pend, lambda = 1 - pend)
}
