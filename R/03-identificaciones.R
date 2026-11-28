

# ACF, PACF y EACF muestrales


# =============================================================
# 03-identificacion.R
# ACF, PACF y EACF muestrales con sus bandas
# =============================================================


correlograma <- function(datos,m = NULL){
  
  # Extraemos la serie numérica
  datos <- datos$y
  
  
  T <- length(datos) #Numero de observaciones
  
  #Valida si m es null o diferente 
  if(is.null(m)){
    m <- min(floor(T/4),24)
  }
  
  #-----------Funcion ACF Manual------------------#
  ACF <- function(datos, h) {
    
    T <- length(datos)
    promedio <- sum(datos) / T
    denominador <- sum((datos - promedio)^2)
    
    producto <- c()
    for (i in 1:T) {
      resultado <- datos[i] - promedio
      producto <- c(producto, resultado)
    }
    
    valor <- c()
    for (i in 1:(T - h)) {
      resultado <- producto[i + h] * producto[i]
      valor <- c(valor, resultado)
    }
    
    acf <- sum(valor) / denominador
    
    return(acf)
  }
  
  # Calculamos la ACF manual desde 1 hasta m
  acf_manual <- c()
  
  for (h in 1:m) {
    resultado <- ACF(datos, h)
    acf_manual <- c(acf_manual, resultado)
  }
  
  # ACF calculada por R
  acf_r <- acf(datos, lag.max = m, plot = FALSE)$acf[-1]
  
  # Diferencias absolutas entre la ACF manual y la de R
  diferencias <- abs(acf_manual - acf_r)
  
  # Maxima diferencia
  diferencia_maxima <- max(diferencias)
  
  #-----------Calculamos la PACF con la funcion de R---------------#
  pacf_r <- pacf(datos, lag.max = m, plot = FALSE)$acf
  
  # Creamos los datos para graficar la ACF
  datos_acf <- data.frame(
    rezago = 1:m,
    valor = acf_manual
  )
  
  # Creamos los datos para graficar la PACF
  datos_pacf <- data.frame(
    rezago = 1:m,
    valor = as.numeric(pacf_r)
  )
  
  # Calculamos la banda de confianza del 95 %
  banda <- qnorm((1 + 0.95) / 2) / sqrt(T)
  
  # Funcion para crear los graficos ACF y PACF
  graficar_correlacion <- function(datos, titulo, etiqueta_y) {
    
    ggplot(datos, aes(x = rezago, y = valor)) +
      
      # Barras verticales
      geom_segment(
        aes(xend = rezago, yend = 0)
      ) +
      
      # Linea central
      geom_hline(yintercept = 0) +
      
      # Banda superior
      geom_hline(
        yintercept = banda,
        linetype = "dashed"
      ) +
      
      # Banda inferior
      geom_hline(
        yintercept = -banda,
        linetype = "dashed"
      ) +
      
      labs(
        title = titulo,
        x = "Rezago",
        y = etiqueta_y
      ) +
      
      theme_minimal()
  }
  
  # Grafico de la ACF
  grafico_acf <- graficar_correlacion(
    datos_acf,
    "Función de Autocorrelación (ACF)",
    "ACF"
  )
  
  # Grafico de la PACF
  grafico_pacf <- graficar_correlacion(
    datos_pacf,
    "Función de Autocorrelación Parcial (PACF)",
    "PACF"
  )
  
  # Unimos los graficos
  panel <- grafico_acf / grafico_pacf
  
  return(list(acf =acf_manual, panel = panel,diferencia_maxima= diferencia_maxima))
  
}

# Bandas de Bartlett para la ACF: bajo H0 de MA(k-1), Var(r_k) ≈ (1 + 2 Σ_{j<k} r_j²)/N
bandas_bartlett <- function(rho, N, nivel = 0.95) {
  m     <- length(rho)
  var_k <- 1 + 2 * c(0, cumsum(rho^2)[-m])
  qnorm((1 + nivel) / 2) * sqrt(var_k / N)
}

dibujar_correlograma <- function(valor, banda, etiqueta, titulo) {
  m     <- length(valor)
  rango <- range(c(valor, banda, -banda, 0))
  plot(seq_len(m), valor, type = "h", lwd = 2, ylim = rango,
       xlab = "Rezago", ylab = etiqueta, main = titulo)
  abline(h = 0)
  lines(seq_len(m), banda, lty = 2, col = "blue")
  lines(seq_len(m), -banda, lty = 2, col = "blue")
}

# ACF (bandas de Bartlett) y PACF (bandas de ruido blanco) muestrales
acf_pacf_muestral <- function(x, m = 24, nivel = 0.95, titulo = "") {
  x   <- as.numeric(x)
  N   <- length(x)
  rho <- stats::acf(x, lag.max = m, plot = FALSE)$acf[-1]
  phi <- stats::pacf(x, lag.max = m, plot = FALSE)$acf[, 1, 1]
  
  banda_acf  <- bandas_bartlett(rho, N, nivel)
  banda_pacf <- rep(qnorm((1 + nivel) / 2) / sqrt(N), m)
  
  op <- par(mfrow = c(1, 2))
  on.exit(par(op))
  dibujar_correlograma(rho, banda_acf,  "ACF",  paste("ACF",  titulo))
  dibujar_correlograma(phi, banda_pacf, "PACF", paste("PACF", titulo))
  
  invisible(list(acf = rho, pacf = phi, banda_acf = banda_acf,
                 banda_pacf = banda_pacf, N = N))
}