# ==============================================================================
# EBEC — VALIDACIÓN METODOLÓGICA COMPLETA
# Excess Brier Error Concentration
#
# Objetivo:
#   1. Validar EBEC bajo H0
#   2. Evaluar error tipo I
#   3. Evaluar cobertura del IC
#   4. Evaluar potencia bajo H1
#   5. Analizar sensibilidad a n, prevalencia y q
#   6. Comparar EBEC con Brier score
#
# SIN PAQUETES EXTERNOS
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACIÓN
# ==============================================================================

set.seed(20261002)

DIR_SALIDA <- "EBEC_VALIDACION"

if (!dir.exists(DIR_SALIDA)) {
  dir.create(DIR_SALIDA)
}

# Número de simulaciones
N_SIM_H0 <- 2000
N_SIM_H1 <- 1000

# Para el cálculo interno de cada p-value Monte Carlo
N_MC <- 1000

# Bootstrap
N_BOOT <- 1000

# Nivel alfa
ALPHA <- 0.05


# ==============================================================================
# 1. FUNCIONES BÁSICAS
# ==============================================================================


# ------------------------------------------------------------------------------
# 1.1. Concentración del error Brier
# ------------------------------------------------------------------------------

ebec_concentration <- function(y, p, q = 0.05) {

  stopifnot(length(y) == length(p))

  n <- length(y)

  if (any(!y %in% c(0, 1))) {
    stop("y debe contener únicamente 0 y 1.")
  }

  if (any(!is.finite(p)) || any(p <= 0) || any(p >= 1)) {
    stop("p debe contener probabilidades estrictamente entre 0 y 1.")
  }

  if (q <= 0 || q >= 1) {
    stop("q debe estar entre 0 y 1.")
  }

  loss <- (y - p)^2

  k <- ceiling(q * n)

  if (sum(loss) == 0) {
    return(0)
  }

  loss_sorted <- sort(loss, decreasing = TRUE)

  sum(loss_sorted[seq_len(k)]) / sum(loss)
}


# ------------------------------------------------------------------------------
# 1.2. EBEC
# ------------------------------------------------------------------------------

calculate_ebec <- function(
    y,
    p,
    q = 0.05,
    B = 1000,
    alpha = 0.05,
    bootstrap = TRUE,
    seed = NULL
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n <- length(y)

  # --------------------------------------------------------------------------
  # Observed concentration
  # --------------------------------------------------------------------------

  observed <- ebec_concentration(
    y = y,
    p = p,
    q = q
  )

  # --------------------------------------------------------------------------
  # H0
  #
  # y ~ Bernoulli(p)
  #
  # p permanece fijo.
  # --------------------------------------------------------------------------

  null_concentration <- numeric(B)

  for (b in seq_len(B)) {

    y_null <- rbinom(
      n = n,
      size = 1,
      prob = p
    )

    null_concentration[b] <- ebec_concentration(
      y = y_null,
      p = p,
      q = q
    )
  }

  # --------------------------------------------------------------------------
  # Expected concentration
  # --------------------------------------------------------------------------

  expected <- mean(null_concentration)

  # --------------------------------------------------------------------------
  # EBEC
  # --------------------------------------------------------------------------

  excess <- observed - expected

  # --------------------------------------------------------------------------
  # Normalized EBEC
  #
  # EBEC_norm = EBEC / (1-q)
  # --------------------------------------------------------------------------

  normalized <- excess / (1 - q)

  # --------------------------------------------------------------------------
  # Monte Carlo p-value
  #
  # Test unilateral:
  # concentración observada > concentración esperada
  # --------------------------------------------------------------------------

  p_value <- (
    1 + sum(null_concentration >= observed)
  ) / (B + 1)

  # --------------------------------------------------------------------------
  # Bootstrap CI
  #
  # Se mantiene fija la expectativa H0 y se remuestrean pares (y,p).
  # --------------------------------------------------------------------------

  if (bootstrap) {

    boot_excess <- numeric(N_BOOT)

    for (b in seq_len(N_BOOT)) {

      ind <- sample.int(
        n = n,
        size = n,
        replace = TRUE
      )

      boot_concentration <- ebec_concentration(
        y = y[ind],
        p = p[ind],
        q = q
      )

      boot_excess[b] <- boot_concentration - expected
    }

    ci <- quantile(
      boot_excess,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE
    )

  } else {

    ci <- c(NA, NA)
  }

  # --------------------------------------------------------------------------
  # Resultado
  # --------------------------------------------------------------------------

  list(
    n = n,
    q = q,
    k = ceiling(q * n),
    observed = observed,
    expected = expected,
    excess = excess,
    normalized = normalized,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p_value = p_value,
    null_distribution = null_concentration
  )
}


# ==============================================================================
# 2. GENERADOR DE PROBABILIDADES PREDICHAS
# ==============================================================================

generate_predictions <- function(
    n,
    prevalence = 0.50
) {

  # Beta con media aproximadamente igual a la prevalencia.
  #
  # La dispersión se mantiene moderada para representar probabilidades
  # heterogéneas entre individuos.

  concentration <- 8

  alpha <- max(prevalence * concentration, 0.1)
  beta  <- max((1 - prevalence) * concentration, 0.1)

  p <- rbeta(
    n = n,
    shape1 = alpha,
    shape2 = beta
  )

  # Evitamos probabilidades exactamente 0/1.
  p <- pmin(
    pmax(p, 0.001),
    0.999
  )

  p
}


# ==============================================================================
# 3. SIMULACIÓN BAJO H0
# ==============================================================================

simulate_h0 <- function(
    n,
    prevalence,
    q,
    n_sim = 2000,
    B = 1000,
    alpha = 0.05
) {

  results <- vector(
    "list",
    n_sim
  )

  for (s in seq_len(n_sim)) {

    # --------------------------------------------------------------------------
    # Predicciones
    # --------------------------------------------------------------------------

    p <- generate_predictions(
      n = n,
      prevalence = prevalence
    )

    # --------------------------------------------------------------------------
    # H0 verdadera:
    #
    # y ~ Bernoulli(p)
    # --------------------------------------------------------------------------

    y <- rbinom(
      n = n,
      size = 1,
      prob = p
    )

    # --------------------------------------------------------------------------
    # EBEC
    # --------------------------------------------------------------------------

    ebec <- calculate_ebec(
      y = y,
      p = p,
      q = q,
      B = B,
      alpha = alpha,
      bootstrap = FALSE
    )

    # --------------------------------------------------------------------------
    # Brier
    # --------------------------------------------------------------------------

    brier <- mean(
      (y - p)^2
    )

    results[[s]] <- data.frame(
      simulation = s,
      n = n,
      prevalence = prevalence,
      q = q,
      observed = ebec$observed,
      expected = ebec$expected,
      ebec = ebec$excess,
      ebec_normalized = ebec$normalized,
      p_value = ebec$p_value,
      brier = brier,
      significant = ebec$p_value < alpha
    )
  }

  do.call(
    rbind,
    results
  )
}


# ==============================================================================
# 4. RESUMEN BAJO H0
# ==============================================================================

summarize_h0 <- function(
    data,
    alpha = 0.05
) {

  data.frame(

    n = unique(data$n),

    prevalence = unique(data$prevalence),

    q = unique(data$q),

    mean_EBEC = mean(data$ebec),

    sd_EBEC = sd(data$ebec),

    median_EBEC = median(data$ebec),

    q025_EBEC = quantile(
      data$ebec,
      0.025
    ),

    q975_EBEC = quantile(
      data$ebec,
      0.975
    ),

    type_I_error = mean(
      data$p_value < alpha
    ),

    mean_Brier = mean(
      data$brier
    ),

    sd_Brier = sd(
      data$brier
    )
  )
}


# ==============================================================================
# 5. SIMULACIÓN BAJO H1
#
# Introducimos concentración deliberada del error.
#
# Una proporción contamination_fraction de observaciones recibe:
#
#     p_true = (1-strength)*p + strength*(1-p)
#
# strength = 0  -> H0
# strength = 1  -> probabilidad completamente invertida
#
# Esto genera deliberadamente peor predicción en un subconjunto.
# ==============================================================================

simulate_h1 <- function(
    n,
    prevalence,
    q,
    contamination_fraction,
    strength,
    n_sim = 1000,
    B = 1000,
    alpha = 0.05
) {

  results <- vector(
    "list",
    n_sim
  )

  for (s in seq_len(n_sim)) {

    # --------------------------------------------------------------------------
    # Predicciones
    # --------------------------------------------------------------------------

    p <- generate_predictions(
      n = n,
      prevalence = prevalence
    )

    # --------------------------------------------------------------------------
    # Selección del subconjunto contaminado
    # --------------------------------------------------------------------------

    m <- max(
      1,
      round(contamination_fraction * n)
    )

    contaminated <- sample.int(
      n = n,
      size = m,
      replace = FALSE
    )

    # --------------------------------------------------------------------------
    # Probabilidad verdadera
    # --------------------------------------------------------------------------

    p_true <- p

    p_true[contaminated] <-
      (1 - strength) * p[contaminated] +
      strength * (1 - p[contaminated])

    p_true <- pmin(
      pmax(p_true, 0.001),
      0.999
    )

    # --------------------------------------------------------------------------
    # Outcome
    # --------------------------------------------------------------------------

    y <- rbinom(
      n = n,
      size = 1,
      prob = p_true
    )

    # --------------------------------------------------------------------------
    # EBEC
    # --------------------------------------------------------------------------

    ebec <- calculate_ebec(
      y = y,
      p = p,
      q = q,
      B = B,
      alpha = alpha,
      bootstrap = FALSE
    )

    # --------------------------------------------------------------------------
    # Brier
    # --------------------------------------------------------------------------

    brier <- mean(
      (y - p)^2
    )

    results[[s]] <- data.frame(
      simulation = s,
      n = n,
      prevalence = prevalence,
      q = q,
      contamination_fraction = contamination_fraction,
      strength = strength,
      observed = ebec$observed,
      expected = ebec$expected,
      ebec = ebec$excess,
      ebec_normalized = ebec$normalized,
      p_value = ebec$p_value,
      brier = brier,
      significant = ebec$p_value < alpha
    )
  }

  do.call(
    rbind,
    results
  )
}


# ==============================================================================
# 6. RESUMEN DE POTENCIA
# ==============================================================================

summarize_power <- function(
    data,
    alpha = 0.05
) {

  data.frame(

    n = unique(data$n),

    prevalence = unique(data$prevalence),

    q = unique(data$q),

    contamination_fraction =
      unique(data$contamination_fraction),

    strength =
      unique(data$strength),

    mean_EBEC =
      mean(data$ebec),

    median_EBEC =
      median(data$ebec),

    q025_EBEC =
      quantile(data$ebec, 0.025),

    q975_EBEC =
      quantile(data$ebec, 0.975),

    power =
      mean(data$p_value < alpha),

    mean_Brier =
      mean(data$brier),

    sd_Brier =
      sd(data$brier)
  )
}


# ==============================================================================
# 7. VALIDACIÓN H0
#
# Diseño principal
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("VALIDACIÓN H0\n")
cat("============================================================\n")

H0_RESULTS <- list()
H0_SUMMARY <- list()

counter <- 1

for (n in c(100, 250, 500, 1000)) {

  for (prevalence in c(0.10, 0.30, 0.50)) {

    for (q in c(0.01, 0.05, 0.10, 0.20)) {

      cat(
        "H0:",
        "n =", n,
        "| prevalence =", prevalence,
        "| q =", q,
        "\n"
      )

      tmp <- simulate_h0(
        n = n,
        prevalence = prevalence,
        q = q,
        n_sim = N_SIM_H0,
        B = N_MC,
        alpha = ALPHA
      )

      H0_RESULTS[[counter]] <- tmp

      H0_SUMMARY[[counter]] <- summarize_h0(
        tmp,
        alpha = ALPHA
      )

      counter <- counter + 1
    }
  }
}

H0_RESULTS <- do.call(
  rbind,
  H0_RESULTS
)

H0_SUMMARY <- do.call(
  rbind,
  H0_SUMMARY
)


# ==============================================================================
# 8. GUARDAR H0
# ==============================================================================

write.csv(
  H0_RESULTS,
  file.path(
    DIR_SALIDA,
    "EBEC_H0_all_simulations.csv"
  ),
  row.names = FALSE
)

write.csv(
  H0_SUMMARY,
  file.path(
    DIR_SALIDA,
    "EBEC_H0_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. SIMULACIÓN H1 — POTENCIA
#
# Concentración dentro de una fracción del dataset.
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("VALIDACIÓN H1 / POTENCIA\n")
cat("============================================================\n")

H1_RESULTS <- list()
H1_SUMMARY <- list()

counter <- 1

for (n in c(100, 250, 500, 1000)) {

  for (prevalence in c(0.10, 0.30, 0.50)) {

    for (q in c(0.01, 0.05, 0.10, 0.20)) {

      for (contamination_fraction in c(0.01, 0.05, 0.10)) {

        for (strength in c(0.25, 0.50, 0.75, 1.00)) {

          cat(
            "H1:",
            "n =", n,
            "| prev =", prevalence,
            "| q =", q,
            "| contamination =", contamination_fraction,
            "| strength =", strength,
            "\n"
          )

          tmp <- simulate_h1(
            n = n,
            prevalence = prevalence,
            q = q,
            contamination_fraction = contamination_fraction,
            strength = strength,
            n_sim = N_SIM_H1,
            B = N_MC,
            alpha = ALPHA
          )

          H1_RESULTS[[counter]] <- tmp

          H1_SUMMARY[[counter]] <- summarize_power(
            tmp,
            alpha = ALPHA
          )

          counter <- counter + 1
        }
      }
    }
  }
}

H1_RESULTS <- do.call(
  rbind,
  H1_RESULTS
)

H1_SUMMARY <- do.call(
  rbind,
  H1_SUMMARY
)


# ==============================================================================
# 10. GUARDAR H1
# ==============================================================================

write.csv(
  H1_RESULTS,
  file.path(
    DIR_SALIDA,
    "EBEC_H1_all_simulations.csv"
  ),
  row.names = FALSE
)

write.csv(
  H1_SUMMARY,
  file.path(
    DIR_SALIDA,
    "EBEC_H1_power_summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 11. ERROR TIPO I GLOBAL
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("ERROR TIPO I\n")
cat("============================================================\n")

TYPE_I_SUMMARY <- aggregate(
  type_I_error ~ n + prevalence + q,
  data = H0_SUMMARY,
  FUN = mean
)

print(
  TYPE_I_SUMMARY
)

write.csv(
  TYPE_I_SUMMARY,
  file.path(
    DIR_SALIDA,
    "EBEC_type_I_error.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 12. GRÁFICO — ERROR TIPO I
# ==============================================================================

png(
  file.path(
    DIR_SALIDA,
    "EBEC_type_I_error.png"
  ),
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  TYPE_I_SUMMARY$n,
  TYPE_I_SUMMARY$type_I_error,
  type = "n",
  ylim = c(0, 0.15),
  xlab = "Sample size",
  ylab = "Empirical type I error",
  main = "EBEC: empirical type I error"
)

abline(
  h = ALPHA,
  lty = 2
)

points(
  TYPE_I_SUMMARY$n,
  TYPE_I_SUMMARY$type_I_error,
  pch = 19
)

dev.off()


# ==============================================================================
# 13. POTENCIA — RESUMEN
# ==============================================================================

POWER_TABLE <- H1_SUMMARY[
  order(
    H1_SUMMARY$n,
    H1_SUMMARY$q,
    H1_SUMMARY$contamination_fraction,
    H1_SUMMARY$strength
  ),
]

write.csv(
  POWER_TABLE,
  file.path(
    DIR_SALIDA,
    "EBEC_power_table.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 14. GRÁFICO DE POTENCIA
#
# Escenario central:
# n = 500
# prevalence = .50
# q = .05
# ==============================================================================

POWER_CENTRAL <- subset(
  H1_SUMMARY,
  n == 500 &
    prevalence == 0.50 &
    q == 0.05
)

png(
  file.path(
    DIR_SALIDA,
    "EBEC_power_central.png"
  ),
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  POWER_CENTRAL$strength[
    POWER_CENTRAL$contamination_fraction == 0.05
  ],
  POWER_CENTRAL$power[
    POWER_CENTRAL$contamination_fraction == 0.05
  ],
  type = "b",
  pch = 19,
  ylim = c(0, 1),
  xlab = "Strength of error concentration",
  ylab = "Power",
  main = "EBEC power: n = 500, prevalence = .50, q = .05"
)

abline(
  h = ALPHA,
  lty = 2
)

dev.off()


# ==============================================================================
# 15. RELACIÓN ENTRE EBEC Y BRIER
#
# Comprobamos que EBEC y Brier no son simplemente la misma información.
# ==============================================================================

CORRELATION_DATA <- H1_RESULTS[
  H1_RESULTS$n == 500 &
    H1_RESULTS$prevalence == 0.50 &
    H1_RESULTS$q == 0.05,
]

correlation_ebec_brier <- cor(
  CORRELATION_DATA$ebec,
  CORRELATION_DATA$brier,
  method = "spearman"
)

cat("\n")
cat("============================================================\n")
cat("CORRELACIÓN EBEC - BRIER\n")
cat("============================================================\n")

cat(
  "Spearman rho =",
  round(correlation_ebec_brier, 4),
  "\n"
)


# ==============================================================================
# 16. GRÁFICO EBEC VS BRIER
# ==============================================================================

png(
  file.path(
    DIR_SALIDA,
    "EBEC_vs_Brier.png"
  ),
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  CORRELATION_DATA$brier,
  CORRELATION_DATA$ebec,
  pch = 19,
  xlab = "Brier score",
  ylab = "EBEC",
  main = "EBEC vs. Brier score"
)

abline(
  h = 0,
  lty = 2
)

dev.off()


# ==============================================================================
# 17. ANÁLISIS DE SENSIBILIDAD DE q
# ==============================================================================

SENS_Q <- subset(
  H1_SUMMARY,
  n == 500 &
    prevalence == 0.50 &
    contamination_fraction == 0.05 &
    strength == 0.75
)

write.csv(
  SENS_Q,
  file.path(
    DIR_SALIDA,
    "EBEC_sensitivity_q.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 18. RESUMEN FINAL AUTOMÁTICO
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("RESUMEN EBEC\n")
cat("============================================================\n")

cat("\n")
cat("H0:\n")
cat(
  "Rango de error tipo I observado:",
  round(
    min(H0_SUMMARY$type_I_error),
    4
  ),
  "-",
  round(
    max(H0_SUMMARY$type_I_error),
    4
  ),
  "\n"
)

cat(
  "Error tipo I medio:",
  round(
    mean(H0_SUMMARY$type_I_error),
    4
  ),
  "\n"
)

cat("\n")
cat("H1:\n")
cat(
  "Potencia mínima:",
  round(
    min(H1_SUMMARY$power),
    4
  ),
  "\n"
)

cat(
  "Potencia máxima:",
  round(
    max(H1_SUMMARY$power),
    4
  ),
  "\n"
)

cat("\n")
cat(
  "Correlación Spearman EBEC-Brier:",
  round(
    correlation_ebec_brier,
    4
  ),
  "\n"
)

cat("\n")
cat(
  "Todos los resultados guardados en:",
  DIR_SALIDA,
  "\n"
)

cat("\n")
cat("FIN\n")


# ==============================================================================
# SELECCIÓN EXPERIMENTAL DEL q ÓPTIMO PARA EBEC
# ==============================================================================
#
# Objetivo:
# Determinar qué proporción q de las observaciones con mayor pérdida Brier
# ofrece el mejor comportamiento global de EBEC.
#
# Criterio:
# 1. Bajo H0: error tipo I próximo a alpha = .05.
# 2. Bajo H1: alta potencia.
# 3. El rendimiento debe ser estable a través de n, prevalencia,
#    proporción contaminada e intensidad de contaminación.
#
# IMPORTANTE:
# q NO se selecciona únicamente por tener la mayor potencia.
# Un q que infle el error tipo I no es aceptable aunque tenga mucha potencia.
# ==============================================================================

set.seed(20261002)

# ------------------------------------------------------------------------------
# 1. CONFIGURACIÓN
# ------------------------------------------------------------------------------

Q_VALUES <- c(
  0.005,
  0.010,
  0.020,
  0.030,
  0.050,
  0.075,
  0.100,
  0.150,
  0.200,
  0.250,
  0.300
)

N_H0_Q <- 1000
N_H1_Q <- 500
B_Q    <- 500

ALPHA_Q <- 0.05

# Escenarios H0
N_VALUES_Q <- c(100, 250, 500, 1000)
PREV_VALUES_Q <- c(0.10, 0.30, 0.50)

# Escenarios H1
CONTAMINATION_VALUES_Q <- c(0.05, 0.10, 0.20)
STRENGTH_VALUES_Q <- c(0.25, 0.50, 0.75, 1.00)


# ==============================================================================
# 2. FUNCIÓN PARA CALCULAR EBEC SIN BOOTSTRAP
# ==============================================================================

calculate_ebec_q <- function(y, p, q, B = 500, alpha = 0.05) {

  observed <- ebec_concentration(y, p, q)

  n <- length(y)

  null_concentration <- numeric(B)

  for (b in seq_len(B)) {

    y_null <- rbinom(
      n = n,
      size = 1,
      prob = p
    )

    null_concentration[b] <- ebec_concentration(
      y_null,
      p,
      q
    )
  }

  expected <- mean(null_concentration)

  excess <- observed - expected

  normalized <- excess / (1 - q)

  p_value <- (
    1 +
      sum(null_concentration >= observed)
  ) / (B + 1)

  significant <- p_value < alpha

  list(
    observed = observed,
    expected = expected,
    excess = excess,
    normalized = normalized,
    p_value = p_value,
    significant = significant
  )
}


# ==============================================================================
# 3. H0
# ==============================================================================

cat("\n============================================\n")
cat("INICIANDO SELECCIÓN DE q — H0\n")
cat("============================================\n\n")

results_H0_q <- list()

counter <- 0

total_H0 <- length(N_VALUES_Q) *
  length(PREV_VALUES_Q) *
  length(Q_VALUES)

for (n in N_VALUES_Q) {

  for (prev in PREV_VALUES_Q) {

    for (q in Q_VALUES) {

      counter <- counter + 1

      cat(
        sprintf(
          "H0 %d/%d | n=%d | prev=%.2f | q=%.3f\n",
          counter,
          total_H0,
          n,
          prev,
          q
        )
      )

      type1 <- numeric(N_H0_Q)
      ebec_values <- numeric(N_H0_Q)

      for (sim in seq_len(N_H0_Q)) {

        # Predicciones
        p <- generate_predictions(
          n = n,
          prevalence = prev
        )

        # Resultado compatible con las predicciones
        y <- rbinom(
          n = n,
          size = 1,
          prob = p
        )

        res <- calculate_ebec_q(
          y = y,
          p = p,
          q = q,
          B = B_Q,
          alpha = ALPHA_Q
        )

        type1[sim] <- res$significant
        ebec_values[sim] <- res$observed
      }

      results_H0_q[[length(results_H0_q) + 1]] <- data.frame(
        n = n,
        prevalence = prev,
        q = q,
        type1_error = mean(type1),
        mean_ebec = mean(ebec_values),
        median_ebec = median(ebec_values),
        sd_ebec = sd(ebec_values)
      )
    }
  }
}

results_H0_q <- do.call(
  rbind,
  results_H0_q
)


# ==============================================================================
# 4. RESUMEN H0 POR q
# ==============================================================================

summary_H0_q <- aggregate(
  type1_error ~ q,
  data = results_H0_q,
  FUN = mean
)

names(summary_H0_q)[2] <- "mean_type1"

summary_H0_q$abs_deviation <- abs(
  summary_H0_q$mean_type1 - ALPHA_Q
)

# Máxima desviación entre escenarios
max_type1 <- aggregate(
  type1_error ~ q,
  data = results_H0_q,
  FUN = max
)

names(max_type1)[2] <- "max_type1"

min_type1 <- aggregate(
  type1_error ~ q,
  data = results_H0_q,
  FUN = min
)

names(min_type1)[2] <- "min_type1"

summary_H0_q <- merge(
  summary_H0_q,
  min_type1,
  by = "q"
)

summary_H0_q <- merge(
  summary_H0_q,
  max_type1,
  by = "q"
)

# Criterio práctico:
# consideramos aceptable un q si el error tipo I medio está dentro
# de ±0.01 de .05.

summary_H0_q$H0_acceptable <- (
  abs(summary_H0_q$mean_type1 - ALPHA_Q) <= 0.01
)


# ==============================================================================
# 5. H1
# ==============================================================================

cat("\n============================================\n")
cat("INICIANDO SELECCIÓN DE q — H1\n")
cat("============================================\n\n")

results_H1_q <- list()

counter <- 0

total_H1 <- length(N_VALUES_Q) *
  length(PREV_VALUES_Q) *
  length(Q_VALUES) *
  length(CONTAMINATION_VALUES_Q) *
  length(STRENGTH_VALUES_Q)

for (n in N_VALUES_Q) {

  for (prev in PREV_VALUES_Q) {

    for (q in Q_VALUES) {

      for (contamination in CONTAMINATION_VALUES_Q) {

        for (strength in STRENGTH_VALUES_Q) {

          counter <- counter + 1

          cat(
            sprintf(
              "H1 %d/%d | n=%d | prev=%.2f | q=%.3f | contamination=%.2f | strength=%.2f\n",
              counter,
              total_H1,
              n,
              prev,
              q,
              contamination,
              strength
            )
          )

          power <- numeric(N_H1_Q)

          for (sim in seq_len(N_H1_Q)) {

            # Predicciones originales
            p <- generate_predictions(
              n = n,
              prevalence = prev
            )

            # Seleccionar observaciones contaminadas
            n_contam <- max(
              1,
              ceiling(n * contamination)
            )

            contaminated <- sample(
              seq_len(n),
              size = n_contam,
              replace = FALSE
            )

            # Invertir parcialmente las predicciones
            p_true <- p

            p_true[contaminated] <-
              (1 - strength) * p[contaminated] +
              strength * (1 - p[contaminated])

            # Generar resultados
            y <- rbinom(
              n = n,
              size = 1,
              prob = p_true
            )

            res <- calculate_ebec_q(
              y = y,
              p = p,
              q = q,
              B = B_Q,
              alpha = ALPHA_Q
            )

            power[sim] <- res$significant
          }

          results_H1_q[[length(results_H1_q) + 1]] <- data.frame(
            n = n,
            prevalence = prev,
            q = q,
            contamination = contamination,
            strength = strength,
            power = mean(power)
          )
        }
      }
    }
  }
}

results_H1_q <- do.call(
  rbind,
  results_H1_q
)


# ==============================================================================
# 6. RESUMEN DE POTENCIA POR q
# ==============================================================================

summary_H1_q <- aggregate(
  power ~ q,
  data = results_H1_q,
  FUN = mean
)

names(summary_H1_q)[2] <- "mean_power"

# Potencia mínima entre escenarios
min_power <- aggregate(
  power ~ q,
  data = results_H1_q,
  FUN = min
)

names(min_power)[2] <- "min_power"

# Potencia máxima
max_power <- aggregate(
  power ~ q,
  data = results_H1_q,
  FUN = max
)

names(max_power)[2] <- "max_power"

# Variabilidad de la potencia
sd_power <- aggregate(
  power ~ q,
  data = results_H1_q,
  FUN = sd
)

names(sd_power)[2] <- "sd_power"

summary_H1_q <- merge(
  summary_H1_q,
  min_power,
  by = "q"
)

summary_H1_q <- merge(
  summary_H1_q,
  max_power,
  by = "q"
)

summary_H1_q <- merge(
  summary_H1_q,
  sd_power,
  by = "q"
)


# ==============================================================================
# 7. UNIR H0 + H1
# ==============================================================================

q_comparison <- merge(
  summary_H0_q,
  summary_H1_q,
  by = "q"
)

# ==============================================================================
# 8. CRITERIO DE SELECCIÓN
# ==============================================================================

# Primero:
#   H0 aceptable
#
# Después:
#   máxima potencia media entre los q aceptables.
#
# En caso de empate práctico:
#   menor desviación del error tipo I.

q_comparison$selection_score <- ifelse(
  q_comparison$H0_acceptable,
  q_comparison$mean_power,
  NA
)

q_comparison <- q_comparison[
  order(
    -q_comparison$selection_score,
    q_comparison$abs_deviation
  ),
]

q_recommended <- q_comparison[
  !is.na(q_comparison$selection_score),
][1, ]


# ==============================================================================
# 9. RESULTADOS
# ==============================================================================

cat("\n\n============================================\n")
cat("RESULTADOS — SELECCIÓN DE q\n")
cat("============================================\n\n")

print(
  q_comparison,
  row.names = FALSE
)

cat("\n--------------------------------------------\n")
cat("q RECOMENDADO SEGÚN ESTE CRITERIO\n")
cat("--------------------------------------------\n\n")

print(
  q_recommended,
  row.names = FALSE
)


# ==============================================================================
# 10. TABLA FINAL SIMPLIFICADA
# ==============================================================================

q_final <- q_comparison[, c(
  "q",
  "mean_type1",
  "min_type1",
  "max_type1",
  "H0_acceptable",
  "mean_power",
  "min_power",
  "sd_power"
)]

q_final <- q_final[
  order(q_final$q),
]

print(
  q_final,
  row.names = FALSE
)


# ==============================================================================
# 11. GRÁFICO: ERROR TIPO I Y POTENCIA
# ==============================================================================

par(mfrow = c(1, 2))

plot(
  q_final$q,
  q_final$mean_type1,
  type = "b",
  pch = 19,
  xlab = "q",
  ylab = "Mean Type I error",
  main = "Type I error under H0"
)

abline(
  h = ALPHA_Q,
  lty = 2
)

abline(
  h = ALPHA_Q + 0.01,
  lty = 3
)

abline(
  h = ALPHA_Q - 0.01,
  lty = 3
)


plot(
  q_final$q,
  q_final$mean_power,
  type = "b",
  pch = 19,
  xlab = "q",
  ylab = "Mean power",
  main = "Power under H1"
)

par(mfrow = c(1, 1))


# ==============================================================================
# 12. GUARDAR RESULTADOS
# ==============================================================================

write.csv(
  results_H0_q,
  "EBEC_q_selection_H0_all.csv",
  row.names = FALSE
)

write.csv(
  results_H1_q,
  "EBEC_q_selection_H1_all.csv",
  row.names = FALSE
)

write.csv(
  q_comparison,
  "EBEC_q_selection_comparison.csv",
  row.names = FALSE
)

write.csv(
  q_final,
  "EBEC_q_selection_final.csv",
  row.names = FALSE
)

cat("\nArchivos guardados correctamente.\n")
