#' Excess Brier Error Concentration (EBEC)
#'
#' Tests whether the worst predictions of a binary model accumulate more
#' squared error (Brier loss) than expected if the predicted probabilities
#' were true. Unlike the mean Brier score, the test is not diluted when only a
#' few predictions are badly wrong.
#'
#' @details
#' \strong{Components.} The test looks at several sums of losses at once: the
#' total Brier loss; the loss accumulated by the \code{k} worst cases, for
#' several tail sizes (\code{k = ceiling(q * n)}, by default the worst 1 percent,
#' 5, 10 and 20 percent of the cases); and, by default (\code{classes = "fp"}),
#' the same tail sums for false positives (non-events, loss \code{p^2}). With
#' \code{classes = "both"} false negatives (events, loss \code{(1 - p)^2}) are
#' added as an exploratory option.
#'
#' \strong{Statistic and p-value.} Each component is compared with its
#' expected value if the predicted probabilities were true (obtained by
#' simulation) and standardized (z). The statistic is the largest z. Its
#' Monte Carlo p-value (the proportion of simulated statistics at least as
#' large as the observed one, counting the observed one) already accounts for
#' the choice among components, so it is the only p-value to report. The statistic has no fixed range (it can be any real number) and the
#' test is one-sided: only an excess of loss gives evidence.
#'
#' \strong{Effect size.} For every component the table gives
#' \code{ratio = observed / expected} (1 = as expected) and
#' \code{excess_pct = 100 * (ratio - 1)}, with a bootstrap interval
#' (\code{ci = TRUE}). In a small simulation (n = 400, about 8 percent events) the
#' coverage of the 95 percent interval was about 88-93 percent. These intervals are descriptive: they are not adjusted
#' for looking at several components and the largest excess is a selected
#' maximum, so it tends to be optimistic.
#'
#' \strong{Validity.} The p-value is valid when the predicted probabilities are
#' fixed in advance (held-out or cross-validated), and the outcomes are
#' independent. With predictions from a model fitted on the same data the test
#' is conservative. Type I error was checked at 3-50 percent prevalence; the power
#' advantage over mean-based calibration tests was shown at 10-50 percent prevalence
#' and, at 3 percent, only for false-positive and mixed faults. A false-negative-specific property has not been
#' demonstrated; the default components were chosen from simulations in which
#' the false-negative components did not add power.
#'
#' Printing the result shows the table, how to read it, a sentence ready to
#' report, and a short guide (\code{print(x, explain = FALSE)} hides the guide;
#' \code{print(x, lang = "es")} prints it in Spanish).
#'
#' @param y Binary outcome vector (0/1).
#' @param p Predicted probabilities, between 0 and 1.
#' @param q Fractions of the sample used as tail sizes.
#'   Default \code{c(0.01, 0.05, 0.10, 0.20)}.
#' @param B Number of Monte Carlo draws for the null distribution. Default 999.
#' @param classes Which class-specific components to include: \code{"fp"}
#'   (default; false positives), \code{"none"} (total and tail sums only) or
#'   \code{"both"} (false positives and false negatives). \code{TRUE} means
#'   \code{"both"} and \code{FALSE} means \code{"none"}.
#' @param ci Logical. Compute bootstrap intervals for the percentage excess?
#'   Default \code{TRUE}; use \code{FALSE} to skip it (much faster).
#' @param B_ci Number of bootstrap resamples for the intervals. Default 300.
#' @param B_inner Number of simulated outcome vectors used inside each
#'   bootstrap resample to obtain the expected value. Default 50.
#' @param level Confidence level of the intervals. Default 0.95.
#' @param seed Optional random seed. If \code{NULL} (default) the current
#'   random number state is used.
#'
#' @return An object of class \code{ebec}: a list with the Monte Carlo p-value
#'   (\code{p_value}), the statistic (\code{statistic}, the largest z), the
#'   component that produced it (\code{driver}), and a data frame
#'   \code{components} with, for each component, the share of the sample in
#'   the tail (\code{tail_pct}), \code{observed}, \code{expected}, \code{ratio},
#'   \code{excess_pct}, its interval (\code{ci_low}, \code{ci_high}), \code{z}
#'   and an unadjusted \code{p_component} (descriptive only).
#' @importFrom stats rbinom sd quantile
#' @export
#'
#' @examples
#' set.seed(1)
#' p <- runif(200, 0.2, 0.8)
#' y <- rbinom(200, 1, p)
#' ebec(y, p, B = 199, ci = FALSE)
ebec <- function(y, p, q = c(0.01, 0.05, 0.10, 0.20), B = 999,
                 classes = c("fp", "none", "both"), ci = TRUE,
                 B_ci = 300, B_inner = 50, level = 0.95, seed = NULL) {
  classes <- .scan_mode(classes)
  .scan_check(y, p, q, B)
  stopifnot(is.numeric(level), length(level) == 1, level > 0, level < 1,
            B_ci >= 1, B_inner >= 1)
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  ks <- .scan_ks(q, n)
  obs <- .scan_stats(y, p, ks, classes)
  null <- .scan_null(n, p, ks, classes, B)
  res <- .scan_result(obs, null, n, B, classes, ks)
  if (isTRUE(ci)) {
    bounds <- .scan_ci(y, p, ks, classes, B_ci, B_inner, level)
    res$components$ci_low <- bounds[, 1]
    res$components$ci_high <- bounds[, 2]
  } else {
    res$components$ci_low <- NA_real_
    res$components$ci_high <- NA_real_
  }
  res$level <- level
  res$n_events <- sum(y)
  res$ci <- isTRUE(ci)
  res$components <- res$components[, c("component", "tail_pct", "observed", "expected",
                                       "ratio", "excess_pct", "ci_low", "ci_high",
                                       "z", "p_component")]
  res
}

#' @param x An object of class \code{ebec}.
#' @param digits Number of digits.
#' @param explain Print the short guide (what EBEC is, how to read and report
#'   it)? Default \code{TRUE}.
#' @param lang Language of the explanations: \code{"en"} (default) or
#'   \code{"es"}.
#' @param ... Unused.
#' @rdname ebec
#' @export
print.ebec <- function(x, digits = 3, explain = TRUE, lang = c("en", "es"), ...) {
  lang <- match.arg(lang)
  es <- lang == "es"
  tx <- function(en, sp) if (es) sp else en
  cat(tx("EBEC - Excess Brier Error Concentration\n",
         "EBEC - Concentraci\u00f3n del Exceso de Error de Brier\n"))
  cat(sprintf(tx("n = %d cases, %d events, %d Monte Carlo draws\n",
                 "n = %d casos, %d eventos, %d simulaciones Monte Carlo\n"),
              x$n, x$n_events, x$B))
  cat("\n")
  df <- x$components
  df$tail_pct <- ifelse(is.na(df$tail_pct), "all", sprintf("%.1f%%", df$tail_pct))
  num <- c("observed", "expected", "ratio", "excess_pct", "ci_low", "ci_high", "z", "p_component")
  df[, num] <- lapply(df[, num], function(v) round(v, digits))
  names(df)[names(df) == "tail_pct"] <- tx("cases", "casos")
  print(df, row.names = FALSE)
  cat("\n")
  cat(sprintf(tx("EBEC statistic (largest z) = %.2f   driver: %s\n",
                 "Estad\u00edstico EBEC (z m\u00e1xima) = %.2f   componente: %s\n"),
              x$statistic, x$driver))
  cat(sprintf(tx("Monte Carlo p-value (one-sided, adjusted for the choice of component) = %.4f\n",
                 "p-valor Monte Carlo (unilateral, ajustado por la elecci\u00f3n de componente) = %.4f\n"),
              x$p_value))
  cat("\n")
  cat(.scan_reading(x, lang), "\n", sep = "")
  cat("\n")
  cat(tx("How to report:\n", "C\u00f3mo reportarlo:\n"))
  cat("  ", .scan_report(x, lang), "\n", sep = "")
  if (isTRUE(explain)) {
    cat("\n")
    cat(.scan_guide(lang), sep = "\n")
  }
  invisible(x)
}

# ---- internal helpers -------------------------------------------------------

.scan_check <- function(y, p, q, B) {
  stopifnot(!anyNA(y), !anyNA(p), length(y) == length(p), length(y) > 0,
            all(y %in% c(0, 1)), all(p >= 0 & p <= 1),
            length(q) >= 1, all(q > 0 & q < 1), length(B) == 1, B >= 1)
  invisible(TRUE)
}

.scan_mode <- function(classes) {
  if (is.logical(classes) && length(classes) == 1L && !is.na(classes))
    return(if (classes) "both" else "none")
  match.arg(classes, c("fp", "none", "both"))
}

.scan_ks <- function(q, n) {
  unique(pmin(n, pmax(1L, as.integer(ceiling(q * n)))))
}

# Observed (or simulated) component vector for outcome y and probabilities p.
# f = loss of false negatives (events), g = loss of false positives (non-events).
.scan_stats <- function(y, p, ks, classes) {
  f <- y * (1 - p)^2
  g <- (1 - y) * p^2
  l <- f + g
  classes <- .scan_mode(classes)
  tops <- function(v) cumsum(sort(v, decreasing = TRUE))[ks]
  out <- c(sum(l), tops(l))
  nm <- c("total (Brier sum)", paste0("worst ", ks, " cases"))
  if (classes == "both") {
    out <- c(out, tops(f), tops(g))
    nm <- c(nm, paste0("FN worst ", ks), paste0("FP worst ", ks))
  } else if (classes == "fp") {
    out <- c(out, tops(g))
    nm <- c(nm, paste0("FP worst ", ks))
  }
  names(out) <- nm
  out
}

.scan_null <- function(n, p, ks, classes, B) {
  nc <- length(.scan_stats(rep(0, n), p, ks, classes))
  out <- matrix(NA_real_, nrow = B, ncol = nc)
  for (b in seq_len(B)) out[b, ] <- .scan_stats(rbinom(n, 1, p), p, ks, classes)
  colnames(out) <- names(.scan_stats(rep(0, n), p, ks, classes))
  out
}

# Max-z statistic and Monte Carlo p-value. The observed vector is pooled with
# the null draws, which keeps the B + 1 values exchangeable under H0.
.scan_p <- function(obs, null, cols = seq_along(obs)) {
  allv <- rbind(obs[cols], null[, cols, drop = FALSE])
  mu <- colMeans(allv)
  sdv <- apply(allv, 2, sd)
  keep <- is.finite(sdv) & sdv > 0
  Z <- sweep(sweep(allv, 2, mu, "-"), 2, sdv, "/")
  Z[, !keep] <- -Inf
  M <- apply(Z, 1, max)
  list(M = M, Z = Z, mu = mu, p_value = mean(M >= M[1]), keep = keep)
}

.scan_result <- function(obs, null, n, B, classes, ks = NULL) {
  r <- .scan_p(obs, null)
  zobs <- r$Z[1, ]
  nm <- names(obs)
  tail_k <- suppressWarnings(as.numeric(sub("^(FN |FP )?(worst )?([0-9]+)( cases)?$", "\\3", nm)))
  tail_k[nm == "total (Brier sum)"] <- NA
  comp <- data.frame(
    component = nm,
    tail_pct = 100 * tail_k / n,
    observed = as.numeric(obs),
    expected = colMeans(null),
    stringsAsFactors = FALSE
  )
  comp$ratio <- ifelse(comp$expected > 0, comp$observed / comp$expected, NA_real_)
  comp$excess_pct <- 100 * (comp$ratio - 1)
  comp$z <- as.numeric(zobs)
  comp$p_component <- vapply(seq_along(obs), function(j) {
    (1 + sum(null[, j] >= obs[j])) / (nrow(null) + 1)
  }, numeric(1))
  structure(list(p_value = r$p_value, statistic = r$M[1],
                 driver = comp$component[which.max(zobs)],
                 components = comp, n = n, B = B, classes = classes),
            class = "ebec")
}

# Bootstrap interval for the percentage excess of every component. Cases
# (y, p) are resampled with replacement; for each resample the expected value
# is re-estimated by simulating outcomes from the resampled probabilities.
.scan_ci <- function(y, p, ks, classes, B_ci, B_inner, level) {
  n <- length(y)
  nc <- length(.scan_stats(y, p, ks, classes))
  R <- matrix(NA_real_, B_ci, nc)
  for (b in seq_len(B_ci)) {
    i <- sample.int(n, n, replace = TRUE)
    yb <- y[i]; pb <- p[i]
    ob <- .scan_stats(yb, pb, ks, classes)
    ex <- rowMeans(vapply(seq_len(B_inner), function(j)
      .scan_stats(rbinom(n, 1, pb), pb, ks, classes), numeric(nc)))
    R[b, ] <- ifelse(ex > 0, 100 * (ob / ex - 1), NA_real_)
  }
  a <- (1 - level) / 2
  t(apply(R, 2, quantile, probs = c(a, 1 - a), na.rm = TRUE, names = FALSE))
}

.scan_fmt_pct <- function(x) sprintf("%+.1f%%", x)

# Plain-language reading of the result.
.scan_reading <- function(x, lang) {
  es <- lang == "es"
  sig <- x$p_value < 0.05
  d <- x$components[x$components$component == x$driver, ]
  if (sig) {
    if (es) sprintf("Lectura: p < 0,05. Hay m\u00e1s p\u00e9rdida en los peores casos de la esperable si las probabilidades fueran ciertas. El mayor exceso est\u00e1 en '%s' (raz\u00f3n %.2f).", x$driver, d$ratio)
    else sprintf("Reading: p < 0.05. The worst cases accumulate more loss than expected if the probabilities were true. The largest excess is in '%s' (ratio %.2f).", x$driver, d$ratio)
  } else {
    if (es) "Lectura: p >= 0,05. No hay evidencia de que los peores casos acumulen m\u00e1s p\u00e9rdida de la esperada. No demuestra que el modelo sea correcto, sobre todo con pocos datos."
    else "Reading: p >= 0.05. No evidence that the worst cases accumulate more loss than expected. This does not prove the model is right, especially with few data."
  }
}

# One sentence ready to report.
.scan_report <- function(x, lang) {
  es <- lang == "es"
  d <- x$components[x$components$component == x$driver, ]
  pv <- if (x$p_value < 0.001) "< 0.001" else sprintf("= %.3f", x$p_value)
  if (es) pv <- sub("\\.", ",", pv)
  where <- if (is.na(d$tail_pct)) {
    if (es) "en la p\u00e9rdida total" else "in the total loss"
  } else {
    sprintf(if (es) "en '%s' (%.0f%% de los casos)" else "in '%s' (%.0f%% of cases)",
            d$component, d$tail_pct)
  }
  cin <- if (isTRUE(x$ci) && is.finite(d$ci_low))
    sprintf(if (es) ", IC %.0f%% %.1f%% a %.1f%%" else ", %.0f%% CI %.1f%% to %.1f%%",
            100 * x$level, d$ci_low, d$ci_high) else ""
  if (es) sprintf("EBEC: p %s. Mayor exceso %s: p\u00e9rdida observada %.2f frente a %.2f esperada (%s%s).",
                  pv, where, d$observed, d$expected, .scan_fmt_pct(d$excess_pct), cin)
  else sprintf("EBEC: p %s. Largest excess %s: observed loss %.2f vs %.2f expected (%s%s).",
               pv, where, d$observed, d$expected, .scan_fmt_pct(d$excess_pct), cin)
}

.scan_guide <- function(lang) {
  if (lang == "es") return(c(
    "----------------------------------------------------------------------",
    "GU\u00cdA R\u00c1PIDA",
    "Qu\u00e9 es EBEC: un test que comprueba si los PEORES casos del modelo",
    "  acumulan m\u00e1s error cuadr\u00e1tico (Brier) del esperable si sus",
    "  probabilidades fueran ciertas. Brier medio es la nota media; EBEC mira",
    "  los peores suspensos, que la media esconde.",
    "La tabla: cada fila es una medida de p\u00e9rdida.",
    "  component   p\u00e9rdida total, o suma de las k peores p\u00e9rdidas",
    "              ('FP' = peores falsos positivos: dijo alto y no ocurri\u00f3).",
    "  cases       % de la muestra que son esos k casos (1%, 5%, 10%, 20%).",
    "  observed    esa suma en tus datos.",
    "  expected    su valor medio si las probabilidades fueran ciertas.",
    "  ratio       observed / expected (1 = lo esperado; 1,10 = 10% m\u00e1s).",
    "  excess_pct  exceso en %, = 100 * (ratio - 1), con su IC bootstrap.",
    "  z           exceso estandarizado (cu\u00e1nto supera lo esperado).",
    "  p_component p sin ajustar de esa fila: solo descriptivo, no reportar.",
    "El % de peores casos: NO se elige. EBEC mira 1%, 5%, 10% y 20% a la vez",
    "  y el p-valor global ya corrige esa elecci\u00f3n.",
    "Rango: el estad\u00edstico (z m\u00e1xima) no tiene l\u00edmite (-inf a +inf). El",
    "  p-valor va de 1/(B+1) a 1. El test es unilateral: solo un exceso de",
    "  p\u00e9rdida da evidencia; z <= 0 nunca es significativo.",
    "Valor y p-valor: z alta y p < 0,05 = m\u00e1s error grave del esperable;",
    "  z alta y p >= 0,05 = sospechoso pero compatible con el azar;",
    "  z cerca de 0 o negativa = sin problema (negativa = menos p\u00e9rdida de",
    "  la esperada, no se eval\u00faa).",
    "Reportar: el p-valor global + el exceso % del componente principal con su",
    "  IC (frase arriba). La tabla completa, como material suplementario.",
    "Cautelas:",
    "  * Usa predicciones fuera de muestra o de validaci\u00f3n cruzada; con las del",
    "    mismo modelo ajustado el test es conservador.",
    "  * El componente principal es un m\u00e1ximo seleccionado y su IC no est\u00e1",
    "    corregido por mirar varios: es orientativo. No lo presentes como",
    "    hallazgo confirmado.",
    "  * El IC es bootstrap y descriptivo. En una simulaci\u00f3n peque\u00f1a (n = 400,",
    "    ~8% de eventos) su cobertura fue de ~88-93% para un 95% nominal.",
    "  * Complementa a los tests de calibraci\u00f3n cl\u00e1sicos (p. ej. Spiegelhalter),",
    "    que son mejores para desajustes suaves y globales.",
    "  * Tama\u00f1o comprobado con prevalencia de 3 a 50%. La ventaja de potencia",
    "    se mostr\u00f3 con 10-50%; con 3% solo para fallos de falsos positivos",
    "    y mixtos (ning\u00fan test detect\u00f3 los fallos de falsos negativos).",
    "  * No se ha demostrado ninguna propiedad espec\u00edfica de falsos negativos.",
    "----------------------------------------------------------------------"))
  c(
    "----------------------------------------------------------------------",
    "QUICK GUIDE",
    "What EBEC is: a test of whether the WORST cases of a model accumulate more",
    "  squared (Brier) error than expected if its probabilities were true. The",
    "  mean Brier score is the average grade; EBEC looks at the worst fails,",
    "  which the mean hides.",
    "The table: each row is one loss measure.",
    "  component   total loss, or sum of the k worst losses ('FP' = worst false",
    "              positives: predicted high and it did not happen).",
    "  cases       share of the sample that those k cases represent",
    "              (1%, 5%, 10%, 20%).",
    "  observed    that sum in your data.",
    "  expected    its average if the probabilities were true.",
    "  ratio       observed / expected (1 = as expected; 1.10 = 10% more).",
    "  excess_pct  excess in %, = 100 * (ratio - 1), with its bootstrap CI.",
    "  z           standardized excess (how far above what is expected).",
    "  p_component unadjusted p of that row: descriptive only, do not report.",
    "The % of worst cases: it is NOT chosen. EBEC looks at 1%, 5%, 10% and 20%",
    "  at once and the global p-value already corrects for that choice.",
    "Range: the statistic (largest z) is unbounded (-inf to +inf). The p-value",
    "  goes from 1/(B+1) to 1. The test is one-sided: only an excess of loss is",
    "  evidence; z <= 0 is never significant.",
    "Value and p-value: high z and p < 0.05 = more serious error than expected;",
    "  high z and p >= 0.05 = suspicious but compatible with chance;",
    "  z near 0 or negative = no problem (negative = less loss than expected,",
    "  not evaluated).",
    "Reporting: the global p-value + the % excess of the main component with",
    "  its CI (sentence above). Full table as supplementary material.",
    "Cautions:",
    "  * Use held-out or cross-validated predictions; with predictions from the",
    "    same fitted model the test is conservative.",
    "  * The main component is a selected maximum and its CI is not adjusted",
    "    for looking at several: indicative only. Do not present it as a",
    "    confirmed finding.",
    "  * The CI is a bootstrap and descriptive. In a small simulation (n = 400,",
    "    about 8% events) its coverage was about 88-93% for a nominal 95%.",
    "  * It complements classical calibration tests (e.g. Spiegelhalter), which",
    "    are better for mild, global miscalibration.",
    "  * Size checked at 3-50% prevalence. The power advantage was shown at",
    "    10-50% prevalence; at 3% only for false-positive and mixed faults",
    "    (false-negative faults were not detected by any test examined).",
    "  * No false-negative-specific property has been demonstrated.",
    "----------------------------------------------------------------------")
}
