# ==============================================================================
# EBEC / EBEC_FN / EBEC_FP -- VALIDACION EMPIRICA POR SIMULACION
# ------------------------------------------------------------------------------
# Donde va: <paquete>/validation/EBEC_validation_final.R
# Como se ejecuta: abrir el proyecto EBEC en RStudio y ejecutar
#       source("validation/EBEC_validation_final.R")
# (el directorio de trabajo debe ser la raiz del paquete)
# Prueba rapida (minutos):  Sys.setenv(EBEC_QUICK = "1")  antes del source().
#
# Este script valida las funciones REALES del paquete (R/ebec.R), no copias.
# Diseno fijado de antemano: se informan TODOS los escenarios, sin seleccionar.
# ==============================================================================

# ---- 0. CONFIGURACION (no tocar tras ver resultados) -------------------------
QUICK <- identical(Sys.getenv("EBEC_QUICK"), "1")

SEED       <- 20261018
N_REP_H0   <- if (QUICK) 100 else 2000   # replicas bajo H0 por celda
N_REP_H1   <- if (QUICK) 50  else 300    # replicas bajo H1 por celda
N_REP_CI   <- if (QUICK) 50  else 300    # replicas para cobertura del IC
B_MC       <- if (QUICK) 99  else 499    # simulaciones Monte Carlo por test
B_BOOT     <- if (QUICK) 99  else 299    # remuestras bootstrap (solo bloque IC)
ALPHA      <- 0.05
N_CORES    <- max(1, parallel::detectCores() - 1)
if (QUICK) N_CORES <- min(N_CORES, 2)

N_H0   <- c(100, 250, 500, 1000)
N_H1   <- c(250, 500, 1000)
PREV   <- c(0.10, 0.30, 0.50)
Q_GRID <- c(0.05, 0.10, 0.20)
if (QUICK) { N_H0 <- c(100, 500); N_H1 <- 500; PREV <- c(0.10, 0.50) }

OUT_DIR <- file.path("validation", "results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PKG_FILE <- file.path("R", "ebec.R")
stopifnot(file.exists(PKG_FILE))
sys.source(PKG_FILE, envir = globalenv())

# ---- 1. GENERADORES ----------------------------------------------------------
# p NUNCA usa y. Probabilidades predichas heterogeneas (Beta, media = prev).
gen_p <- function(n, prev) {
  a <- max(prev * 8, 0.1); b <- max((1 - prev) * 8, 0.1)
  pmin(pmax(rbeta(n, a, b), 0.001), 0.999)
}

# Escenarios H1: la probabilidad VERDADERA difiere de la predicha en un subconjunto
H1_SCENARIOS <- list()
for (cc in c(0.05, 0.10)) for (ss in c(0.5, 1))
  H1_SCENARIOS[[length(H1_SCENARIOS) + 1]] <-
    list(name = sprintf("inversion c=%.2f s=%.1f", cc, ss), type = "inv", c = cc, s = ss)
for (cc in c(0.05, 0.10)) for (dd in c(1.5, 3))
  H1_SCENARIOS[[length(H1_SCENARIOS) + 1]] <-
    list(name = sprintf("FN_infraprediccion c=%.2f d=%.1f", cc, dd), type = "fn", c = cc, d = dd)
for (cc in c(0.05, 0.10)) for (dd in c(1.5, 3))
  H1_SCENARIOS[[length(H1_SCENARIOS) + 1]] <-
    list(name = sprintf("FP_sobreprediccion c=%.2f d=%.1f", cc, dd), type = "fp", c = cc, d = dd)
for (kk in c(0.5, 0.75))
  H1_SCENARIOS[[length(H1_SCENARIOS) + 1]] <-
    list(name = sprintf("sobreconfianza global k=%.2f", kk), type = "sharp", k = kk)

gen_h1 <- function(scn, n, prev) {
  p <- gen_p(n, prev); pt <- p
  if (scn$type == "inv") {                       # inversion parcial de un subconjunto
    i <- sample.int(n, max(1, round(scn$c * n)))
    pt[i] <- (1 - scn$s) * p[i] + scn$s * (1 - p[i])
  } else if (scn$type == "fn") {                 # eventos con p baja: la verdadera es mayor
    pool <- which(p <= quantile(p, 0.30))
    i <- pool[sample.int(length(pool), min(length(pool), max(1, round(scn$c * n))))]
    pt[i] <- plogis(qlogis(p[i]) + scn$d)
  } else if (scn$type == "fp") {                 # p alta: la verdadera es menor
    pool <- which(p >= quantile(p, 0.70))
    i <- pool[sample.int(length(pool), min(length(pool), max(1, round(scn$c * n))))]
    pt[i] <- plogis(qlogis(p[i]) - scn$d)
  } else if (scn$type == "sharp") {              # el modelo es demasiado extremo
    pt <- plogis(scn$k * qlogis(p))
  }
  pt <- pmin(pmax(pt, 0.001), 0.999)
  list(y = rbinom(n, 1, pt), p = p)
}

# Comparador clasico: test de Spiegelhalter (unilateral, Brier mayor de lo esperado)
spiegelhalter_p <- function(y, p) {
  z <- (mean((y - p)^2) - mean(p * (1 - p))) /
    sqrt(sum(p * (1 - p) * (1 - 2 * p)^2) / length(p)^2)
  1 - pnorm(z)
}

# Un conjunto de datos -> p-valores y excesos de los 3 tests en cada q + Spiegelhalter
run_tests <- function(y, p, qs, B) {
  out <- vector("list", length(qs) * 3 + 1); j <- 0
  for (q in qs) {
    r <- list(EBEC    = suppressWarnings(ebec(y, p, q = q, B = B, ci = FALSE)),
              EBEC_FN = suppressWarnings(ebec_fn(y, p, q = q, B = B, ci = FALSE)),
              EBEC_FP = suppressWarnings(ebec_fp(y, p, q = q, B = B, ci = FALSE)))
    for (nm in names(r)) {
      j <- j + 1
      out[[j]] <- data.frame(test = nm, q = q, p_value = r[[nm]]$p_value,
                             excess = r[[nm]]$excess, k = r[[nm]]$k)
    }
  }
  j <- j + 1
  out[[j]] <- data.frame(test = "Spiegelhalter", q = NA_real_,
                         p_value = spiegelhalter_p(y, p), excess = NA_real_, k = NA_integer_)
  do.call(rbind, out)
}

# ---- 2. CELDAS ---------------------------------------------------------------
run_cell <- function(cell) {
  set.seed(cell$seed)
  res <- vector("list", cell$n_rep)
  for (r in seq_len(cell$n_rep)) {
    d <- if (cell$hyp == "H0") { p <- gen_p(cell$n, cell$prev); list(y = rbinom(cell$n, 1, p), p = p) }
         else gen_h1(cell$scn, cell$n, cell$prev)
    tmp <- run_tests(d$y, d$p, cell$qs, cell$B)
    tmp$rep <- r
    res[[r]] <- tmp
  }
  res <- do.call(rbind, res)
  res$hyp <- cell$hyp; res$scenario <- cell$scenario
  res$n <- cell$n; res$prev <- cell$prev
  res
}

cells <- list(); id <- 0
for (n in N_H0) for (pv in PREV) {
  id <- id + 1
  cells[[id]] <- list(hyp = "H0", scenario = "H0", scn = NULL, n = n, prev = pv,
                      qs = Q_GRID, n_rep = N_REP_H0, B = B_MC, seed = SEED + id)
}
for (scn in H1_SCENARIOS) for (n in N_H1) for (pv in PREV) {
  id <- id + 1
  cells[[id]] <- list(hyp = "H1", scenario = scn$name, scn = scn, n = n, prev = pv,
                      qs = Q_GRID, n_rep = N_REP_H1, B = B_MC, seed = SEED + id)
}
cat(sprintf("Celdas: %d (H0: %d, H1: %d) | nucleos: %d | modo %s\n",
            length(cells), length(N_H0) * length(PREV),
            length(cells) - length(N_H0) * length(PREV), N_CORES,
            if (QUICK) "RAPIDO" else "COMPLETO"))

t0 <- Sys.time()
if (N_CORES > 1) {
  cl <- parallel::makeCluster(N_CORES)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("gen_p", "gen_h1", "spiegelhalter_p", "run_tests",
                                "run_cell", "H1_SCENARIOS", "PKG_FILE"), envir = environment())
  parallel::clusterEvalQ(cl, sys.source(PKG_FILE, envir = globalenv()))
  out <- parallel::parLapplyLB(cl, cells, run_cell)
} else {
  out <- lapply(cells, run_cell)
}
ALL <- do.call(rbind, out)
cat(sprintf("Simulacion terminada en %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(ALL, file.path(OUT_DIR, "EBEC_all_pvalues.csv"), row.names = FALSE)

# ---- 3. TIPO I BAJO H0 -------------------------------------------------------
wilson <- function(x, n, z = 1.96) {
  ph <- x / n; d <- 1 + z^2 / n
  c(lo = (ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / d,
    hi = (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / d)
}
summ_rej <- function(df) {
  ok <- !is.na(df$p_value); x <- sum(df$p_value[ok] < ALPHA); m <- sum(ok)
  w <- wilson(x, m)
  data.frame(n_rep = m, rate = x / m, ci_lo = w["lo"], ci_hi = w["hi"],
             rate_01 = mean(df$p_value[ok] < 0.01), rate_10 = mean(df$p_value[ok] < 0.10),
             mean_excess = mean(df$excess, na.rm = TRUE), prop_NA = mean(!ok),
             mean_k = mean(df$k, na.rm = TRUE), row.names = NULL)
}
agg <- function(df, by) {
  key <- do.call(paste, c(df[by], sep = "|"))
  res <- lapply(split(df, key), function(d) cbind(d[1, by, drop = FALSE], summ_rej(d)))
  r <- do.call(rbind, res); rownames(r) <- NULL; r
}

H0 <- ALL[ALL$hyp == "H0", ]
T1 <- agg(H0, c("test", "n", "prev", "q"))
T1$pass <- T1$ci_lo <= ALPHA & T1$ci_hi >= ALPHA          # IC de Wilson contiene 0.05
T1$verdict <- ifelse(T1$pass, "OK", ifelse(T1$rate > ALPHA, "ANTICONSERVADOR", "CONSERVADOR"))
write.csv(T1, file.path(OUT_DIR, "EBEC_typeI_by_cell.csv"), row.names = FALSE)

T1_global <- agg(H0, c("test", "q"))
T1_global$cells_pass <- sapply(seq_len(nrow(T1_global)), function(i) {
  s <- T1[T1$test == T1_global$test[i] & ((is.na(T1$q) & is.na(T1_global$q[i])) | T1$q %in% T1_global$q[i]), ]
  sprintf("%d/%d", sum(s$pass), nrow(s))
})
write.csv(T1_global, file.path(OUT_DIR, "EBEC_typeI_global.csv"), row.names = FALSE)

# ---- 4. POTENCIA BAJO H1 -----------------------------------------------------
H1 <- ALL[ALL$hyp == "H1", ]
PW <- agg(H1, c("scenario", "test", "n", "prev", "q"))
write.csv(PW, file.path(OUT_DIR, "EBEC_power_by_cell.csv"), row.names = FALSE)

pw_scn <- agg(H1, c("scenario", "test", "q"))
pw_scn <- pw_scn[, c("scenario", "test", "q", "rate")]
names(pw_scn)[4] <- "power"
pw_wide <- reshape(transform(pw_scn, col = paste0(test, ifelse(is.na(q), "", paste0("_q", q)))),
                   idvar = "scenario", timevar = "col", direction = "wide",
                   v.names = "power", drop = c("test", "q"))
names(pw_wide) <- sub("^power\\.", "", names(pw_wide))
write.csv(pw_wide, file.path(OUT_DIR, "EBEC_power_by_scenario.csv"), row.names = FALSE)

# ---- 5. COBERTURA DEL INTERVALO BOOTSTRAP (bloque separado) ------------------
cat("Cobertura del IC bootstrap...\n")
ci_cell <- function(a) {
  set.seed(a$seed)
  cov <- matrix(NA, a$n_rep, 3); wid <- matrix(NA, a$n_rep, 3)
  for (r in seq_len(a$n_rep)) {
    p <- gen_p(a$n, a$prev); y <- rbinom(a$n, 1, p)
    rr <- suppressWarnings(list(ebec(y, p, q = a$q, B = a$B, ci = TRUE),
                                ebec_fn(y, p, q = a$q, B = a$B, ci = TRUE),
                                ebec_fp(y, p, q = a$q, B = a$B, ci = TRUE)))
    for (j in 1:3) { cov[r, j] <- rr[[j]]$ci_lower <= 0 & rr[[j]]$ci_upper >= 0
                     wid[r, j] <- rr[[j]]$ci_upper - rr[[j]]$ci_lower }
  }
  data.frame(test = c("EBEC", "EBEC_FN", "EBEC_FP"), n = a$n, prev = a$prev, q = a$q,
             coverage_of_0 = colMeans(cov, na.rm = TRUE), mean_width = colMeans(wid, na.rm = TRUE))
}
ci_args <- list(); k <- 0
for (n in N_H1[N_H1 <= 500]) for (q in c(0.10, 0.20)) {
  k <- k + 1
  ci_args[[k]] <- list(n = n, prev = 0.30, q = q, n_rep = N_REP_CI, B = B_BOOT, seed = SEED + 5000 + k)
}
if (N_CORES > 1) {
  parallel::clusterExport(cl, c("ci_cell"), envir = environment())
  CI <- do.call(rbind, parallel::parLapplyLB(cl, ci_args, ci_cell))
} else CI <- do.call(rbind, lapply(ci_args, ci_cell))
write.csv(CI, file.path(OUT_DIR, "EBEC_ci_coverage.csv"), row.names = FALSE)

# ---- 6. INFORME AUTOMATICO ---------------------------------------------------
f <- function(x) formatC(x, digits = 3, format = "f")
rep_lines <- c(
  "# EBEC / EBEC_FN / EBEC_FP - informe de validacion por simulacion", "",
  sprintf("Generado: %s | semilla base %d | modo %s", format(Sys.time(), "%Y-%m-%d %H:%M"), SEED,
          if (QUICK) "RAPIDO (NO usar en el paper)" else "COMPLETO"),
  sprintf("Replicas: H0 = %d por celda, H1 = %d por celda; B Monte Carlo = %d; alpha = %.2f",
          N_REP_H0, N_REP_H1, B_MC, ALPHA), "",
  "## 1. Error tipo I (H0 verdadera: y ~ Bernoulli(p), p fijo)", "",
  "Criterio fijado de antemano: la celda pasa si el IC de Wilson al 95% de la tasa de rechazo contiene 0.05.", "",
  "| test | q | tasa media | min celda | max celda | celdas OK |", "|---|---|---|---|---|---|")
for (i in seq_len(nrow(T1_global))) {
  s <- T1[T1$test == T1_global$test[i] & ((is.na(T1$q) & is.na(T1_global$q[i])) | T1$q %in% T1_global$q[i]), ]
  rep_lines <- c(rep_lines, sprintf("| %s | %s | %s | %s | %s | %s |", T1_global$test[i],
        ifelse(is.na(T1_global$q[i]), "-", T1_global$q[i]), f(T1_global$rate[i]),
        f(min(s$rate)), f(max(s$rate)), T1_global$cells_pass[i]))
}
bad <- T1[!T1$pass, c("test", "n", "prev", "q", "rate", "verdict")]
rep_lines <- c(rep_lines, "", sprintf("Celdas que NO pasan: %d de %d.", nrow(bad), nrow(T1)))
if (nrow(bad) > 0) {
  bad$rate <- f(bad$rate)
  rep_lines <- c(rep_lines, "", "| test | n | prev | q | tasa | veredicto |", "|---|---|---|---|---|---|",
                 apply(bad, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
rep_lines <- c(rep_lines, "",
  "Nota: con n pequeno y q pequeno el subconjunto de cola tiene pocas observaciones y el p-valor es muy discreto (tests conservadores).", "",
  "## 2. Potencia (media sobre n y prevalencia; por escenario y q)", "",
  "Detalle completo en EBEC_power_by_cell.csv. Se informan TODOS los escenarios.", "")
pw_tab <- agg(H1, c("scenario", "test", "q"))
for (sc in unique(pw_tab$scenario)) {
  s <- pw_tab[pw_tab$scenario == sc, ]
  rep_lines <- c(rep_lines, sprintf("**%s**", sc), "", "| test | q | potencia media |", "|---|---|---|",
                 sprintf("| %s | %s | %s |", s$test, ifelse(is.na(s$q), "-", s$q), f(s$rate)), "")
}
rep_lines <- c(rep_lines, "## 3. Cobertura del IC bootstrap (la verdad es exceso = 0 bajo H0)", "",
  "| test | n | q | cobertura de 0 | anchura media |", "|---|---|---|---|---|",
  sprintf("| %s | %d | %.2f | %s | %s |", CI$test, CI$n, CI$q, f(CI$coverage_of_0), f(CI$mean_width)), "",
  "Cobertura esperada: 0.95. Si es claramente menor, el IC no esta validado y no debe usarse en el paper.")
writeLines(rep_lines, file.path(OUT_DIR, "EBEC_validation_report.md"))
cat(paste(rep_lines, collapse = "\n"), "\n")
cat("\nArchivos en:", normalizePath(OUT_DIR), "\n")
