# ==============================================================================
# EBEC - ANALISIS DE RESULTADOS (NO SIMULA NADA, tarda segundos)
# Lee EBEC_all_pvalues.csv ya generado y produce tablas + informe.
# Uso:  source("scripts/EBEC_analyze.R")
# ==============================================================================
ALPHA <- 0.05
cand <- list.files(".", pattern = "^EBEC_all_pvalues\\.csv$", recursive = TRUE, full.names = TRUE)
cand <- cand[!grepl("results_quick", cand)]
if (length(cand) == 0) stop("No encuentro EBEC_all_pvalues.csv. Ejecuta list.files(recursive=TRUE, pattern='pvalues') y dime donde esta.")
ALL_FILE <- cand[which.max(file.mtime(cand))]
OUT_DIR <- dirname(ALL_FILE)
cat("Leyendo:", normalizePath(ALL_FILE), "\n")
ALL <- read.csv(ALL_FILE, stringsAsFactors = FALSE)
cat("Filas:", nrow(ALL), "\n"); print(table(ALL$hyp)); print(table(ALL$test))
QUICK <- FALSE; SEED <- NA; N_REP_H0 <- max(table(paste(ALL$n, ALL$prev, ALL$q, ALL$test)[ALL$hyp == "H0"]))
N_REP_H1 <- NA; B_MC <- NA
CI_FILE <- file.path(OUT_DIR, "EBEC_ci_coverage.csv")
CI <- if (file.exists(CI_FILE)) read.csv(CI_FILE) else data.frame(test = "no calculado", n = NA_real_, q = NA_real_, coverage_of_0 = NA_real_, mean_width = NA_real_)

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

# ---- 6. INFORME AUTOMATICO ---------------------------------------------------
f <- function(x) formatC(x, digits = 3, format = "f")
rep_lines <- c(
  "# EBEC / EBEC_FN / EBEC_FP - informe de validacion por simulacion", "",
  sprintf("Generado: %s | analisis de %s", format(Sys.time(), "%Y-%m-%d %H:%M"), basename(ALL_FILE)),
  sprintf("Replicas H0 por celda y test: %d; alpha = %.2f", N_REP_H0, ALPHA), "",
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
  sprintf("| %s | %s | %s | %s | %s |", CI$test, CI$n, CI$q, f(CI$coverage_of_0), f(CI$mean_width)), "",
  "Cobertura esperada: 0.95. Si es claramente menor, el IC no esta validado y no debe usarse en el paper.")
writeLines(rep_lines, file.path(OUT_DIR, "EBEC_validation_report.md"))
cat(paste(rep_lines, collapse = "\n"), "\n")
cat("\nArchivos en:", normalizePath(OUT_DIR), "\n")
