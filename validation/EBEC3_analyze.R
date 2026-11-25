# ==============================================================================
# EBEC-scan -- ANALISIS (no simula). Lee las celdas guardadas y aplica C1-C5.
# Uso:  Sys.setenv(EBEC_MODE = "pilot"); source("validation/EBEC3_analyze.R")
# ==============================================================================
MODE <- Sys.getenv("EBEC_MODE", "dry")
ALPHA <- 0.05
OUT_DIR <- file.path("validation", paste0("results_v3_", MODE))
CELL_DIR <- file.path(OUT_DIR, "cells")
files <- list.files(CELL_DIR, pattern = "^cell_[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("No hay celdas en ", CELL_DIR)
LST <- lapply(files, readRDS)
is_ins <- vapply(LST, function(d) "scan_plain" %in% names(d), logical(1))
D <- do.call(rbind, LST[!is_ins])                       # H0 y H1
INS <- if (any(is_ins)) do.call(rbind, LST[is_ins]) else NULL
cat("Celdas leidas:", length(files), "| filas H0/H1:", nrow(D), "| in-sample:", if (is.null(INS)) 0 else nrow(INS), "\n")

wilson <- function(x, n, z = 1.96) { ph <- x / n; d <- 1 + z^2 / n
  c(lo = (ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / d,
    hi = (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / d) }
TESTS <- c("scan_fp", "scan_both", "scan_tot", "spiegelhalter", "hl", "cox")
rej <- function(p) mean(p < ALPHA, na.rm = TRUE)
f3 <- function(x) ifelse(is.na(x), "NA", formatC(as.numeric(x), digits = 3, format = "f"))
tag <- function(x) if (isTRUE(x)) "PASA" else if (identical(x, FALSE)) "FALLA" else "NO EVALUABLE"
L <- character()
say <- function(...) { s <- sprintf(...); L <<- c(L, s); cat(s, "\n", sep = "") }

# ---- C1: tipo I --------------------------------------------------------------
H0 <- D[D$hyp == "H0", ]
say("# EBEC-scan: informe de validacion v3 (scan_fp) (modo %s)", MODE); say("")
say("## C1. Error tipo I bajo H0")
cellk <- unique(H0[, c("id", "n", "prev")])
rows <- list()
for (i in seq_len(nrow(cellk))) { d <- H0[H0$id == cellk$id[i], ]
  for (t in TESTS) { x <- sum(d[[t]] < ALPHA, na.rm = TRUE); m <- sum(!is.na(d[[t]])); w <- wilson(x, m)
    rows[[length(rows) + 1]] <- data.frame(test = t, n = cellk$n[i], prev = cellk$prev[i], rate = x / m, lo = w[1], hi = w[2], reps = m) } }
T1 <- do.call(rbind, rows); rownames(T1) <- NULL
write.csv(T1, file.path(OUT_DIR, "typeI_by_cell.csv"), row.names = FALSE)
G1 <- do.call(rbind, lapply(TESTS, function(t) { s <- T1[T1$test == t, ]; N <- sum(s$reps)
  data.frame(test = t, pooled_rate = sum(s$rate * s$reps) / N, min_cell = min(s$rate), max_cell = max(s$rate),
             cells = nrow(s), cells_liberal = sum(s$lo > 0.05), cells_conservative = sum(s$hi < 0.05), N = N) }))
write.csv(G1, file.path(OUT_DIR, "typeI_global.csv"), row.names = FALSE)
say(""); say("| test | tasa agregada | min celda | max celda | celdas | liberales | conservadoras |"); say("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(G1))) say("| %s | %s | %s | %s | %d | %d | %d |", G1$test[i], f3(G1$pooled_rate[i]), f3(G1$min_cell[i]), f3(G1$max_cell[i]), G1$cells[i], G1$cells_liberal[i], G1$cells_conservative[i])
g <- G1[G1$test == "scan_fp", ]; se <- sqrt(0.05 * 0.95 / g$N)
C1 <- g$pooled_rate <= 0.05 + 3 * se && (g$N < 1000 || g$pooled_rate >= 0.03) && g$cells_liberal <= 0.10 * g$cells
say(""); say("**C1 (tamano) = %s**  (tasa %.4f; limite superior %.4f; celdas liberales %d de %d)", tag(C1), g$pooled_rate, 0.05 + 3 * se, g$cells_liberal, g$cells)

# ---- potencia ----------------------------------------------------------------
H1 <- D[D$hyp == "H1", ]
C2 <- C3 <- NA
if (nrow(H1)) {
  PW <- do.call(rbind, lapply(split(H1, list(H1$scenario, H1$n, H1$prev), drop = TRUE), function(d)
    data.frame(scenario = d$scenario[1], group = d$group[1], target = d$target[1], n = d$n[1], prev = d$prev[1],
               t(sapply(TESTS, function(t) rej(d[[t]]))), reps = nrow(d))))
  rownames(PW) <- NULL
  write.csv(PW, file.path(OUT_DIR, "power_by_cell.csv"), row.names = FALSE)
  PS <- do.call(rbind, lapply(split(PW, PW$scenario), function(d)
    data.frame(scenario = d$scenario[1], group = d$group[1], target = d$target[1], t(colMeans(d[, TESTS])), cells = nrow(d))))
  rownames(PS) <- NULL
  write.csv(PS, file.path(OUT_DIR, "power_by_scenario.csv"), row.names = FALSE)
  say(""); say("## Potencia por escenario (media de las celdas ejecutadas)"); say("")
  say("| escenario | grupo | scan_fp | scan_both | scan_tot | Spiegelhalter | HL | Cox |"); say("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(PS))) say("| %s | %s | %s | %s | %s | %s | %s | %s |", PS$scenario[i], PS$group[i], f3(PS$scan_fp[i]), f3(PS$scan_both[i]), f3(PS$scan_tot[i]), f3(PS$spiegelhalter[i]), f3(PS$hl[i]), f3(PS$cox[i]))
  # C2
  r <- PS[PS$group == "rare", ]
  if (nrow(r)) { k2 <- sum(r$scan_fp - r$spiegelhalter >= 0.15); need <- ceiling(nrow(r) * 8 / 12); C2 <- k2 >= need
    say(""); say("**C2 (errores raros confiados) = %s**  (%d de %d escenarios con ventaja >= 0.15 sobre Spiegelhalter; se exigen %d)", tag(C2), k2, nrow(r), need) }
  # C3 (enmienda 2): solo FP es afirmable; FN se informa sin criterio
  c3 <- function(tg) { s <- PS[PS$target == tg, ]; if (!nrow(s)) return(NA); sum(s$scan_fp - s$scan_tot >= 0.05) / nrow(s) }
  fpk <- c3("FP"); fnk <- c3("FN")
  C3 <- if (is.na(fpk)) NA else (fpk >= 0.5)
  say(""); say("**C3 (componente FP) = %s**  (escenarios FP con scan_fp - scan_tot >= 0.05: %s)", tag(C3), f3(fpk))
  say("(Informativo, sin criterio: escenarios FN con scan_fp - scan_tot >= 0.05: %s)", f3(fnk))
  # C6: coste de incluir el componente FP donde no hace falta
  nf <- PS[!PS$target %in% "FP", ]
  C6 <- if (nrow(nf)) mean(nf$scan_fp - nf$scan_tot) >= -0.03 else NA
  say(""); say("**C6 (coste fuera de FP, informativo) = %s**  (media scan_fp - scan_tot en escenarios no FP = %s; se tolera >= -0.03)", tag(C6), f3(if (nrow(nf)) mean(nf$scan_fp - nf$scan_tot) else NA))
  say("Referencia: media scan_fp - scan_both en todos los escenarios = %s", f3(mean(PS$scan_fp - PS$scan_both)))
  # C4
  gl <- PS[PS$group == "global", ]
  if (nrow(gl)) { say(""); say("**C4 (informativo) global:** scan_fp - Spiegelhalter = %s", paste(sprintf("%s: %s", gl$scenario, f3(gl$scan_fp - gl$spiegelhalter)), collapse = "; ")) }
}

# ---- C5: in-sample -----------------------------------------------------------
if (!is.null(INS) && nrow(INS)) {
  say(""); say("## C5. Predicciones de un modelo ajustado (tipo I)"); say(""); say("| n | ajustadas en la misma muestra | validacion cruzada 5-fold | repeticiones |"); say("|---|---|---|---|")
  for (n in sort(unique(INS$n))) { d <- INS[INS$n == n, ]; say("| %d | %s | %s | %d |", n, f3(rej(d$scan_plain)), f3(rej(d$scan_cv)), nrow(d)) }
}

# ---- veredicto ---------------------------------------------------------------
say(""); say("## Veredicto segun PREREGISTRATION.md")
verdict <- local({
  if (!C1) {
    "C1 FALLA: el test no controla el error tipo I. Parar y arreglar antes de nada."
  } else if (isTRUE(C2) && isTRUE(C3)) {
    "C1, C2 y C3 PASAN: test valido con ventaja demostrada en errores raros y con componente de falsos positivos. FN: sin afirmacion."
  } else if (isTRUE(C2)) {
    "C1 y C2 PASAN, C3 FALLA: ventaja en errores raros; NO afirmar utilidad del componente FP."
  } else if (identical(C2, FALSE)) {
    "C1 PASA, C2 FALLA: test valido pero sin ventaja de potencia demostrada; presentar como diagnostico descriptivo."
  } else {
    "C1 PASA; C2/C3 no evaluables en este modo (faltan escenarios)."
  }
})
say("%s", verdict)
if (MODE == "dry") say("(Modo dry: PRUEBA TECNICA del codigo, sin significado estadistico.)") else if (MODE == "pilot") say("(Modo pilot: orientativo, NO usar en el paper.)")
writeLines(L, file.path(OUT_DIR, "REPORT.md"))
cat("\nInforme guardado en:", normalizePath(file.path(OUT_DIR, "REPORT.md")), "\n")
