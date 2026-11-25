# ==============================================================================
# EBEC-scan -- VALIDACION v3 (default FP; modos dry / pilot / full, con checkpoints)
# Donde va: <paquete>/validation/EBEC3_validation.R
# Ejecutar con el directorio de trabajo en la RAIZ del paquete:
#   Sys.setenv(EBEC_MODE = "dry");   source("validation/EBEC3_validation.R")   # 1-3 min
#   Sys.setenv(EBEC_MODE = "pilot"); source("validation/EBEC3_validation.R")   # ~15-30 min
#   Sys.setenv(EBEC_MODE = "full");  source("validation/EBEC3_validation.R")   # horas
# Criterios fijados en validation/PREREGISTRATION.md, Enmienda 2 (no tocar tras ver resultados).
# Se puede cortar y reanudar: cada celda terminada se guarda en disco.
# ==============================================================================
MODE <- Sys.getenv("EBEC_MODE", "dry")
if (!MODE %in% c("dry", "pilot", "full", "lowprev")) stop("EBEC_MODE debe ser dry, pilot, full o lowprev")
SEED <- if (MODE == "lowprev") 20261021 else 20261020
ALPHA <- 0.05
Q_GRID <- c(0.01, 0.05, 0.10, 0.20)

N_CORES <- max(1L, parallel::detectCores() - 1L)
if (nzchar(Sys.getenv("EBEC_CORES"))) N_CORES <- as.integer(Sys.getenv("EBEC_CORES"))
if (MODE == "dry") N_CORES <- 2L

CFG <- switch(MODE,
  dry   = list(reps_h0 = 20,   reps_h1 = 10,  reps_ins = 10,  B = 49,  B_ins = 29,
               n_h0 = 250, prev_h0 = 0.30, n_h1 = 250, prev_h1 = 0.30, n_ins = 250, small = TRUE),
  pilot = list(reps_h0 = 200,  reps_h1 = 100, reps_ins = 100, B = 199, B_ins = 99,
               n_h0 = c(250, 500), prev_h0 = c(0.30, 0.50), n_h1 = 500, prev_h1 = 0.30, n_ins = 250, small = FALSE),
  full  = list(reps_h0 = 2000, reps_h1 = 300, reps_ins = 300, B = 499, B_ins = 199,
               n_h0 = c(100, 250, 500, 1000), prev_h0 = c(0.10, 0.30, 0.50),
               n_h1 = c(250, 500, 1000), prev_h1 = c(0.10, 0.30, 0.50), n_ins = c(250, 500), small = FALSE),
  lowprev = list(reps_h0 = 1000, reps_h1 = 200, reps_ins = 0, B = 499, B_ins = 199,
               n_h0 = c(500, 1000, 2000), prev_h0 = 0.03, n_h1 = c(500, 1000, 2000), prev_h1 = 0.03,
               n_ins = numeric(0), small = FALSE))
SCALE <- as.numeric(Sys.getenv("EBEC_SCALE", "1"))
if (SCALE < 1) { CFG$reps_h0 <- max(2, round(CFG$reps_h0 * SCALE)); CFG$reps_h1 <- max(2, round(CFG$reps_h1 * SCALE)); CFG$B <- 49 }

OUT_DIR  <- file.path("validation", paste0("results_v3_", MODE))
CELL_DIR <- file.path(OUT_DIR, "cells")
dir.create(CELL_DIR, recursive = TRUE, showWarnings = FALSE)
PROG <- file.path(OUT_DIR, "progress.log")
PKG_FILES <- c(file.path("R", "ebec.R"))   # contiene ebec() y los helpers .scan_*
stopifnot(all(file.exists(PKG_FILES)))
for (f in PKG_FILES) sys.source(f, envir = globalenv())

# ---- generadores (p NUNCA usa y) ---------------------------------------------
gen_p <- function(n, prev) {
  a <- max(prev * 8, 0.1); b <- max((1 - prev) * 8, 0.1)
  pmin(pmax(rbeta(n, a, b), 0.001), 0.999)
}
gen_h0 <- function(sc, n, prev) { p <- gen_p(n, prev); list(y = rbinom(n, 1, p), p = p) }
gen_h1 <- function(sc, n, prev) {
  p <- gen_p(n, prev); pt <- p
  if (sc$kind == "inv") {
    i <- sample.int(n, max(1, round(sc$c * n)))
    pt[i] <- (1 - sc$s) * p[i] + sc$s * (1 - p[i])
  } else if (sc$kind == "fn") {
    pool <- which(p <= quantile(p, 0.30))
    i <- pool[sample.int(length(pool), min(length(pool), max(1, round(sc$c * n))))]
    pt[i] <- plogis(qlogis(p[i]) + sc$d)
  } else if (sc$kind == "fp") {
    pool <- which(p >= quantile(p, 0.70))
    i <- pool[sample.int(length(pool), min(length(pool), max(1, round(sc$c * n))))]
    pt[i] <- plogis(qlogis(p[i]) - sc$d)
  } else if (sc$kind == "sharp") {
    pt <- plogis(sc$k * qlogis(p))
  } else if (sc$kind == "rare") {
    # El modelo declara probabilidades extremas en una pequena fraccion de unidades;
    # la probabilidad verdadera es la original. Los extremos NO dependen de y.
    m <- max(1, ceiling(sc$c * n)); i <- sample.int(n, m)
    ps <- p
    ps[i] <- if (sc$type == "FN") sc$mag else if (sc$type == "FP") 1 - sc$mag
             else ifelse(runif(m) < 0.5, sc$mag, 1 - sc$mag)
    y <- rbinom(n, 1, p)
    return(list(y = y, p = pmin(pmax(ps, 0.001), 0.999)))
  }
  pt <- pmin(pmax(pt, 0.001), 0.999)
  list(y = rbinom(n, 1, pt), p = p)
}

SC <- list()
add <- function(name, group, target, ...) SC[[length(SC) + 1]] <<- c(list(name = name, group = group, target = target), list(...))
for (cc in c(0.05, 0.10)) for (ss in c(0.5, 1))
  add(sprintf("Inversion c=%.2f s=%.1f", cc, ss), "inv", "any", kind = "inv", c = cc, s = ss)
for (cc in c(0.05, 0.10)) for (dd in c(1.5, 3))
  add(sprintf("FN under-prediction c=%.2f d=%.1f", cc, dd), "fn_under", "FN", kind = "fn", c = cc, d = dd)
for (cc in c(0.05, 0.10)) for (dd in c(1.5, 3))
  add(sprintf("FP over-prediction c=%.2f d=%.1f", cc, dd), "fp_over", "FP", kind = "fp", c = cc, d = dd)
for (kk in c(0.5, 0.75))
  add(sprintf("Global overconfidence k=%.2f", kk), "global", "any", kind = "sharp", k = kk)
for (ty in c("mixed", "FN", "FP")) for (cc in c(0.01, 0.02)) for (mg in c(0.02, 0.10))
  add(sprintf("Rare confident %s c=%.2f m=%.2f", ty, cc, mg), "rare",
      if (ty == "mixed") "any" else ty, kind = "rare", type = ty, c = cc, mag = mg)
if (isTRUE(CFG$small)) SC <- SC[c(1, 5, length(SC) - 4)]   # dry: 3 escenarios
if (MODE == "lowprev") {
  keep <- c("FN under-prediction c=0.10 d=3.0", "FP over-prediction c=0.10 d=3.0",
            "Global overconfidence k=0.50", "Inversion c=0.10 s=1.0")
  SC <- Filter(function(s) s$group == "rare" || s$name %in% keep, SC)
}

# ---- comparadores ------------------------------------------------------------
spieg_p <- function(y, p) {
  z <- (mean((y - p)^2) - mean(p * (1 - p))) / sqrt(sum(p * (1 - p) * (1 - 2 * p)^2) / length(p)^2)
  1 - pnorm(z)
}
hl_p <- function(y, p, g = 10) {
  o <- order(p); ch <- 0
  for (idx in split(o, cut(seq_along(o), g, labels = FALSE))) {
    e <- sum(p[idx]); ob <- sum(y[idx]); pb <- e / length(idx)
    if (pb <= 0 || pb >= 1) next
    ch <- ch + (ob - e)^2 / (length(idx) * pb * (1 - pb))
  }
  1 - pchisq(ch, g)
}
cox_p <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  tryCatch({
    f <- suppressWarnings(glm(y ~ lp, family = binomial()))
    ll1 <- as.numeric(logLik(f)); ll0 <- sum(y * log(pmax(p, 1e-12)) + (1 - y) * log(pmax(1 - p, 1e-12)))
    1 - pchisq(max(0, 2 * (ll1 - ll0)), 2)
  }, error = function(e) NA_real_)
}

run_one <- function(y, p, B) {
  n <- length(y); ks <- .scan_ks(Q_GRID, n)
  obs <- .scan_stats(y, p, ks, "both"); null <- .scan_null(n, p, ks, "both", B)
  nm <- names(obs)
  both <- .scan_p(obs, null)
  fp   <- .scan_p(obs, null, which(!grepl("^FN ", nm)))
  tot  <- .scan_p(obs, null, which(!grepl("^(FN|FP) ", nm)))
  c(scan_fp = fp$p_value, scan_both = both$p_value, scan_tot = tot$p_value,
    spiegelhalter = spieg_p(y, p), hl = hl_p(y, p), cox = cox_p(y, p))
}
driver_of <- function(y, p, B) NA_character_

# ---- celdas ------------------------------------------------------------------
cells <- list(); id <- 0
for (n in CFG$n_h0) for (pv in CFG$prev_h0) {
  id <- id + 1
  cells[[id]] <- list(id = id, hyp = "H0", scn = list(name = "H0", group = "H0", target = "none"),
                      n = n, prev = pv, reps = CFG$reps_h0, B = CFG$B)
}
for (sc in SC) for (n in CFG$n_h1) for (pv in CFG$prev_h1) {
  id <- id + 1
  cells[[id]] <- list(id = id, hyp = "H1", scn = sc, n = n, prev = pv, reps = CFG$reps_h1, B = CFG$B)
}
for (n in CFG$n_ins) {
  id <- id + 1
  cells[[id]] <- list(id = id, hyp = "INS", scn = list(name = "In-sample logistic fit", group = "ins", target = "none"),
                      n = n, prev = NA_real_, reps = CFG$reps_ins, B = CFG$B_ins)
}

run_cell <- function(cell) {
  fn <- file.path(CELL_DIR, sprintf("cell_%04d.rds", cell$id))
  if (file.exists(fn)) return("skipped")
  t0 <- Sys.time()
  res <- tryCatch({
    set.seed(SEED + cell$id)
    R <- cell$reps
    if (cell$hyp == "INS") {
      M <- matrix(NA_real_, R, 2, dimnames = list(NULL, c("scan_plain", "scan_cv")))
      for (r in seq_len(R)) {
        x1 <- rnorm(cell$n); x2 <- rnorm(cell$n)
        y <- rbinom(cell$n, 1, plogis(-0.5 + 0.8 * x1 - 0.6 * x2))
        dd <- data.frame(y = y, x1 = x1, x2 = x2)
        fit <- suppressWarnings(glm(y ~ x1 + x2, data = dd, family = binomial()))
        M[r, 1] <- ebec(y, as.numeric(fitted(fit)), q = Q_GRID, B = cell$B, ci = FALSE)$p_value
        folds <- sample(rep(1:5, length.out = cell$n)); pcv <- numeric(cell$n)
        for (k in 1:5) {
          fk <- suppressWarnings(glm(y ~ x1 + x2, data = dd[folds != k, ], family = binomial()))
          pcv[folds == k] <- predict(fk, newdata = dd[folds == k, ], type = "response")
        }
        M[r, 2] <- ebec(y, pmin(pmax(pcv, 0.001), 0.999), q = Q_GRID, B = cell$B, ci = FALSE)$p_value
      }
    } else {
      M <- matrix(NA_real_, R, 6, dimnames = list(NULL, c("scan_fp", "scan_both", "scan_tot", "spiegelhalter", "hl", "cox")))
      for (r in seq_len(R)) {
        d <- if (cell$hyp == "H0") gen_h0(cell$scn, cell$n, cell$prev) else gen_h1(cell$scn, cell$n, cell$prev)
        M[r, ] <- run_one(d$y, d$p, cell$B)
      }
    }
    out <- data.frame(id = cell$id, hyp = cell$hyp, scenario = cell$scn$name, group = cell$scn$group,
                      target = cell$scn$target, n = cell$n, prev = cell$prev, rep = seq_len(R), M,
                      row.names = NULL)
    saveRDS(out, paste0(fn, ".tmp")); file.rename(paste0(fn, ".tmp"), fn)
    "ok"
  }, error = function(e) {
    writeLines(conditionMessage(e), sub("\\.rds$", ".err", fn)); "ERROR"
  })
  cat(sprintf("%s cell %d (%s, n=%s) %s %.1fs\n", format(Sys.time(), "%H:%M:%S"), cell$id,
              cell$scn$name, cell$n, res, as.numeric(difftime(Sys.time(), t0, units = "secs"))),
      file = PROG, append = TRUE)
  res
}

# ---- estimacion de tiempo (antes de empezar) ---------------------------------
cat(sprintf("Modo %s | %d celdas | %d nucleos | salida: %s\n", MODE, length(cells), N_CORES, OUT_DIR))
pending <- Filter(function(cl) !file.exists(file.path(CELL_DIR, sprintf("cell_%04d.rds", cl$id))), cells)
cat(sprintf("Celdas pendientes: %d (ya hechas: %d)\n", length(pending), length(cells) - length(pending)))
if (length(pending) > 0) {
  probe <- pending[[which.max(vapply(pending, function(cl) cl$n * cl$reps, numeric(1)))]]
  probe$reps <- 3; probe$id <- 9999
  t_probe <- system.time({ set.seed(1); suppressWarnings(for (r in 1:3) {
    d <- if (probe$hyp == "H0") gen_h0(NULL, probe$n, probe$prev) else if (probe$hyp == "H1") gen_h1(probe$scn, probe$n, probe$prev) else list(y = rbinom(probe$n, 1, .3), p = rep(.3, probe$n))
    run_one(d$y, d$p, probe$B) }) })[["elapsed"]] / 3
  est <- sum(vapply(pending, function(cl) cl$reps * t_probe * (cl$n / probe$n), numeric(1))) / N_CORES / 60
  cat(sprintf("Tiempo estimado (cota prudente): ~%.0f minutos\n", est))
  if (interactive() && MODE %in% c("full", "lowprev")) {
    ans <- readline("Continuar? [s/N] "); if (!tolower(ans) %in% c("s", "y")) stop("Cancelado por el usuario.")
  }
}

# ---- ejecucion ---------------------------------------------------------------
t0 <- Sys.time()
if (length(pending) > 0) {
  if (N_CORES > 1) {
    cl <- parallel::makeCluster(N_CORES)
    tryCatch({
      parallel::clusterExport(cl, c("MODE", "SEED", "ALPHA", "Q_GRID", "CFG", "CELL_DIR", "PROG", "PKG_FILES",
        "gen_p", "gen_h0", "gen_h1", "SC", "spieg_p", "hl_p", "cox_p", "run_one", "driver_of", "run_cell"), envir = globalenv())
      invisible(parallel::clusterEvalQ(cl, for (f in PKG_FILES) sys.source(f, envir = globalenv())))
      status <- parallel::parLapplyLB(cl, pending, run_cell)
    }, finally = try(parallel::stopCluster(cl), silent = TRUE))
  } else status <- lapply(pending, run_cell)
  cat(sprintf("Simulacion terminada en %.1f min. Celdas con error: %d\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(unlist(status) == "ERROR")))
}
errs <- list.files(CELL_DIR, pattern = "\\.err$", full.names = TRUE)
if (length(errs)) { cat("ERRORES en celdas:\n"); for (e in errs) cat(" ", basename(e), ":", readLines(e)[1], "\n") }

# ---- analisis automatico (si falla, los datos ya estan guardados) -------------
tryCatch(sys.source(file.path("validation", "EBEC3_analyze.R"), envir = globalenv()),
         error = function(e) cat("El analisis fallo (los datos estan guardados). Ejecuta source('validation/EBEC3_analyze.R'). Error:", conditionMessage(e), "\n"))
