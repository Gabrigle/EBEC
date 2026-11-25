#' Tail loss share index (legacy; descriptive)
#'
#' Legacy index of the share of Brier loss in the worst q fraction. It is kept
#' for description only: for detection use \code{\link{ebec}}. Tests whether the worst q fraction of predictions concentrates more
#' squared (Brier) loss than expected under model self-consistency.
#'
#' @param y Binary outcome vector (0/1).
#' @param p Predicted probability vector, between 0 and 1.
#' @param q Fraction of worst-loss cases to examine. Default 0.05.
#' @param B Number of Monte Carlo replicates. Default 1000.
#' @param seed Optional random seed. If NULL (default), the current random
#'   number state is used.
#' @param min_k Minimum recommended tail size before flagging instability.
#'   Default 10.
#' @param ci Logical. If TRUE (default) a bootstrap interval for the excess is
#'   computed; use FALSE to skip it (much faster, e.g. in simulations).
#'
#' @return An object of class \code{ebec_share}.
#' @importFrom stats rbinom quantile
#' @examples
#' set.seed(1)
#' p <- runif(200, 0.2, 0.8)
#' y <- rbinom(200, 1, p)
#' ebec_share(y, p, q = 0.05, B = 200, seed = 1)
#' @export
ebec_share <- function(y, p, q = 0.05, B = 1000, seed = NULL, min_k = 10, ci = TRUE) {
  stopifnot(length(y) == length(p))
  stopifnot(all(y %in% c(0, 1)))
  stopifnot(all(p >= 0 & p <= 1))
  stopifnot(q > 0, q < 1, B >= 1)

  n <- length(y)
  k <- ceiling(q * n)
  small_sample_flag <- k < min_k

  # q = 0.05 is a conventional default: justify it for each analysis.
  if (missing(q)) {
    warning(
      "EBEC: you are using q = 0.05 (default value). Choose q according to ",
      "your context (e.g. the proportion of predictions considered 'high risk').",
      call. = FALSE
    )
  }
  if (small_sample_flag) {
    warning(sprintf(
      "EBEC: tail size k = %d (q = %.3f, n = %d) is below the recommended minimum of %d. Results may be unstable; consider increasing q.",
      k, q, n, min_k), call. = FALSE)
  }

  share <- function(yy, pp) {
    loss <- sort((yy - pp)^2, decreasing = TRUE)
    total <- sum(loss)
    if (total == 0) return(NA_real_)
    sum(loss[seq_len(k)]) / total
  }

  observed <- share(y, p)

  if (!is.null(seed)) set.seed(seed)
  null <- replicate(B, share(rbinom(n, 1, p), p))
  null <- null[!is.na(null)]
  expected <- mean(null)
  excess <- observed - expected
  p_value <- (1 + sum(null >= observed)) / (length(null) + 1)

  flags <- .ebec_flags(null, p_value, "EBEC")

  if (ci) {
    if (!is.null(seed)) set.seed(seed + 1)
    boot_excess <- replicate(B, {
      idx <- sample(seq_len(n), n, replace = TRUE)
      share(y[idx], p[idx]) - expected
    })
    ci_v <- unname(quantile(boot_excess[!is.na(boot_excess)], c(0.025, 0.975)))
  } else {
    ci_v <- c(NA_real_, NA_real_)
  }

  structure(
    list(observed = observed, expected = expected, excess = excess,
         excess_norm = excess / (1 - q),
         ci_lower = ci_v[1], ci_upper = ci_v[2],
         ci_lower_norm = ci_v[1] / (1 - q), ci_upper_norm = ci_v[2] / (1 - q),
         p_value = p_value, q = q, B = B, n = n, k = k,
         small_sample_flag = small_sample_flag,
         low_resolution_flag = flags$low_resolution_flag,
         mc_se_flag = flags$mc_se_flag, mc_se = flags$mc_se),
    class = "ebec_share")
}

#' Excess Brier Error Concentration for False Negatives (EBEC_FN)
#'
#' Among the observed events (y = 1), tests whether the worst q fraction of
#' false-negative errors, measured as (1 - p)^2, concentrates more of the
#' total false-negative loss than expected if the predicted probabilities were
#' the true ones. The null distribution simulates the whole outcome vector
#' from Bernoulli(p) and re-selects the events in every simulation, so the
#' selection of the y = 1 subset is reproduced under the null.
#'
#' @inheritParams ebec_share
#' @param q Fraction of worst cases examined within the y = 1 subset.
#'   Default 0.05.
#' @return An object of class \code{ebec_share_fn}.
#' @importFrom stats rbinom quantile
#' @examples
#' set.seed(1)
#' p <- runif(200, 0.2, 0.8)
#' y <- rbinom(200, 1, p)
#' ebec_share_fn(y, p, q = 0.10, B = 200, seed = 1)
#' @export
ebec_share_fn <- function(y, p, q = 0.05, B = 1000, seed = NULL, min_k = 10, ci = TRUE) {
  .ebec_cond(y, p, q, B, seed, min_k, type = "FN", q_missing = missing(q), ci = ci)
}

#' Excess Brier Error Concentration for False Positives (EBEC_FP)
#'
#' Among the observed non-events (y = 0), tests whether the worst q fraction of
#' false-positive errors, measured as p^2, concentrates more of the total
#' false-positive loss than expected if the predicted probabilities were the
#' true ones. The null distribution simulates the whole outcome vector from
#' Bernoulli(p) and re-selects the non-events in every simulation.
#'
#' @inheritParams ebec_share
#' @param q Fraction of worst cases examined within the y = 0 subset.
#'   Default 0.05.
#' @return An object of class \code{ebec_share_fp}.
#' @importFrom stats rbinom quantile
#' @examples
#' set.seed(1)
#' p <- runif(200, 0.2, 0.8)
#' y <- rbinom(200, 1, p)
#' ebec_share_fp(y, p, q = 0.10, B = 200, seed = 1)
#' @export
ebec_share_fp <- function(y, p, q = 0.05, B = 1000, seed = NULL, min_k = 10, ci = TRUE) {
  .ebec_cond(y, p, q, B, seed, min_k, type = "FP", q_missing = missing(q), ci = ci)
}

# Shared engine for ebec_fn() and ebec_fp()
.ebec_cond <- function(y, p, q, B, seed, min_k, type, q_missing, ci = TRUE) {
  stopifnot(length(y) == length(p))
  stopifnot(all(y %in% c(0, 1)))
  stopifnot(all(p >= 0 & p <= 1))
  stopifnot(q > 0, q < 1, B >= 1)

  target <- if (type == "FN") 1 else 0
  loss_fun <- if (type == "FN") function(pp) (1 - pp)^2 else function(pp) pp^2
  nm <- paste0("EBEC_", type)
  n <- length(y)

  # share of the loss of the worst ceiling(q * m) cases among the m selected
  share <- function(pp_sel) {
    m <- length(pp_sel)
    if (m == 0) return(NA_real_)
    loss <- sort(loss_fun(pp_sel), decreasing = TRUE)
    total <- sum(loss)
    if (total == 0) return(NA_real_)
    sum(loss[seq_len(ceiling(q * m))]) / total
  }

  sel <- which(y == target)
  n_sel <- length(sel)
  if (n_sel == 0) {
    warning(nm, ": no ", if (type == "FN") "events (y = 1)" else "non-events (y = 0)",
            " found. Returning NA.", call. = FALSE)
    return(structure(list(observed = NA_real_, expected = NA_real_, excess = NA_real_,
                          excess_norm = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                          ci_lower_norm = NA_real_, ci_upper_norm = NA_real_,
                          p_value = NA_real_, q = q, B = B, n = 0L, k = 0L,
                          small_sample_flag = FALSE, low_resolution_flag = FALSE,
                          mc_se_flag = FALSE, mc_se = NA_real_),
                     class = paste0("ebec_share_", tolower(type))))
  }
  k <- ceiling(q * n_sel)
  small_sample_flag <- k < min_k

  if (q_missing) {
    warning(nm, ": you are using q = 0.05 (default value). Choose q according to your context.",
            call. = FALSE)
  }
  if (small_sample_flag) {
    warning(sprintf(
      "%s: tail size k = %d (q = %.3f, subset size = %d) is below the recommended minimum of %d. Results may be unstable; consider increasing q.",
      nm, k, q, n_sel, min_k), call. = FALSE)
  }

  observed <- share(p[sel])

  # Null: simulate the whole outcome vector, then re-select the subset.
  if (!is.null(seed)) set.seed(seed)
  null <- replicate(B, {
    ys <- rbinom(n, 1, p)
    share(p[ys == target])
  })
  null <- null[!is.na(null)]
  expected <- mean(null)
  excess <- observed - expected
  p_value <- (1 + sum(null >= observed)) / (length(null) + 1)

  flags <- .ebec_flags(null, p_value, nm)

  if (ci) {
    if (!is.null(seed)) set.seed(seed + 1)
    boot_excess <- replicate(B, {
      idx <- sample(seq_len(n_sel), n_sel, replace = TRUE)
      share(p[sel][idx]) - expected
    })
    ci_v <- unname(quantile(boot_excess[!is.na(boot_excess)], c(0.025, 0.975)))
  } else {
    ci_v <- c(NA_real_, NA_real_)
  }

  structure(
    list(observed = observed, expected = expected, excess = excess,
         excess_norm = excess / (1 - q),
         ci_lower = ci_v[1], ci_upper = ci_v[2],
         ci_lower_norm = ci_v[1] / (1 - q), ci_upper_norm = ci_v[2] / (1 - q),
         p_value = p_value, q = q, B = B, n = n_sel, k = k,
         small_sample_flag = small_sample_flag,
         low_resolution_flag = flags$low_resolution_flag,
         mc_se_flag = flags$mc_se_flag, mc_se = flags$mc_se),
    class = paste0("ebec_share_", tolower(type)))
}

.ebec_flags <- function(null, p_value, nm) {
  n_unique_null <- length(unique(round(null, 6)))
  low_resolution_flag <- n_unique_null < 20
  mc_se <- sqrt(p_value * (1 - p_value) / length(null))
  mc_se_flag <- abs(p_value - 0.05) < 2 * mc_se
  if (low_resolution_flag) {
    warning(sprintf(
      "%s: the null distribution has only %d distinct values - p-value resolution may be too coarse. Consider increasing q or n.",
      nm, n_unique_null), call. = FALSE)
  }
  if (mc_se_flag) {
    warning(sprintf(
      "%s: the p-value (%.4f) is within 2 Monte Carlo SE (%.4f) of 0.05, so the decision at alpha = 0.05 could change. Increase B.",
      nm, p_value, mc_se), call. = FALSE)
  }
  list(low_resolution_flag = low_resolution_flag, mc_se_flag = mc_se_flag, mc_se = mc_se)
}

.print_ebec <- function(x, title, n_label, tag, digits = 4) {
  w <- 50
  cat(title, "\n", sep = "")
  cat(strrep("-", w), "\n", sep = "")
  cat(sprintf("  %s = %d   worst %.1f%% tail  ->  k = %d observations\n",
              n_label, x$n, 100 * x$q, x$k))
  cat(strrep("-", w), "\n", sep = "")
  cat(sprintf("  %-30s %8.*f\n", "Observed concentration:", digits, x$observed))
  cat(sprintf("  %-30s %8.*f\n", "Expected concentration (H0):", digits, x$expected))
  cat(sprintf("  %-30s %8.*f\n", paste0(tag, ":"), digits, x$excess))
  cat(sprintf("  %-30s [%.*f, %.*f]\n", "  95% CI (boot.):",
              digits, x$ci_lower, digits, x$ci_upper))
  cat(sprintf("  %-30s %8.*f\n", paste0(tag, " normalized:"), digits, x$excess_norm))
  cat(sprintf("  %-30s [%.*f, %.*f]\n", "  95% CI (normalized):",
              digits, x$ci_lower_norm, digits, x$ci_upper_norm))
  cat(sprintf("  %-30s %8.*f\n", "p-value (MC):", digits, x$p_value))
  cat(strrep("-", w), "\n", sep = "")
  cat(sprintf("  Monte Carlo: %d null simulations\n", x$B))
  if (isTRUE(x$small_sample_flag)) cat("  [!] Small tail size (k < 10) - interpret with caution.\n")
  if (isTRUE(x$low_resolution_flag)) cat("  [!] Low null-distribution resolution - p-value may be coarse.\n")
  if (isTRUE(x$mc_se_flag)) cat(sprintf("  [!] p-value within 2 MC SE (%.4f) of 0.05; increase B.\n", x$mc_se))
  invisible(x)
}

#' @export
print.ebec_share <- function(x, digits = 4, ...) {
  .print_ebec(x, "Tail loss share index (legacy, see ebec())", "n", "EBEC", digits)
}

#' @export
print.ebec_share_fn <- function(x, digits = 4, ...) {
  .print_ebec(x, "Tail loss share index, false negatives (legacy, descriptive)",
              "n_events", "EBEC_FN", digits)
}

#' @export
print.ebec_share_fp <- function(x, digits = 4, ...) {
  .print_ebec(x, "Tail loss share index, false positives (legacy, descriptive)",
              "n_nonevents", "EBEC_FP", digits)
}
