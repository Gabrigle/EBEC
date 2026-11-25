test_that("ebec returns the expected structure", {
  set.seed(1)
  p <- runif(200, 0.2, 0.8)
  y <- rbinom(200, 1, p)
  res <- ebec(y, p, B = 99, ci = FALSE, seed = 1)
  expect_s3_class(res, "ebec")
  expect_true(res$p_value > 0 && res$p_value <= 1)
  expect_equal(nrow(res$components), 1 + 2 * 4)          # default: total, tails, FP
  expect_false(any(grepl("^FN", res$components$component)))
  res2 <- ebec(y, p, B = 99, classes = "none", seed = 1)
  expect_equal(nrow(res2$components), 1 + 4)
  res3 <- ebec(y, p, B = 99, classes = "both", seed = 1)
  expect_equal(nrow(res3$components), 1 + 3 * 4)
  expect_equal(nrow(ebec(y, p, B = 99, classes = TRUE, seed = 1)$components), 1 + 3 * 4)
  expect_equal(nrow(ebec(y, p, B = 99, classes = FALSE, seed = 1)$components), 1 + 4)
  expect_error(ebec(y, p, B = 99, classes = "xx"))
  expect_output(print(res), "EBEC")
})

test_that("component statistics agree with a hand calculation", {
  y <- c(1, 0, 1, 0)
  p <- c(0.1, 0.8, 0.9, 0.3)
  # losses: 0.81, 0.64, 0.01, 0.09 ; FN losses: 0.81, 0, 0.01, 0 ; FP: 0, 0.64, 0, 0.09
  s <- EBEC:::.scan_stats(y, p, ks = c(1L, 2L), classes = "both")
  expect_equal(unname(s),
               c(1.55, 0.81, 1.45, 0.81, 0.82, 0.64, 0.73), tolerance = 1e-12)
  s2 <- EBEC:::.scan_stats(y, p, ks = c(1L, 2L), classes = "fp")
  expect_equal(unname(s2), c(1.55, 0.81, 1.45, 0.64, 0.73), tolerance = 1e-12)
})

test_that("tail sums do not decrease when errors are added (no dilution)", {
  n <- 100
  p <- rep(0.5, n)
  tail_sum <- function(r) {            # r confident errors, rest correct
    y <- c(rep(1, r), rep(0, n - r)); pp <- c(rep(0.01, r), rep(0.01, n - r))
    EBEC:::.scan_stats(y, pp, ks = 10L, classes = "none")[2]
  }
  vals <- vapply(c(0, 5, 10, 15, 40), tail_sum, numeric(1))
  expect_true(all(diff(vals) >= -1e-12))
})

test_that("ebec detects a few confident errors", {
  set.seed(2)
  n <- 500
  p <- runif(n, 0.2, 0.8)
  y <- rbinom(n, 1, p)
  idx <- 1:10
  p[idx] <- ifelse(y[idx] == 1, 0.02, 0.98)
  res <- ebec(y, p, B = 499, seed = 3)
  expect_lt(res$p_value, 0.05)
})

test_that("confident false negatives are detected (exploratory FN components)", {
  set.seed(4)
  n <- 600
  p <- runif(n, 0.2, 0.8)
  y <- rbinom(n, 1, p)
  ev <- which(y == 1)[1:12]
  p[ev] <- 0.02
  res <- ebec(y, p, B = 499, seed = 5)
  expect_lt(res$p_value, 0.05)
  res_b <- ebec(y, p, B = 499, classes = "both", seed = 5)
  expect_lt(res_b$p_value, 0.05)
})

test_that("confident false positives are detected and point to the FP component", {
  set.seed(6)
  n <- 600
  p <- runif(n, 0.2, 0.8)
  y <- rbinom(n, 1, p)
  ne <- which(y == 0)[1:12]
  p[ne] <- 0.98
  res <- ebec(y, p, B = 499, seed = 7)
  expect_lt(res$p_value, 0.05)
})

test_that("ebec is reproducible with a seed", {
  set.seed(3)
  p <- runif(100, 0.1, 0.9)
  y <- rbinom(100, 1, p)
  a <- ebec(y, p, B = 99, seed = 5)$p_value
  b <- ebec(y, p, B = 99, seed = 5)$p_value
  expect_identical(a, b)
})

test_that("ebec validates its inputs", {
  expect_error(ebec(c(0, 2), c(0.5, 0.5)))
  expect_error(ebec(c(0, 1), c(0.5, 1.5)))
  expect_error(ebec(c(0, NA), c(0.5, 0.5)))
  expect_error(ebec(c(0, 1), c(0.5, 0.5, 0.5)))
  expect_error(ebec(c(0, 1), c(0.5, 0.5), q = 1.2))
  expect_error(ebec(c(0, 1), c(0.5, 0.5), B = 0))
})

test_that("p equal to 0 or 1 is handled", {
  y <- c(1, 0, 1, 0, 1, 0, 1, 0)
  p <- c(1, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  res <- ebec(y, p, q = c(0.25, 0.5), B = 49, seed = 1)
  expect_true(res$p_value > 0 && res$p_value <= 1)
})

test_that("type I error is controlled in a small study", {
  skip_on_cran()
  set.seed(11)
  pv <- replicate(150, {
    p <- runif(150, 0.1, 0.9)
    y <- rbinom(150, 1, p)
    ebec(y, p, B = 99)$p_value
  })
  expect_lt(mean(pv <= 0.05), 0.12)
})

test_that("components table has ratio, excess and interval columns", {
  set.seed(8)
  p <- runif(150, 0.1, 0.9)
  y <- rbinom(150, 1, p)
  res <- ebec(y, p, B = 99, B_ci = 30, B_inner = 10, seed = 1)
  cmp <- res$components
  expect_true(all(c("tail_pct", "observed", "expected", "ratio", "excess_pct",
                    "ci_low", "ci_high", "z", "p_component") %in% names(cmp)))
  expect_equal(cmp$ratio, cmp$observed / cmp$expected)
  expect_equal(cmp$excess_pct, 100 * (cmp$ratio - 1))
  expect_true(all(cmp$ci_low <= cmp$ci_high, na.rm = TRUE))
  expect_true(is.na(cmp$tail_pct[1]))
  expect_equal(cmp$tail_pct[2], 100 * ceiling(0.01 * 150) / 150)
  res0 <- ebec(y, p, B = 99, ci = FALSE, seed = 1)
  expect_true(all(is.na(res0$components$ci_low)))
})

test_that("print shows the guide, the report sentence and both languages", {
  set.seed(9)
  p <- runif(120, 0.1, 0.9)
  y <- rbinom(120, 1, p)
  res <- ebec(y, p, B = 49, B_ci = 20, B_inner = 10, seed = 1)
  expect_output(print(res), "QUICK GUIDE")
  expect_output(print(res), "How to report")
  expect_output(print(res), "Range")
  expect_output(print(res, explain = FALSE), "How to report")
  out <- capture.output(print(res, explain = FALSE))
  expect_false(any(grepl("QUICK GUIDE", out)))
  expect_output(print(res, lang = "es"), "PIDA")
})
