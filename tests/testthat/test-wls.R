# Tests for the weighted least-squares helpers in R/wls.R.
# stats::lm.wfit() is used here ONLY as an independent reference; it is
# never called from R/.

test_that("wls_solve matches lm.wfit and the normal equations", {
  set.seed(1)
  n <- 50
  X <- cbind("(Intercept)" = 1, x1 = rnorm(n), x2 = runif(n))
  z <- rnorm(n)
  w <- rexp(n)

  beta <- wls_solve(X, w, z)
  expect_equal(unname(beta), unname(lm.wfit(X, z, w)$coefficients),
               tolerance = 1e-10)
  # Same as solving (X'WX) beta = X'Wz directly.
  expect_equal(unname(beta),
               unname(drop(solve(crossprod(X, w * X), crossprod(X, w * z)))),
               tolerance = 1e-10)
  expect_named(beta, colnames(X))
})

test_that("wls_solve tolerates zero weights", {
  # Zero weights occur in the ZINB count M-step (a_i = 1 - tau_i = 0).
  set.seed(2)
  X <- cbind(1, rnorm(20))
  z <- rnorm(20)
  w <- c(rep(0, 5), rep(1, 15))
  expect_equal(unname(wls_solve(X, w, z)),
               unname(lm.wfit(X[6:20, ], z[6:20], w[6:20])$coefficients),
               tolerance = 1e-10)
})

test_that("wls_solve errors informatively on a rank-deficient design", {
  X <- cbind(1, 1:10, 2 * (1:10))          # third column = 2 * second
  expect_error(wls_solve(X, rep(1, 10), rnorm(10)), "rank deficient")
})

test_that("xtwx_inverse equals solve(X'WX)", {
  set.seed(3)
  X <- cbind(a = 1, b = rnorm(30))
  w <- rexp(30)
  V <- xtwx_inverse(X, w)
  expect_equal(V, solve(crossprod(X, w * X)), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(dimnames(V), list(c("a", "b"), c("a", "b")))
})
