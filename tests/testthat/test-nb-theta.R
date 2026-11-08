# Tests for nb_theta_derivs() and nb_newton_theta() (Sec 3.1(b)).
# MASS::theta.ml() is used ONLY as a reference here.

# Central finite difference f'(x) ~ [f(x+h) - f(x-h)] / (2h): an
# independent check of the analytic derivatives, using base R only.
num_deriv <- function(f, x, h = 1e-5 * max(1, abs(x))) {
  (f(x + h) - f(x - h)) / (2 * h)
}

test_that("U(theta) and H(theta) match finite differences", {
  set.seed(4)
  y  <- rnbinom(100, size = 1.5, mu = 3)
  mu <- runif(100, 1, 5)
  a  <- runif(100)
  for (theta in c(0.4, 1.5, 8)) {
    d <- nb_theta_derivs(y, mu, theta, a)
    # U = d loglik / d theta
    expect_equal(d$U, num_deriv(function(t) nb_loglik(y, mu, t, a), theta),
                 tolerance = 1e-6)
    # H = d U / d theta
    expect_equal(d$H,
                 num_deriv(function(t) nb_theta_derivs(y, mu, t, a)$U, theta),
                 tolerance = 1e-6)
  }
})

test_that("Newton on log(theta) finds the root of U and matches theta.ml", {
  skip_if_not_installed("MASS")
  set.seed(5)
  y  <- rnbinom(500, size = 2, mu = 4)
  mu <- rep(mean(y), 500)
  ref <- MASS::theta.ml(y, mu, limit = 100, eps = 1e-12)

  for (start in c(0.05, 1, 50)) {          # from far below and far above
    fit <- nb_newton_theta(y, mu, theta_start = start)
    expect_true(fit$converged)
    expect_lt(abs(fit$U), 1e-6)            # score is zero at the MLE
    expect_equal(fit$theta, as.numeric(ref), tolerance = 1e-8)
    # SE(theta) = (-H)^(-1/2) (Sec 3.1(d)) equals theta.ml's SE.
    expect_equal(1 / sqrt(-fit$H), attr(ref, "SE"), tolerance = 1e-6)
  }
})

test_that("prior weights enter the theta update", {
  skip_if_not_installed("MASS")
  set.seed(6)
  y  <- rnbinom(300, size = 1, mu = 2)
  mu <- rep(2, 300)
  a  <- runif(300)
  fit <- nb_newton_theta(y, mu, 1, weights = a)
  ref <- MASS::theta.ml(y, mu, weights = a, limit = 100, eps = 1e-12)
  expect_equal(fit$theta, as.numeric(ref), tolerance = 1e-8)
})

test_that("Poisson-like data push theta upward without failing", {
  set.seed(7)
  y  <- rpois(300, 3)                       # no overdispersion: MLE theta -> Inf
  mu <- rep(mean(y), 300)
  fit <- nb_newton_theta(y, mu, 1, maxit = 50)
  expect_gt(fit$theta, 100)
  expect_true(is.finite(fit$loglik))
})
