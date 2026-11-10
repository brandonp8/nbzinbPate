# Regression tests for a bug first seen on macOS CI (4 failing tests there,
# all passing on Windows/Linux):
#
# Near the optimum a correct IRLS/Newton step changes the log-likelihood by
# less than its floating-point rounding error. The step-halving safeguards
# used to compare `ll_new < ll_old` exactly, so a rounding wobble of -1e-13
# (which macOS's math library produced, Windows's did not) made them halve
# a correct step 30 times; the parameters then stopped moving and the
# "change < tol" rule reported a false convergence.
#
# To reproduce this on any platform, we replace the log-likelihood
# functions with versions that add +/- 1e-12 of random noise -- the size of
# real rounding error -- and check that every algorithm still converges to
# the true optimum.

test_that("ll_decreased ignores rounding-level changes only", {
  expect_false(ll_decreased(-1561 - 1e-12, -1561))   # rounding noise
  expect_false(ll_decreased(-1561 + 1e-3, -1561))    # an increase
  expect_true(ll_decreased(-1562, -1561))            # a real decrease
  expect_true(ll_decreased(NaN, -1561))
  expect_true(ll_decreased(-Inf, -1561))
})

noisy <- function(f) {
  force(f)
  function(...) f(...) + stats::runif(1, -1e-12, 1e-12)
}

test_that("NB IRLS converges to a zero score despite rounding noise", {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  X <- model.matrix(~ fem + mar + kid5 + phd + ment, d)
  local_mocked_bindings(nb_or_pois_loglik = noisy(nb_or_pois_loglik))
  set.seed(61)
  for (rep in 1:5) {
    fit <- nb_irls_beta(X, d$art, theta = 2.264)
    expect_true(fit$converged)
    score <- crossprod(X, 2.264 * (d$art - fit$mu) / (2.264 + fit$mu))
    expect_lt(max(abs(score)), 1e-6)
  }
})

test_that("Newton for theta reaches the root despite rounding noise", {
  set.seed(62)
  y  <- rnbinom(500, size = 2, mu = 4)
  mu <- rep(mean(y), 500)
  local_mocked_bindings(nb_loglik = noisy(nb_loglik))
  for (rep in 1:5) {
    fit <- nb_newton_theta(y, mu, theta_start = 0.5)
    expect_true(fit$converged)
    expect_lt(abs(fit$U), 1e-6)
  }
})

test_that("logistic IRLS reaches a zero score despite rounding noise", {
  set.seed(63)
  Z   <- cbind(1, rnorm(300), rbinom(300, 1, 0.4))
  tau <- runif(300) * rbinom(300, 1, 0.5)
  local_mocked_bindings(logit_frac_loglik = noisy(logit_frac_loglik))
  for (rep in 1:5) {
    fit <- logit_irls(Z, tau)
    expect_true(fit$converged)
    expect_lt(max(abs(crossprod(Z, tau - fit$pi))), 1e-8)
  }
})

test_that("theta is capped at the Poisson limit and nbreg still converges", {
  # Underdispersed counts: the NB likelihood increases as theta -> Inf.
  set.seed(64)
  d <- data.frame(x = rnorm(300))
  d$y <- rbinom(300, size = 4, prob = 0.5)
  expect_warning(fit <- nbreg(y ~ x, data = d), "Poisson")
  expect_true(fit$converged)
  # theta ends somewhere in the Poisson limit (at or below the 1e8 cap):
  # out there U(theta) is pure rounding noise, so Newton stops once its
  # steps fall below the tolerance, not necessarily exactly at the cap.
  expect_gt(fit$theta, 1e5)
  expect_lte(fit$theta, 1e8)
  expect_lt(fit$iter, 20)
})
