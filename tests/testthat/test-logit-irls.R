# Tests for logit_irls() (Sec 3.2 zero-part M-step).
# glm() is used ONLY as a reference here.

make_Z <- function(n, seed) {
  set.seed(seed)
  cbind("(Intercept)" = 1, z1 = rnorm(n), z2 = rbinom(n, 1, 0.4))
}

test_that("0/1 response reproduces logistic glm()", {
  Z <- make_Z(400, 21)
  tau <- rbinom(400, 1, plogis(-0.5 + 0.8 * Z[, 2] - 0.6 * Z[, 3]))
  fit <- logit_irls(Z, tau)
  ref <- glm(tau ~ Z - 1, family = binomial,
             control = glm.control(epsilon = 1e-14, maxit = 100))
  expect_true(fit$converged)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(fit$loglik, as.numeric(logLik(ref)), tolerance = 1e-10)
})

test_that("fractional response reproduces quasibinomial glm()", {
  # glm(family = quasibinomial) maximizes the same Q(gamma) for responses
  # in [0, 1] (it only rescales the standard errors), so its estimates are
  # an independent reference for the fractional case.
  Z <- make_Z(400, 22)
  tau <- plogis(-1 + 0.7 * Z[, 2] + 0.3 * Z[, 3] + rnorm(400, sd = 0.5))
  tau[sample(400, 150)] <- 0                  # EM gives tau = 0 when y > 0
  fit <- logit_irls(Z, tau)
  ref <- glm(tau ~ Z - 1, family = quasibinomial,
             control = glm.control(epsilon = 1e-14, maxit = 100))
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
})

test_that("the score Z'(tau - pi) is zero at the solution", {
  Z <- make_Z(300, 23)
  tau <- runif(300) * rbinom(300, 1, 0.5)
  fit <- logit_irls(Z, tau)
  expect_lt(max(abs(crossprod(Z, tau - fit$pi))), 1e-8)
})

test_that("canonical link: Newton converges in few iterations", {
  Z <- make_Z(500, 24)
  tau <- runif(500)
  fit <- logit_irls(Z, tau)
  expect_true(fit$converged)
  expect_lte(fit$iter, 10)
})

test_that("a warm start converges to the same answer", {
  Z <- make_Z(300, 25)
  tau <- runif(300)
  cold <- logit_irls(Z, tau)
  warm <- logit_irls(Z, tau, gamma_start = cold$coefficients + 0.2)
  expect_equal(warm$coefficients, cold$coefficients, tolerance = 1e-9)
})

test_that("intercept-only: pi-hat equals the mean of tau", {
  tau <- c(0, 0, 0.3, 0.9, 1, 0.2)
  fit <- logit_irls(matrix(1, 6, 1), tau)
  expect_equal(plogis(fit$coefficients[[1]]), mean(tau), tolerance = 1e-10)
})
