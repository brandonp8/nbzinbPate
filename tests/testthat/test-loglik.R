# Tests for the log-likelihood helpers in R/loglik.R.
# stats::dnbinom() is the independent reference: in R, dnbinom(size = theta,
# mu = mu) is exactly the NB2 pmf of Sec 2.1.

test_that("nb_logdens matches dnbinom(log = TRUE)", {
  y  <- c(0, 1, 2, 5, 10, 40)
  mu <- c(0.5, 1, 3, 4, 12, 30)
  for (theta in c(0.3, 1, 2.5, 50)) {
    expect_equal(nb_logdens(y, mu, theta),
                 dnbinom(y, size = theta, mu = mu, log = TRUE),
                 tolerance = 1e-12)
  }
})

test_that("nb_logdens stays accurate for huge theta (Poisson limit)", {
  # As theta -> Inf the NB2 pmf tends to the Poisson pmf.
  # Note: for astronomically large theta (~1e10) lgamma(y + theta) -
  # lgamma(theta) subtracts two numbers of size ~1e11, so double precision
  # leaves an absolute error ~1e-5 (MASS::glm.nb uses the same formula and
  # has the same limit). theta = 1e6 is already deep in the Poisson regime.
  y  <- c(0, 1, 3, 7)
  mu <- c(0.7, 2, 3, 6)
  expect_equal(nb_logdens(y, mu, 1e6), dpois(y, mu, log = TRUE),
               tolerance = 1e-5)
  expect_equal(nb_logdens(y, mu, 1e6),
               dnbinom(y, size = 1e6, mu = mu, log = TRUE), tolerance = 1e-8)
})

test_that("nb_loglik applies prior weights", {
  y  <- c(0, 2, 4)
  mu <- c(1, 2, 3)
  a  <- c(0.2, 1, 0.5)
  expect_equal(nb_loglik(y, mu, 2, a),
               sum(a * dnbinom(y, size = 2, mu = mu, log = TRUE)))
  # Default weights are all 1.
  expect_equal(nb_loglik(y, mu, 2),
               sum(dnbinom(y, size = 2, mu = mu, log = TRUE)))
})

test_that("zinb_loglik matches the mixture formula of Sec 2.2", {
  y     <- c(0, 0, 1, 3, 0, 8)
  mu    <- c(1, 2.5, 1.5, 3, 0.4, 6)
  pi    <- c(0.1, 0.3, 0.2, 0.05, 0.6, 0.4)
  theta <- 1.7
  f  <- dnbinom(y, size = theta, mu = mu)
  # P(Y=0) = pi + (1-pi) f(0);  P(Y=y) = (1-pi) f(y) for y > 0.
  ref <- sum(log(ifelse(y == 0, pi + (1 - pi) * f, (1 - pi) * f)))
  # zinb_loglik takes the zero-model linear predictor eta = logit(pi).
  expect_equal(zinb_loglik(y, mu, qlogis(pi), theta), ref, tolerance = 1e-12)
})

test_that("zinb_loglik reduces to the NB log-likelihood when pi -> 0", {
  y  <- c(0, 1, 4, 0, 2)
  mu <- c(1, 2, 3, 0.5, 2)
  # eta = -800 means pi = 1/(1 + e^800), which is 0 to machine precision.
  expect_equal(zinb_loglik(y, mu, rep(-800, 5), 3), nb_loglik(y, mu, 3))
})

test_that("zinb_loglik stays finite and correct at extreme eta", {
  y  <- c(0, 0, 3, 5)
  mu <- c(1, 2, 2, 4)
  theta <- 2
  # pi rounds to exactly 1 for eta >= 37, so the naive log(1 - pi) is -Inf.
  # The stable version must instead give log(1 - pi) = -eta (to 1e-16).
  eta <- c(-800, 50, -800, 40)
  ll <- zinb_loglik(y, mu, eta, theta)
  expect_true(is.finite(ll))
  ref <- log(dnbinom(0, size = theta, mu = 1)) +   # pi ~ 0: NB zero
    0 +                                            # pi ~ 1: log(pi + ...) ~ 0
    dnbinom(3, size = theta, mu = 2, log = TRUE) + # pi ~ 0
    (-40 + dnbinom(5, size = theta, mu = 4, log = TRUE))  # log(1-pi) ~ -40
  expect_equal(ll, ref, tolerance = 1e-10)
})
