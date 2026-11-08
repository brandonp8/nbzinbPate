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
  expect_equal(zinb_loglik(y, mu, pi, theta), ref, tolerance = 1e-12)
})

test_that("zinb_loglik reduces to the NB log-likelihood when pi = 0", {
  y  <- c(0, 1, 4, 0, 2)
  mu <- c(1, 2, 3, 0.5, 2)
  expect_equal(zinb_loglik(y, mu, rep(0, 5), 3), nb_loglik(y, mu, 3))
})
