# Tests for nb_irls_beta() (Sec 3.1(a): IRLS for beta with theta fixed).
# glm() and MASS::negative.binomial() are used ONLY as references here.

bio_design <- function() {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  list(X = model.matrix(~ fem + mar + kid5 + phd + ment, data = d),
       y = d$art)
}

test_that("theta = Inf reproduces Poisson glm()", {
  b <- bio_design()
  fit <- nb_irls_beta(b$X, b$y, theta = Inf)
  ref <- glm(b$y ~ b$X - 1, family = poisson,
             control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_true(fit$converged)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(fit$loglik, as.numeric(logLik(ref)), tolerance = 1e-8)
})

test_that("fixed theta reproduces glm() with the NB family at that theta", {
  skip_if_not_installed("MASS")
  b <- bio_design()
  theta <- 2.264                       # roughly glm.nb's theta-hat on these data
  fit <- nb_irls_beta(b$X, b$y, theta = theta)
  ref <- glm(b$y ~ b$X - 1, family = MASS::negative.binomial(theta),
             control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_true(fit$converged)

  # Exact characterization of the MLE for fixed theta: the score
  # U(beta) = X' [a_i theta (y_i - mu_i) / (theta + mu_i)] is zero.
  score <- crossprod(b$X, theta * (b$y - fit$mu) / (theta + fit$mu))
  expect_lt(max(abs(score)), 1e-6)

  # glm() stops on a relative-deviance criterion that, even with
  # epsilon = 1e-12, leaves its score around 1e-3 here, so it is only
  # accurate to ~1e-6 in beta. Compare at the project's Table 4 tolerance.
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-5)
  # Working weights at convergence match glm's final IRLS weights.
  expect_equal(unname(fit$w), unname(ref$weights), tolerance = 1e-5)
})

test_that("offset and prior weights are handled", {
  skip_if_not_installed("MASS")
  set.seed(10)
  n <- 200
  x <- rnorm(n)
  expo <- runif(n, 1, 5)
  y <- rnbinom(n, size = 2, mu = expo * exp(0.3 + 0.5 * x))
  a <- runif(n)                        # fractional prior weights, as in EM
  X <- cbind("(Intercept)" = 1, x = x)

  fit <- nb_irls_beta(X, y, theta = 2, offset = log(expo), weights = a)
  ref <- suppressWarnings(              # glm warns about non-integer weights*y
    glm(y ~ x, offset = log(expo), weights = a,
        family = MASS::negative.binomial(2),
        control = glm.control(epsilon = 1e-12, maxit = 100)))
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-8)
})

test_that("a warm start converges to the same answer in fewer iterations", {
  b <- bio_design()
  cold <- nb_irls_beta(b$X, b$y, theta = 2)
  warm <- nb_irls_beta(b$X, b$y, theta = 2,
                       beta_start = cold$coefficients + 0.01)
  expect_equal(warm$coefficients, cold$coefficients, tolerance = 1e-8)
  expect_lte(warm$iter, cold$iter)
})
