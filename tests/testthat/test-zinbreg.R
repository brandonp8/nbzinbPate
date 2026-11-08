# Tests for zinbreg() and the EM engine (Sec 3.2).
# pscl::zeroinfl() is used ONLY as a reference here.

max_abs_diff <- function(a, b) max(abs(unname(a) - unname(b)))

ref_zeroinfl <- function(formula, data) {
  pscl::zeroinfl(formula, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12))
}

# Table 4 checks for coefficients, log(theta), and log-likelihood.
# (Standard errors are added in step 7.)
expect_matches_zeroinfl <- function(fit, ref) {
  expect_true(fit$converged)
  expect_lt(max_abs_diff(fit$coefficients$count, ref$coefficients$count), 1e-3)
  expect_lt(max_abs_diff(fit$coefficients$zero, ref$coefficients$zero), 1e-3)
  expect_lt(abs(log(fit$theta) - log(ref$theta)), 1e-3)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-3)
}

test_that("zinbreg matches zeroinfl on bioChemists (Table 4)", {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  fit <- zinbreg(f, data = d)
  expect_matches_zeroinfl(fit, ref_zeroinfl(f, d))
})

test_that("zinbreg matches zeroinfl on NMES1988 (Table 4)", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  fit <- zinbreg(f, data = NMES1988)
  expect_matches_zeroinfl(fit, ref_zeroinfl(f, NMES1988))
})

test_that("EM never decreases the observed-data log-likelihood", {
  skip_if_not_installed("pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = pscl::bioChemists)
  # Allow only floating-point rounding (relative 1e-12).
  expect_true(all(diff(fit$loglik_path) > -1e-12 * abs(fit$loglik)))
  expect_equal(fit$loglik, max(fit$loglik_path))
})

test_that("a formula without | uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  a <- zinbreg(art ~ fem + ment, data = d)
  b <- zinbreg(art ~ fem + ment | fem + ment, data = d)
  expect_equal(a$coefficients, b$coefficients)
  expect_named(a$coefficients$zero, c("(Intercept)", "femWomen", "ment"))
})

test_that("different covariates in the two parts", {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  fit <- zinbreg(art ~ fem + kid5 + ment | ment, data = d)
  expect_matches_zeroinfl(fit, ref_zeroinfl(art ~ fem + kid5 + ment | ment, d))
  expect_named(fit$coefficients$zero, c("(Intercept)", "ment"))
  expect_equal(fit$df, 4L + 2L + 1L)
  expect_false(fit$boundary)
})

test_that("no zero inflation: boundary is detected, loglik equals NB", {
  # With an intercept-only zero model these data show no zero inflation:
  # the likelihood increases as pi -> 0 (gamma -> -Inf), so the MLE is on
  # the boundary. (zeroinfl stops at gamma ~ -14.5 only because BFGS stops
  # seeing improvement.) EM should stop, flag it, warn, and reach the NB
  # log-likelihood, which is the supremum.
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  expect_warning(fit <- zinbreg(art ~ fem + kid5 + ment | 1, data = d),
                 "no zero inflation")
  expect_true(fit$converged)
  expect_true(fit$boundary)
  nb <- nbreg(art ~ fem + kid5 + ment, data = d)
  expect_lt(abs(fit$loglik - nb$loglik), 1e-6)
  expect_lt(max_abs_diff(fit$coefficients$count, nb$coefficients), 1e-4)
  # At least as high as zeroinfl's (which stopped early on the way to -Inf).
  ref <- ref_zeroinfl(art ~ fem + kid5 + ment | 1, d)
  expect_gte(fit$loglik, as.numeric(logLik(ref)) - 1e-8)
})

test_that("zinbreg handles an offset (argument or formula term)", {
  skip_if_not_installed("pscl")
  set.seed(31)
  n <- 2000
  d <- data.frame(x = rnorm(n), z = rnorm(n), expo = runif(n, 1, 3))
  d$y <- ifelse(rbinom(n, 1, plogis(-0.5 + 1 * d$z)) == 1, 0,
                rnbinom(n, size = 2, mu = d$expo * exp(0.5 + 0.7 * d$x)))

  fit_arg  <- zinbreg(y ~ x | z, data = d, offset = log(expo))
  fit_term <- zinbreg(y ~ x + offset(log(expo)) | z, data = d)
  ref <- ref_zeroinfl(y ~ x + offset(log(expo)) | z, d)

  expect_matches_zeroinfl(fit_arg, ref)
  expect_equal(fit_arg$coefficients, fit_term$coefficients, tolerance = 1e-10)
  # (Recovery of the true parameters is assessed properly, over R = 100
  # replications, by the simulation study in inst/benchmark/.)
})

test_that("E-step tau equals the Sec 3.2 formula", {
  y   <- c(0, 0, 2, 0, 5)
  mu  <- c(1, 3, 2, 0.5, 4)
  eta <- c(-1, 0.5, 0, 2, -3)
  theta <- 1.3
  pi <- plogis(eta)
  f0 <- dnbinom(0, size = theta, mu = mu)
  ref <- ifelse(y == 0, pi / (pi + (1 - pi) * f0), 0)
  expect_equal(zinb_estep(y, mu, eta, theta), ref, tolerance = 1e-12)
})

test_that("fitted values are (1 - pi) mu and the object has its parts", {
  skip_if_not_installed("pscl")
  fit <- zinbreg(art ~ fem + ment | ment, data = pscl::bioChemists)
  expect_s3_class(fit, "zinbreg")
  expect_equal(fit$fitted.values, (1 - fit$pi) * fit$mu)
  expect_equal(fit$df, 3L + 2L + 1L)
  expect_equal(fit$nobs, nrow(pscl::bioChemists))
})

test_that("informative errors for bad input", {
  d <- data.frame(x = 1:6, z = c(1, 0, 1, 0, 1, 0), expo = 1:6,
                  y = c(0, 1, 0, 3, 2, 0))
  expect_error(zinbreg(~ x, data = d), "two-sided")
  expect_error(zinbreg(I(y - 1) ~ x | z, data = d), "non-negative integer")
  expect_error(zinbreg(I(y + 1) ~ x | z, data = d), "no zeros")
  expect_error(zinbreg(y ~ x | z + offset(log(expo)), data = d),
               "count part")
  expect_error(zinbreg(y ~ x | z, data = d, control = nbreg_control()),
               "zinbreg_control")
  expect_error(zinbreg_control(par_tol = 0), "positive number")
})

test_that("warning when EM hits maxit", {
  skip_if_not_installed("pscl")
  expect_warning(zinbreg(art ~ fem + ment, data = pscl::bioChemists,
                         control = zinbreg_control(maxit = 3)),
                 "did not converge")
})
