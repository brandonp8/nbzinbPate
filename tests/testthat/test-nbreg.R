# Tests for nbreg(): agreement with MASS::glm.nb() at the project's Table 4
# tolerances (max absolute difference), plus basic behaviour.
# MASS::glm.nb() is used ONLY as a reference here.

# Max absolute difference, the measure used in Table 4.
max_abs_diff <- function(a, b) max(abs(unname(a) - unname(b)))

# glm.nb run with tight tolerances, as the project suggests.
ref_glm_nb <- function(formula, data) {
  MASS::glm.nb(formula, data = data,
               control = glm.control(epsilon = 1e-10, maxit = 100))
}

test_that("nbreg matches glm.nb on bioChemists (Table 4)", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  f <- art ~ fem + mar + kid5 + phd + ment
  fit <- nbreg(f, data = d)
  ref <- ref_glm_nb(f, d)

  expect_true(fit$converged)
  expect_lt(max_abs_diff(fit$coefficients, coef(ref)), 1e-5)
  expect_lt(max_abs_diff(sqrt(diag(fit$vcov)), sqrt(diag(vcov(ref)))), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$SE.theta - ref$SE.theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
})

test_that("nbreg matches glm.nb on NMES1988 (Table 4)", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  fit <- nbreg(f, data = NMES1988)
  ref <- ref_glm_nb(f, NMES1988)

  expect_true(fit$converged)
  expect_lt(max_abs_diff(fit$coefficients, coef(ref)), 1e-5)
  expect_lt(max_abs_diff(sqrt(diag(fit$vcov)), sqrt(diag(vcov(ref)))), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
})

test_that("an offset (argument or formula term) is handled like glm.nb", {
  skip_if_not_installed("MASS")
  set.seed(11)
  n <- 400
  d <- data.frame(x = rnorm(n), expo = runif(n, 1, 4))
  d$y <- rnbinom(n, size = 1.5, mu = d$expo * exp(0.2 + 0.6 * d$x))

  fit_arg  <- nbreg(y ~ x, data = d, offset = log(expo))
  fit_term <- nbreg(y ~ x + offset(log(expo)), data = d)
  ref <- MASS::glm.nb(y ~ x + offset(log(expo)), data = d,
                      control = glm.control(epsilon = 1e-10, maxit = 100))

  expect_lt(max_abs_diff(fit_arg$coefficients, coef(ref)), 1e-5)
  expect_equal(fit_arg$coefficients, fit_term$coefficients, tolerance = 1e-10)
  expect_lt(abs(fit_arg$theta - ref$theta), 1e-4)
})

test_that("a single-covariate and an intercept-only model work", {
  skip_if_not_installed("MASS")
  set.seed(12)
  d <- data.frame(x = runif(250))
  d$y <- rnbinom(250, size = 3, mu = exp(1 - d$x))

  fit1 <- nbreg(y ~ x, data = d)
  ref1 <- ref_glm_nb(y ~ x, d)
  expect_lt(max_abs_diff(fit1$coefficients, coef(ref1)), 1e-5)

  fit0 <- nbreg(y ~ 1, data = d)
  # Intercept-only: the MLE of mu is the sample mean.
  expect_equal(unname(exp(fit0$coefficients)), mean(d$y), tolerance = 1e-8)
})

test_that("nbreg returns an 'nbreg' object with the documented parts", {
  set.seed(13)
  d <- data.frame(x = rnorm(100))
  d$y <- rnbinom(100, size = 2, mu = exp(d$x))
  fit <- nbreg(y ~ x, data = d)
  expect_s3_class(fit, "nbreg")
  expect_named(fit$coefficients, c("(Intercept)", "x"))
  expect_equal(dim(fit$vcov), c(2L, 2L))
  expect_equal(fit$df, 3L)
  expect_equal(fit$nobs, 100L)
  expect_equal(length(fit$fitted.values), 100L)
})

test_that("informative errors for bad input", {
  d <- data.frame(x = 1:5, y = c(0, 1, 2, 3, 4))
  expect_error(nbreg(cbind(y, y) ~ x, data = d), "numeric vector")
  expect_error(nbreg(I(y - 1) ~ x, data = d), "non-negative integer")
  expect_error(nbreg(I(y + 0.5) ~ x, data = d), "non-negative integer")
  expect_error(nbreg(y ~ x, data = d, control = list(tol = 1)), "nbreg_control")
  expect_error(nbreg_control(tol = -1), "positive number")
  expect_error(nbreg_control(maxit = 2.5), "whole number")
  expect_error(nbreg_control(trace = "yes"), "TRUE or FALSE")
})

test_that("warnings for non-convergence and for the Poisson limit", {
  set.seed(14)
  d <- data.frame(x = rnorm(300))
  d$y <- rnbinom(300, size = 2, mu = exp(0.5 + d$x))
  expect_warning(nbreg(y ~ x, data = d, control = nbreg_control(maxit = 1)),
                 "did not converge")

  # Underdispersed counts (Binomial(4, 0.5): variance/mean = 0.5 < 1), so
  # the likelihood keeps increasing as theta -> Inf (the Poisson limit).
  # (Plain Poisson samples are not reliable here: by chance their sample
  # variance often exceeds the mean, giving a legitimately finite theta-hat.)
  d$y <- rbinom(300, size = 4, prob = 0.5)
  expect_warning(fit <- nbreg(y ~ x, data = d), "Poisson")
  expect_gt(fit$theta, 1e5)
})

test_that("trace prints progress", {
  set.seed(15)
  d <- data.frame(x = rnorm(50))
  d$y <- rnbinom(50, size = 2, mu = exp(d$x))
  expect_output(nbreg(y ~ x, data = d, control = nbreg_control(trace = TRUE)),
                "nbreg iter")
})
