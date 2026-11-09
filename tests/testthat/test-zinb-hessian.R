# Tests for the analytic ZINB score and Hessian (R/zinb_hessian.R).
# numDeriv is used ONLY here, as an independent numerical check of our
# analytic derivatives; pscl::zeroinfl() only as a reference for the SEs.

# Log-likelihood as a function of the stacked parameter vector
# (beta, gamma, log theta), for numDeriv.
loglik_par <- function(par, X, Z, y, offset) {
  p <- ncol(X); q <- ncol(Z)
  beta  <- par[seq_len(p)]
  gamma <- par[p + seq_len(q)]
  theta <- exp(par[p + q + 1])
  zinb_loglik(y, exp(drop(X %*% beta) + offset), drop(Z %*% gamma), theta)
}

bio_setup <- function() {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  list(X = model.matrix(~ fem + mar + kid5 + phd + ment, d),
       Z = model.matrix(~ kid5 + ment, d),
       y = d$art, off = rep(0, nrow(d)))
}

test_that("analytic score and Hessian match numDeriv at an arbitrary point", {
  skip_if_not_installed("numDeriv")
  s <- bio_setup()
  set.seed(41)
  # An arbitrary (non-optimal) point, so all terms are exercised.
  beta  <- c(0.3, -0.2, 0.1, -0.15, 0.02, 0.03) + rnorm(6, sd = 0.05)
  gamma <- c(-1, 0.4, -0.3)
  theta <- 1.8
  an <- zinb_score_hessian(s$X, s$Z, s$y, s$off, beta, gamma, theta)
  par <- c(beta, gamma, log(theta))

  num_score <- numDeriv::grad(loglik_par, par, X = s$X, Z = s$Z, y = s$y,
                              offset = s$off)
  num_hess  <- numDeriv::hessian(loglik_par, par, X = s$X, Z = s$Z, y = s$y,
                                 offset = s$off)
  expect_equal(unname(an$score), num_score, tolerance = 1e-7)
  expect_equal(unname(an$hessian), num_hess, tolerance = 1e-6)
  expect_true(isSymmetric(unname(an$hessian)))
})

test_that("Hessian is correct when theta is small and when it is large", {
  skip_if_not_installed("numDeriv")
  s <- bio_setup()
  beta  <- c(0.3, -0.2, 0.1, -0.15, 0.02, 0.03)
  gamma <- c(-0.5, 0.2, -0.4)
  for (theta in c(0.3, 40)) {
    an  <- zinb_score_hessian(s$X, s$Z, s$y, s$off, beta, gamma, theta)
    num <- numDeriv::hessian(loglik_par, c(beta, gamma, log(theta)),
                             X = s$X, Z = s$Z, y = s$y, offset = s$off)
    expect_equal(unname(an$hessian), num, tolerance = 1e-6)
  }
})

test_that("score is ~0 at the zinbreg estimates", {
  skip_if_not_installed("pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment | kid5 + ment,
                 data = pscl::bioChemists)
  expect_lt(max(abs(fit$score)), 1e-3)
})

test_that("standard errors match zeroinfl on bioChemists and NMES1988 (Table 4)", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  cases <- list(
    list(art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment,
         pscl::bioChemists),
    list(visits ~ hospital + health + chronic + gender + school + insurance |
           hospital + chronic + gender + school + insurance, NMES1988))
  for (cs in cases) {
    fit <- zinbreg(cs[[1]], data = cs[[2]])
    ref <- pscl::zeroinfl(cs[[1]], data = cs[[2]], dist = "negbin",
                          control = pscl::zeroinfl.control(reltol = 1e-12))
    se  <- sqrt(diag(fit$vcov))
    # zeroinfl's vcov() covers beta and gamma (not log theta), same names.
    se_ref <- sqrt(diag(vcov(ref)))
    expect_lt(max(abs(se[names(se_ref)] - se_ref)), 1e-3)
    expect_lt(abs(fit$SE.logtheta - ref$SE.logtheta), 1e-3)
  }
})

test_that("SE(theta) is the delta-method transform of SE(log theta)", {
  skip_if_not_installed("pscl")
  fit <- zinbreg(art ~ fem + ment | ment, data = pscl::bioChemists)
  expect_equal(fit$SE.theta, fit$theta * fit$SE.logtheta)
  expect_equal(unname(sqrt(fit$vcov["Log(theta)", "Log(theta)"])),
               fit$SE.logtheta)
})

test_that("at the no-zero-inflation boundary the SEs are NA", {
  skip_if_not_installed("pscl")
  fit <- suppressWarnings(zinbreg(art ~ fem + kid5 + ment | 1,
                                  data = pscl::bioChemists))
  expect_true(fit$boundary)
  expect_true(all(is.na(fit$vcov)))
})
