# Tests for the S3 methods of "nbreg" and "zinbreg" objects.
# MASS::glm.nb() and pscl::zeroinfl() are used ONLY as references here.

# ---- nbreg -------------------------------------------------------------------

nb_pair <- function() {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  f <- art ~ fem + mar + kid5 + phd + ment
  list(d = d,
       fit = nbreg(f, data = d),
       ref = MASS::glm.nb(f, data = d,
                          control = glm.control(epsilon = 1e-10, maxit = 100)))
}

test_that("nbreg accessors, logLik, AIC, BIC agree with glm.nb", {
  p <- nb_pair()
  expect_equal(coef(p$fit), coef(p$ref), tolerance = 1e-6)
  expect_equal(vcov(p$fit), vcov(p$ref), tolerance = 1e-5)
  expect_equal(nobs(p$fit), nobs(p$ref))
  ll <- logLik(p$fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), attr(logLik(p$ref), "df"))
  expect_equal(AIC(p$fit), AIC(p$ref), tolerance = 1e-8)
  expect_equal(BIC(p$fit), BIC(p$ref), tolerance = 1e-8)
  expect_equal(fitted(p$fit), fitted(p$ref), tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("nbreg residuals match glm.nb (response and Pearson)", {
  p <- nb_pair()
  expect_equal(residuals(p$fit), residuals(p$ref, type = "response"),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(residuals(p$fit, type = "pearson"),
               residuals(p$ref, type = "pearson"),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("predict.nbreg: no newdata, newdata with factors, both types", {
  p <- nb_pair()
  expect_equal(predict(p$fit), fitted(p$fit))
  expect_equal(predict(p$fit, type = "link"), log(fitted(p$fit)))
  # newdata whose factors use only some levels: dummies must still line up.
  nd <- data.frame(fem = factor("Women", levels = levels(p$d$fem)),
                   mar = factor(c("Single", "Married"),
                                levels = levels(p$d$mar)),
                   kid5 = c(0, 2), phd = c(3, 4), ment = c(5, 20))
  expect_equal(predict(p$fit, newdata = nd),
               predict(p$ref, newdata = nd, type = "response"),
               tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(predict(p$fit, newdata = nd, type = "link"),
               predict(p$ref, newdata = nd, type = "link"),
               tolerance = 1e-6, ignore_attr = TRUE)
  # Plain character columns are matched to the fitted factor levels too.
  nd_chr <- transform(nd, fem = "Women",
                      mar = c("Single", "Married"))
  expect_equal(predict(p$fit, newdata = nd_chr), predict(p$fit, newdata = nd))
})

test_that("predict.nbreg uses the offset in newdata", {
  set.seed(51)
  d <- data.frame(x = rnorm(200), expo = runif(200, 1, 4))
  d$y <- rnbinom(200, size = 2, mu = d$expo * exp(0.3 + 0.5 * d$x))
  nd <- data.frame(x = c(0, 1), expo = c(1, 3))
  f_arg  <- nbreg(y ~ x, data = d, offset = log(expo))
  f_term <- nbreg(y ~ x + offset(log(expo)), data = d)
  manual <- nd$expo * exp(drop(cbind(1, nd$x) %*% coef(f_arg)))
  expect_equal(predict(f_arg, newdata = nd), manual, ignore_attr = TRUE)
  expect_equal(predict(f_term, newdata = nd), manual, tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("nbreg print and summary produce the expected output", {
  p <- nb_pair()
  expect_output(print(p$fit), "Negative binomial")
  s <- summary(p$fit)
  expect_s3_class(s, "summary.nbreg")
  # Coefficient table equals glm.nb's (estimates, SEs, z, p).
  expect_equal(s$coefficients, summary(p$ref)$coefficients, tolerance = 1e-5,
               ignore_attr = TRUE)
  expect_output(print(s), "Pr\\(>\\|z\\|\\)")
  expect_output(print(s), "Theta")
})

# ---- zinbreg -------------------------------------------------------------------

zinb_pair <- function() {
  skip_if_not_installed("pscl")
  d <- pscl::bioChemists
  f <- art ~ fem + mar + kid5 + phd + ment | kid5 + ment
  list(d = d,
       fit = zinbreg(f, data = d),
       ref = pscl::zeroinfl(f, data = d, dist = "negbin",
                            control = pscl::zeroinfl.control(reltol = 1e-12)))
}

test_that("zinbreg accessors, logLik, AIC, BIC agree with zeroinfl", {
  p <- zinb_pair()
  expect_equal(coef(p$fit), coef(p$ref), tolerance = 1e-4)
  expect_named(coef(p$fit), names(coef(p$ref)))
  expect_equal(coef(p$fit, "count"), coef(p$ref, "count"), tolerance = 1e-4)
  expect_equal(coef(p$fit, "zero"), coef(p$ref, "zero"), tolerance = 1e-4)
  expect_equal(vcov(p$fit), vcov(p$ref), tolerance = 1e-3)
  expect_equal(dimnames(vcov(p$fit, "zero")),
               list(names(coef(p$fit, "zero")), names(coef(p$fit, "zero"))))
  expect_equal(nobs(p$fit), nrow(p$d))          # zeroinfl has no nobs() method
  expect_equal(attr(logLik(p$fit), "df"), attr(logLik(p$ref), "df"))
  expect_equal(AIC(p$fit), AIC(p$ref), tolerance = 1e-8)
  expect_equal(BIC(p$fit), BIC(p$ref), tolerance = 1e-8)
  expect_equal(fitted(p$fit), fitted(p$ref), tolerance = 1e-5,
               ignore_attr = TRUE)
})

test_that("zinbreg residuals match zeroinfl (response and Pearson)", {
  p <- zinb_pair()
  expect_equal(residuals(p$fit), residuals(p$ref, type = "response"),
               tolerance = 1e-5, ignore_attr = TRUE)
  expect_equal(residuals(p$fit, type = "pearson"),
               residuals(p$ref, type = "pearson"),
               tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("predict.zinbreg types match zeroinfl, with and without newdata", {
  # Tolerance 1e-4 (relative): gamma-hat agrees with zeroinfl to ~1e-5
  # (Table 4 allows 1e-3), which moves the predicted pi by ~2e-5 relative.
  p <- zinb_pair()
  for (ty in c("response", "count", "zero")) {
    expect_equal(predict(p$fit, type = ty), predict(p$ref, type = ty),
                 tolerance = 1e-4, ignore_attr = TRUE)
  }
  nd <- data.frame(fem = factor(c("Men", "Women"), levels = levels(p$d$fem)),
                   mar = factor("Married", levels = levels(p$d$mar)),
                   kid5 = c(0, 1), phd = c(2.5, 4), ment = c(3, 30))
  for (ty in c("response", "count", "zero")) {
    expect_equal(predict(p$fit, newdata = nd, type = ty),
                 predict(p$ref, newdata = nd, type = ty),
                 tolerance = 1e-4, ignore_attr = TRUE)
  }
  # response = (1 - pi) mu
  expect_equal(predict(p$fit, nd),
               (1 - predict(p$fit, nd, type = "zero")) *
                 predict(p$fit, nd, type = "count"))
})

test_that("predict.zinbreg applies the count offset in newdata", {
  set.seed(52)
  n <- 600
  d <- data.frame(x = rnorm(n), z = rnorm(n), expo = runif(n, 1, 3))
  d$y <- ifelse(rbinom(n, 1, plogis(-0.5 + d$z)) == 1, 0,
                rnbinom(n, size = 2, mu = d$expo * exp(0.5 + 0.5 * d$x)))
  nd <- data.frame(x = c(0, 1), z = c(0, 0), expo = c(1, 2))
  fit <- zinbreg(y ~ x | z, data = d, offset = log(expo))
  manual <- nd$expo * exp(drop(cbind(1, nd$x) %*% coef(fit, "count")))
  expect_equal(predict(fit, nd, type = "count"), manual, ignore_attr = TRUE)
})

test_that("zinbreg print and summary produce the expected output", {
  p <- zinb_pair()
  expect_output(print(p$fit), "Zero-inflation model coefficients")
  s <- summary(p$fit)
  expect_s3_class(s, "summary.zinbreg")
  ref_s <- summary(p$ref)$coefficients
  # Both tables (count incl. Log(theta), and zero) equal zeroinfl's.
  expect_equal(rownames(s$coefficients$count), rownames(ref_s$count))
  expect_equal(s$coefficients$count, ref_s$count, tolerance = 1e-3,
               ignore_attr = TRUE)
  expect_equal(s$coefficients$zero, ref_s$zero, tolerance = 1e-3,
               ignore_attr = TRUE)
  expect_output(print(s), "Log\\(theta\\)")
  expect_output(print(s), "EM iterations")
})
