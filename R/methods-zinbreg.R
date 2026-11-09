# S3 methods for class "zinbreg" (objects returned by zinbreg()).


#' Methods for zinbreg objects
#'
#' Standard accessor, summary, and prediction methods for fits from
#' [zinbreg()].
#'
#' @param object,x An object of class `"zinbreg"` (or `"summary.zinbreg"`).
#' @param model For `coef()` and `vcov()`: `"full"` (count and zero
#'   coefficients, named `count_*` and `zero_*` as in `pscl::zeroinfl()`),
#'   `"count"`, or `"zero"`. \eqn{\log\theta} is not included (as in
#'   `zeroinfl()`); its standard error is in `object$SE.logtheta`.
#' @param type For `residuals()`: `"response"` (\eqn{y_i - \hat E(Y_i)}) or
#'   `"pearson"` (response residual divided by \eqn{\sqrt{\widehat{\mathrm{Var}}(Y_i)}}).
#'   For `predict()`: `"response"` (\eqn{(1-\hat\pi)\hat\mu}, the mean of
#'   \eqn{Y}), `"count"` (\eqn{\hat\mu}, the mean of the NB component), or
#'   `"zero"` (\eqn{\hat\pi}, the probability of a structural zero).
#' @param newdata Optional data frame for prediction.
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (ignored or passed on to [stats::printCoefmat()]).
#' @return `coef()`, `vcov()`: estimates and covariance matrix. `logLik()`:
#'   a `"logLik"` object with `df` = p + q + 1, so `AIC()`/`BIC()` work.
#'   `fitted()`: \eqn{(1-\hat\pi_i)\hat\mu_i}. `residuals()`, `predict()`:
#'   see `type`. `nobs()`: number of observations. `summary()`: object of
#'   class `"summary.zinbreg"` with count and zero coefficient tables.
#' @examples
#' set.seed(3)
#' n <- 500
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' d$y <- ifelse(rbinom(n, 1, plogis(-1 + d$z)) == 1, 0,
#'               rnbinom(n, size = 1.5, mu = exp(1 + 0.5 * d$x)))
#' fit <- zinbreg(y ~ x | z, data = d)
#' summary(fit)
#' AIC(fit)
#' predict(fit, newdata = data.frame(x = 0, z = c(-1, 0, 1)), type = "zero")
#' @name zinbreg-methods
NULL


#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         count = cf$count,
         zero  = cf$zero,
         full  = c(stats::setNames(cf$count, paste0("count_", names(cf$count))),
                   stats::setNames(cf$zero,  paste0("zero_",  names(cf$zero)))))
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  V <- object$vcov
  keep <- switch(model,
                 full  = grep("^(count|zero)_", rownames(V)),
                 count = grep("^count_", rownames(V)),
                 zero  = grep("^zero_", rownames(V)))
  V <- V[keep, keep, drop = FALSE]
  if (model != "full") {
    # Drop the prefix so names match coef(object, model = model).
    rownames(V) <- colnames(V) <- sub("^(count|zero)_", "", rownames(V))
  }
  V
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) object$nobs

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  mu <- object$mu
  pi <- object$pi
  r  <- object$y - (1 - pi) * mu                   # y - E(Y), E(Y) = (1 - pi) mu
  if (type == "pearson") {
    # ZINB variance: E(Y^2) = (1-pi)(mu + mu^2 + mu^2/theta), so
    # Var(Y) = E(Y^2) - E(Y)^2 = (1-pi) mu [1 + mu (pi + 1/theta)].
    r <- r / sqrt((1 - pi) * mu * (1 + mu * (pi + 1 / object$theta)))
  }
  r
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata,
                            type = c("response", "count", "zero"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    mu <- object$mu
    pi <- object$pi
  } else {
    tx <- object$terms$count
    tz <- object$terms$zero
    # Count part: eta_c = X beta + offset (offset() terms + offset argument).
    mf_x <- new_model_frame(tx, newdata, object$xlevels$count)
    X    <- model.matrix(delete.response(tx), mf_x,
                         contrasts.arg = object$contrasts$count)
    off  <- count_offset(mf_x, tx) +
      call_offset_newdata(object$call, newdata, environment(tx))
    mu   <- exp(drop(X %*% object$coefficients$count) + off)
    # Zero part: pi = logit^{-1}(Z gamma).
    mf_z <- new_model_frame(tz, newdata, object$xlevels$zero)
    Z    <- model.matrix(tz, mf_z, contrasts.arg = object$contrasts$zero)
    pi   <- plogis(drop(Z %*% object$coefficients$zero))
  }
  switch(type,
         response = (1 - pi) * mu,
         count    = mu,
         zero     = pi)
}

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Zero-inflated negative binomial regression, fitted by EM\n\n")
  cat("Count model coefficients (NB2 with log link):\n")
  print.default(format(x$coefficients$count, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits),
                print.gap = 2L, quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Log-likelihood:", format(x$loglik, digits = digits + 2L),
      "on", x$df, "Df\n")
  if (!x$converged) cat("WARNING: EM did not converge.\n")
  if (x$boundary) cat("WARNING: no zero inflation (zero model on the boundary).\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  b  <- object$coefficients$count
  g  <- object$coefficients$zero
  # Count table: beta rows plus a Log(theta) row, as zeroinfl() reports.
  count_tab <- coef_table(c(b, `Log(theta)` = log(object$theta)),
                          c(se[paste0("count_", names(b))], object$SE.logtheta))
  zero_tab  <- coef_table(g, se[paste0("zero_", names(g))])
  rownames(count_tab) <- c(names(b), "Log(theta)")
  rownames(zero_tab)  <- names(g)
  structure(list(call         = object$call,
                 coefficients = list(count = count_tab, zero = zero_tab),
                 theta        = object$theta,
                 SE.theta     = object$SE.theta,
                 loglik       = object$loglik,
                 df           = object$df,
                 aic          = AIC(object),
                 bic          = BIC(object),
                 nobs         = object$nobs,
                 iter         = object$iter,
                 converged    = object$converged,
                 boundary     = object$boundary,
                 pearson      = residuals(object, type = "pearson")),
            class = "summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  print(structure(quantile(x$pearson), names = c("Min", "1Q", "Median",
                                                 "3Q", "Max")),
        digits = digits)
  cat("\nCount model coefficients (NB2 with log link):\n")
  printCoefmat(x$coefficients$count, digits = digits, ...)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  printCoefmat(x$coefficients$zero, digits = digits, ...)
  cat("\nTheta:", format(x$theta, digits = digits),
      " (SE", format(x$SE.theta, digits = digits), ")\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      x$df, "Df;  AIC:", format(x$aic, digits = digits + 2L),
      "  BIC:", format(x$bic, digits = digits + 2L), "\n")
  cat("Number of observations:", x$nobs, "  EM iterations:", x$iter, "\n")
  if (!x$converged) cat("WARNING: EM did not converge.\n")
  if (x$boundary) cat("WARNING: no zero inflation (zero model on the boundary).\n")
  invisible(x)
}
