# S3 methods for class "nbreg" (objects returned by nbreg()).


#' Methods for nbreg objects
#'
#' Standard accessor, summary, and prediction methods for fits from
#' [nbreg()].
#'
#' @param object,x An object of class `"nbreg"` (or `"summary.nbreg"` for
#'   the print method of a summary).
#' @param type For `residuals()`: `"response"` (\eqn{y_i - \hat\mu_i}) or
#'   `"pearson"` (\eqn{(y_i-\hat\mu_i)/\sqrt{\hat\mu_i + \hat\mu_i^2/\hat\theta}}).
#'   For `predict()`: `"response"` (\eqn{\hat\mu}) or `"link"`
#'   (\eqn{\hat\eta = \log\hat\mu}).
#' @param newdata Optional data frame in which to look for covariates (and
#'   the offset variable) for prediction. If omitted, fitted values are used.
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (ignored or passed on to [stats::printCoefmat()]).
#' @return `coef()`: estimated coefficients. `vcov()`: their covariance
#'   matrix. `logLik()`: a `"logLik"` object with `df` = number of
#'   coefficients + 1 (for \eqn{\theta}), so [stats::AIC()] and [stats::BIC()]
#'   work. `fitted()`: \eqn{\hat\mu}. `residuals()`: see `type`. `nobs()`:
#'   number of observations. `predict()`: numeric vector. `summary()`: an
#'   object of class `"summary.nbreg"` with a coefficient table.
#' @examples
#' set.seed(2)
#' d <- data.frame(x = rnorm(200), g = factor(sample(c("a", "b"), 200, TRUE)))
#' d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.4 * d$x + 0.3 * (d$g == "b")))
#' fit <- nbreg(y ~ x + g, data = d)
#' summary(fit)
#' AIC(fit); BIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1), g = "b"))
#' @name nbreg-methods
NULL


#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) object$nobs

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  # A "logLik" object carries df (number of estimated parameters) and nobs;
  # AIC() = -2 logLik + 2 df and BIC() = -2 logLik + log(nobs) df use them.
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  r  <- object$y - mu                                # response residual y - mu
  if (type == "pearson") {
    # Divide by the NB2 standard deviation sqrt(V(mu)), V(mu) = mu + mu^2/theta.
    r <- r / sqrt(mu + mu^2 / object$theta)
  }
  r
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    # Rebuild the design matrix for newdata with the fitted terms, factor
    # levels and contrasts, then eta = X beta + offset.
    mf  <- new_model_frame(object$terms, newdata, object$xlevels)
    X   <- model.matrix(delete.response(object$terms), mf,
                        contrasts.arg = object$contrasts)
    off <- model.offset(mf)                          # offset() terms, or NULL
    if (is.null(off)) off <- 0
    off <- off + call_offset_newdata(object$call, newdata,
                                     environment(object$terms))
    eta <- drop(X %*% object$coefficients) + off
  }
  if (type == "link") eta else exp(eta)              # mu = exp(eta)
}

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Negative binomial (NB2) regression, fitted by IRLS + Newton-Raphson\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Log-likelihood:", format(x$loglik, digits = digits + 2L),
      "on", x$df, "Df\n")
  if (!x$converged) cat("WARNING: the algorithm did not converge.\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  structure(list(call         = object$call,
                 coefficients = coef_table(object$coefficients, se),
                 theta        = object$theta,
                 SE.theta     = object$SE.theta,
                 loglik       = object$loglik,
                 df           = object$df,
                 aic          = AIC(object),
                 bic          = BIC(object),
                 nobs         = object$nobs,
                 iter         = object$iter,
                 converged    = object$converged,
                 pearson      = residuals(object, type = "pearson")),
            class = "summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  print(structure(quantile(x$pearson), names = c("Min", "1Q", "Median",
                                                 "3Q", "Max")),
        digits = digits)
  cat("\nCoefficients (NB2 with log link):\n")
  printCoefmat(x$coefficients, digits = digits, ...)
  cat("\nTheta:", format(x$theta, digits = digits),
      " (SE", format(x$SE.theta, digits = digits), ")\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2L), "on",
      x$df, "Df;  AIC:", format(x$aic, digits = digits + 2L),
      "  BIC:", format(x$bic, digits = digits + 2L), "\n")
  cat("Number of observations:", x$nobs,
      "  Outer iterations (IRLS + Newton):", x$iter, "\n")
  if (!x$converged) cat("WARNING: the algorithm did not converge.\n")
  invisible(x)
}
