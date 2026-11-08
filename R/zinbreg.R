# User-facing formula interface for ZINB regression.
# The EM algorithm itself is in zinbreg_fit() (R/zinbreg_fit.R); this file
# parses the two-part formula y ~ x | z, builds the two design matrices,
# checks the input, and packages the result as an S3 object "zinbreg".


#' Zero-inflated negative binomial (ZINB) regression
#'
#' Fits the ZINB model of Lambert (1992) with an NB2 count component:
#' \deqn{P(Y_i = 0) = \pi_i + (1-\pi_i) f_{NB}(0;\mu_i,\theta), \quad
#'       P(Y_i = y) = (1-\pi_i) f_{NB}(y;\mu_i,\theta),\ y \ge 1,}
#' with \eqn{\log\mu_i = x_i^\top\beta + o_i} (count model) and
#' \eqn{\mathrm{logit}\,\pi_i = z_i^\top\gamma} (zero-inflation model).
#' Estimation is by the EM algorithm: the E-step computes the posterior
#' probability \eqn{\tau_i} that each zero is a structural zero, and the
#' M-step fits a fractional-response logistic regression for \eqn{\gamma}
#' and a weighted NB regression (weights \eqn{1-\tau_i}) for
#' \eqn{(\beta,\theta)}, both by IRLS. No general-purpose optimizer is used.
#' The parameterization is the same as
#' `pscl::zeroinfl(..., dist = "negbin")`.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before `|`, zero-model covariates after it. Without `|`,
#'   the same covariates are used in both parts (as in `pscl::zeroinfl()`).
#'   `offset()` terms are allowed in the count part only.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional count-model offset (e.g. `log(exposure)`),
#'   evaluated in `data`.
#' @param control A list from [zinbreg_control()].
#' @return An object of class `"zinbreg"`: a list with components
#'   \describe{
#'     \item{coefficients}{List with `count` (\eqn{\beta}) and `zero` (\eqn{\gamma}).}
#'     \item{theta}{\eqn{\hat\theta}.}
#'     \item{loglik, loglik_path}{Maximized observed-data log-likelihood,
#'       and its value at every EM iteration (non-decreasing).}
#'     \item{df}{Number of estimated parameters (\eqn{p + q + 1}).}
#'     \item{fitted.values}{\eqn{E(Y_i) = (1-\hat\pi_i)\hat\mu_i}.}
#'     \item{mu, pi}{Count-model means \eqn{\hat\mu_i} and zero probabilities
#'       \eqn{\hat\pi_i}.}
#'     \item{iter, converged}{Number of EM iterations and convergence flag.}
#'     \item{boundary}{`TRUE` if the zero-inflation component vanished
#'       (all \eqn{\hat\pi_i < 10^{-8}}); a warning is then issued.}
#'     \item{y, offset, nobs}{Response, count offset, number of observations.}
#'     \item{call, terms, xlevels, contrasts, model}{Information used by
#'       `predict()` and other methods.}
#'   }
#' @references Lambert, D. (1992). Zero-inflated Poisson regression, with an
#'   application to defects in manufacturing. *Technometrics*, 34(1), 1-14.
#' @seealso [zinbreg_control()], [nbreg()]
#' @examples
#' # Simulated ZINB data: one covariate in each part
#' set.seed(1)
#' n <- 500
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' structural_zero <- rbinom(n, 1, plogis(-1 + d$z))
#' d$y <- ifelse(structural_zero == 1, 0,
#'               rnbinom(n, size = 1.5, mu = exp(1 + 0.5 * d$x)))
#' fit <- zinbreg(y ~ x | z, data = d)
#' fit$coefficients
#' fit$theta
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control()) {
  cl <- match.call()

  # ---- Split the two-part formula y ~ count_rhs | zero_rhs -------------------
  parts <- split_zinb_formula(formula)

  # ---- One model frame holding every variable of both parts ------------------
  # model.frame() is called with the combined formula y ~ count_rhs +
  # zero_rhs, so that missing values are dropped from both parts at once
  # and both design matrices have the same rows.
  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf$formula <- parts$combined
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Separate terms objects for the two parts. `data` is passed so that a
  # `.` in the formula expands to the columns of `data`.
  dat <- if (missing(data)) NULL else data
  mt_x <- terms(parts$count, data = dat)
  mt_z <- terms(parts$zero, data = dat)
  if (length(attr(mt_z, "offset")) > 0L) {
    stop("offset() terms are only supported in the count part of the formula.",
         call. = FALSE)
  }

  y <- model.response(mf, "numeric")
  X <- model.matrix(mt_x, mf)                    # count-model design
  Z <- model.matrix(mt_z, mf)                    # zero-model design
  off <- count_offset(mf, mt_x)

  check_count_response(y)
  if (!any(y == 0)) {
    stop("the response has no zeros, so there is no zero inflation to model.",
         call. = FALSE)
  }
  if (!inherits(control, "zinbreg_control")) {
    stop("'control' must be created by zinbreg_control().", call. = FALSE)
  }
  if (length(off) != length(y) || any(!is.finite(off))) {
    stop("'offset' must be finite and have one value per observation.",
         call. = FALSE)
  }

  # ---- Fit by EM -------------------------------------------------------------
  fit <- zinbreg_fit(X, Z, y, offset = off, control = control)

  if (!fit$converged) {
    warning("zinbreg did not converge in ", control$maxit, " EM iterations; ",
            "results may be unreliable. Try increasing 'maxit' in ",
            "zinbreg_control().", call. = FALSE)
  }
  if (fit$boundary) {
    warning("the zero-inflation probabilities converge to 0 (all pi-hat < ",
            "1e-8): the data show no zero inflation, the MLE of the zero ",
            "model is on the boundary (coefficients -> -Inf), and the NB ",
            "model nbreg() is adequate. Zero-model estimates are not ",
            "meaningful.", call. = FALSE)
  }
  if (fit$theta > control$theta_warn) {
    warning("theta-hat = ", format(fit$theta, digits = 3), " is very large: ",
            "the count part shows little or no overdispersion, so a ",
            "zero-inflated Poisson model may be adequate.", call. = FALSE)
  }

  # ---- Package as a "zinbreg" object -----------------------------------------
  structure(
    list(coefficients  = fit$coefficients,
         theta         = fit$theta,
         loglik        = fit$loglik,
         loglik_path   = fit$loglik_path,
         df            = ncol(X) + ncol(Z) + 1L,   # beta, gamma, theta
         fitted.values = (1 - fit$pi) * fit$mu,    # E(Y) = (1 - pi) mu
         mu            = fit$mu,
         pi            = fit$pi,
         iter          = fit$iter,
         converged     = fit$converged,
         boundary      = fit$boundary,
         y             = y,
         offset        = off,
         nobs          = length(y),
         call          = cl,
         terms         = list(count = mt_x, zero = mt_z,
                              full = attr(mf, "terms")),
         xlevels       = list(count = .getXlevels(mt_x, mf),
                              zero  = .getXlevels(mt_z, mf)),
         contrasts     = list(count = attr(X, "contrasts"),
                              zero  = attr(Z, "contrasts")),
         model         = mf,
         control       = control),
    class = "zinbreg")
}


#' Split a ZINB formula `y ~ a + b | c + d` into its parts (base R only)
#'
#' A formula is a call: `y ~ a + b | c + d` is `~`(y, `|`(a + b, c + d)),
#' because `|` binds more loosely than `+`. So element [[2]] is the
#' response and element [[3]] the right-hand side; if the right-hand side
#' is a call to `|`, its two arguments are the count and zero parts.
#'
#' @param formula The user's formula.
#' @return A list with `count` (`y ~ a + b`), `zero` (`~ c + d`), and
#'   `combined` (`y ~ a + b + (c + d)`, used to build one model frame).
#' @noRd
split_zinb_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula such as y ~ x1 + x2 | z1.",
         call. = FALSE)
  }
  env <- environment(formula)
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    count_rhs <- rhs[[2L]]                     # before the bar
    zero_rhs  <- rhs[[3L]]                     # after the bar
  } else {
    count_rhs <- rhs                           # no bar: same covariates
    zero_rhs  <- rhs                           # in both parts
  }

  # Build the formulas from pieces with bquote(); .() splices in a piece.
  count    <- eval(bquote(.(lhs) ~ .(count_rhs)))
  zero     <- eval(bquote(~ .(zero_rhs)))
  combined <- eval(bquote(.(lhs) ~ .(count_rhs) + .(zero_rhs)))
  environment(count) <- environment(zero) <- environment(combined) <- env
  list(count = count, zero = zero, combined = combined)
}


#' Count-model offset: the `offset` argument plus offset() terms of the
#' count formula
#'
#' In the model frame, each `offset(...)` term is a column named by its
#' deparsed expression, and the `offset` argument is the column
#' `"(offset)"`. `attr(terms, "offset")` gives the positions of the
#' offset terms among the formula's variables.
#' @noRd
count_offset <- function(mf, mt_x) {
  off <- rep(0, nrow(mf))
  if (!is.null(mf[["(offset)"]])) off <- off + mf[["(offset)"]]
  vars <- attr(mt_x, "variables")              # call: list(y, x1, offset(...), ...)
  for (i in attr(mt_x, "offset")) {
    off <- off + mf[[deparse(vars[[i + 1L]])]] # +1 skips the `list` head
  }
  off
}
