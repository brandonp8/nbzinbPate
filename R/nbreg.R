# User-facing formula interface for NB2 regression.
# The estimation itself happens in nbreg_fit() (R/nbreg_fit.R); this file
# only turns (formula, data, offset) into (X, y, offset), checks the input,
# and packages the result as an S3 object of class "nbreg".


#' Negative binomial (NB2) regression
#'
#' Fits the NB2 regression model
#' \deqn{Y_i \sim \mathrm{NB}(\mu_i, \theta), \quad \log\mu_i = x_i^\top\beta + o_i,
#'   \quad \mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta,}
#' by maximum likelihood. The coefficients \eqn{\beta} are updated by
#' iteratively reweighted least squares (IRLS) with \eqn{\theta} held fixed,
#' and \eqn{\log\theta} by Newton-Raphson with \eqn{\beta} held fixed; the
#' two steps alternate until convergence. No general-purpose optimizer is
#' used. The parameterization is the same as [MASS::glm.nb()].
#'
#' Standard errors for \eqn{\beta} are \eqn{(X^\top\hat W X)^{-1/2}} (diagonal),
#' treating \eqn{\theta} as known, as `glm.nb()` does; the standard error of
#' \eqn{\hat\theta} is \eqn{\{-H(\hat\theta)\}^{-1/2}}, where \eqn{H} is the
#' second derivative of the log-likelihood in \eqn{\theta}.
#'
#' @param formula A model formula such as `y ~ x1 + x2`. The response must
#'   be a vector of non-negative integer counts. `offset()` terms are allowed.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset vector (e.g. `log(exposure)`), evaluated in
#'   `data`; it is added to any `offset()` terms in the formula.
#' @param control A list from [nbreg_control()].
#' @return An object of class `"nbreg"`: a list with components
#'   \describe{
#'     \item{coefficients}{Estimated \eqn{\beta}.}
#'     \item{vcov}{Estimated covariance matrix of \eqn{\hat\beta}.}
#'     \item{theta, SE.theta}{\eqn{\hat\theta} and its standard error.}
#'     \item{loglik}{Maximized log-likelihood.}
#'     \item{df}{Number of estimated parameters (\eqn{p + 1}).}
#'     \item{fitted.values, linear.predictors}{\eqn{\hat\mu_i} and \eqn{\hat\eta_i}.}
#'     \item{iter, converged}{Number of outer iterations and convergence flag.}
#'     \item{y, offset, nobs}{The response, offset, and number of observations.}
#'     \item{call, terms, xlevels, contrasts, model}{Information used by
#'       `predict()` and other methods.}
#'   }
#' @seealso [nbreg_control()], [zinbreg()]
#' @examples
#' # Simulated data with an exposure offset
#' set.seed(1)
#' n <- 300
#' d <- data.frame(x = rnorm(n), expo = runif(n, 1, 3))
#' d$y <- rnbinom(n, size = 2, mu = d$expo * exp(0.5 + 0.8 * d$x))
#' fit <- nbreg(y ~ x, data = d, offset = log(expo))
#' fit$coefficients
#' fit$theta
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control()) {
  cl <- match.call()

  # ---- Build the model frame -----------------------------------------------------
  # Standard R idiom (used by lm/glm): keep only the arguments model.frame()
  # understands and evaluate a call to stats::model.frame() in the caller's
  # environment. This evaluates `offset` inside `data` (so offset =
  # log(expo) works), drops rows with missing values, and records the terms.
  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  y  <- model.response(mf, "numeric")
  X  <- model.matrix(mt, mf)
  # model.offset() adds up offset() terms in the formula and the `offset`
  # argument; it is NULL when there is neither.
  off <- model.offset(mf)
  if (is.null(off)) off <- rep(0, length(y))

  check_count_response(y)
  if (!inherits(control, "nbreg_control")) {
    stop("'control' must be created by nbreg_control().", call. = FALSE)
  }
  if (length(off) != length(y) || any(!is.finite(off))) {
    stop("'offset' must be finite and have one value per observation.",
         call. = FALSE)
  }

  # ---- Fit -------------------------------------------------------------------------
  fit <- nbreg_fit(X, y, offset = off, control = control)

  # ---- Package as an "nbreg" object ---------------------------------------------
  structure(
    list(coefficients      = fit$coefficients,
         vcov              = fit$vcov,
         theta             = fit$theta,
         SE.theta          = fit$SE.theta,
         loglik            = fit$loglik,
         df                = ncol(X) + 1L,       # p coefficients + theta
         fitted.values     = fit$mu,
         linear.predictors = fit$eta,
         working.weights   = fit$w,
         iter              = fit$iter,
         converged         = fit$converged,
         y                 = y,
         offset            = off,
         nobs              = length(y),
         call              = cl,
         terms             = mt,
         xlevels           = .getXlevels(mt, mf),
         contrasts         = attr(X, "contrasts"),
         model             = mf,
         control           = control),
    class = "nbreg")
}


#' Stop unless y is a non-empty vector of non-negative integer counts
#' @noRd
check_count_response <- function(y) {
  if (is.null(y) || length(y) == 0L) {
    stop("the model has no response or no observations.", call. = FALSE)
  }
  if (!is.numeric(y) || is.matrix(y)) {
    stop("the response must be a numeric vector of counts.", call. = FALSE)
  }
  if (any(y < 0) || any(y != round(y))) {
    stop("the response must contain non-negative integer counts only.",
         call. = FALSE)
  }
}
