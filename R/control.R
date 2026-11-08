# Control functions: collect tolerances, iteration limits, and the trace
# option in one validated list, in the style of glm.control().


#' Control settings for [nbreg()]
#'
#' @param tol Outer-loop tolerance (Sec 3.1(c)). The fit has converged when
#'   the log-likelihood changes by less than `tol` **and** no coefficient or
#'   \eqn{\log\theta} changes by more than `tol` between outer iterations.
#' @param maxit Maximum number of outer iterations (each one = an IRLS fit
#'   for \eqn{\beta} followed by a Newton fit for \eqn{\theta}).
#' @param irls_tol,irls_maxit Tolerance (largest change in \eqn{\beta}) and
#'   iteration limit for the inner IRLS of Sec 3.1(a).
#' @param theta_tol,theta_maxit Tolerance (change in \eqn{\log\theta}) and
#'   iteration limit for the inner Newton-Raphson of Sec 3.1(b).
#' @param theta_warn A warning is issued if \eqn{\hat\theta} exceeds this
#'   value: the data then show essentially no overdispersion (the Poisson
#'   limit), and a Poisson model may be adequate.
#' @param trace Logical; if `TRUE`, print the log-likelihood and
#'   \eqn{\theta} at every outer iteration.
#' @return A list of class `"nbreg_control"` with the components above.
#' @examples
#' nbreg_control(tol = 1e-10, trace = TRUE)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100,
                          irls_tol = 1e-10, irls_maxit = 100,
                          theta_tol = 1e-10, theta_maxit = 100,
                          theta_warn = 1e5, trace = FALSE) {
  check_positive_number(tol, "tol")
  check_positive_number(irls_tol, "irls_tol")
  check_positive_number(theta_tol, "theta_tol")
  check_positive_number(theta_warn, "theta_warn")
  check_count(maxit, "maxit")
  check_count(irls_maxit, "irls_maxit")
  check_count(theta_maxit, "theta_maxit")
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    stop("'trace' must be TRUE or FALSE.", call. = FALSE)
  }
  structure(list(tol = tol, maxit = maxit,
                 irls_tol = irls_tol, irls_maxit = irls_maxit,
                 theta_tol = theta_tol, theta_maxit = theta_maxit,
                 theta_warn = theta_warn, trace = trace),
            class = "nbreg_control")
}


# ---- small validators shared by the control functions ------------------------

#' Stop unless x is a single positive finite number
#' @noRd
check_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("'", name, "' must be a single positive number.", call. = FALSE)
  }
}

#' Stop unless x is a single positive whole number
#' @noRd
check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
      x != round(x)) {
    stop("'", name, "' must be a single positive whole number.", call. = FALSE)
  }
}


#' Control settings for [zinbreg()]
#'
#' @param tol EM tolerance (Sec 3.2). EM has converged when the
#'   observed-data log-likelihood changes by less than `tol` **and** no
#'   parameter (\eqn{\beta}, \eqn{\gamma}, \eqn{\log\theta}) changes by more
#'   than `par_tol` between EM iterations.
#' @param par_tol Tolerance on the largest parameter change; see `tol`.
#' @param maxit Maximum number of EM iterations.
#' @param mstep_tol,mstep_maxit Tolerance and iteration limit for the
#'   M-step fits. Each M-step is a *full* fit, warm-started at the current
#'   values: the weighted NB fit uses [nbreg_control()] with
#'   `tol = irls_tol = theta_tol = mstep_tol` and these iteration limits,
#'   and the logistic IRLS uses the same tolerance and limit.
#' @param theta_warn Warn if \eqn{\hat\theta} exceeds this value.
#' @param trace Logical; if `TRUE`, print the log-likelihood at every EM
#'   iteration.
#' @return A list of class `"zinbreg_control"`.
#' @examples
#' zinbreg_control(tol = 1e-10, maxit = 2000)
#' @export
zinbreg_control <- function(tol = 1e-8, par_tol = 1e-6, maxit = 1000,
                            mstep_tol = 1e-10, mstep_maxit = 100,
                            theta_warn = 1e5, trace = FALSE) {
  check_positive_number(tol, "tol")
  check_positive_number(par_tol, "par_tol")
  check_positive_number(mstep_tol, "mstep_tol")
  check_positive_number(theta_warn, "theta_warn")
  check_count(maxit, "maxit")
  check_count(mstep_maxit, "mstep_maxit")
  if (!is.logical(trace) || length(trace) != 1L || is.na(trace)) {
    stop("'trace' must be TRUE or FALSE.", call. = FALSE)
  }
  structure(list(tol = tol, par_tol = par_tol, maxit = maxit,
                 mstep_tol = mstep_tol, mstep_maxit = mstep_maxit,
                 theta_warn = theta_warn, trace = trace),
            class = "zinbreg_control")
}
