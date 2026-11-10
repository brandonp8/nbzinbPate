# Helpers shared by the fitting code and the S3 methods.


#' Did the log-likelihood really decrease (beyond rounding error)?
#'
#' Used by the step-halving safeguards in nb_irls_beta(), nb_newton_theta()
#' and logit_irls(). Near the optimum a correct step changes the
#' log-likelihood by less than its own floating-point rounding error (for
#' |l| ~ 1500, doubles resolve changes only down to ~1e-13). Comparing
#' `ll_new < ll_old` exactly then fails at random -- on some platforms
#' rounding gives -1e-13 -- and halving a correct step 30 times stalls the
#' iteration and fakes convergence. So a decrease counts only if it exceeds
#' a relative 1e-10 (about 1.6e-7 when |l| ~ 1500): far above rounding noise
#' (~1e-12), far below the decrease caused by a real overshoot from a poor
#' start (order 1 or more). A non-finite new value always counts.
#'
#' @param ll_new,ll_old New and previous log-likelihood (scalars).
#' @return `TRUE` if the step should be shortened.
#' @noRd
ll_decreased <- function(ll_new, ll_old) {
  !is.finite(ll_new) || ll_new < ll_old - 1e-10 * (1 + abs(ll_old))
}


#' Wald coefficient table: estimate, SE, z value, two-sided p-value
#'
#' \eqn{z_j = \hat\beta_j / SE_j} and \eqn{p_j = 2\{1 - \Phi(|z_j|)\}}, computed
#' as `2 * pnorm(-|z|)` (more accurate than `2 * (1 - pnorm(|z|))` for
#' large |z|, where 1 - pnorm(|z|) would round to 0).
#' @param est Named vector of estimates.
#' @param se Vector of standard errors (same order).
#' @return Matrix with columns as in summary.glm().
#' @noRd
coef_table <- function(est, se) {
  z <- est / se
  cbind(Estimate     = est,
        `Std. Error` = se,
        `z value`    = z,
        `Pr(>|z|)`   = 2 * pnorm(-abs(z)))
}


#' Model frame for new data, consistent with the fitted model
#'
#' Rebuilds a model frame from `newdata` using the fitted terms (with the
#' response removed, since new data need not contain it) and the factor
#' levels seen at fitting time (`xlev`), so that a factor in `newdata`
#' with fewer levels still gets the same dummy columns.
#' @noRd
new_model_frame <- function(terms, newdata, xlev) {
  model.frame(delete.response(terms), newdata, na.action = na.pass,
              xlev = xlev)
}


#' Offset given through the `offset =` argument, re-evaluated in newdata
#'
#' At fitting time `offset = log(expo)` was evaluated inside `data`; for
#' prediction the same expression (stored in the matched call) is evaluated
#' inside `newdata`, falling back to the formula's environment for
#' variables not in `newdata`. Returns 0 if no `offset` argument was used.
#' @param call The fitted model's matched call.
#' @param newdata The new data frame.
#' @param env The environment of the model formula.
#' @noRd
call_offset_newdata <- function(call, newdata, env) {
  expr <- call$offset
  if (is.null(expr)) return(0)
  eval(expr, newdata, env)
}
