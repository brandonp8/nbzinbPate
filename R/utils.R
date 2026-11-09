# Helpers shared by the nbreg and zinbreg S3 methods.


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
