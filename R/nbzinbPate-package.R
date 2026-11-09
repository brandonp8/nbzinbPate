# Package-level documentation.
#
# The special string "_PACKAGE" tells roxygen2 to build the help page
# ?nbzinbPate from the Title and Description fields of DESCRIPTION.
# Package-wide @importFrom tags (functions we borrow from base package
# 'stats') are collected here, so NAMESPACE has one obvious source.
# None of these are fitting routines: the model.* functions only build
# design matrices, plogis/qlogis are the logistic CDF and its inverse,
# pnorm gives p-values, and coef/vcov/... are the generics our S3 methods
# extend.

#' @keywords internal
#' @importFrom stats model.frame model.response model.matrix model.offset
#' @importFrom stats .getXlevels plogis qlogis terms delete.response na.pass
#' @importFrom stats pnorm quantile printCoefmat AIC BIC
#' @importFrom stats coef vcov nobs logLik fitted residuals predict
"_PACKAGE"
