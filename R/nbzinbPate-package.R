# Package-level documentation.
#
# The special string "_PACKAGE" tells roxygen2 to build the help page
# ?nbzinbPate from the Title and Description fields of DESCRIPTION.
# Package-wide @importFrom tags (functions we borrow from base package
# 'stats') are collected here, so NAMESPACE has one obvious source.
# None of these are fitting routines: the model.* functions only build
# design matrices, and plogis/qlogis are the logistic CDF and its inverse.

#' @keywords internal
#' @importFrom stats model.frame model.response model.matrix model.offset
#' @importFrom stats .getXlevels plogis qlogis
"_PACKAGE"
