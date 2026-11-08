# Weighted least squares: the linear-algebra core of every IRLS step.
#
# Every IRLS update in this package (NB count model, Sec 3.1(a); logistic
# zero model, Sec 3.2) has the same form
#
#     beta_new = (X' W X)^{-1} X' W z,      W = diag(w_1, ..., w_n),
#
# i.e. a weighted least-squares regression of the working response z on X.
# The project forbids lm(), lm.wfit() and glm.fit() for this, so we solve
# it directly with base R's QR decomposition.


#' Solve a weighted least-squares problem by QR decomposition
#'
#' Computes \eqn{\hat\beta = (X^\top W X)^{-1} X^\top W z}.
#'
#' **Why QR instead of `solve(crossprod(X, w * X), ...)`?** Define
#' \eqn{\tilde X = W^{1/2} X} and \eqn{\tilde z = W^{1/2} z}. Then
#' \eqn{X^\top W X = \tilde X^\top \tilde X} and
#' \eqn{X^\top W z = \tilde X^\top \tilde z}, so \eqn{\hat\beta} is the
#' ordinary least-squares solution of \eqn{\tilde z} on \eqn{\tilde X}.
#' QR solves that without ever forming \eqn{X^\top W X}, whose condition
#' number is the *square* of that of \eqn{\tilde X}; this is the same
#' approach glm.fit() uses internally and gives more accurate answers when
#' covariates are nearly collinear.
#'
#' @param X n x p design matrix.
#' @param w Length-n vector of non-negative working weights.
#' @param z Length-n working response.
#' @return Length-p coefficient vector (named by `colnames(X)`).
#' @noRd
wls_solve <- function(X, w, z) {
  sw <- sqrt(w)                    # W^{1/2} (diagonal, so elementwise sqrt)
  Xt <- sw * X                     # X~ = W^{1/2} X  (row i multiplied by sqrt(w_i))
  zt <- sw * z                     # z~ = W^{1/2} z

  qrX <- qr(Xt)                    # X~ = Q R
  if (qrX$rank < ncol(X)) {
    # A rank-deficient design (e.g. a duplicated or constant column) has no
    # unique solution; stop with a clear message instead of returning junk.
    stop("design matrix is rank deficient (rank ", qrX$rank, " < ",
         ncol(X), " columns); remove redundant covariates.", call. = FALSE)
  }

  # qr.coef() solves R beta = Q' z~ by back-substitution, which is the
  # least-squares solution of z~ on X~, i.e. (X'WX)^{-1} X'Wz.
  beta <- qr.coef(qrX, zt)
  names(beta) <- colnames(X)
  beta
}


#' Inverse of the weighted cross-product matrix \eqn{(X^\top W X)^{-1}}
#'
#' Used for standard errors: at convergence of the NB fit,
#' \eqn{\widehat{\mathrm{Var}}(\hat\beta) = (X^\top \hat W X)^{-1}}
#' (Sec 3.1(d)).
#'
#' `chol()` factors the symmetric positive-definite matrix as
#' \eqn{X^\top W X = R^\top R}, and `chol2inv(R)` returns
#' \eqn{(R^\top R)^{-1}}; this is cheaper and more stable than `solve()`.
#'
#' @inheritParams wls_solve
#' @return p x p matrix with row/column names from `colnames(X)`.
#' @noRd
xtwx_inverse <- function(X, w) {
  XtWX <- crossprod(X, w * X)      # X' W X  (w * X scales row i by w_i)
  V <- chol2inv(chol(XtWX))
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}
