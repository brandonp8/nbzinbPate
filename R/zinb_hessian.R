# Analytic score and Hessian of the ZINB observed-data log-likelihood with
# respect to (beta, gamma, phi = log theta), for standard errors (Sec 3.2,
# "Standard errors": EM gives no SEs, so we invert the observed
# information -H at the MLE).
#
# ---------------------------------------------------------------------------
# Derivation (also in the report). Per observation, write
#   eta_c = x'beta + o, mu = exp(eta_c),  eta_z = z'gamma, pi = plogis(eta_z),
#   phi = log(theta),  lf = log f_NB(y; mu, theta),
#   tau = pi / [pi + (1 - pi) f_NB(0)] if y = 0, tau = 0 if y > 0 (E-step).
# beta enters only through eta_c and gamma only through eta_z, so we
# differentiate l_i with respect to the scalars (eta_z, eta_c, phi) and use
# d eta_c / d beta = x_i, d eta_z / d gamma = z_i at the end.
#
# For y = 0, l = log P with P = pi + (1 - pi) f0. Since tau = pi / P and
# 1 - tau = (1 - pi) f0 / P, with a, b in {eta_c, phi}:
#   dl/d eta_z = pi(1-pi)(1-f0)/P = (1-pi) tau - pi (1-tau) = tau - pi
#   dl/da      = (1-pi) f0 (d lf/da) / P             = (1-tau) d lf/da
# For y > 0, l = log(1-pi) + lf: the same formulas with tau = 0.
# Differentiating again, using d tau/d eta_z = tau(1-tau) and
# d(1-tau)/db = tau(1-tau) d lf/db (from log(1-tau) = log(1-pi) + lf - l):
#   d2l/d eta_z^2    = tau(1-tau) - pi(1-pi)
#   d2l/d eta_z da   = -tau(1-tau) d lf/da
#   d2l/da db        = (1-tau) d2lf/da db + tau(1-tau) (d lf/da)(d lf/db)
# (This is Louis' (1982) formula: the observed information equals the
# expected complete-data information minus Var(complete-data score);
# tau(1-tau) = Var(D_i | y_i) is the "missing information".)
#
# NB pieces (Sec 3.1 with the chain rule, d mu/d eta_c = mu,
# d theta/d phi = theta):
#   d lf/d eta_c        = theta (y - mu) / (theta + mu)
#   d2lf/d eta_c^2      = -theta mu (theta + y) / (theta + mu)^2
#   d lf/d phi          = theta U_i
#   d2lf/d phi^2        = theta U_i + theta^2 H_i
#   d2lf/d eta_c d phi  = theta mu (y - mu) / (theta + mu)^2
# where U_i, H_i are the per-observation terms of U(theta), H(theta).
# ---------------------------------------------------------------------------


#' Per-observation derivatives of log f_NB(y; mu, theta) in (eta_c, phi)
#'
#' See the derivation at the top of this file.
#' @return List of vectors `c`, `cc`, `p`, `pp`, `cp` (c = eta_c, p = phi).
#' @noRd
nb_lf_derivs <- function(y, mu, theta) {
  tm <- theta + mu
  # Per-observation U_i and H_i of Sec 3.1(b) (same regrouping as in
  # nb_theta_derivs(), without the sum).
  U_i <- digamma(y + theta) - digamma(theta) - log1p(mu / theta) + (mu - y) / tm
  H_i <- trigamma(y + theta) - trigamma(theta) + mu / (theta * tm) +
    (y - mu) / tm^2
  list(c  = theta * (y - mu) / tm,                  # d lf / d eta_c
       cc = -theta * mu * (theta + y) / tm^2,       # d2 lf / d eta_c^2
       p  = theta * U_i,                            # d lf / d phi
       pp = theta * U_i + theta^2 * H_i,            # d2 lf / d phi^2
       cp = theta * mu * (y - mu) / tm^2)           # d2 lf / d eta_c d phi
}


#' Analytic score vector and Hessian of the ZINB log-likelihood
#'
#' Parameters are ordered \eqn{(\beta, \gamma, \log\theta)}, as in
#' `pscl::zeroinfl()`. Formulas: see the derivation at the top of this file.
#'
#' @param X,Z Count- and zero-model design matrices.
#' @param y Counts.
#' @param offset Count-model offset.
#' @param beta,gamma,theta Parameter values.
#' @return List with `score` (length p + q + 1) and `hessian`
#'   ((p+q+1) x (p+q+1)), both named.
#' @noRd
zinb_score_hessian <- function(X, Z, y, offset, beta, gamma, theta) {
  eta_c <- drop(X %*% beta) + offset
  mu    <- exp(eta_c)
  eta_z <- drop(Z %*% gamma)
  pi    <- plogis(eta_z)
  tau   <- zinb_estep(y, mu, eta_z, theta)   # posterior P(D_i = 1 | y_i)
  v     <- tau * (1 - tau)                   # Var(D_i | y_i): missing information
  d     <- nb_lf_derivs(y, mu, theta)

  # ---- Per-observation first derivatives ------------------------------------
  s_z <- tau - pi                            # dl / d eta_z
  s_c <- (1 - tau) * d$c                     # dl / d eta_c
  s_p <- (1 - tau) * d$p                     # dl / d phi

  # ---- Per-observation second derivatives ------------------------------------
  h_zz <- v - pi * (1 - pi)                  # d2l / d eta_z^2
  h_zc <- -v * d$c                           # d2l / d eta_z d eta_c
  h_zp <- -v * d$p                           # d2l / d eta_z d phi
  h_cc <- (1 - tau) * d$cc + v * d$c^2       # d2l / d eta_c^2
  h_cp <- (1 - tau) * d$cp + v * d$c * d$p   # d2l / d eta_c d phi
  h_pp <- (1 - tau) * d$pp + v * d$p^2       # d2l / d phi^2

  # ---- Chain rule to (beta, gamma, phi) ---------------------------------------
  # Score:   d l / d beta = X' s_c,   d l / d gamma = Z' s_z,   d l / d phi = sum(s_p)
  score <- c(drop(crossprod(X, s_c)), drop(crossprod(Z, s_z)), sum(s_p))
  # Hessian blocks: X' diag(h) X, X' diag(h) Z, X' h, etc.
  # (h * Z scales row i of Z by h_i, so crossprod(X, h * Z) = X' diag(h) Z.)
  H_bb <- crossprod(X, h_cc * X)
  H_bg <- crossprod(X, h_zc * Z)
  H_gg <- crossprod(Z, h_zz * Z)
  H_bp <- crossprod(X, h_cp)
  H_gp <- crossprod(Z, h_zp)
  H_pp <- sum(h_pp)
  hessian <- rbind(cbind(H_bb,    H_bg,    H_bp),
                   cbind(t(H_bg), H_gg,    H_gp),
                   cbind(t(H_bp), t(H_gp), H_pp))

  nm <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
          "Log(theta)")
  names(score) <- nm
  dimnames(hessian) <- list(nm, nm)
  list(score = score, hessian = hessian)
}


#' Covariance matrix of (beta, gamma, log theta) from the analytic Hessian
#'
#' \eqn{\widehat{\mathrm{Var}} = \{-H(\hat\beta,\hat\gamma,\log\hat\theta)\}^{-1}},
#' the inverse observed information. `chol()` succeeds only if \eqn{-H} is
#' positive definite, which it is at a proper interior maximum. If it fails
#' (for example at the no-zero-inflation boundary, where the information
#' about \eqn{\gamma} is ~0), the covariance is returned as `NA` with a
#' warning rather than as meaningless numbers.
#' @noRd
zinb_vcov <- function(hessian) {
  V <- tryCatch(chol2inv(chol(-hessian)), error = function(e) NULL)
  if (is.null(V)) {
    warning("the observed information matrix is not positive definite; ",
            "standard errors are not available.", call. = FALSE)
    V <- matrix(NA_real_, nrow(hessian), ncol(hessian))
  }
  dimnames(V) <- dimnames(hessian)
  V
}
