# Logistic regression with a fractional response, fitted by IRLS.
#
# This is the zero-part M-step of the ZINB EM algorithm (Sec 3.2 of the
# project PDF). Given the E-step probabilities tau_i in [0, 1], it maximizes
#
#     Q(gamma) = sum_i [ tau_i log(pi_i) + (1 - tau_i) log(1 - pi_i) ],
#     logit(pi_i) = z_i' gamma,
#
# i.e. a Bernoulli log-likelihood whose "response" is a probability rather
# than a 0/1 outcome. glm() is not allowed in R/, so we code the IRLS.


#' Fractional-response logistic log-likelihood Q(gamma)
#'
#' \eqn{Q = \sum_i [\tau_i\log\pi_i + (1-\tau_i)\log(1-\pi_i)]}, evaluated
#' from the linear predictor \eqn{\eta_i = z_i^\top\gamma}.
#' `plogis(eta, log.p = TRUE)` returns \eqn{\log\pi} and
#' `plogis(-eta, log.p = TRUE)` returns \eqn{\log(1-\pi)} without ever
#' forming \eqn{\pi} itself, so neither underflows to `log(0) = -Inf` when
#' \eqn{|\eta|} is large.
#' @noRd
logit_frac_loglik <- function(tau, eta) {
  sum(tau * plogis(eta, log.p = TRUE) +
        (1 - tau) * plogis(-eta, log.p = TRUE))
}


#' IRLS for logistic regression with a fractional response (Sec 3.2)
#'
#' Iterates, with \eqn{\pi_i = \mathrm{logit}^{-1}(z_i^\top\gamma)},
#' \deqn{w_i = \pi_i(1-\pi_i), \qquad
#'       \tilde z_i = z_i^\top\gamma + \frac{\tau_i - \pi_i}{\pi_i(1-\pi_i)}, \qquad
#'       \gamma^{(t+1)} = (Z^\top W Z)^{-1} Z^\top W \tilde z.}
#'
#' Why this is Newton-Raphson: the score of \eqn{Q} is \eqn{Z^\top(\tau-\pi)}
#' (since \eqn{d\log\pi/d\eta = 1-\pi} and \eqn{d\log(1-\pi)/d\eta = -\pi})
#' and the Hessian is \eqn{-Z^\top W Z}. The Newton step
#' \eqn{\gamma + (Z^\top W Z)^{-1} Z^\top(\tau-\pi)} rearranges to the WLS
#' form above. The logit is the canonical link, so the Hessian does not
#' involve \eqn{\tau}: IRLS here is exact Newton and converges quadratically.
#'
#' @param Z n x q design matrix of the zero-inflation model.
#' @param tau Fractional responses in \[0, 1\] (the E-step \eqn{\tau_i}).
#' @param gamma_start Optional warm start (EM passes the previous gamma).
#'   If `NULL`, start from \eqn{\pi_i = (\tau_i + 0.5)/2}, the same start
#'   glm() uses for binomial data, which keeps \eqn{\pi} away from 0 and 1.
#' @param tol Convergence tolerance on \eqn{\max_j|\Delta\gamma_j|}.
#' @param maxit Maximum number of iterations.
#' @return A list with `coefficients`, `eta`, `pi`, `loglik` (value of
#'   \eqn{Q}), `iter`, `converged`.
#' @noRd
logit_irls <- function(Z, tau, gamma_start = NULL, tol = 1e-10, maxit = 100) {

  # ---- Starting values ------------------------------------------------------
  if (is.null(gamma_start)) {
    pi    <- (tau + 0.5) / 2
    eta   <- qlogis(pi)                       # eta = logit(pi)
    gamma <- NULL
  } else {
    gamma <- gamma_start
    eta   <- drop(Z %*% gamma)
    pi    <- plogis(eta)                      # pi = 1 / (1 + exp(-eta))
  }
  ll_old <- logit_frac_loglik(tau, eta)

  converged <- FALSE
  for (iter in seq_len(maxit)) {

    # ---- One IRLS (= Newton) step -------------------------------------------
    # w_i = pi_i (1 - pi_i). If |eta_i| is huge, pi_i(1 - pi_i) can
    # underflow to 0 and the division below would give 0/0; flooring w at
    # machine epsilon prevents that. Such an observation then gets a
    # negligible weight, which is what its true weight (~0) means anyway.
    w <- pmax(pi * (1 - pi), .Machine$double.eps)
    z <- eta + (tau - pi) / w                 # z_i = z_i'gamma + (tau_i - pi_i)/w_i
    gamma_new <- wls_solve(Z, w, z)           # (Z'WZ)^{-1} Z'W z

    eta_new <- drop(Z %*% gamma_new)
    ll_new  <- logit_frac_loglik(tau, eta_new)

    # ---- Step-halving safeguard (same as nb_irls_beta) -----------------------
    # Q is concave, so a full Newton step rarely overshoots, but it can from
    # a poor start; if Q went down by more than rounding error (see
    # ll_decreased()), halve the step. If that fails 30 times, stop and
    # report non-convergence rather than faking convergence.
    if (!is.null(gamma)) {
      halvings <- 0
      while (ll_decreased(ll_new, ll_old) && halvings < 30) {
        gamma_new <- (gamma + gamma_new) / 2
        eta_new   <- drop(Z %*% gamma_new)
        ll_new    <- logit_frac_loglik(tau, eta_new)
        halvings  <- halvings + 1
      }
      if (ll_decreased(ll_new, ll_old)) break
    }

    # ---- Convergence: largest change in gamma (no check on the first
    #      pass from the default start, as there is no previous gamma) -------
    delta <- if (is.null(gamma)) Inf else max(abs(gamma_new - gamma))

    gamma  <- gamma_new
    eta    <- eta_new
    pi     <- plogis(eta)
    ll_old <- ll_new

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  list(coefficients = gamma,
       eta          = eta,
       pi           = pi,
       loglik       = ll_old,
       iter         = iter,
       converged    = converged)
}
