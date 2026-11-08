# Log-likelihood functions for the NB2 and ZINB models.
#
# These are used (1) to monitor convergence of the NB outer loop and the EM
# algorithm, (2) to report logLik / AIC / BIC, and (3) in the tests, where
# they are checked against stats::dnbinom().
#
# All references "Sec x.y" are to project-docs/STA651_Final_Project.pdf.


#' Log of the NB2 probability mass function, one value per observation
#'
#' Implements the log of the NB2 pmf from Sec 2.1:
#' \deqn{\log f_{NB}(y; \mu, \theta) = \log\Gamma(y+\theta) - \log\Gamma(\theta)
#'   - \log y! + \theta\log\theta + y\log\mu - (\theta+y)\log(\theta+\mu).}
#'
#' @param y Non-negative integer counts (vector).
#' @param mu Means \eqn{\mu_i > 0} (vector, same length as `y`).
#' @param theta Dispersion \eqn{\theta > 0} (scalar).
#' @return Vector of \eqn{\log f_{NB}(y_i; \mu_i, \theta)}.
#' @noRd
nb_logdens <- function(y, mu, theta) {
  # The three Gamma-function terms: log Gamma(y+theta) - log Gamma(theta) - log y!
  # (log y! is written lgamma(y + 1) because y! = Gamma(y + 1)).
  gamma_terms <- lgamma(y + theta) - lgamma(theta) - lgamma(y + 1)

  # The remaining terms  theta*log(theta) + y*log(mu) - (theta+y)*log(theta+mu)
  # are regrouped, purely algebraically, as
  #     theta * log(theta / (theta + mu))  +  y * log(mu / (theta + mu))
  #   = -theta * log1p(mu / theta)         +  y * (log(mu) - log(theta + mu)).
  # Why: when theta is very large (near the Poisson limit) theta*log(theta)
  # and theta*log(theta+mu) are two huge, almost equal numbers, and
  # subtracting them loses most significant digits. log1p(x) = log(1 + x)
  # is computed accurately even when x is tiny, so the regrouped form stays
  # accurate for any theta.
  gamma_terms - theta * log1p(mu / theta) + y * (log(mu) - log(theta + mu))
}


#' NB2 log-likelihood with optional prior weights
#'
#' Implements Sec 2.1:
#' \deqn{\ell(\beta,\theta) = \sum_i a_i \log f_{NB}(y_i; \mu_i, \theta).}
#' The prior weights \eqn{a_i} are all 1 for an ordinary fit; the ZINB
#' count-part M-step uses \eqn{a_i = 1 - \tau_i} (Sec 3.2).
#'
#' @inheritParams nb_logdens
#' @param weights Prior weights \eqn{a_i \ge 0}; default all 1.
#' @return Scalar log-likelihood.
#' @noRd
nb_loglik <- function(y, mu, theta, weights = rep(1, length(y))) {
  sum(weights * nb_logdens(y, mu, theta))
}


#' ZINB observed-data log-likelihood
#'
#' Implements Sec 2.2:
#' \deqn{\ell = \sum_{y_i=0} \log[\pi_i + (1-\pi_i) f_{NB}(0;\mu_i,\theta)]
#'   + \sum_{y_i>0} [\log(1-\pi_i) + \log f_{NB}(y_i;\mu_i,\theta)].}
#' This is the quantity EM must never decrease, so it is evaluated at every
#' EM iteration.
#'
#' @inheritParams nb_logdens
#' @param pi Structural-zero probabilities \eqn{\pi_i \in (0,1)}.
#' @return Scalar log-likelihood.
#' @noRd
zinb_loglik <- function(y, mu, pi, theta) {
  is_zero <- (y == 0)

  # f_NB(0; mu, theta) = (theta / (theta + mu))^theta  (Sec 2.2), written as
  # exp(-theta * log1p(mu/theta)) for the same accuracy reason as above.
  f0 <- exp(-theta * log1p(mu[is_zero] / theta))

  # Zero observations: log[ pi + (1 - pi) * f_NB(0) ]
  ll_zero <- log(pi[is_zero] + (1 - pi[is_zero]) * f0)

  # Positive observations: log(1 - pi) + log f_NB(y).
  # log1p(-pi) = log(1 - pi), accurate when pi is close to 0.
  ll_pos <- log1p(-pi[!is_zero]) +
    nb_logdens(y[!is_zero], mu[!is_zero], theta)

  sum(ll_zero) + sum(ll_pos)
}
