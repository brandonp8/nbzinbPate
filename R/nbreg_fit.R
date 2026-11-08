# Matrix-level engine for NB2 regression (Sec 3.1 of the project PDF).
#
# This file works only with matrices and vectors (X, y, offset, weights);
# the user-facing formula interface nbreg() lives in R/nbreg.R. Keeping the
# engine formula-free lets the ZINB count-part M-step (Sec 3.2) call it
# directly with prior weights a_i = 1 - tau_i.
#
# Pieces, matching Sec 3.1:
#   (a) nb_irls_beta()    IRLS for beta with theta held fixed
#   (b) nb_newton_theta() Newton-Raphson for log(theta), beta held fixed
#   (c) nbreg_fit()       outer loop alternating (a) and (b)
#   (d) standard errors, computed at the end of nbreg_fit()


#' IRLS working weights for the NB2 / Poisson log-link GLM
#'
#' Sec 3.1(a):
#' \deqn{w_i = a_i \frac{(\partial\mu_i/\partial\eta_i)^2}{V(\mu_i)}
#'           = a_i \frac{\theta\mu_i}{\theta + \mu_i}.}
#' Derivation: with the log link, \eqn{\mu = e^\eta} so
#' \eqn{\partial\mu/\partial\eta = \mu}; with \eqn{V(\mu) = \mu + \mu^2/\theta},
#' \eqn{\mu^2 / (\mu + \mu^2/\theta) = \theta\mu/(\theta+\mu)}.
#'
#' `theta = Inf` gives the Poisson limit \eqn{V(\mu) = \mu}, so
#' \eqn{w_i = a_i\mu_i}. This is how the Poisson starting fit of Sec 3.1(c)
#' reuses the same IRLS code. (Plugging Inf into the formula directly would
#' give Inf/Inf = NaN, hence the special case.)
#'
#' @param mu Current means.
#' @param theta Dispersion (scalar, may be `Inf`).
#' @param weights Prior weights \eqn{a_i}.
#' @noRd
nb_working_weights <- function(mu, theta, weights) {
  if (is.infinite(theta)) {
    weights * mu
  } else {
    weights * theta * mu / (theta + mu)
  }
}


#' Log-likelihood used to monitor IRLS: NB2, or Poisson when theta = Inf
#'
#' Poisson: \eqn{\sum_i a_i [y_i\log\mu_i - \mu_i - \log y_i!]}.
#' @noRd
nb_or_pois_loglik <- function(y, mu, theta, weights) {
  if (is.infinite(theta)) {
    sum(weights * (y * log(mu) - mu - lgamma(y + 1)))
  } else {
    nb_loglik(y, mu, theta, weights)
  }
}


#' IRLS for beta with theta held fixed (Sec 3.1(a))
#'
#' Repeats the Fisher-scoring / weighted-least-squares update
#' \deqn{z_i = (\eta_i - o_i) + \frac{y_i - \mu_i}{\mu_i}, \qquad
#'       \beta^{(t+1)} = (X^\top W X)^{-1} X^\top W z}
#' until \eqn{\max_j |\beta_j^{(t+1)} - \beta_j^{(t)}| < } `tol`.
#'
#' Why this is Fisher scoring: for the log link the score is
#' \eqn{X^\top W z^*} with \eqn{z^*_i = (y_i-\mu_i)/\mu_i} and the expected
#' information is \eqn{X^\top W X}. The scoring step
#' \eqn{\beta + (X^\top W X)^{-1} X^\top W z^*} equals
#' \eqn{(X^\top W X)^{-1} X^\top W (X\beta + z^*)}, and
#' \eqn{X\beta = \eta - o}, which gives the working response above.
#'
#' @param X n x p design matrix.
#' @param y Counts.
#' @param theta Fixed dispersion (`Inf` = Poisson).
#' @param offset Offset \eqn{o_i} (length n).
#' @param weights Prior weights \eqn{a_i} (length n).
#' @param beta_start Optional starting coefficients (a warm start, used by
#'   the outer loop and by EM). If `NULL`, start from
#'   \eqn{\mu_i = y_i + 0.1}, the same start glm() uses for Poisson.
#' @param tol Convergence tolerance on the largest absolute change in beta.
#' @param maxit Maximum number of IRLS iterations.
#' @return A list with `coefficients`, `eta`, `mu`, `w` (working weights at
#'   the final `mu`), `loglik`, `iter`, `converged`.
#' @noRd
nb_irls_beta <- function(X, y, theta, offset = rep(0, length(y)),
                         weights = rep(1, length(y)), beta_start = NULL,
                         tol = 1e-10, maxit = 100) {

  # ---- Starting values ------------------------------------------------------
  if (is.null(beta_start)) {
    # No beta yet: start from fitted means mu = y + 0.1 (the +0.1 avoids
    # log(0) for zero counts). The first WLS step then produces a beta.
    mu   <- y + 0.1
    eta  <- log(mu)
    beta <- NULL
  } else {
    beta <- beta_start
    eta  <- drop(X %*% beta) + offset            # eta = X beta + o
    mu   <- exp(eta)                             # log link: mu = exp(eta)
  }
  ll_old <- nb_or_pois_loglik(y, mu, theta, weights)

  converged <- FALSE
  for (iter in seq_len(maxit)) {

    # ---- One IRLS step (Sec 3.1(a)) -----------------------------------------
    w <- nb_working_weights(mu, theta, weights)  # w_i = a_i theta mu_i/(theta+mu_i)
    z <- (eta - offset) + (y - mu) / mu          # z_i = (eta_i - o_i) + (y_i - mu_i)/mu_i
    beta_new <- wls_solve(X, w, z)               # (X'WX)^{-1} X'Wz

    eta_new <- drop(X %*% beta_new) + offset
    mu_new  <- exp(eta_new)
    ll_new  <- nb_or_pois_loglik(y, mu_new, theta, weights)

    # ---- Step-halving safeguard ----------------------------------------------
    # Fisher scoring almost always increases the log-likelihood near the
    # optimum, but from a poor start (or with a non-canonical variance) a
    # full step can overshoot. If the log-likelihood went down (or became
    # non-finite), move back halfway toward the previous beta, up to 30
    # times. This only shortens the step along the same direction; it is
    # not an optimizer.
    if (!is.null(beta)) {
      halvings <- 0
      while ((!is.finite(ll_new) || ll_new < ll_old) && halvings < 30) {
        beta_new <- (beta + beta_new) / 2
        eta_new  <- drop(X %*% beta_new) + offset
        mu_new   <- exp(eta_new)
        ll_new   <- nb_or_pois_loglik(y, mu_new, theta, weights)
        halvings <- halvings + 1
      }
    }

    # ---- Convergence check ---------------------------------------------------
    # Stop when the largest change in any coefficient, max_j |beta_j^(t+1) -
    # beta_j^(t)|, is below tol (Sec 3.1(a) allows "the change in beta").
    # Why beta and not the log-likelihood: the log-likelihood is flat
    # (quadratic) at its maximum, so a log-likelihood change of 1e-7 can
    # still leave beta wrong in the 6th digit. Also, the log link is not
    # the canonical link for the NB, so Fisher scoring converges only
    # linearly and needs a criterion on beta itself to reach 1e-8 accuracy.
    # (No check on the first pass when started from mu = y + 0.1, because
    # there is no previous beta to compare with.)
    delta <- if (is.null(beta)) Inf else max(abs(beta_new - beta))

    beta   <- beta_new
    eta    <- eta_new
    mu     <- mu_new
    ll_old <- ll_new

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  list(coefficients = beta,
       eta          = eta,
       mu           = mu,
       w            = nb_working_weights(mu, theta, weights),  # for SEs (Sec 3.1(d))
       loglik       = ll_old,
       iter         = iter,
       converged    = converged)
}


#' Score and second derivative of the NB2 log-likelihood in theta
#'
#' With beta (hence mu) held fixed, Sec 3.1(b):
#' \deqn{U(\theta) = \sum_i a_i\Big[\psi(y_i+\theta) - \psi(\theta) + \log\theta
#'   + 1 - \log(\theta+\mu_i) - \frac{y_i+\theta}{\theta+\mu_i}\Big],}
#' \deqn{H(\theta) = \sum_i a_i\Big[\psi'(y_i+\theta) - \psi'(\theta)
#'   + \frac{1}{\theta} - \frac{2}{\theta+\mu_i}
#'   + \frac{y_i+\theta}{(\theta+\mu_i)^2}\Big],}
#' where \eqn{\psi} = `digamma()` and \eqn{\psi'} = `trigamma()`.
#'
#' Each formula is coded in an algebraically identical regrouped form that
#' avoids subtracting nearly equal numbers when theta is large:
#' * \eqn{\log\theta + 1 - \log(\theta+\mu) - (y+\theta)/(\theta+\mu)
#'   = -\log(1+\mu/\theta) + (\mu - y)/(\theta+\mu)}, using
#'   \eqn{1 - (y+\theta)/(\theta+\mu) = (\mu-y)/(\theta+\mu)};
#' * \eqn{1/\theta - 2/(\theta+\mu) + (y+\theta)/(\theta+\mu)^2
#'   = \mu/(\theta(\theta+\mu)) + (y-\mu)/(\theta+\mu)^2}, using
#'   \eqn{1/\theta - 1/(\theta+\mu) = \mu/(\theta(\theta+\mu))} and
#'   \eqn{-1/(\theta+\mu) + (y+\theta)/(\theta+\mu)^2 = (y-\mu)/(\theta+\mu)^2}.
#'
#' @param y Counts.
#' @param mu Fitted means (fixed).
#' @param theta Current dispersion.
#' @param weights Prior weights \eqn{a_i}.
#' @return `list(U = , H = )`.
#' @noRd
nb_theta_derivs <- function(y, mu, theta, weights) {
  tm <- theta + mu                                        # theta + mu_i

  # U(theta): psi(y+theta) - psi(theta) + [log theta + 1 - log(theta+mu)
  #           - (y+theta)/(theta+mu)]   (bracket in regrouped form)
  U <- sum(weights * (digamma(y + theta) - digamma(theta) -
                        log1p(mu / theta) + (mu - y) / tm))

  # H(theta): psi'(y+theta) - psi'(theta) + [1/theta - 2/(theta+mu)
  #           + (y+theta)/(theta+mu)^2]  (bracket in regrouped form)
  H <- sum(weights * (trigamma(y + theta) - trigamma(theta) +
                        mu / (theta * tm) + (y - mu) / tm^2))

  list(U = U, H = H)
}


#' Newton-Raphson for log(theta) with beta held fixed (Sec 3.1(b))
#'
#' Works on \eqn{\phi = \log\theta}, which keeps \eqn{\theta = e^\phi > 0}
#' automatically. By the chain rule (\eqn{d\theta/d\phi = \theta}):
#' \deqn{U_\phi = \theta U(\theta), \qquad
#'       H_\phi = \frac{d}{d\phi}[\theta U(\theta)] = \theta^2 H(\theta) + \theta U(\theta),}
#' and each Newton step is \eqn{\phi \leftarrow \phi - U_\phi / H_\phi}.
#'
#' Two safeguards keep the iteration on track without turning it into a
#' general optimizer:
#' 1. If \eqn{H_\phi \ge 0} the log-likelihood is not locally concave in
#'    \eqn{\phi} and the Newton step would point the wrong way; then we
#'    simply move \eqn{\phi} one unit in the uphill direction
#'    (\eqn{\mathrm{sign}(U_\phi)}), i.e. multiply or divide theta by e.
#' 2. Step-halving: if a step lowers the log-likelihood, halve it (up to
#'    30 times), exactly as in nb_irls_beta().
#'
#' @inheritParams nb_theta_derivs
#' @param theta_start Starting value for theta.
#' @param tol Convergence tolerance on \eqn{|\Delta\phi|}.
#' @param maxit Maximum number of Newton iterations.
#' @return A list with `theta`, `loglik`, `U`, `H` (in the theta scale, at
#'   the final theta; \eqn{-H} is the observed information used for
#'   \eqn{SE(\hat\theta)} in Sec 3.1(d)), `iter`, `converged`.
#' @noRd
nb_newton_theta <- function(y, mu, theta_start, weights = rep(1, length(y)),
                            tol = 1e-10, maxit = 100) {
  phi <- log(theta_start)                        # phi = log(theta)
  ll  <- nb_loglik(y, mu, exp(phi), weights)

  converged <- FALSE
  for (iter in seq_len(maxit)) {
    theta <- exp(phi)
    d     <- nb_theta_derivs(y, mu, theta, weights)
    U_phi <- theta * d$U                         # U_phi = theta U
    H_phi <- theta^2 * d$H + theta * d$U         # H_phi = theta^2 H + theta U

    # Newton step phi <- phi - U_phi / H_phi (safeguard 1 if H_phi >= 0)
    step <- if (H_phi < 0) -U_phi / H_phi else sign(U_phi)

    # Safeguard 2: step-halving if the log-likelihood does not increase.
    phi_new  <- phi + step
    ll_new   <- nb_loglik(y, mu, exp(phi_new), weights)
    halvings <- 0
    while ((!is.finite(ll_new) || ll_new < ll) && halvings < 30) {
      step     <- step / 2
      phi_new  <- phi + step
      ll_new   <- nb_loglik(y, mu, exp(phi_new), weights)
      halvings <- halvings + 1
    }
    if (!is.finite(ll_new) || ll_new < ll) {
      # Even a tiny step does not improve: we are at the maximum to
      # machine precision, so keep the current phi.
      phi_new <- phi
      ll_new  <- ll
    }

    delta <- abs(phi_new - phi)                  # |change in log(theta)|
    phi   <- phi_new
    ll    <- ll_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  theta <- exp(phi)
  d     <- nb_theta_derivs(y, mu, theta, weights)  # U, H at the final theta
  list(theta     = theta,
       loglik    = ll,
       U         = d$U,
       H         = d$H,
       iter      = iter,
       converged = converged)
}


#' Fit an NB2 regression from matrices: the outer loop of Sec 3.1(c)-(d)
#'
#' Alternates IRLS for \eqn{\beta} (theta fixed) and Newton-Raphson for
#' \eqn{\log\theta} (beta fixed) until neither the log-likelihood nor the
#' parameters change, then computes standard errors.
#'
#' This is the engine behind [nbreg()] and the count-part M-step of the
#' ZINB EM algorithm (Sec 3.2), which calls it with prior weights
#' \eqn{a_i = 1-\tau_i} and warm starts.
#'
#' @param X n x p design matrix.
#' @param y Counts.
#' @param offset Offset (length n).
#' @param weights Prior weights \eqn{a_i} (length n).
#' @param control Output of [nbreg_control()].
#' @param beta_start,theta_start Optional warm starts. If `NULL`, beta comes
#'   from a Poisson IRLS fit and theta from a moment estimate.
#' @param warn If `TRUE`, warn about non-convergence or a huge theta.
#'   (EM turns this off for its inner fits and checks things itself.)
#' @return A list with `coefficients`, `vcov`, `theta`, `SE.theta`,
#'   `loglik`, `mu`, `eta`, `w`, `iter`, `converged`.
#' @noRd
nbreg_fit <- function(X, y, offset = rep(0, length(y)),
                      weights = rep(1, length(y)), control = nbreg_control(),
                      beta_start = NULL, theta_start = NULL, warn = TRUE) {

  # ---- Starting values (Sec 3.1(c)) -------------------------------------------
  if (is.null(beta_start)) {
    # Poisson IRLS fit = the same IRLS with theta = Inf (V(mu) = mu).
    pois <- nb_irls_beta(X, y, theta = Inf, offset = offset, weights = weights,
                         tol = control$irls_tol, maxit = control$irls_maxit)
    beta <- pois$coefficients
  } else {
    beta <- beta_start
  }
  mu <- exp(drop(X %*% beta) + offset)

  if (is.null(theta_start)) {
    # Moment estimate. Under NB2, E[(y/mu - 1)^2] = Var(y)/mu^2
    # = 1/mu + 1/theta, which is >= 1/theta, so 1/theta is roughly the
    # average of (y/mu - 1)^2. Setting 1/theta equal to that weighted
    # average gives theta0 = sum(a) / sum(a (y/mu - 1)^2) (the same start
    # MASS uses). It is always positive.
    theta <- sum(weights) / sum(weights * (y / mu - 1)^2)
  } else {
    theta <- theta_start
  }
  ll <- nb_loglik(y, mu, theta, weights)

  if (control$trace) {
    cat(sprintf("nbreg start: loglik = %.8f, theta = %.6g\n", ll, theta))
  }

  # ---- Outer loop: alternate (a) and (b) (Sec 3.1(c)) ----------------------------
  converged <- FALSE
  for (iter in seq_len(control$maxit)) {

    # (a) IRLS for beta, theta held fixed; warm start at the current beta.
    irls <- nb_irls_beta(X, y, theta = theta, offset = offset,
                         weights = weights, beta_start = beta,
                         tol = control$irls_tol, maxit = control$irls_maxit)

    # (b) Newton-Raphson for log(theta), beta (so mu) held fixed;
    #     warm start at the current theta.
    newt <- nb_newton_theta(y, irls$mu, theta_start = theta, weights = weights,
                            tol = control$theta_tol,
                            maxit = control$theta_maxit)

    # Size of this outer iteration's move in all parameters.
    delta_par <- max(abs(irls$coefficients - beta),
                     abs(log(newt$theta) - log(theta)))
    delta_ll  <- abs(newt$loglik - ll)

    beta  <- irls$coefficients
    mu    <- irls$mu
    theta <- newt$theta
    ll    <- newt$loglik                # = nb_loglik(y, mu, theta, weights)

    if (control$trace) {
      cat(sprintf("nbreg iter %3d: loglik = %.8f, theta = %.6g\n",
                  iter, ll, theta))
    }

    # Converged when the log-likelihood changes by < tol (Sec 3.1(c)) and,
    # for the reason explained in nb_irls_beta(), the parameters have also
    # stopped moving.
    if (delta_ll < control$tol && delta_par < control$tol) {
      converged <- TRUE
      break
    }
  }

  # ---- Standard errors (Sec 3.1(d)) ---------------------------------------------
  eta <- drop(X %*% beta) + offset
  # Working weights at the final (mu, theta): w_i = a_i theta mu_i/(theta+mu_i).
  w <- nb_working_weights(mu, theta, weights)
  # Var(beta-hat) = (X' W-hat X)^{-1}, treating theta as known (as glm.nb).
  V <- xtwx_inverse(X, w)
  # SE(theta-hat) = {-H(theta-hat)}^{-1/2}, H from Sec 3.1(b).
  H <- nb_theta_derivs(y, mu, theta, weights)$H
  se_theta <- if (H < 0) 1 / sqrt(-H) else NaN

  # ---- Edge-case warnings --------------------------------------------------------
  if (warn) {
    if (!converged) {
      warning("nbreg did not converge in ", control$maxit,
              " outer iterations; results may be unreliable. ",
              "Try increasing 'maxit' in nbreg_control().", call. = FALSE)
    }
    if (theta > control$theta_warn) {
      warning("theta-hat = ", format(theta, digits = 3), " is very large: ",
              "the data show little or no overdispersion (the Poisson ",
              "limit), so a Poisson model may be adequate.", call. = FALSE)
    }
  }

  list(coefficients = beta,
       vcov         = V,
       theta        = theta,
       SE.theta     = se_theta,
       loglik       = ll,
       mu           = mu,
       eta          = eta,
       w            = w,
       iter         = iter,
       converged    = converged)
}
