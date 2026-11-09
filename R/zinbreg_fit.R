# Matrix-level engine for ZINB regression: the EM algorithm of Sec 3.2.
#
# Why EM: the ZINB log-likelihood (Sec 2.2) has log[pi + (1-pi) f(0)] for
# every zero, a log of a SUM, so it is not an exponential-family GLM
# log-likelihood and IRLS cannot be applied to it directly. EM introduces
# the latent indicator D_i = 1 if observation i is a structural zero. If
# the D_i were observed, the complete-data log-likelihood would split into
#   sum_i [D_i log pi_i + (1-D_i) log(1-pi_i)]       (a logistic regression)
# + sum_i (1-D_i) log f_NB(y_i; mu_i, theta)          (a weighted NB fit),
# each of which IS fitted by IRLS. EM replaces D_i by its conditional
# expectation tau_i (E-step) and maximizes both parts (M-step).
#
# The formula interface zinbreg() lives in R/zinbreg.R.


#' E-step: posterior probability of being a structural zero (Sec 3.2)
#'
#' \deqn{\tau_i = P(D_i = 1 \mid y_i) =
#'   \frac{\pi_i}{\pi_i + (1-\pi_i) f_{NB}(0;\mu_i,\theta)} \text{ if } y_i = 0,
#'   \qquad \tau_i = 0 \text{ if } y_i > 0.}
#'
#' Computed in an algebraically identical but stable form. Dividing top and
#' bottom by \eqn{\pi_i}:
#' \deqn{\tau_i = \frac{1}{1 + \frac{1-\pi_i}{\pi_i} f_{NB}(0)}
#'   = \frac{1}{1 + \exp\{-\eta^z_i + \log f_{NB}(0)\}}
#'   = \mathrm{logit}^{-1}\{\eta^z_i - \log f_{NB}(0)\},}
#' because \eqn{(1-\pi_i)/\pi_i = e^{-\eta^z_i}}. `plogis()` evaluates this
#' without forming \eqn{\pi_i} (which could round to 0 or 1).
#'
#' @param y Counts.
#' @param mu Count-model means \eqn{\mu_i}.
#' @param eta_zero Zero-model linear predictors \eqn{\eta^z_i = z_i^\top\gamma}.
#' @param theta NB dispersion.
#' @return Vector of \eqn{\tau_i}.
#' @noRd
zinb_estep <- function(y, mu, eta_zero, theta) {
  tau <- numeric(length(y))                    # tau_i = 0 when y_i > 0
  z0  <- (y == 0)
  tau[z0] <- plogis(eta_zero[z0] - nb_log_f0(mu[z0], theta))
  tau
}


#' Fit a ZINB regression from matrices by EM (Sec 3.2)
#'
#' Each EM iteration:
#' 1. **E-step**: \eqn{\tau_i} from `zinb_estep()`.
#' 2. **M-step, zero part**: maximize
#'    \eqn{\sum_i[\tau_i\log\pi_i + (1-\tau_i)\log(1-\pi_i)]} over \eqn{\gamma}
#'    by fractional logistic IRLS (`logit_irls()`).
#' 3. **M-step, count part**: maximize
#'    \eqn{\sum_i (1-\tau_i)\log f_{NB}(y_i;\mu_i,\theta)} over
#'    \eqn{(\beta,\theta)}: exactly the NB fit of Sec 3.1 with prior weights
#'    \eqn{a_i = 1-\tau_i}, so `nbreg_fit()` is reused.
#'
#' Both M-steps are *full* fits (run to convergence), warm-started at the
#' current values, so each is cheap after the first few EM iterations.
#'
#' Why EM never decreases the observed-data log-likelihood \eqn{\ell}: each
#' M-step maximizes the expected complete-data log-likelihood
#' \eqn{Q(\cdot \mid \text{current})}, and \eqn{\ell(\text{new}) - \ell(\text{old})
#' \ge Q(\text{new}) - Q(\text{old}) \ge 0} (Dempster, Laird and Rubin, 1977).
#' The code checks this at every iteration; a decrease beyond rounding error
#' would indicate a bug, so it triggers a warning.
#'
#' @param X n x p count-model design matrix.
#' @param Z n x q zero-model design matrix.
#' @param y Counts.
#' @param offset Count-model offset.
#' @param control Output of [zinbreg_control()].
#' @return A list with `coefficients` (list of `count` and `zero`),
#'   `theta`, `loglik`, `loglik_path` (log-likelihood at every iteration),
#'   `mu`, `eta_count`, `eta_zero`, `pi`, `tau`, `iter`, `converged`, and
#'   `boundary` (`TRUE` if every \eqn{\hat\pi_i < 10^{-8}}: no zero inflation).
#' @noRd
zinbreg_fit <- function(X, Z, y, offset = rep(0, length(y)),
                        control = zinbreg_control()) {

  # Settings for the full M-step fits (see zinbreg_control()).
  nb_ctrl <- nbreg_control(tol = control$mstep_tol,
                           maxit = control$mstep_maxit,
                           irls_tol = control$mstep_tol,
                           irls_maxit = control$mstep_maxit,
                           theta_tol = control$mstep_tol,
                           theta_maxit = control$mstep_maxit)

  # ---- Starting values (Sec 3.2) ------------------------------------------------
  # (beta, theta) from the ordinary NB fit of Sec 3.1.
  nb    <- nbreg_fit(X, y, offset = offset, control = nb_ctrl, warn = FALSE)
  beta  <- nb$coefficients
  theta <- nb$theta
  # gamma from our own logistic IRLS of the indicator I(y_i = 0) on z_i.
  gamma <- logit_irls(Z, as.numeric(y == 0), tol = control$mstep_tol,
                      maxit = control$mstep_maxit)$coefficients

  eta_count <- drop(X %*% beta) + offset       # log(mu) = x'beta + o
  mu        <- exp(eta_count)
  eta_zero  <- drop(Z %*% gamma)               # logit(pi) = z'gamma
  ll        <- zinb_loglik(y, mu, eta_zero, theta)
  ll_path   <- ll

  if (control$trace) {
    cat(sprintf("zinbreg start: loglik = %.8f, theta = %.6g\n", ll, theta))
  }

  converged <- FALSE
  for (iter in seq_len(control$maxit)) {

    # ---- E-step ------------------------------------------------------------------
    tau <- zinb_estep(y, mu, eta_zero, theta)

    # ---- M-step, zero part: logistic IRLS with fractional response tau ----------
    zero_fit  <- logit_irls(Z, tau, gamma_start = gamma,
                            tol = control$mstep_tol,
                            maxit = control$mstep_maxit)

    # ---- M-step, count part: NB fit with prior weights a_i = 1 - tau_i ----------
    count_fit <- nbreg_fit(X, y, offset = offset, weights = 1 - tau,
                           control = nb_ctrl, beta_start = beta,
                           theta_start = theta, warn = FALSE)

    # Largest change in any parameter (beta, gamma, log theta) this iteration.
    delta_par <- max(abs(count_fit$coefficients - beta),
                     abs(zero_fit$coefficients - gamma),
                     abs(log(count_fit$theta) - log(theta)))

    beta      <- count_fit$coefficients
    theta     <- count_fit$theta
    gamma     <- zero_fit$coefficients
    eta_count <- count_fit$eta
    mu        <- count_fit$mu
    eta_zero  <- zero_fit$eta

    # ---- Observed-data log-likelihood and the monotonicity check ---------------
    ll_new  <- zinb_loglik(y, mu, eta_zero, theta)
    ll_path <- c(ll_path, ll_new)
    # Allow a decrease only at the level of floating-point rounding.
    if (ll_new < ll - 1e-8 * max(1, abs(ll))) {
      warning(sprintf(paste0("EM log-likelihood decreased at iteration %d ",
                             "(%.10f -> %.10f); this should not happen."),
                      iter, ll, ll_new), call. = FALSE)
    }
    delta_ll <- abs(ll_new - ll)
    ll       <- ll_new

    if (control$trace) {
      cat(sprintf("zinbreg iter %4d: loglik = %.10f, theta = %.6g\n",
                  iter, ll, theta))
    }

    # ---- Boundary check: is the zero-inflation component vanishing? -------------
    # If the data show no zero inflation, the likelihood keeps increasing
    # as every pi_i -> 0, i.e. some gamma_j -> -Inf: the MLE is on the
    # boundary of the parameter space and gamma never settles. We call the
    # zero component degenerate once every pi_i < 1e-8, i.e. every
    # eta_zero_i < logit(1e-8) (about -18.4). The log-likelihood is then
    # that of the NB model to ~8 digits.
    boundary <- max(eta_zero) < qlogis(1e-8)

    # ---- Convergence (Sec 3.2) -----------------------------------------------------
    # The observed-data log-likelihood must change by < tol (as the project
    # specifies); as in nbreg, we also require the parameters to have
    # settled, because EM converges linearly and the log-likelihood is flat
    # at its maximum -- unless we are at the boundary, where gamma
    # cannot settle (it drifts toward -Inf) but the log-likelihood has.
    if (delta_ll < control$tol && (delta_par < control$par_tol || boundary)) {
      converged <- TRUE
      break
    }
  }

  list(coefficients = list(count = beta, zero = gamma),
       theta        = theta,
       loglik       = ll,
       loglik_path  = ll_path,
       mu           = mu,
       eta_count    = eta_count,
       eta_zero     = eta_zero,
       pi           = plogis(eta_zero),
       tau          = zinb_estep(y, mu, eta_zero, theta),
       iter         = iter,
       converged    = converged,
       boundary     = boundary)
}
