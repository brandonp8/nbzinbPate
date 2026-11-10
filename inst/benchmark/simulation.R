# Simulation study (project Sec 5.1): does zinbreg() recover known
# parameters, and are its standard errors and Wald intervals accurate?
#
# Design: R = 100 data sets of size n = 500 from the ZINB model
#   count part:  log(mu_i)   = b0 + b1 x1_i + b2 x2_i,  x1 ~ N(0,1), x2 ~ Bernoulli(0.5)
#   zero part:   logit(pi_i) = g0 + g1 z1_i,            z1 ~ N(0,1), independent of x
#   dispersion:  theta
#   Y_i = 0 with probability pi_i (structural zero), else Y_i ~ NB(mu_i, theta).
# Covariates are redrawn for every data set. The seed makes it reproducible.
#
# For each parameter (beta, gamma, log theta) we report the bias of the
# estimates, the Monte Carlo SD of the estimates vs the mean estimated SE,
# and the empirical coverage of 95% Wald intervals est +/- 1.96 SE.
# Coverage for theta is assessed on the log scale, where its SE is computed.
#
# Run from the package root, with the package installed:
#   Rscript inst/benchmark/simulation.R
# Results are saved to inst/benchmark/results/simulation.rds.

library(nbzinbPate)

set.seed(651)

R_reps <- 100
n      <- 500
beta   <- c("(Intercept)" = 1.0, x1 = 0.5, x2 = -0.3)
gamma  <- c("(Intercept)" = -0.5, z1 = 0.8)
theta  <- 1.5

truth <- c(setNames(beta,  paste0("count_", names(beta))),
           setNames(gamma, paste0("zero_",  names(gamma))),
           "Log(theta)" = log(theta))

out_dir <- file.path("inst", "benchmark", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Simulate one data set from the ZINB model ------------------------------------
simulate_zinb <- function(n) {
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  z1 <- rnorm(n)
  mu <- exp(beta[1] + beta[2] * x1 + beta[3] * x2)     # count mean
  pi <- plogis(gamma[1] + gamma[2] * z1)                # structural-zero prob.
  structural_zero <- rbinom(n, 1, pi)                   # latent D_i
  y <- ifelse(structural_zero == 1, 0,
              rnbinom(n, size = theta, mu = mu))        # NB2(mu, theta)
  data.frame(y, x1, x2, z1)
}


# ---- Fit one data set, recording estimates, SEs and any warnings ---------------
fit_one <- function(d) {
  warns <- character(0)
  fit <- withCallingHandlers(
    zinbreg(y ~ x1 + x2 | z1, data = d),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  est <- c(coef(fit), "Log(theta)" = log(fit$theta))
  se  <- sqrt(diag(fit$vcov))[names(est)]   # vcov includes Log(theta)
  list(est = est, se = se, converged = fit$converged, iter = fit$iter,
       boundary = fit$boundary, warnings = warns,
       prop_zero = mean(d$y == 0))
}


# ---- Run the replications ---------------------------------------------------------------
time <- system.time(
  reps <- lapply(seq_len(R_reps), function(r) fit_one(simulate_zinb(n)))
)[["elapsed"]]

est <- t(sapply(reps, `[[`, "est"))           # R x 6 matrix of estimates
se  <- t(sapply(reps, `[[`, "se"))            # R x 6 matrix of SEs
stopifnot(identical(colnames(est), names(truth)))


# ---- Summaries ---------------------------------------------------------------------------
lower <- est - qnorm(0.975) * se
upper <- est + qnorm(0.975) * se
covered <- sweep(lower, 2, truth, "<=") & sweep(upper, 2, truth, ">=")

summary_table <- data.frame(
  parameter  = names(truth),
  true       = unname(truth),
  mean_est   = colMeans(est),
  bias       = colMeans(est) - truth,
  mc_sd      = apply(est, 2, sd),          # Monte Carlo SD of the estimates
  mean_se    = colMeans(se),               # average estimated SE
  se_ratio   = colMeans(se) / apply(est, 2, sd),
  coverage   = colMeans(covered),          # empirical 95% Wald coverage
  row.names  = NULL)

diagnostics <- list(
  R = R_reps, n = n,
  n_converged   = sum(sapply(reps, `[[`, "converged")),
  n_boundary    = sum(sapply(reps, `[[`, "boundary")),
  n_with_warnings = sum(lengths(lapply(reps, `[[`, "warnings")) > 0),
  warnings      = unique(unlist(lapply(reps, `[[`, "warnings"))),
  em_iterations = summary(sapply(reps, `[[`, "iter")),
  prop_zero     = summary(sapply(reps, `[[`, "prop_zero")),
  # Monte Carlo SE of a coverage estimate near 0.95 with R replications.
  coverage_mcse = sqrt(0.95 * 0.05 / R_reps),
  time_seconds  = time)

results <- list(design = list(n = n, R = R_reps, beta = beta, gamma = gamma,
                              theta = theta, seed = 651),
                summary = summary_table, diagnostics = diagnostics,
                estimates = est, std_errors = se,
                session = list(R = R.version.string,
                               nbzinbPate = as.character(packageVersion("nbzinbPate")),
                               date = Sys.Date()))
saveRDS(results, file.path(out_dir, "simulation.rds"))


# ---- Print -------------------------------------------------------------------------------
print(summary_table, digits = 3)
cat(sprintf("\nConverged: %d/%d   boundary: %d   with warnings: %d\n",
            diagnostics$n_converged, R_reps, diagnostics$n_boundary,
            diagnostics$n_with_warnings))
cat(sprintf("Coverage Monte Carlo SE: %.3f   EM iterations (median): %.0f   time: %.0fs\n",
            diagnostics$coverage_mcse, median(sapply(reps, `[[`, "iter")), time))
cat(sprintf("Proportion of zeros per data set: %.2f to %.2f\n",
            min(sapply(reps, `[[`, "prop_zero")), max(sapply(reps, `[[`, "prop_zero"))))
cat("Saved to", file.path(out_dir, "simulation.rds"), "\n")
