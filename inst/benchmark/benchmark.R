# Benchmark: nbzinbPate vs the reference implementations (project Sec 5).
#
# For each data set (pscl::bioChemists, AER::NMES1988) and each model
# (NB, ZINB) this script fits our function and the reference function,
# puts the two sets of estimates side by side, computes their maximum
# absolute difference, compares that with the Table 4 tolerances, and
# times both implementations.
#
# The reference functions MASS::glm.nb() and pscl::zeroinfl() are used
# here only for comparison (allowed in inst/benchmark/), with the tight
# tolerances the project recommends.
#
# Run from the package root, with the package installed:
#   Rscript inst/benchmark/benchmark.R
# Results are saved to inst/benchmark/results/benchmark.rds, which the
# report reads (so the report does not have to refit everything).

library(nbzinbPate)

out_dir <- file.path("inst", "benchmark", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Data sets and model formulas ------------------------------------------------

data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

datasets <- list(
  bioChemists = list(
    data = bioChemists,
    nb   = art ~ fem + mar + kid5 + phd + ment,
    zinb = art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ),
  NMES1988 = list(
    data = NMES1988,
    # Count covariates required by Sec 5.2; zero covariates = the minimum set.
    nb   = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance
  )
)

# Table 4: required agreement (maximum absolute difference).
tolerance <- list(
  NB   = c(Coefficients = 1e-5, `Standard errors` = 1e-5,
           theta = 1e-4, `Log-likelihood` = 1e-4),
  ZINB = c(Coefficients = 1e-3, `Standard errors` = 1e-3,
           `log(theta)` = 1e-3, `Log-likelihood` = 1e-3)
)


# ---- Helpers ---------------------------------------------------------------------------

# Median elapsed time (seconds) over `reps` runs of expr.
time_it <- function(expr_fun, reps) {
  median(replicate(reps, system.time(expr_fun())[["elapsed"]]))
}

# Side-by-side table: parameter, ours, reference, |difference|.
side_by_side <- function(ours, ref, quantity) {
  stopifnot(identical(names(ours), names(ref)))
  data.frame(quantity  = quantity,
             parameter = names(ours),
             ours      = unname(ours),
             reference = unname(ref),
             abs_diff  = abs(unname(ours) - unname(ref)),
             row.names = NULL)
}

# Agreement summary: max |difference| per quantity vs the Table 4 tolerance.
agreement <- function(tab, tol) {
  qs <- names(tol)
  max_diff <- vapply(qs, function(q) max(tab$abs_diff[tab$quantity == q]),
                     numeric(1))
  data.frame(quantity = qs, max_abs_diff = max_diff, tolerance = unname(tol),
             pass = max_diff <= tol, row.names = NULL)
}


# ---- NB: nbreg() vs MASS::glm.nb() --------------------------------------------------

bench_nb <- function(d, f, reps = 5) {
  ours <- nbreg(f, data = d)
  ref  <- MASS::glm.nb(f, data = d,
                       control = glm.control(epsilon = 1e-10, maxit = 100))

  tab <- rbind(
    side_by_side(coef(ours), coef(ref), "Coefficients"),
    side_by_side(sqrt(diag(vcov(ours))), sqrt(diag(vcov(ref))),
                 "Standard errors"),
    side_by_side(c(theta = ours$theta), c(theta = ref$theta), "theta"),
    side_by_side(c(`SE(theta)` = ours$SE.theta),
                 c(`SE(theta)` = ref$SE.theta), "SE(theta)"),
    side_by_side(c(logLik = ours$loglik),
                 c(logLik = as.numeric(logLik(ref))), "Log-likelihood"))

  time <- c(
    ours      = time_it(function() nbreg(f, data = d), reps),
    reference = time_it(function()
      MASS::glm.nb(f, data = d,
                   control = glm.control(epsilon = 1e-10, maxit = 100)), reps))

  list(table = tab, agreement = agreement(tab, tolerance$NB), time = time,
       iterations = ours$iter)
}


# ---- ZINB: zinbreg() vs pscl::zeroinfl() ----------------------------------------------

bench_zinb <- function(d, f, reps = 3) {
  ours <- zinbreg(f, data = d)
  ref  <- pscl::zeroinfl(f, data = d, dist = "negbin",
                         control = pscl::zeroinfl.control(reltol = 1e-12))

  # Both packages name coefficients count_* and zero_* in coef() and vcov().
  tab <- rbind(
    side_by_side(coef(ours), coef(ref), "Coefficients"),
    side_by_side(sqrt(diag(vcov(ours))), sqrt(diag(vcov(ref))),
                 "Standard errors"),
    side_by_side(c(`log(theta)` = log(ours$theta)),
                 c(`log(theta)` = log(ref$theta)), "log(theta)"),
    side_by_side(c(`SE(log theta)` = ours$SE.logtheta),
                 c(`SE(log theta)` = ref$SE.logtheta), "SE(log theta)"),
    side_by_side(c(logLik = ours$loglik),
                 c(logLik = as.numeric(logLik(ref))), "Log-likelihood"))

  time <- c(
    ours      = time_it(function() zinbreg(f, data = d), reps),
    reference = time_it(function()
      pscl::zeroinfl(f, data = d, dist = "negbin",
                     control = pscl::zeroinfl.control(reltol = 1e-12)), reps))

  list(table = tab, agreement = agreement(tab, tolerance$ZINB), time = time,
       iterations = ours$iter,
       # Sec 5 tip: report it if our log-likelihood is higher than zeroinfl's.
       loglik_higher = ours$loglik > as.numeric(logLik(ref)))
}


# ---- Run everything --------------------------------------------------------------------

results <- lapply(datasets, function(ds) {
  list(NB   = bench_nb(ds$data, ds$nb),
       ZINB = bench_zinb(ds$data, ds$zinb))
})

results$session <- list(R = R.version.string,
                        nbzinbPate = as.character(packageVersion("nbzinbPate")),
                        MASS = as.character(packageVersion("MASS")),
                        pscl = as.character(packageVersion("pscl")),
                        date = Sys.Date())

saveRDS(results, file.path(out_dir, "benchmark.rds"))


# ---- Print a summary -------------------------------------------------------------------

for (nm in names(datasets)) {
  for (m in c("NB", "ZINB")) {
    r <- results[[nm]][[m]]
    cat("\n==", nm, "-", m, "==\n")
    print(r$agreement, digits = 3)
    cat(sprintf("Time (median): ours %.3fs, reference %.3fs; iterations: %d\n",
                r$time[["ours"]], r$time[["reference"]], r$iterations))
    if (!is.null(r$loglik_higher) && r$loglik_higher) {
      cat("Note: our log-likelihood is higher than zeroinfl's.\n")
    }
  }
}
cat("\nSaved to", file.path(out_dir, "benchmark.rds"), "\n")
