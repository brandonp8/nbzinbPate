# Enforces the project rule (Sec 3.3, Table 2): the package's fitting code
# must not call optim, glm, lm, MASS/pscl fitting functions, or any other
# optimizer.
#
# How: instead of searching the source text (which would also match
# comments, and is unavailable under R CMD check, where tests run against
# the installed package), we inspect the parsed code of every function in
# the package namespace. all.names(body(f)) returns every name used in the
# function body -- including the names of called functions and both parts
# of pkg::fun -- but never anything inside a comment.

prohibited_functions <- c(
  "optim", "optimize", "optimise", "nlm", "nlminb", "uniroot", "optimHess",
  "glm", "glm.fit", "lm", "lm.fit", "lm.wfit",
  "glm.nb", "theta.ml", "theta.md", "theta.mm", "zeroinfl", "hurdle",
  "optimx", "maxLik", "mle", "mle2", "glmmTMB", "gamlss", "vglm"
)
prohibited_packages <- c("MASS", "pscl", "countreg", "glmmTMB", "gamlss",
                         "VGAM", "optimx", "maxLik", "stats4", "bbmle",
                         "numDeriv")   # numDeriv: we use the analytic Hessian

test_that("no package function calls a prohibited optimizer or fitter", {
  ns  <- asNamespace("nbzinbPate")
  fns <- Filter(function(nm) is.function(get(nm, envir = ns)), ls(ns, all.names = TRUE))
  expect_gt(length(fns), 10)                  # sanity: we really found them

  offenders <- character(0)
  for (nm in fns) {
    used <- unique(c(all.names(body(get(nm, envir = ns))),
                     all.names(as.call(c(as.name("list"),
                                         formals(get(nm, envir = ns)))))))
    bad <- intersect(used, c(prohibited_functions, prohibited_packages))
    if (length(bad) > 0) {
      offenders <- c(offenders, paste0(nm, ": ", paste(bad, collapse = ", ")))
    }
  }
  expect_identical(offenders, character(0))
})

test_that("the check itself would catch a prohibited call", {
  # Guard against the test above passing vacuously.
  f <- function(x) stats::optim(1, function(p) (p - x)^2)$par
  expect_true("optim" %in% all.names(body(f)))
  g <- function(d) MASS::glm.nb(y ~ x, data = d)
  expect_true(all(c("MASS", "glm.nb") %in% all.names(body(g))))
  # ...and that comments are invisible to it.
  h <- function() {
    # optim() is not used here
    1
  }
  expect_false("optim" %in% all.names(body(h)))
})
