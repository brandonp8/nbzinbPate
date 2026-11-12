# Zero-inflated negative binomial (ZINB) regression

Fits the ZINB model of Lambert (1992) with an NB2 count component:
\$\$P(Y_i = 0) = \pi_i + (1-\pi_i) f\_{NB}(0;\mu_i,\theta), \quad P(Y_i
= y) = (1-\pi_i) f\_{NB}(y;\mu_i,\theta),\\ y \ge 1,\$\$ with
\\\log\mu_i = x_i^\top\beta + o_i\\ (count model) and
\\\mathrm{logit}\\\pi_i = z_i^\top\gamma\\ (zero-inflation model).
Estimation is by the EM algorithm: the E-step computes the posterior
probability \\\tau_i\\ that each zero is a structural zero, and the
M-step fits a fractional-response logistic regression for \\\gamma\\ and
a weighted NB regression (weights \\1-\tau_i\\) for \\(\beta,\theta)\\,
both by IRLS. No general-purpose optimizer is used. The parameterization
is the same as `pscl::zeroinfl(..., dist = "negbin")`.

## Usage

``` r
zinbreg(formula, data, offset, control = zinbreg_control())
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before `|`, zero-model covariates after it. Without `|`, the same
  covariates are used in both parts (as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)).
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms are allowed in
  the count part only.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional count-model offset (e.g. `log(exposure)`), evaluated in
  `data`.

- control:

  A list from
  [`zinbreg_control()`](https://brandonp8.github.io/nbzinbPate/reference/zinbreg_control.md).

## Value

An object of class `"zinbreg"`: a list with components

- coefficients:

  List with `count` (\\\beta\\) and `zero` (\\\gamma\\).

- vcov:

  Estimated covariance matrix of \\(\hat\beta, \hat\gamma,
  \log\hat\theta)\\: the inverse of the observed information, from the
  analytic Hessian of the observed-data log-likelihood. Row/column names
  are `count_*`, `zero_*`, and `Log(theta)`, as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html).

- theta, SE.theta, SE.logtheta:

  \\\hat\theta\\, its standard error (delta method), and the standard
  error of \\\log\hat\theta\\.

- score:

  Analytic score vector at the estimates (should be ~0).

- loglik, loglik_path:

  Maximized observed-data log-likelihood, and its value at every EM
  iteration (non-decreasing).

- df:

  Number of estimated parameters (\\p + q + 1\\).

- fitted.values:

  \\E(Y_i) = (1-\hat\pi_i)\hat\mu_i\\.

- mu, pi:

  Count-model means \\\hat\mu_i\\ and zero probabilities \\\hat\pi_i\\.

- iter, converged:

  Number of EM iterations and convergence flag.

- boundary:

  `TRUE` if the zero-inflation component vanished (all \\\hat\pi_i \<
  10^{-8}\\); a warning is then issued.

- y, offset, nobs:

  Response, count offset, number of observations.

- call, terms, xlevels, contrasts, model:

  Information used by
  [`predict()`](https://rdrr.io/r/stats/predict.html) and other methods.

## Details

EM does not produce standard errors, so after convergence the observed
information \\-H\\ of the observed-data log-likelihood with respect to
\\(\beta, \gamma, \log\theta)\\ is computed from analytic second
derivatives and inverted.

## References

Lambert, D. (1992). Zero-inflated Poisson regression, with an
application to defects in manufacturing. *Technometrics*, 34(1), 1-14.

## See also

[`zinbreg_control()`](https://brandonp8.github.io/nbzinbPate/reference/zinbreg_control.md),
[`nbreg()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg.md)

## Examples

``` r
# Simulated ZINB data: one covariate in each part
set.seed(1)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
structural_zero <- rbinom(n, 1, plogis(-1 + d$z))
d$y <- ifelse(structural_zero == 1, 0,
              rnbinom(n, size = 1.5, mu = exp(1 + 0.5 * d$x)))
fit <- zinbreg(y ~ x | z, data = d)
fit$coefficients
#> $count
#> (Intercept)           x 
#>   0.9979526   0.4374179 
#> 
#> $zero
#> (Intercept)           z 
#>   -1.653210    1.253677 
#> 
fit$theta
#> [1] 1.157176
```
