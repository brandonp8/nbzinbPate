# Negative binomial (NB2) regression

Fits the NB2 regression model \$\$Y_i \sim \mathrm{NB}(\mu_i, \theta),
\quad \log\mu_i = x_i^\top\beta + o_i, \quad \mathrm{Var}(Y_i) = \mu_i +
\mu_i^2/\theta,\$\$ by maximum likelihood. The coefficients \\\beta\\
are updated by iteratively reweighted least squares (IRLS) with
\\\theta\\ held fixed, and \\\log\theta\\ by Newton-Raphson with
\\\beta\\ held fixed; the two steps alternate until convergence. No
general-purpose optimizer is used. The parameterization is the same as
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

## Usage

``` r
nbreg(formula, data, offset, control = nbreg_control())
```

## Arguments

- formula:

  A model formula such as `y ~ x1 + x2`. The response must be a vector
  of non-negative integer counts.
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms are allowed.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset vector (e.g. `log(exposure)`), evaluated in `data`; it
  is added to any [`offset()`](https://rdrr.io/r/stats/offset.html)
  terms in the formula.

- control:

  A list from
  [`nbreg_control()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg_control.md).

## Value

An object of class `"nbreg"`: a list with components

- coefficients:

  Estimated \\\beta\\.

- vcov:

  Estimated covariance matrix of \\\hat\beta\\.

- theta, SE.theta:

  \\\hat\theta\\ and its standard error.

- loglik:

  Maximized log-likelihood.

- df:

  Number of estimated parameters (\\p + 1\\).

- fitted.values, linear.predictors:

  \\\hat\mu_i\\ and \\\hat\eta_i\\.

- iter, converged:

  Number of outer iterations and convergence flag.

- y, offset, nobs:

  The response, offset, and number of observations.

- call, terms, xlevels, contrasts, model:

  Information used by
  [`predict()`](https://rdrr.io/r/stats/predict.html) and other methods.

## Details

Standard errors for \\\beta\\ are \\(X^\top\hat W X)^{-1/2}\\
(diagonal), treating \\\theta\\ as known, as `glm.nb()` does; the
standard error of \\\hat\theta\\ is \\\\-H(\hat\theta)\\^{-1/2}\\, where
\\H\\ is the second derivative of the log-likelihood in \\\theta\\.

## See also

[`nbreg_control()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg_control.md),
[`zinbreg()`](https://brandonp8.github.io/nbzinbPate/reference/zinbreg.md)

## Examples

``` r
# Simulated data with an exposure offset
set.seed(1)
n <- 300
d <- data.frame(x = rnorm(n), expo = runif(n, 1, 3))
d$y <- rnbinom(n, size = 2, mu = d$expo * exp(0.5 + 0.8 * d$x))
fit <- nbreg(y ~ x, data = d, offset = log(expo))
fit$coefficients
#> (Intercept)           x 
#>   0.3999129   0.8652156 
fit$theta
#> [1] 2.119073
```
