# Methods for nbreg objects

Standard accessor, summary, and prediction methods for fits from
[`nbreg()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'nbreg'
predict(object, newdata, type = c("response", "link"), ...)

# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x:

  An object of class `"nbreg"` (or `"summary.nbreg"` for the print
  method of a summary).

- ...:

  Further arguments (ignored or passed on to
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html)).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (\\y_i - \hat\mu_i\\) or `"pearson"`
  (\\(y_i-\hat\mu_i)/\sqrt{\hat\mu_i + \hat\mu_i^2/\hat\theta}\\). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (\\\hat\mu\\) or `"link"` (\\\hat\eta = \log\hat\mu\\).

- newdata:

  Optional data frame in which to look for covariates (and the offset
  variable) for prediction. If omitted, fitted values are used.

- digits:

  Number of significant digits to print.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): estimated coefficients.
[`vcov()`](https://rdrr.io/r/stats/vcov.html): their covariance matrix.
[`logLik()`](https://rdrr.io/r/stats/logLik.html): a `"logLik"` object
with `df` = number of coefficients + 1 (for \\\theta\\), so
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html): \\\hat\mu\\.
[`residuals()`](https://rdrr.io/r/stats/residuals.html): see `type`.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): number of observations.
[`predict()`](https://rdrr.io/r/stats/predict.html): numeric vector.
[`summary()`](https://rdrr.io/r/base/summary.html): an object of class
`"summary.nbreg"` with a coefficient table.

## Examples

``` r
set.seed(2)
d <- data.frame(x = rnorm(200), g = factor(sample(c("a", "b"), 200, TRUE)))
d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.4 * d$x + 0.3 * (d$g == "b")))
fit <- nbreg(y ~ x + g, data = d)
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x + g, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1946 -0.7139 -0.2793  0.4788  3.6406 
#> 
#> Coefficients (NB2 with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.58657    0.09615   6.101 1.06e-09 ***
#> x            0.28609    0.06268   4.564 5.01e-06 ***
#> gb           0.28533    0.13330   2.140   0.0323 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.474  (SE 0.5774 )
#> Log-likelihood: -378.056 on 4 Df;  AIC: 764.112   BIC: 777.306 
#> Number of observations: 200   Outer iterations (IRLS + Newton): 4 
AIC(fit); BIC(fit)
#> [1] 764.1122
#> [1] 777.3055
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -1.0288890 -1.0359851 -0.3374709 -0.2136776 -1.0963069  1.7307515 
predict(fit, newdata = data.frame(x = c(-1, 0, 1), g = "b"))
#>        1        2        3 
#> 1.796449 2.391451 3.183522 
```
