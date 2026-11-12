# Methods for zinbreg objects

Standard accessor, summary, and prediction methods for fits from
[`zinbreg()`](https://brandonp8.github.io/nbzinbPate/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'zinbreg'
predict(object, newdata, type = c("response", "count", "zero"), ...)

# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x:

  An object of class `"zinbreg"` (or `"summary.zinbreg"`).

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (count and
  zero coefficients, named `count_*` and `zero_*` as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)),
  `"count"`, or `"zero"`. \\\log\theta\\ is not included (as in
  `zeroinfl()`); its standard error is in `object$SE.logtheta`.

- ...:

  Further arguments (ignored or passed on to
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html)).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (\\y_i - \hat E(Y_i)\\) or `"pearson"` (response residual
  divided by \\\sqrt{\widehat{\mathrm{Var}}(Y_i)}\\). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (\\(1-\hat\pi)\hat\mu\\, the mean of \\Y\\), `"count"` (\\\hat\mu\\,
  the mean of the NB component), or `"zero"` (\\\hat\pi\\, the
  probability of a structural zero).

- newdata:

  Optional data frame for prediction.

- digits:

  Number of significant digits to print.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html): estimates and covariance
matrix. [`logLik()`](https://rdrr.io/r/stats/logLik.html): a `"logLik"`
object with `df` = p + q + 1, so
[`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
work. [`fitted()`](https://rdrr.io/r/stats/fitted.values.html):
\\(1-\hat\pi_i)\hat\mu_i\\.
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`predict()`](https://rdrr.io/r/stats/predict.html): see `type`.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): number of observations.
[`summary()`](https://rdrr.io/r/base/summary.html): object of class
`"summary.zinbreg"` with count and zero coefficient tables.

## Examples

``` r
set.seed(3)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- ifelse(rbinom(n, 1, plogis(-1 + d$z)) == 1, 0,
              rnbinom(n, size = 1.5, mu = exp(1 + 0.5 * d$x)))
fit <- zinbreg(y ~ x | z, data = d)
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1146 -0.6764 -0.3953  0.3179  6.0266 
#> 
#> Count model coefficients (NB2 with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.10088    0.06628  16.608  < 2e-16 ***
#> x            0.51751    0.05361   9.652  < 2e-16 ***
#> Log(theta)   0.65011    0.18921   3.436 0.000591 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -0.9439     0.2215  -4.260 2.04e-05 ***
#> z             1.2761     0.2235   5.709 1.14e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 1.916  (SE 0.3625 )
#> Log-likelihood: -927.701 on 5 Df;  AIC: 1865.4   BIC: 1886.47 
#> Number of observations: 500   EM iterations: 37 
AIC(fit)
#> [1] 1865.402
predict(fit, newdata = data.frame(x = 0, z = c(-1, 0, 1)), type = "zero")
#>         1         2         3 
#> 0.0979755 0.2801227 0.5822989 
```
