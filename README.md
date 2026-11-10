
<!-- README.md is generated from README.Rmd. Please edit README.Rmd, then run
     devtools::build_readme() to regenerate README.md. -->

# nbzinbPate

<!-- badges: start -->

[![R-CMD-check](https://github.com/brandonp8/nbzinbPate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/brandonp8/nbzinbPate/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**nbzinbPate** fits

- **negative binomial (NB2) regression** with `nbreg()`, comparable to
  `MASS::glm.nb()`, and
- **zero-inflated negative binomial (ZINB) regression** with
  `zinbreg()`, comparable to `pscl::zeroinfl(..., dist = "negbin")`,

with all of the estimation implemented from scratch. No general-purpose
optimizer (`optim()`, `nlm()`, …) and no GLM fitter (`glm()`, `lm()`, …)
is used:

- **NB**: iteratively reweighted least squares (IRLS) for the
  coefficients with $\theta$ fixed, alternated with Newton-Raphson for
  $\log\theta$ with the coefficients fixed.
- **ZINB**: an EM algorithm. The E-step computes the posterior
  probability that each zero is a structural zero. The M-step fits a
  fractional-response logistic regression for the zero model and a
  weighted NB regression for the count model, both by IRLS.
- **Standard errors**: $(X^\top \hat W X)^{-1}$ and
  $\{-H(\hat\theta)\}^{-1/2}$ for NB. For ZINB, the inverse observed
  information, computed from an **analytic** Hessian of the
  observed-data log-likelihood with respect to
  $(\beta, \gamma, \log\theta)$.

This package was written for the STA 651 final project (Fall 2026).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("brandonp8/nbzinbPate")
```

The package only depends on base R (`stats`). The reference packages
`MASS`, `pscl`, `AER` and `numDeriv` are needed only to run the tests.

## Example

The `bioChemists` data (from **pscl**) record the number of articles
published by 915 biochemistry PhD students.

``` r
library(nbzinbPate)
data("bioChemists", package = "pscl")

# Negative binomial regression
m_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(m_nb)
#> 
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.2677 -0.8754 -0.2613  0.4983  6.6573 
#> 
#> Coefficients (NB2 with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .  
#> femWomen    -0.216418   0.072636  -2.979 0.002887 ** 
#> marMarried   0.150489   0.082097   1.833 0.066791 .  
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326    
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.264  (SE 0.2712 )
#> Log-likelihood: -1560.96 on 7 Df;  AIC: 3135.92   BIC: 3169.65 
#> Number of observations: 915   Outer iterations (IRLS + Newton): 5
```

``` r
# Zero-inflated NB: count model before "|", zero-inflation model after it
m_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | kid5 + ment,
                  data = bioChemists)
summary(m_zinb)
#> 
#> Call:
#> zinbreg(formula = art ~ fem + mar + kid5 + phd + ment | kid5 + 
#>     ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3041 -0.7666 -0.2675  0.4677  6.3591 
#> 
#> Count model coefficients (NB2 with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.401517   0.141742   2.833  0.00462 ** 
#> femWomen    -0.211890   0.071913  -2.946  0.00321 ** 
#> marMarried   0.139455   0.081129   1.719  0.08563 .  
#> kid5        -0.161571   0.054769  -2.950  0.00318 ** 
#> phd          0.002164   0.035587   0.061  0.95152    
#> ment         0.024335   0.003520   6.912 4.77e-12 ***
#> Log(theta)   1.003303   0.142698   7.031 2.05e-12 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.8756     0.4028  -2.174   0.0297 *
#> kid5          0.1483     0.3634   0.408   0.6831  
#> ment         -0.6089     0.2445  -2.491   0.0127 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta: 2.727  (SE 0.3892 )
#> Log-likelihood: -1553.19 on 10 Df;  AIC: 3126.39   BIC: 3174.58 
#> Number of observations: 915   EM iterations: 133
```

The usual methods work: `coef()`, `vcov()`, `logLik()`, `AIC()`,
`BIC()`, `fitted()`, `residuals()`, `nobs()`, and `predict()`, including
`predict(m_zinb, newdata, type = "response" / "count" / "zero")`.

``` r
AIC(m_nb, m_zinb)
#>        df      AIC
#> m_nb    7 3135.917
#> m_zinb 10 3126.387
predict(m_zinb, newdata = data.frame(fem = "Women", mar = "Married",
                                     kid5 = 0, phd = 3, ment = c(0, 10, 40)),
        type = "zero")
#>            1            2            3 
#> 2.940806e-01 9.440864e-04 1.102945e-11
```

## Agreement with the reference implementations

On `bioChemists` and `AER::NMES1988`, the estimates match
`MASS::glm.nb()` and `pscl::zeroinfl()` far more closely than the
project tolerances require.

- **NB vs `glm.nb()`**: the coefficients, their standard errors,
  $\theta$ and the log-likelihood agree to about 1e-10. SE($\theta$)
  agrees to about 1e-6, because `glm.nb()` evaluates it at the
  second-to-last Newton step for $\theta$.
- **ZINB vs `zeroinfl()`**: the coefficients, standard errors and
  $\log\theta$ agree to about 1e-5 or better, and the log-likelihood to
  about 1e-9.

The package tests (`tests/testthat/`) check this agreement
automatically.

## R CMD check

`devtools::check()` returns 0 errors, 0 warnings, and 0 notes.
