# nbzinbPate: Negative Binomial and Zero-Inflated Negative Binomial Regression

Fits negative binomial (NB2) regression by iteratively reweighted least
squares (IRLS) for the regression coefficients, alternated with
Newton-Raphson for the log dispersion parameter, and zero-inflated
negative binomial (ZINB) regression by an EM algorithm whose M-steps are
weighted IRLS fits. Standard errors for the ZINB model come from an
analytic Hessian of the observed-data log-likelihood. No general-purpose
optimizer is used.

## See also

Useful links:

- <https://github.com/brandonp8/nbzinbPate>

- Report bugs at <https://github.com/brandonp8/nbzinbPate/issues>

## Author

**Maintainer**: Brandon Pate <patebd2@miamioh.edu>

Authors:

- Brandon Pate <patebd2@miamioh.edu>
