# Control settings for [`zinbreg()`](https://brandonp8.github.io/nbzinbPate/reference/zinbreg.md)

Control settings for
[`zinbreg()`](https://brandonp8.github.io/nbzinbPate/reference/zinbreg.md)

## Usage

``` r
zinbreg_control(
  tol = 1e-08,
  par_tol = 1e-06,
  maxit = 1000,
  mstep_tol = 1e-10,
  mstep_maxit = 100,
  theta_warn = 1e+05,
  trace = FALSE
)
```

## Arguments

- tol:

  EM tolerance (Sec 3.2). EM has converged when the observed-data
  log-likelihood changes by less than `tol` **and** no parameter
  (\\\beta\\, \\\gamma\\, \\\log\theta\\) changes by more than `par_tol`
  between EM iterations.

- par_tol:

  Tolerance on the largest parameter change; see `tol`.

- maxit:

  Maximum number of EM iterations.

- mstep_tol, mstep_maxit:

  Tolerance and iteration limit for the M-step fits. Each M-step is a
  *full* fit, warm-started at the current values: the weighted NB fit
  uses
  [`nbreg_control()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg_control.md)
  with `tol = irls_tol = theta_tol = mstep_tol` and these iteration
  limits, and the logistic IRLS uses the same tolerance and limit.

- theta_warn:

  Warn if \\\hat\theta\\ exceeds this value.

- trace:

  Logical; if `TRUE`, print the log-likelihood at every EM iteration.

## Value

A list of class `"zinbreg_control"`.

## Examples

``` r
zinbreg_control(tol = 1e-10, maxit = 2000)
#> $tol
#> [1] 1e-10
#> 
#> $par_tol
#> [1] 1e-06
#> 
#> $maxit
#> [1] 2000
#> 
#> $mstep_tol
#> [1] 1e-10
#> 
#> $mstep_maxit
#> [1] 100
#> 
#> $theta_warn
#> [1] 1e+05
#> 
#> $trace
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "zinbreg_control"
```
