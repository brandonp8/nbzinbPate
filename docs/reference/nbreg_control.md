# Control settings for [`nbreg()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg.md)

Control settings for
[`nbreg()`](https://brandonp8.github.io/nbzinbPate/reference/nbreg.md)

## Usage

``` r
nbreg_control(
  tol = 1e-08,
  maxit = 100,
  irls_tol = 1e-10,
  irls_maxit = 100,
  theta_tol = 1e-10,
  theta_maxit = 100,
  theta_warn = 1e+05,
  trace = FALSE
)
```

## Arguments

- tol:

  Outer-loop tolerance (Sec 3.1(c)). The fit has converged when the
  log-likelihood changes by less than `tol` **and** no coefficient or
  \\\log\theta\\ changes by more than `tol` between outer iterations.

- maxit:

  Maximum number of outer iterations (each one = an IRLS fit for
  \\\beta\\ followed by a Newton fit for \\\theta\\).

- irls_tol, irls_maxit:

  Tolerance (largest change in \\\beta\\) and iteration limit for the
  inner IRLS of Sec 3.1(a).

- theta_tol, theta_maxit:

  Tolerance (change in \\\log\theta\\) and iteration limit for the inner
  Newton-Raphson of Sec 3.1(b).

- theta_warn:

  A warning is issued if \\\hat\theta\\ exceeds this value: the data
  then show essentially no overdispersion (the Poisson limit), and a
  Poisson model may be adequate.

- trace:

  Logical; if `TRUE`, print the log-likelihood and \\\theta\\ at every
  outer iteration.

## Value

A list of class `"nbreg_control"` with the components above.

## Examples

``` r
nbreg_control(tol = 1e-10, trace = TRUE)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 100
#> 
#> $irls_tol
#> [1] 1e-10
#> 
#> $irls_maxit
#> [1] 100
#> 
#> $theta_tol
#> [1] 1e-10
#> 
#> $theta_maxit
#> [1] 100
#> 
#> $theta_warn
#> [1] 1e+05
#> 
#> $trace
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "nbreg_control"
```
