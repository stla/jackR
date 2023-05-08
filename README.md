The ‘jack’ package: Jack polynomials
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/jackR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/jackR/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-valgrind](https://github.com/stla/jackR/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/jackR/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->

``` r
library(jack)
library(microbenchmark)
```

Schur polynomials have applications in combinatorics and zonal
polynomials have applications in multivariate statistics. They are
particular cases of [Jack
polynomials](https://en.wikipedia.org/wiki/Jack_function). This package
allows to evaluate these polynomials. It can also compute their symbolic
form.

The functions `JackPol`, `ZonalPol`, `ZonalQPol` and `SchurPol`
respectively return the Jack polynomial, the zonal polynomial, the
quaternionic zonal polynomial, and the Schur polynomial.

Each of these polynomials corresponds is given by a positive integer,
the number of variables, and an integer partition, the `lambda`
argument; the Jack polynomial has one more parameter, the `alpha`
argument, a positive number.

To get an exact symbolic polynomial with `JackPol`, you have to supply a
`bigq` rational number for the parameter `alpha`:

``` r
jpol <- JackPol(2, lambda = c(3, 1), alpha = gmp::as.bigq("2/5"))
jpol
## 98/25*x^(3, 1) + 98/25*x^(1, 3) + 28/5*x^(2, 2)
```

This is a `qspray` object, from the
[**qspray**](https://github.com/stla/qspray) package. Here is how you
can evaluate this polynomial:

``` r
qspray::evalQspray(jpol, c("2", "3/2"))
## Big Rational ('bigq') :
## [1] 1239/10
```

By default, `ZonalPol`, `ZonalQPol` and `SchurPol` return exact symbolic
polynomials.

``` r
zpol <- ZonalPol(2, lambda = c(3, 1))
zpol
## 24/7*x^(3, 1) + 24/7*x^(1, 3) + 16/7*x^(2, 2)
```

It is also possible to convert a `qspray` polynomial to a function whose
evaluation is performed by the **Ryacas** package:

``` r
zyacas <- as.function(zpol)
```

You can provide the values of the variables of this function as numbers
or character strings:

``` r
zyacas(2, "3/2")
## [1] "594/7"
```

You can even pass a variable name to this function:

``` r
zyacas("x", "x")
## [1] "(64*x^4)/7"
```

If you want to substitute a variable with a complex number, use a
character string which represents this number, with `I` denoting the
imaginary unit:

``` r
zyacas("2 + 2*I", "2/3")
## [1] "Complex((-2176)/63,2944/63)"
```

## Jack polynomials with Julia

As of version 2.0.0, the Jack polynomials can be calculated with Julia.
The speed is amazing:

``` r
julia <- Jack_julia()
## Starting Julia ...
x <- c(1/2, 2/3, 1, 2/3, -1, -2, 1)
lambda <- c(5, 3, 2, 2, 1)
alpha <- 3
print(
  microbenchmark(
        R = Jack(x, lambda, alpha),
    Julia = julia$Jack(x, lambda, alpha),
    times = 6L,
    unit  = "seconds"
  ),
  signif = 6L
)
## Unit: seconds
##   expr       min        lq     mean    median       uq     max neval
##      R 6.8420000 6.8521900 6.897750 6.8773600 6.952530 6.98506     6
##  Julia 0.0034501 0.0919164 0.217018 0.0945618 0.112051 0.90557     6
```

`Jack_julia()` returns a list of functions. `ZonalPol`, `ZonalQPol` and
`SchurPol` always return an exact expression of the polynomial,
i.e. with rational coefficients (integers for `SchurPol`). If you want
an exact expression with `JackPol`, you have to give a rational number
for the argument `alpha`, as a character string:

``` r
JP <- julia$JackPol(m = 2, lambda = c(3, 1), alpha = "2/5")
JP
## 98/25*x^(1, 3) + 28/5*x^(2, 2) + 98/25*x^(3, 1)
```

Again, Julia is faster:

``` r
n <- 5
lambda <- c(4, 3, 3)
alpha <- "2/3"
alphaq <- gmp::as.bigq(alpha)
print(
  microbenchmark(
        R = JackPol(n, lambda, alphaq),
    Julia = julia$JackPol(n, lambda, alpha),
    times = 6L
  ),
signif = 2L)
## Unit: milliseconds
##   expr  min   lq mean median   uq  max neval
##      R 1000 1000 1000   1000 1000 1100     6
##  Julia  450  510  540    520  520  720     6
```

## ‘Rcpp’ implementation of the polynomials

As of version 5.0.0, a ‘Rcpp’ implementation of the polynomials is
provided by the package. It is faster than Julia (though I didn’t
compare in pure Julia - the Julia execution time is slowed down by the
‘JuliaConnectoR’ package):

``` r
n <- 5
lambda <- c(4, 3, 3, 2)
print(
  microbenchmark(
     Rcpp = SchurPolCPP(n, lambda),
    Julia = julia$SchurPol(n, lambda),
    times = 6L
  ), 
signif = 2L)
## Unit: milliseconds
##   expr   min    lq  mean median    uq    max neval
##   Rcpp   6.5   6.5   6.6    6.7   6.7    6.8     6
##  Julia 530.0 530.0 650.0  540.0 550.0 1200.0     6
```

``` r
n <- 5
lambda <- c(4, 3, 3, 2)
alpha <- "2/3"
print(
  microbenchmark(
     Rcpp = JackPolCPP(n, lambda, alpha),
    Julia = julia$JackPol(n, lambda, alpha),
    times = 6L
  ), 
signif = 2L)
## Unit: milliseconds
##   expr min  lq mean median  uq max neval
##   Rcpp  23  23   24     24  24  24     6
##  Julia 370 380  410    380 460 500     6
```

As of version 5.1.0, there’s also a ‘Rcpp’ implementation of the
evaluation of the polynomials.

``` r
x <- c("1/2", "2/3", "1", "2/3", "1", "5/4")
lambda <- c(5, 3, 2, 2, 1)
alpha <- "3"
print(
  microbenchmark(
        R = Jack(gmp::as.bigq(x), lambda, gmp::as.bigq(alpha)),
     Rcpp = JackCPP(x, lambda, alpha),    
    times = 10L,
    unit  = "seconds"
  ),
  signif = 6L
)
## Unit: seconds
##  expr       min        lq      mean    median       uq       max neval
##     R 15.302300 15.649700 16.026900 16.276300 16.37570 16.442200    10
##  Rcpp  0.146324  0.149152  0.152302  0.153073  0.15542  0.156996    10
```

``` r
JuliaConnectoR::stopJulia()
```
