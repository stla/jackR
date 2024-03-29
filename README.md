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

As of version 2.0.0, it was possible to calculate the Jack polynomials
with Julia. This feature has been removed in version 5.3.0. Use the
Julia package **JackPolynomials.jl** instead.

## ‘Rcpp’ implementation of the polynomials

As of version 5.0.0, a ‘Rcpp’ implementation of the polynomials is
provided by the package.

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
    times = 6L,
    unit  = "seconds"
  ),
  signif = 2L
)
## Unit: seconds
##  expr    min    lq  mean median    uq   max neval
##     R 110.00 130.0 140.0  130.0 160.0 160.0     6
##  Rcpp   0.98   1.2   1.3    1.2   1.5   1.6     6
```
