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
##   expr       min        lq     mean    median       uq      max neval cld
##      R 6.4281400 6.4819200 6.513880 6.5269900 6.553240 6.566030     6   b
##  Julia 0.0465011 0.0482123 0.213404 0.0702411 0.104002 0.941224     6  a
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
microbenchmark(
      R = JackPol(n, lambda, alphaq),
  Julia = julia$JackPol(n, lambda, alpha),
  times = 6L
)
## Unit: milliseconds
##   expr      min       lq     mean   median       uq       max neval cld
##      R 950.0712 965.2956 990.0872 978.1418 982.1807 1086.6921     6   b
##  Julia 460.0946 464.8219 555.1198 545.3837 604.2341  710.8007     6  a
```

## ‘Rcpp’ implementation of the Schur polynomials

As of version 4.1.0, a ‘Rcpp’ implementation of the Schur polynomials is
provided by the package. It is faster than Julia (though I didn’t
compare in pure Julia - the Julia execution time is slowed down by the
‘JuliaConnectoR’ package):

``` r
n <- 5
lambda <- c(4, 3, 3, 2)
microbenchmark(
   Rcpp = SchurPolCPP(n, lambda),
  Julia = julia$SchurPol(n, lambda),
  times = 6L
)
## Unit: milliseconds
##   expr       min        lq      mean    median        uq        max neval cld
##   Rcpp   5.92909   6.01620   6.11295   6.13169   6.21107    6.25796     6  a
##  Julia 529.34127 532.54833 635.09685 534.60691 538.85875 1140.61894     6   b
```

``` r
JuliaConnectoR::stopJulia()
```
