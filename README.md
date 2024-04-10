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

Each of these polynomials is given by a positive integer, the number of
variables, and an integer partition, the `lambda` argument; the Jack
polynomial has one more parameter, the `alpha` argument, a number called
the *Jack parameter*.

To get an exact symbolic polynomial with `JackPol`, you have to supply a
`bigq` rational number for the Jack parameter `alpha`:

``` r
jpol <- JackPol(2, lambda = c(3, 1), alpha = gmp::as.bigq("2/5"))
jpol
## 98/25*x^3y + 28/5*x^2y^2 + 98/25*xy^3
```

This is a `qspray` object, from the
[**qspray**](https://github.com/stla/qspray) package. Here is how you
can evaluate this polynomial:

``` r
evalQspray(jpol, c("2", "3/2"))
## Big Rational ('bigq') :
## [1] 1239/10
```

By default, `ZonalPol`, `ZonalQPol` and `SchurPol` return exact symbolic
polynomials.

``` r
zpol <- ZonalPol(2, lambda = c(3, 1))
zpol
## 24/7*x^3y + 16/7*x^2y^2 + 24/7*xy^3
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
    times = 5L,
    unit  = "seconds"
  ),
  signif = 2L
)
## Unit: seconds
##  expr   min    lq mean median   uq   max neval cld
##     R 60.00 61.00 78.0  62.00 90.0 120.0     5  a 
##  Rcpp  0.56  0.57  0.9   0.87  1.1   1.4     5   b
```

## Skew Schur polynomials

As of version 6.0.0, the package is able to compute the skew Schur
polynomials, with the function `SkewSchurPol`.

## Symbolic Jack parameter

As of version 6.0.0, it is possible to get a Jack polynomial with a
symbolic Jack parameter in its coefficients, with the `JackSymPol`
function:

``` r
( J <- JackSymPol(2, lambda = c(3, 1)) )
## { [2*a^2 + 4*a + 2] } * X^3Y  +  { [4*a + 4] } * X^2Y^2  +  { [2*a^2 + 4*a + 2] } * XY^3
```

This is a `symbolicQspray` object, from the
[**symbolicQspray**](https://github.com/stla/symbolicQspray) package.

A `symbolicQspray` object corresponds to a multivariate polynomial whose
coefficients are fractions of polynomials with rational coefficients.
The Jack polynomials fit into this category: from their definition,
their coefficients are fractions of polynomials in the Jack parameter.
However you can see in the above output that for this example, the
coefficients are *polynomials* in the Jack parameter (`a`): there’s no
fraction. Actually this is always true for any Jack polynomial. This
fact is established and it is not obvious.

Note that you can change the letters used to denote the variables. By
default, the Jack parameter is denoted by `a` and the variables are
denoted by `X`, `Y`, `Z` if there are at most three variables, otherwise
they are denoted by `X1`, `X2`, … Here is how to change these symbols:

``` r
showSymbolicQsprayOption(J, "a") <- "alpha"
showSymbolicQsprayOption(J, "X") <- "x"
J
## { [2*alpha^2 + 4*alpha + 2] } * x1^3.x2  +  { [4*alpha + 4] } * x1^2.x2^2  +  { [2*alpha^2 + 4*alpha + 2] } * x1.x2^3
```

If you want to have the variables denoted by `x` and `y`, do:

``` r
showSymbolicQsprayOption(J, "showMonomial") <- showMonomialXYZ(c("x", "y"))
J
## { [2*alpha^2 + 4*alpha + 2] } * x^3y  +  { [4*alpha + 4] } * x^2y^2  +  { [2*alpha^2 + 4*alpha + 2] } * xy^3
```

## Compact expression of Jack polynomials

The expression of a Jack polynomial in the canonical basis can be long.
Since these polynomials are symmetric, one can get a considerably
shorter expression by writing the polynomial as a linear combination of
the monomial symmetric polynomials. This is what the function
`compactSymmetricQspray` does:

``` r
( J <- JackPolCPP(3, lambda = c(4, 3, 1), alpha = "2") )
## 3888*x^4y^3z + 2592*x^4y^2z^2 + 3888*x^4yz^3 + 3888*x^3y^4z + 4752*x^3y^3z^2 + 4752*x^3y^2z^3 + 3888*x^3yz^4 + 2592*x^2y^4z^2 + 4752*x^2y^3z^3 + 2592*x^2y^2z^4 + 3888*xy^4z^3 + 3888*xy^3z^4
cat(compactSymmetricQspray(J))
## (3888) * M[4, 3, 1] + (2592) * M[4, 2, 2] + (4752) * M[3, 3, 2]
```

The function `compactSymmetricQspray` is from the **qspray** package but
it is also possible to apply it to a `symbolicQspray` object, like a
Jack polynomial with symbolic Jack parameter.
