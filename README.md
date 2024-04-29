The ‘jack’ package: Jack polynomials
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/jackR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/jackR/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check-valgrind](https://github.com/stla/jackR/actions/workflows/R-CMD-check-valgrind.yaml/badge.svg)](https://github.com/stla/jackR/actions/workflows/R-CMD-check-valgrind.yaml)
<!-- badges: end -->

``` r
library(jack)
```

Schur polynomials have applications in combinatorics and zonal
polynomials have applications in multivariate statistics. They are
particular cases of [Jack
polynomials](https://en.wikipedia.org/wiki/Jack_function "Jack polynomials on Wikipedia").
This package allows to evaluate these polynomials and also to compute
them in symbolic form.

## Breaking change in version 6.0.0

In version 6.0.0, each function whose name ended with the suffix `CPP`
(`JackCPP`, `JackPolCPP`, etc.) has been renamed by removing this
suffix, and the functions `Jack`, `JackPol`, etc. have been renamed by
adding the suffix `R` to their name: `JackR`, `JackPolR`, etc. The
reason of these changes is that a name like `Jack` is more appealing
than `JackCPP` and it is more sensible to assign the more appealing
names to the functions implemented with **Rcpp** since they are highly
more efficient. The interest of the functions `JackR`, `JackPolR`, etc.
is meager.

## Getting the polynomials

The functions `JackPol`, `ZonalPol`, `ZonalQPol` and `SchurPol`
respectively return the Jack polynomial, the zonal polynomial, the
quaternionic zonal polynomial, and the Schur polynomial.

Each of these polynomials is given by a positive integer, the number of
variables (the `n` argument), and an integer partition (the `lambda`
argument); the Jack polynomial has a parameter in addition, the `alpha`
argument, a number called the *Jack parameter*.

To get a Jack polynomial with `JackPol`, you have to supply the Jack
parameter as a `bigq` rational number or as a character string
representing a fraction, e.g. `"2/5"`:

``` r
jpol <- JackPol(2, lambda = c(3, 1), alpha = "2/5")
jpol
## 98/25*x^3.y + 28/5*x^2.y^2 + 98/25*x.y^3
```

This is a `qspray` object, from the [**qspray**
package](https://github.com/stla/qspray "the 'qspray' package on Github").
Here is how you can evaluate this polynomial:

``` r
evalQspray(jpol, c("2", "3/2"))
## Big Rational ('bigq') :
## [1] 1239/10
```

It is also possible to convert a `qspray` polynomial to a function whose
evaluation is performed by the **Ryacas** package:

``` r
jyacas <- as.function(jpol)
```

You can provide the values of the variables of this function as numbers
or character strings:

``` r
jyacas(2, "3/2")
## [1] "1239/10"
```

You can even pass a variable name to this function:

``` r
jyacas("x", "x")
## [1] "(336*x^4)/25"
```

If you want to substitute a variable with a complex number, use a
character string which represents this number, with `I` denoting the
imaginary unit:

``` r
jyacas("2 + 2*I", "2/3")
## [1] "Complex((-26656)/675,43232/675)"
```

## Direct evaluation of the polynomials

If you just have to evaluate a Jack polynomial, you don’t need to resort
to a `qspray` polynomial: you can use the functions `Jack`, `Zonal`,
`ZonalQ` or `Schur`, which directly evaluate the polynomial; this is
much more efficient than computing the `qspray` polynomial and then
applying `evalQspray`.

``` r
Jack(c("2", "3/2"), lambda = c(3, 1), alpha = "2/5")
## Big Rational ('bigq') :
## [1] 1239/10
```

However, if you have to evaluate a Jack polynomial for several values,
it could be better to resort to the `qspray` polynomial.

## Skew Schur polynomials

As of version 6.0.0, the package is able to compute the skew Schur
polynomials, with the function `SkewSchurPol`.

## Symbolic Jack parameter

As of version 6.0.0, it is possible to get a Jack polynomial with a
symbolic Jack parameter in its coefficients, thanks to the
[**symbolicQspray**
package](https://github.com/stla/symbolicQspray "the 'symbolicQspray' package on Github").

``` r
( J <- JackSymPol(2, lambda = c(3, 1)) )
## { [ 2*a^2 + 4*a + 2 ] } * X^3.Y  +  { [ 4*a + 4 ] } * X^2.Y^2  +  { [ 2*a^2 + 4*a + 2 ] } * X.Y^3
```

This is a `symbolicQspray` object, from the **symbolicQspray** package.

A `symbolicQspray` object corresponds to a multivariate polynomial whose
coefficients are fractions of polynomials with rational coefficients.
The variables of these fractions of polynomials can be seen as some
parameters. The Jack polynomials fit into this category: from their
definition, their coefficients are fractions of polynomials in the Jack
parameter. However you can see in the above output that for this
example, the coefficients are *polynomials* in the Jack parameter (`a`):
there’s no fraction. Actually this fact is always true for any Jack
polynomial (for any *J*-Jack polynomial, I should say). This is an
established fact and it is not obvious (it is a consequence of the [Knop
& Sahi
formula](https://en.wikipedia.org/wiki/Jack_function#Combinatorial_formula "Knop & Sahi formula")).

You can substitute a value to the Jack parameter with the help of the
`substituteParameters` function:

``` r
( J5 <- substituteParameters(J, 5) )
## 72*X^3.Y + 24*X^2.Y^2 + 72*X.Y^3
J5 == JackPol(2, lambda = c(3, 1), alpha = "5")
## [1] TRUE
```

Note that you can change the letters used to denote the variables. By
default, the Jack parameter is denoted by `a` and the variables are
denoted by `X`, `Y`, `Z` if there are no more than three variables,
otherwise they are denoted by `X1`, `X2`, … Here is how to change these
symbols:

``` r
showSymbolicQsprayOption(J, "a") <- "alpha"
showSymbolicQsprayOption(J, "X") <- "x"
J
## { [ 2*alpha^2 + 4*alpha + 2 ] } * x1^3.x2  +  { [ 4*alpha + 4 ] } * x1^2.x2^2  +  { [ 2*alpha^2 + 4*alpha + 2 ] } * x1.x2^3
```

If you want to have the variables denoted by `x` and `y`, do:

``` r
showSymbolicQsprayOption(J, "showMonomial") <- showMonomialXYZ(c("x", "y"))
J
## { [ 2*alpha^2 + 4*alpha + 2 ] } * x^3.y  +  { [ 4*alpha + 4 ] } * x^2.y^2  +  { [ 2*alpha^2 + 4*alpha + 2 ] } * x.y^3
```

## Compact expression of Jack polynomials

The expression of a Jack polynomial in the canonical basis can be long.
Since these polynomials are symmetric, one can get a considerably
shorter expression by writing the polynomial as a linear combination of
the monomial symmetric polynomials. This is what the function
`compactSymmetricQspray` does:

``` r
( J <- JackPol(3, lambda = c(4, 3, 1), alpha = "2") )
## 3888*x^4.y^3.z + 2592*x^4.y^2.z^2 + 3888*x^4.y.z^3 + 3888*x^3.y^4.z + 4752*x^3.y^3.z^2 + 4752*x^3.y^2.z^3 + 3888*x^3.y.z^4 + 2592*x^2.y^4.z^2 + 4752*x^2.y^3.z^3 + 2592*x^2.y^2.z^4 + 3888*x.y^4.z^3 + 3888*x.y^3.z^4
cat(compactSymmetricQspray(J))
## 3888*M[4, 3, 1] + 2592*M[4, 2, 2] + 4752*M[3, 3, 2]
```

The function `compactSymmetricQspray` is also applicable to a
`symbolicQspray` object, like a Jack polynomial with symbolic Jack
parameter.

It is easy to figure out what is a monomial symmetric polynomial:
`M[i, j, k]` is the sum of all monomials `x^p.y^q.z^r` where `(p, q, r)`
is a permutation of `(i, j, k)`.

<!-- -------------------- links -------------------- -->
