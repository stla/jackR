The ‘jack’ package: Jack polynomials
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/jackR/workflows/R-CMD-check/badge.svg)](https://github.com/stla/jackR/actions)
<!-- badges: end -->

``` r
library(jack)
library(microbenchmark)
```

As of version 2.0.0, the Jack polynomials can be calculated with Julia.
The speed is amazing:

``` r
julia <- Jack_julia()
## Starting Julia ...
x <- c(1/2, 2/3, 1, 2/3, -1, -2, 1)
lambda <- c(5, 3, 2, 2, 1)
alpha <- 3
microbenchmark(
      R = Jack(x, lambda, alpha),
  Julia = julia$Jack(x, lambda, alpha),
  times = 5
)
## Unit: milliseconds
##   expr          min           lq       mean     median         uq       max
##      R 15767.690900 16131.731601 16546.8662 16381.9088 17123.5399 17329.460
##  Julia     8.395101     9.733001   458.8062    11.0415    90.6334  2174.228
##  neval
##      5
##      5
```

`Jack_julia()` returns a list of functions. `ZonalPol`, `ZonalQPol` and
`SchurPol` always return an exact expression of the polynomial,
i.e. with rational coefficients (integers for `SchurPol`). If you want
an exact expression with `JackPol`, you have to give a rational number
for the argument `alpha`, as a character string:

``` r
JP <- julia$JackPol(m = 2, lambda = c(3,1), alpha = "2/5")
JP
## mvp object algebraically equal to
## 3.92 x1 x2^3  +  5.6 x1^2 x2^2  +  3.92 x1^3 x2
## 
## Exact expression:
## 98/25 * x1 * x2^3  +  28/5 * x1^2 * x2^2  +  98/25 * x1^3 * x2
```

To evaluate a polynomial, use `as.function`:

``` r
jp <- as.function(JP)
```

You can provide the values of the variables of this function as numbers
or character strings:

``` r
jp(2, "3/2")
## [1] "1239/10"
```

The evaluation is performed by the **Ryacas** package. If you want to
substitute a variable with a complex number, use a character string
which represents this number, with `I` denoting the imaginary unit:

``` r
jp("2 + 2*I", "2/3")
## [1] "Complex((-26656)/675,43232/675)"
```

Two functions are provided to print the polynomials with an exact
expression:

``` r
prettyForm(JP)
## 
##   3               2          2               3
## x2  * 98 * x1   x2  * 28 * x1    x2 * 98 * x1 
## ------------- + -------------- + -------------
##      25               5               25
```

``` r
toLaTeX(JP)
## $\frac{x_{2}^{3} 98 x_{1}}{25}  + \frac{x_{2}^{2} 28 x_{1}^{2}}{5}  + \frac{x_{2} 98 x_{1}^{3}}{25} $
```
