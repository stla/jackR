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
##   expr        min         lq      mean     median         uq      max neval
##      R 14.6556000 14.8779000 15.213200 15.1759000 15.6298000 15.76430     6
##  Julia  0.0091368  0.0106161  0.369453  0.0114276  0.0138549  2.16026     6
```

`Jack_julia()` returns a list of functions. `ZonalPol`, `ZonalQPol` and
`SchurPol` always return an exact expression of the polynomial,
i.e. with rational coefficients (integers for `SchurPol`). If you want
an exact expression with `JackPol`, you have to give a rational number
for the argument `alpha`, as a character string:

``` r
JP <- julia$JackPol(m = 2, lambda = c(3, 1), alpha = "2/5")
JP
## mvp object algebraically equal to
## 3.92 x_1 x_2^3  +  5.6 x_1^2 x_2^2  +  3.92 x_1^3 x_2
## 
## Exact expression:
## 98/25 * x1^3 * x2  +  28/5 * x1^2 * x2^2  +  98/25 * x1 * x2^3
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
##             3     2          2     3          
## x2 * 98 * x1    x2  * 28 * x1    x2  * 98 * x1
## ------------- + -------------- + -------------
##      25               5               25
```

``` r
toLaTeX(JP)
## $\frac{x_{2} 98 x_{1}^{3}}{25}  + \frac{x_{2}^{2} 28 x_{1}^{2}}{5}  + \frac{x_{2}^{3} 98 x_{1}}{25} $
```

You can also use the functions `JackPol`, `ZonalPol`, `ZonalQPol` and
`SchurPol` without passing by `Jack_julia()`. They are implemented in R.
To get an exact symbolic polynomial with `JackPol`, you have to supply a
`bigq` rational number for the parameter `alpha`:

``` r
jpol <- JackPol(2, lambda = c(3, 1), alpha = gmp::as.bigq("2/5"))
jpol
## gmpoly object algebraically equal to
## 98/25 x^(1,3) + 28/5 x^(2,2) + 98/25 x^(3,1)
```

This is a `gmpoly` object, from the
[gmpoly](https://github.com/stla/gmpoly) package.

``` r
gmpoly::gmpolyEval(jpol, c(gmp::as.bigq("2"), gmp::as.bigq("3/2")))
## Big Rational ('bigq') :
## [1] 1239/10
```

By default, `ZonalPol`, `ZonalQPol` and `SchurPol` return exact symbolic
polynomials.

``` r
zpol <- ZonalPol(2, lambda = c(3, 1))
zpol
## gmpoly object algebraically equal to
## 24/7 x^(1,3) + 16/7 x^(2,2) + 24/7 x^(3,1)
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
##   expr       min       lq     mean    median        uq      max neval
##      R 5944.5445 6126.592 6518.818 6473.0548 6699.5025 7396.161     6
##  Julia  877.5917  936.721 1031.260  948.7592  993.8618 1481.870     6
```

As of version 3.0.0, one can also get a `gmpoly` polynomial with Julia,
by setting the argument `poly` to `"gmpoly"` in the `XXXPol` functions
in the list returned by `Jack_julia`:

``` r
julia$JackPol(2, lambda = c(3, 1), alpha = "2/5", poly = "gmpoly")
## gmpoly object algebraically equal to
## 98/25 x^(1,3) + 28/5 x^(2,2) + 98/25 x^(3,1)
```

``` r
n <- 5
lambda <- c(4, 3, 3)
alpha <- "2/3"
microbenchmark(
     Julia_mvp = julia$JackPol(n, lambda, alpha),
  Julia_gmpoly = julia$JackPol(n, lambda, alpha, poly = "gmpoly"),
  times = 6L
)
## Unit: milliseconds
##          expr      min       lq     mean   median       uq      max neval
##     Julia_mvp 867.1429 882.4804 906.8851 886.4314 901.5080 1017.316     6
##  Julia_gmpoly 847.6498 855.8250 925.9224 938.1641 974.2659 1001.465     6
```
