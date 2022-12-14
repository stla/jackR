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
##   expr       min        lq     mean    median       uq      max neval cld
##      R 6.9101300 6.9297900 6.995080 7.0125200 7.020270 7.085240     6   b
##  Julia 0.0467379 0.0524001 0.210701 0.0719558 0.107837 0.913322     6  a
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

You can even pass a variable name to this function:

``` r
jp("x", "x")
## [1] "(336*x^4)/25"
```

This evaluation is performed by the **Ryacas** package. If you want to
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
##        3               2     2               3
## 98 * x1  * x2   28 * x1  * x2    98 * x1 * x2 
## ------------- + -------------- + -------------
##      25               5               25
```

``` r
toLaTeX(JP)
## $\frac{98 x_{1}^{3} x_{2}}{25}  + \frac{28 x_{1}^{2} x_{2}^{2}}{5}  + \frac{98 x_{1} x_{2}^{3}}{25} $
```

You can also use the functions `JackPol`, `ZonalPol`, `ZonalQPol` and
`SchurPol` without passing by `Jack_julia()`. They are implemented in R.
To get an exact symbolic polynomial with `JackPol`, you have to supply a
`bigq` rational number for the parameter `alpha`:

``` r
jpol <- JackPol(2, lambda = c(3, 1), alpha = gmp::as.bigq("2/5"))
jpol
## 98/25*x^(3, 1) + 98/25*x^(1, 3) + 28/5*x^(2, 2)
```

This is a `qspray` object, from the
[**qspray**](https://github.com/stla/qspray) package.

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
##   expr       min        lq      mean    median        uq      max neval cld
##      R 1017.5010 1028.0519 1059.4227 1043.2110 1067.5932 1156.968     6   b
##  Julia  783.0187  799.3045  882.7926  869.3771  884.4896 1091.189     6  a
```

As of version 3.0.0, one can also get a `gmpoly` polynomial with Julia,
by setting the argument `poly` to `"gmpoly"` in the `XXXPol` functions
in the list returned by `Jack_julia`:

``` r
julia$JackPol(2, lambda = c(3, 1), alpha = "2/5", poly = "qspray")
## 98/25*x^(1, 3) + 28/5*x^(2, 2) + 98/25*x^(3, 1)
```

``` r
n <- 5
lambda <- c(4, 3, 3)
alpha <- "2/3"
microbenchmark(
     Julia_mvp = julia$JackPol(n, lambda, alpha),
  Julia_gmpoly = julia$JackPol(n, lambda, alpha, poly = "qspray"),
  times = 6L
)
## Unit: milliseconds
##          expr      min       lq     mean   median       uq      max neval cld
##     Julia_mvp 827.5496 855.1814 873.4309 882.5008 895.3674 897.4856     6   b
##  Julia_gmpoly 418.1035 434.2616 491.9679 492.6028 552.8773 561.3593     6  a
```
