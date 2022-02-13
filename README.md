The ‘jack’ package: Jack polynomials
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/jackR/workflows/R-CMD-check/badge.svg)](https://github.com/stla/jackR/actions)
<!-- badges: end -->

``` r
library(jack)
library(microbenchmark)
```

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
##   expr        min         lq       mean     median         uq       max neval
##      R 14883.8315 15254.1332 15897.9037 16104.9934 16440.5175 16806.043     5
##  Julia     6.7366     6.7952   377.3046    10.6128    27.4705  1834.908     5
```
