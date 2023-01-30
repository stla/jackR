# Version 5.0.1 (2023-01-30)

Changed C++17 to C++20.


# Version 5.0.0 (2023-01-27)

- Now there is a 'Rcpp' implementation of the polynomials: functions 
`SchurPolCPP`, `JackPolCPP`, `ZonalPolCPP` and `ZonalQPolCPP`. It is faster 
than the Julia implementation.


# Version 4.0.0 (2022-12-19)

- The package does not longer depend on the 'gmpoly' package. This dependency 
has been replaced with the 'qspray' package.


# Version 3.0.0 (2022-02-21)

- Now one can use a `bigq` number for `alpha` in `JackPol`, thanks to the 
'gmpoly' package, and one can use `exact=TRUE` with `algorithm=DK` for 
`ZonalPol`, `ZonalQPol` and `SchurPol`.

- Now one can get a `gmpoly` polynomial with Julia. 


# Version 2.0.0 (2022-02-14)

- New function `Jack_julia`, to evaluate the polynomials with Julia.


# Version 1.1.1 (2019-09-16)

- Fixed a test of empty partition

- Added more checks of parameters validity

- Added more unit tests

- Improved documentation


# Version 1.1.0 (2019-09-09)

- Some functions didn't handle the empty partition.


