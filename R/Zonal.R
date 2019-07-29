#' Zonal polynomial
#'
#' Evaluates the zonal polynomials.
#'
#' @param x numeric vector
#' @param lambda integer partition
#'
#' @return A number.
#' @export
#'
#' @examples x <- c(3,1)
#' Zonal(x, c(1,1)) + Zonal(x, 2) # sum(x)^2
#' Zonal(x, 3) + Zonal(x, c(2,1)) + Zonal(x, c(1,1,1)) # sum(x)^3
Zonal <- function(x, lambda){
  jack <- Jack(x, lambda, alpha= 2)
  jlambda <- sum(logHookLengths(lambda, alpha = 2))
  n <- sum(lambda)
  exp(n*log(2) + lfactorial(n) - jlambda) * jack
}

ZonalQ <- function(x, lambda){
  jack <- JackQ(x, lambda, alpha= as.bigq(2))
  jlambda <- prod(hookLengths_gmp(lambda, alpha = as.bigq(2)))
  as.bigq(2L)^n * as.bigq(factorialZ(n)) / jlambda * jack
}
