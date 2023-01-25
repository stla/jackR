#' Schur polynomial - C++ implementation
#'
#' Returns the Schur polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A \code{qspray} multivariate polynomial.
#'
#' @export
#' @importFrom qspray qsprayMaker
#'
#' @examples
#' SchurPolCPP(3, lambda = c(3, 1))
SchurPolCPP <- function(n, lambda) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  x <- SchurPolRcpp(as.integer(n), as.integer(lambda))
  qsprayMaker(x[["exponents"]], x[["coeffs"]])
}

#' Jack polynomial - C++ implementation
#'
#' Returns the Jack polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#' @param alpha positive rational number, given as a string such as
#'   \code{"2/3"} or as a \code{bigq} number
#'
#' @return A \code{qspray} multivariate polynomial.
#'
#' @export
#' @importFrom qspray qsprayMaker
#' @importFrom gmp as.bigq
#'
#' @examples
#' JackPolCPP(3, lambda = c(3, 1), alpha = "2/5")
JackPolCPP <- function(n, lambda, alpha) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  alpha <- as.bigq(alpha)
  stopifnot(alpha >= 0)
  alpha <- as.character(alpha)
  x <- JackPolRcpp(as.integer(n), as.integer(lambda), alpha)
  qsprayMaker(x[["exponents"]], x[["coeffs"]])
}

#' Zonal polynomial - C++ implementation
#'
#' Returns the Zonal polynomial.
#'
#' @param m number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A \code{qspray} multivariate polynomial.
#'
#' @export
#' @importFrom gmp as.bigq factorialZ
#' @examples
#' ZonalPolCPP(3, lambda = c(3, 1))
ZonalPolCPP <- function(m, lambda){
  twoq <- as.bigq("2")
  jack <- JackPolCPP(m, lambda, alpha = "2")
  jlambda <- prod(hookLengths_gmp(lambda, alpha = twoq))
  n <- sum(lambda)
  (twoq^n * factorialZ(n) / jlambda) * jack
}

#' Quaternionic Zonal polynomial - C++ implementation
#'
#' Returns the quaternionic Zonal polynomial.
#'
#' @param m number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A \code{qspray} multivariate polynomial.
#'
#' @export
#' @importFrom gmp as.bigq factorialZ
#' @examples
#' ZonalQPolCPP(3, lambda = c(3, 1))
ZonalQPolCPP <- function(m, lambda){
  onehalfq <- as.bigq("1/2")
  jack <- JackPolCPP(m, lambda, alpha = onehalfq)
  jlambda <- prod(hookLengths_gmp(lambda, alpha = onehalfq))
  n <- sum(lambda)
  as.qspray(onehalfq^n * factorialZ(n) / jlambda) * jack
}
