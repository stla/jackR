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

#' Evaluation of Schur polynomial - C++ implementation
#'
#' Evaluates the Schur polynomial.
#'
#' @param x variables, a vector of \code{bigq} numbers, or a vector that can
#'   be coerced as such (e.g. \code{c("2", "5/3")})
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A \code{bigq} number.
#'
#' @export
#' @importFrom gmp as.bigq
#'
#' @examples
#' SchurCPP(c("1", "3/2", "-2/3"), lambda = c(3, 1))
SchurCPP <- function(x, lambda) {
  stopifnot(isPartition(lambda))
  x <- as.character(as.bigq(x))
  if(anyNA(x)) {
    stop("Found missing values in `x`.")
  }
  res <- SchurEvalRcpp(x, as.integer(lambda))
  as.bigq(res)
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

#' Evaluation of Jack polynomial - C++ implementation
#'
#' Evaluates the Jack polynomial.
#'
#' @param x variables, a vector of \code{bigq} numbers, or a vector that can
#'   be coerced as such (e.g. \code{c("2", "5/3")})
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#' @param alpha positive rational number, given as a string such as
#'   \code{"2/3"} or as a \code{bigq} number
#'
#' @return A \code{bigq} number.
#'
#' @export
#' @importFrom gmp as.bigq
#'
#' @examples
#' JackCPP(c("1", "3/2", "-2/3"), lambda = c(3, 1), alpha = "1/4")
JackCPP <- function(x, lambda, alpha) {
  stopifnot(isPartition(lambda))
  alpha <- as.bigq(alpha)
  stopifnot(alpha >= 0)
  alpha <- as.character(alpha)
  x <- as.character(as.bigq(x))
  if(anyNA(x)) {
    stop("Found missing values in `x`.")
  }
  res <- JackEvalRcpp(x, as.integer(lambda), alpha)
  as.bigq(res)
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

#' Evaluation of zonal polynomial - C++ implementation
#'
#' Evaluates the zonal polynomial.
#'
#' @param x variables, a vector of \code{bigq} numbers, or a vector that can
#'   be coerced as such (e.g. \code{c("2", "5/3")})
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A \code{bigq} number.
#'
#' @export
#' @importFrom gmp as.bigq factorialZ
#'
#' @examples
#' ZonalCPP(c("1", "3/2", "-2/3"), lambda = c(3, 1))
ZonalCPP <- function(x, lambda){
  twoq <- as.bigq("2")
  jack <- JackCPP(x, lambda, alpha = "2")
  jlambda <- prod(hookLengths_gmp(lambda, alpha = twoq))
  n <- sum(lambda)
  (twoq^n * factorialZ(n) / jlambda) * jack
}

#' Quaternionic zonal polynomial - C++ implementation
#'
#' Returns the quaternionic zonal polynomial.
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

#' Evaluation of zonal quaternionic polynomial - C++ implementation
#'
#' Evaluates the zonal quaternionic polynomial.
#'
#' @param x variables, a vector of \code{bigq} numbers, or a vector that can
#'   be coerced as such (e.g. \code{c("2", "5/3")})
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A \code{bigq} number.
#'
#' @export
#' @importFrom gmp as.bigq factorialZ
#'
#' @examples
#' ZonalQCPP(c("1", "3/2", "-2/3"), lambda = c(3, 1))
ZonalQCPP <- function(x, lambda){
  onehalfq <- as.bigq("1/2")
  jack <- JackCPP(x, lambda, alpha = "1/2")
  jlambda <- prod(hookLengths_gmp(lambda, alpha = onehalfq))
  n <- sum(lambda)
  (onehalfq^n * factorialZ(n) / jlambda) * jack
}
