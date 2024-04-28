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
#' @importFrom qspray qspray_from_list
#'
#' @examples
#' SchurPol(3, lambda = c(3, 1))
SchurPol <- function(n, lambda) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  qspray_from_list(SchurPolRcpp(as.integer(n), as.integer(lambda)))
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
#' Schur(c("1", "3/2", "-2/3"), lambda = c(3, 1))
Schur <- function(x, lambda) {
  stopifnot(isPartition(lambda))
  if(is.numeric(x)) {
    if(anyNA(x)) {
      stop("Found missing values in `x`.")
    }
    SchurEvalRcpp_double(as.double(x), as.integer(lambda))
  } else {
    x <- as.character(as.bigq(x))
    if(anyNA(x)) {
      stop("Found missing values in `x`.")
    }
    res <- SchurEvalRcpp_gmpq(x, as.integer(lambda))
    as.bigq(res)
  }
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
#' @param which which Jack polynomial, \code{"J"}, \code{"P"} or \code{"Q"}
#'
#' @return A \code{qspray} multivariate polynomial.
#'
#' @export
#' @importFrom gmp as.bigq
#' @importFrom qspray qspray_from_list
#'
#' @examples
#' JackPol(3, lambda = c(3, 1), alpha = "2/5")
JackPol <- function(n, lambda, alpha, which = "J") {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  alpha <- as.bigq(alpha)
  which <- match.arg(which, c("J", "P", "Q"))
  alpha <- as.character(alpha)
  JackPolynomial <-
    qspray_from_list(JackPolRcpp(as.integer(n), as.integer(lambda), alpha))
  if(which != "J") {
    K <- switch(
      which,
      "P" = 1L / prod(hookLengths_gmp(lambda, alpha)[1L, ]),
      "Q" = 1L / prod(hookLengths_gmp(lambda, alpha)[2L, ])
    )
    JackPolynomial <- K * JackPolynomial
  }
  JackPolynomial
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
#' Jack(c("1", "3/2", "-2/3"), lambda = c(3, 1), alpha = "1/4")
Jack <- function(x, lambda, alpha) {
  stopifnot(isPartition(lambda))
  if(is.numeric(x) && is.numeric(alpha)) {
    stopifnot(alpha >= 0)
    x <- as.double(x)
    if(anyNA(x)) {
      stop("Found missing values in `x`.")
    }
    JackEvalRcpp_double(x, as.integer(lambda), alpha)
  } else {
    alpha <- as.bigq(alpha)
    stopifnot(alpha >= 0)
    alpha <- as.character(alpha)
    x <- as.character(as.bigq(x))
    if(anyNA(x)) {
      stop("Found missing values in `x`.")
    }
    res <- JackEvalRcpp_gmpq(x, as.integer(lambda), alpha)
    as.bigq(res)
  }
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
#' ZonalPol(3, lambda = c(3, 1))
ZonalPol <- function(m, lambda){
  twoq <- as.bigq("2")
  jack <- JackPol(m, lambda, alpha = "2")
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
#' @importFrom gmp as.bigq factorialZ asNumeric
#'
#' @examples
#' Zonal(c("1", "3/2", "-2/3"), lambda = c(3, 1))
Zonal <- function(x, lambda){
  if(is.numeric(x)) {
    jack <- Jack(x, lambda, alpha = 2)
    jlambda <- asNumeric(prod(hookLengths_gmp(lambda, alpha = as.bigq("2"))))
    n <- sum(lambda)
    (factorial(n) * 2^n / jlambda) * jack
  } else {
    twoq <- as.bigq("2")
    jack <- Jack(x, lambda, alpha = "2")
    jlambda <- prod(hookLengths_gmp(lambda, alpha = twoq))
    n <- sum(lambda)
    (twoq^n * factorialZ(n) / jlambda) * jack
  }
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
#' ZonalQPol(3, lambda = c(3, 1))
ZonalQPol <- function(m, lambda){
  onehalfq <- as.bigq("1/2")
  jack <- JackPol(m, lambda, alpha = onehalfq)
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
#' @importFrom gmp as.bigq factorialZ asNumeric
#'
#' @examples
#' ZonalQ(c("1", "3/2", "-2/3"), lambda = c(3, 1))
ZonalQ <- function(x, lambda){
  if(is.numeric(x)) {
    jack <- Jack(x, lambda, alpha = 0.5)
    jlambda <- asNumeric(prod(hookLengths_gmp(lambda, alpha = as.bigq("1/2"))))
    n <- sum(lambda)
    (factorial(n) / 2^n / jlambda) * jack
  } else {
    onehalfq <- as.bigq("1/2")
    jack <- Jack(x, lambda, alpha = "1/2")
    jlambda <- prod(hookLengths_gmp(lambda, alpha = onehalfq))
    n <- sum(lambda)
    (onehalfq^n * factorialZ(n) / jlambda) * jack
  }
}
