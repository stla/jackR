#' @importFrom qspray MSFpoly
#' @importFrom DescTools Permn
#' @noRd 
msPowersMatrix <- function(n, lambda) {
  kappa <- integer(n)
  kappa[seq_along(lambda)] <- lambda
  Permn(kappa)
} 

msPowers <- function(n, lambda) {
  M <- msPowersMatrix(n, lambda)
  lapply(seq_len(nrow(M)), function(i) {
    removeTrailingZeros(M[i, ])
  })
} 

#' @importFrom spray spray one zero
#' @noRd 
MSFspray <- function(n, lambda) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  ellLambda <- length(lambda)
  if(ellLambda == 0L) {
    return(one)
  }
  if(ellLambda > n) {
    return(zero)
  }
  Mpowers <- msPowersMatrix(n, lambda)
  spray(Mpowers, rep(1, nrow(Mpowers)))
}

#' Evaluation of monomial symmetric functions
#'
#' Evaluates a monomial symmetric function.
#'
#' @param x a numeric vector or a \code{\link[gmp]{bigq}} vector
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#'
#' @return A number if \code{x} is numeric, a \code{bigq} rational number
#'   if \code{x} is a \code{bigq} vector.
#' 
#' @importFrom gmp as.bigq is.bigq
#' @importFrom utils head
#' @noRd
MSF <- function(x, lambda){
  stopifnot(isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  gmp <- is.bigq(x)
  n <- length(x)
  ellLambda <- length(lambda)
  if(ellLambda == 0L) {
    return(if(gmp) as.bigq(1L) else 1L)
  }
  if(ellLambda > n) {
    return(if(gmp) as.bigq(0L) else 0L)
  }
  powers <- msPowers(n, lambda)
  if(gmp) {
    out <- as.bigq(0L)
    for(exponents in powers){
      m <- length(exponents)
      factors <- as.bigq(integer(m))
      for(j in seq_len(m)){
        factors[j] <- x[j]^exponents[j]
      }
      out <- out + prod(factors)
    }
  } else {
    out <- 0L
    for(exponents in powers){
      m <- length(exponents)
      out <- out + prod(head(x, m)^exponents)
    }
  }
  out
}

#' @importFrom qspray qzero qone
#' @importFrom methods new
#' @noRd
msPolynomialUnsafe <- function(n, lambda) {
  ellLambda <- length(lambda)
  if(ellLambda == 0L) {
    return(qone())
  }
  if(ellLambda > n) {
    return(qzero())
  }
  powers <- msPowers(n, lambda)
  # kappa <- integer(n)
  # kappa[seq_len(ellLambda)] <- lambda
  # perms <- Permn(kappa)
  # powers <- apply(perms, 1L, function(perm) {
  #   removeTrailingZeros(perm)
  # }, simplify = FALSE)
  new(
    "qspray", powers = powers, coeffs = rep("1", length(powers))
  )
}

#' @title Monomial symmetric polynomial
#' @description Returns a monomial symmetric polynomial.
#'
#' @param n integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   nonnegative integers
#'
#' @return A \code{qspray} object.
#' @export
#'
#' @examples
#' library(jack)
#' msPolynomial(3, c(3, 1))
msPolynomial <- function(n, lambda) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  msPolynomialUnsafe(n, lambda)
}