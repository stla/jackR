#' @title Power sum polynomial
#' @description Returns a power sum polynomial.
#'
#' @param n integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   positive integers
#'
#' @return A \code{qspray} object.
#' @export
#' @importFrom qspray qone qzero
#' @importFrom methods new
#'
#' @examples
#' library(jack)
#' psPolynomial(3, c(3, 1))
psPolynomial <- function(n, lambda) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  if(length(lambda) == 0L) {
    qone()
  } else if(n == 0L){
    qzero()
  } else {
    n <- as.integer(n)
    coeffs <- rep("1", n)
    n_ <- seq_len(n)
    out <- qone()
    for(k in lambda) {
      powers <- lapply(n_, function(i) {
        c(rep(0L, i-1L), k)
      })
      pk <- new(
        "qspray", powers = powers, coeffs = coeffs
      )
      out <- out * pk
    }
    out    
  } 
}

#' @title Elementary symmetric polynomial
#' @description Returns an elementary symmetric polynomial.
#'
#' @param n integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#'   nonnegative integers
#'
#' @return A \code{qspray} object.
#' @export
#' @importFrom qspray qone qzero qsprayMaker
#' @importFrom methods new
#' @importFrom DescTools Permn
#'
#' @examples
#' library(jack)
#' esPolynomial(3, c(3, 1))
esPolynomial <- function(n, lambda) {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  if(length(lambda) == 0L) {
    qone()
  } else if(lambda[1L] > n){
    qzero()
  } else {
    out <- qone()
    kappa0 <- integer(n)
    for(k in seq_along(lambda)) {
      kappa <- kappa0
      lambda_k <- lambda[k]
      kappa[seq_len(lambda_k)] <- rep(1L, lambda_k)
      perms <- Permn(kappa)
      powers <- Rows(perms)
      ek <- qsprayMaker(
        powers = powers, coeffs = rep("1", length(powers))
      )
      out <- out * ek
    }
    out    
  } 
}