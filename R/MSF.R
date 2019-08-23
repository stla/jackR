#' Evaluation of monomial symmetric functions
#'
#' Evaluates a monomial symmetric function.
#'
#' @param x a numeric vector or a \code{\link[gmp]{bigq}} vector
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#'
#' @return A number if \code{x} is numeric, a \code{bigq} rational number
#' if \code{x} is a \code{bigq} vector.
#' @importFrom gmp as.bigq is.bigq
#' @importFrom DescTools Permn
#' @export
#'
#' @examples x <- c(1, 2, 5/2)
#' lambda <- c(3, 1)
#' MSF(x, lambda)
#' library(gmp)
#' x <- c(as.bigq(1), as.bigq(2), as.bigq(5,2))
#' MSF(x, lambda)
MSF <- function(x, lambda){
  stopifnot(isPartition(lambda))
  gmp <- is.bigq(x)
  m <- length(x)
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(if(gmp) as.bigq(0L) else 0)
  kappa <- numeric(m)
  kappa[seq_along(lambda)] <- lambda
  perms <- DescTools::Permn(kappa)
  if(gmp){
    out <- as.bigq(0L)
    for(i in 1L:nrow(perms)){
      pows <- as.bigq(rep(0L, m))
      for(j in 1L:m){
        pows[j] <- x[j]^perms[i,j]
      }
      out <- out + prod(pows)
    }
  }else{
    out <- 0
    for(i in 1L:nrow(perms)){
      out <- out + prod(x^perms[i,])
    }
  }
  out
}

#' Monomial symmetric function
#'
#' Returns a monomial symmetric function as a polynomial.
#'
#' @param m integer, the number of variables
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}).
#' @importFrom mvp constant product
#' @importFrom DescTools Permn
#' @export
#'
#' @examples
#' MSFpoly(3, c(3,1))
MSFpoly <- function(m, lambda){
  stopifnot(m > 0L, isPositiveInteger(m), isPartition(lambda))
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(mvp::constant(0))
  kappa <- numeric(m)
  kappa[seq_along(lambda)] <- lambda
  perms <- DescTools::Permn(kappa)
  out <- mvp::constant(0)
  vars <- paste0("x_", 1L:m)
  for(i in 1L:nrow(perms)){
    out <- out + mvp::product(perms[i,], symbols = vars)
  }
  out
}

