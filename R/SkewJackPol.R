#' Skew Jack polynomial
#' @description Computes a skew Jack polynomial.
#'
#' @param n positive integer, the number of variables
#' @param lambda outer integer partition of the skew partition
#' @param mu inner integer partition of the skew partition; it must be a
#'   subpartition of \code{lambda}
#' @param alpha the Jack parameter, an integer or a \code{bigq} number, positive
#' @param which which Jack polynomial, \code{"J"}, \code{"P"} or \code{"Q"}
#'
#' @return A \code{qspray} polynomial.
#' @export
#' @importFrom gmp is.bigq as.bigq
#' @importFrom partitions parts
#' @importFrom qspray PSPexpression HallInnerProduct
#'
#' @details
#' Gröbner bases are used in the algorithm, and this is very slow.
#'
#' @examples
#' SkewJackPol(3, c(3,1), c(2), 2)
SkewJackPol <- function(n, lambda, mu, alpha, which = "J") {
  stopifnot(isPositiveInteger(n))
  stopifnot(isPartition(lambda), isPartition(mu))
  mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  if(any(lambda - mu < 0L)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  stopifnot(isPositiveInteger(alpha) || is.bigq(alpha))
  alpha <- as.bigq(alpha)
  stopifnot(alpha >= 0L)
  Jlambda <- PSPexpression(JackPolCPP(n, lambda, alpha, which))
  Jmu     <- JackPolCPP(n, mu, alpha, which)
  nus <- parts(sum(lambda) - sum(mu))
  terms <- apply(nus, 2L, function(nu) {
    Jnu <- JackPolCPP(n, nu, alpha, which)
    coeff <- HallInnerProduct(Jlambda, Jmu * Jnu, alpha) /
      HallInnerProduct(Jnu, Jnu, alpha)
    coeff * Jnu
  }, simplify = FALSE)
  Reduce(`+`, terms)
}
