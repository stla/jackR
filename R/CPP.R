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
