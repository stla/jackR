#' Jack polynomial with symbolic Jack parameter
#'
#' Returns the Jack polynomial with symbolic Jack parameter.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#' @param which which Jack polynomial, \code{"J"}, \code{"P"} or \code{"Q"}
#'
#' @return A \code{symbolicQspray} object.
#'
#' @export
#' @importFrom symbolicQspray symbolicQspray_from_list
#'
#' @examples
#' JackSymPol(3, lambda = c(3, 1))
JackSymPol <- function(n, lambda, which = "J") {
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  which <- match.arg(which, c("J", "P", "Q"))
  JackPolynomial <- symbolicQspray_from_list(
    JackSymPolRcpp(as.integer(n), as.integer(lambda))
  )
  if(which != "J") {
    invK <- switch(
      which,
      "P" = symbolicJackPcoefficientInverse(lambda),
      "Q" = symbolicJackQcoefficientInverse(lambda)
    )
    # K <- new(
    #   "ratioOfQsprays",
    #   numerator   = as.qspray(1L),
    #   denominator = invK
    # )
    JackPolynomial <- JackPolynomial / invK
  }
  JackPolynomial
}
