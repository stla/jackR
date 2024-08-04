#' @title Green X-polynomials
#' @description Computes the Green X-polynomials for a given partition.
#' 
#' @param rho integer partition
#' 
#' @return The Green X-polynomials, usually denoted by \eqn{X_\rho^\lambda}, 
#'   for all integer partitions \eqn{\lambda} of the same weight as the 
#'   integer partition \eqn{\rho}. 
#'   They are returned in a list. Each element of this list is itself a list, 
#'   with two elements. The first one, called \code{lambda}, represents the 
#'   partition \eqn{\lambda}. The second one, called \code{polynomial}, 
#'   represents the Green X-polynomial \eqn{X_\rho^\lambda}. This is a 
#'   univariate \code{qspray} polynomial whose variable is denoted by \code{t}.
#'   The names of the list encode the partitions \eqn{\lambda}. 
GreenXpolynomials <- function(rho) {
  stopifnot(isPartition(rho))
  rho <- as.integer(removeTrailingZeros(rho))
  if(length(rho) == 0L) {
    out <- list(
      list("lambda" = integer(0L)),
      list("polynomial" = qzero())
    )
    names(out) <- "[]"
    return(out)
  }
}