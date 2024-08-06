#' @title Green X-polynomials
#' @description Computes the Green X-polynomials for a given integer partition.
#' 
#' @param rho an integer partition
#' 
#' @return The Green X-polynomials, usually denoted by \eqn{X_\rho^\lambda}, 
#'   for all integer partitions \eqn{\lambda} of the same weight as the 
#'   given integer partition \eqn{\rho}. 
#'   They are returned in a list. Each element of this list is itself a list, 
#'   with two elements. The first one, called \code{lambda}, represents the 
#'   partition \eqn{\lambda}. The second one, called \code{polynomial}, 
#'   represents the Green X-polynomial \eqn{X_\rho^\lambda}. This is a 
#'   univariate \code{qspray} polynomial whose variable is denoted by 
#'   \code{t}, and all its coefficients are some integers.
#'   The names of the list encode the partitions \eqn{\lambda}. 
#' @importFrom qspray qone PSFpoly showQsprayOption<- showQsprayXYZ
#' @export 
#' @note The Green X-polynomials are a variant of the "true" Green polynomials, 
#'   that we called the Green Q-polynomials (\code{\link{GreenQpolynomials}}).
#' @seealso \code{\link{GreenQpolynomials}}.
#' @examples 
#' GreenXpolynomials(c(2, 1)) 
GreenXpolynomials <- function(rho) {
  stopifnot(isPartition(rho))
  rho <- as.integer(removeTrailingZeros(rho))
  if(length(rho) == 0L) {
    out <- list(
      list("lambda" = integer(0L)),
      list("polynomial" = qone())
    )
    names(out) <- "[]"
    return(out)
  }
  n <- sum(rho)
  psPoly <- PSFpoly(n, rho)
  hlpCombo <- HLPcombination(psPoly)
  lapply(hlpCombo, function(lst) {
    lambda <- lst[["lambda"]]
    qspray <- lst[["coeff"]]
    showQsprayOption(qspray, "showQspray") <- showQsprayXYZ("t")
    list(
      "lambda"     = lambda,
      "polynomial" = qspray
    )
  })
}

#' @title Green Q-polynomials
#' @description Computes the Green Q-polynomials for a given integer partition.
#' 
#' @param rho an integer partition
#' 
#' @return The Green Q-polynomials, usually denoted by \eqn{Q_\rho^\lambda}, 
#'   for all integer partitions \eqn{\lambda} of the same weight as the 
#'   given integer partition \eqn{\rho}. 
#'   They are returned in a list. Each element of this list is itself a list, 
#'   with two elements. The first one, called \code{lambda}, represents the 
#'   partition \eqn{\lambda}. The second one, called \code{polynomial}, 
#'   represents the Green Q-polynomial \eqn{Q_\rho^\lambda}. This is a 
#'   univariate \code{qspray} polynomial whose variable is denoted by 
#'   \code{q}.
#'   The names of the list encode the partitions \eqn{\lambda}. 
#' @importFrom qspray qone PSFpoly showQsprayOption<- showQsprayXYZ
#' @export 
#' @note The Green Q-polynomials are the "true" Green polynomials. 
#'   The Green X-polynomials (\code{\link{GreenXpolynomials}}) are a 
#'   variant of the Green Q-polynomials.
#' @seealso \code{\link{GreenXpolynomials}}.
#' @examples 
#' GreenQpolynomials(c(2, 1)) 
GreenQpolynomials <- function(rho) {
  stopifnot(isPartition(rho))
  rho <- as.integer(removeTrailingZeros(rho))
  if(length(rho) == 0L) {
    out <- list(
      list("lambda" = integer(0L)),
      list("polynomial" = qone())
    )
    names(out) <- "[]"
    return(out)
  }
  n <- sum(rho)
  psPoly <- PSFpoly(n, rho)
  hlpCombo <- HLPcombination(psPoly)
  q <- qlone(1L)
  lapply(hlpCombo, function(lst) {
    lambda <- lst[["lambda"]]
    qspray <- lst[["coeff"]]
    rOQ <- q^(.n(lambda)) * .substitute_invt(qspray)
    qspray <- rOQ@numerator
    showQsprayOption(qspray, "showQspray") <- showQsprayXYZ("q")
    list(
      "lambda"     = lambda,
      "polynomial" = qspray
    )
  })
}