#' Evaluation of Schur polynomials
#'
#' Evaluates the Schur polynomials.
#'
#' @param x numeric vector or \link[gmp]{bigq} vector
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param algorithm the algorithm used, either \code{"DK"} (Demmel-Koev)
#' or \code{"naive"}
#'
#' @return A number or a \code{bigq} rational number.
#' @export
#'
#' @seealso \code{\link{SchurPol}}
#'
#' @examples x <- c(2,3,4)
#' Schur(x, c(2,1,1))
#' prod(x) * sum(x)
Schur <- function(x, lambda, algorithm = "DK"){
  algorithm <- match.arg(algorithm, c("DK", "naive"))
  if(algorithm == "DK"){
    SchurEval(x, lambda)
  }else{
    SchurEvalNaive(x, lambda)
  }
}

