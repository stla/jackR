#' Evaluation of zonal polynomials
#'
#' Evaluates the zonal polynomials.
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
#' @seealso \code{\link{ZonalPol}}
#'
#' @examples lambda <- c(2,2)
#' Zonal(c(1,1), lambda)
#' Zonal(c(gmp::as.bigq(1),gmp::as.bigq(1)), lambda)
#' ##
#' x <- c(3,1)
#' Zonal(x, c(1,1)) + Zonal(x, 2) # sum(x)^2
#' Zonal(x, 3) + Zonal(x, c(2,1)) + Zonal(x, c(1,1,1)) # sum(x)^3
Zonal <- function(x, lambda, algorithm = "DK"){
  algorithm <- match.arg(algorithm, c("DK", "naive"))
  if(algorithm == "DK"){
    ZonalEval(x, lambda)
  }else{
    ZonalEvalNaive(x, lambda)
  }
}
