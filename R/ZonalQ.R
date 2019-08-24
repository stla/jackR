#' Evaluation of quaternionic zonal polynomials
#'
#' Evaluates a quaternionic zonal polynomial.
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
#' @seealso \code{\link{ZonalQPol}}
#'
#' @references F. Li, Y. Xue. \emph{Zonal polynomials and hypergeometric
#' functions of quaternion matrix argument}.
#' Comm. Statist. Theory Methods, 38 (8), 1184-1206, 2009
#'
#' @examples lambda <- c(2,2)
#' ZonalQ(c(3,1), lambda)
#' ZonalQ(c(gmp::as.bigq(3),gmp::as.bigq(1)), lambda)
#' ##
#' x <- c(3,1)
#' ZonalQ(x, c(1,1)) + ZonalQ(x, 2) # sum(x)^2
#' ZonalQ(x, 3) + ZonalQ(x, c(2,1)) + ZonalQ(x, c(1,1,1)) # sum(x)^3
ZonalQ <- function(x, lambda, algorithm = "DK"){
  algorithm <- match.arg(algorithm, c("DK", "naive"))
  if(algorithm == "DK"){
    ZonalQEval(x, lambda)
  }else{
    ZonalQEvalNaive(x, lambda)
  }
}
