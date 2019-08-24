#' Evaluation of Jack polynomials
#'
#' Evaluates a Jack polynomial.
#'
#' @param x numeric vector or \link[gmp]{bigq} vector
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param alpha positive number or \code{bigq} rational number
#' @param algorithm the algorithm used, either \code{"DK"} (Demmel-Koev)
#' or \code{"naive"}
#'
#' @return A numeric scalar or a \code{bigq} rational number.
#' @export
#' @importFrom partitions conjugate
#' @importFrom gmp factorialZ is.bigq
#'
#' @seealso \code{\link{JackPol}}
#'
#' @references \itemize{
#' \item I.G. Macdonald.
#' \emph{Symmetric Functions and Hall Polynomials}.
#' Oxford Mathematical Monographs.
#' The Clarendon Press Oxford University Press,
#' New York, second edition, 1995.
#' \item J. Demmel & P. Koev.
#' \emph{Accurate and efficient evaluation of Schur and Jack functions}.
#' Mathematics of computations, vol. 75, n. 253, 223-229, 2005.
#' \item \emph{Jack polynomials}.
#' \url{https://www.math.upenn.edu/~peal/polynomials/jack.htm}
#' }
#'
#' @examples lambda <- c(2,1,1)
#' Jack(c(1/2, 2/3, 1), lambda, alpha = 3)
#' # exact value:
#' Jack(c(gmp::as.bigq(1,2), gmp::as.bigq(2,3), gmp::as.bigq(1)), lambda,
#'      alpha = gmp::as.bigq(3))
Jack <- function(x, lambda, alpha, algorithm = "DK"){
  if(alpha == 0){
    stopifnot(isPartition(lambda))
    lambdaPrime <- conjugate(lambda)
    if(is.bigq(x)){
      f <- as.bigq(prod(factorialZ(lambdaPrime[lambdaPrime>0L])))
    }else{
      f <- prod(factorial(lambdaPrime[lambdaPrime>0L]))
    }
    return(f * ESF(x, lambdaPrime))
  }
  algorithm <- match.arg(algorithm, c("DK", "naive"))
  if(algorithm == "DK"){
    JackEval(x, lambda, alpha)
  }else{
    JackEvalNaive(x, lambda, alpha)
  }
}

