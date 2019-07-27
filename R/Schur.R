#' Schur polynomial
#'
#' Evaluates the Schur polynomials.
#'
#' @param x numeric vector
#' @param lambda integer partition
#'
#' @return A number.
#' @export
#'
#' @examples x <- c(2,3,4)
#' Schur(x, c(2,1,1))
#' prod(x) * sum(x)
Schur <- function(x, lambda){
  i <- rep(seq_along(lambda), times = lambda)
  j <- unlist(sapply(lambda, seq_len, simplify = FALSE))
  lambdaPrime <- dualPartition(lambda)
  hookslengths <- lambdaPrime[j] - i + lambda[i] - j + 1
  Jack(x, lambda, alpha = 1) / prod(hookslengths)
}
