#' @title Factorial Schur polynomial
#' @description Computes a factorial Schur polynomial.
#'
#' @param n number of variables
#' @param lambda integer partition
#' @param a vector of \code{bigq} numbers, or vector of elements coercible
#'   to \code{bigq} numbers; this vector corresponds to the sequence denoted by
#'   \eqn{a} in the
#'   \href{https://www.kurims.kyoto-u.ac.jp/EMIS/journals/SLC/opapers/s28macdonald.pdf}{reference paper},
#'   section \strong{6th Variation} (in this paper \eqn{a} is a doubly
#'   infinite sequence, but in the case of a non-skew partition, the
#'   non-positive indices of this sequence are not involved); the length of
#'   this vector must be large enough (an error will be thrown if it is too
#'   small) but it is not easy to know the minimal possible length
#'
#' @return A \code{qspray} polynomial.
#' @export
#' @importFrom syt all_ssytx
#' @importFrom qspray qlone
#'
#' @references I.G. Macdonald.
#' \emph{Schur functions: theme and variations}.
#' Publ. IRMA Strasbourg, 1992.
#'
#' @examples
#' # for a=c(0, 0, ...), the factorial Schur polynomial is the Schur polynomial
#' n <- 3
#' lambda <- c(2, 2, 2)
#' a <- c(0, 0, 0, 0)
#' factorialSchurPoly <- factorialSchurPol(n, lambda, a)
#' schurPoly <- SchurPol(n, lambda)
#' factorialSchurPoly == schurPoly # should be TRUE
factorialSchurPol <- function(n, lambda, a) {
  stopifnot(isPositiveInteger(n))
  stopifnot(isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  tableaux <- all_ssytx(lambda, n)
  l <- length(lambda)
  i_ <- 1L:l
  qlones <- lapply(1L:n, qlone)
  toAdd <- lapply(tableaux, function(tableau) {
    factors <- lapply(i_, function(i) {
      toMultiply <- lapply(1L:lambda[i], function(j) {
        entry <- tableau[[i]][j]
        k <- entry + j - i
        qlones[[entry]] + a[k]
      })
      Reduce(`*`, toMultiply)
    })
    Reduce(`*`, factors)
  })
  Reduce(`+`, toAdd)
}

# factorialSchurPol n lambda y
#   | otherwise =
#       sumOfSprays sprays
#   where
#     tableaux = semiStandardYoungTableaux n (toPartition lambda)
#     lones = [lone i | i <- [1 .. n]]
#     idx tableau i j =
#       let row = tableau !! (i-1)
#           a = row !! (j-1)
#       in (a, a + j - i)
#     factor tableau i j =
#       let (a, k) = idx tableau i j in lones !! (a-1) <+ y !! (k-1)
#     i_ = [1 .. length lambda]
#     ij_ = [(i, j) | i <- i_, j <- [1 .. lambda !! (i-1)]]
#     factors tableau = [factor tableau i j | (i, j) <- ij_]
#     spray tableau = productOfSprays (factors tableau)
#     sprays = map spray tableaux
