#' Kostka numbers
#'
#' The Kostka numbers for partitions of a given weight.
#'
#' @param n positive integer, the weight of the partitions
#' @param alpha the Jack parameter, a \code{bigq} number or an object coercible
#'   to a \code{bigq} number; setting \code{alpha=NULL} is equivalent to set
#'   \code{alpha=1}
#'
#' @return A matrix of character strings representing integers or fractions.
#' @export
#' @importFrom gmp as.bigq c_bigq
#'
#' @examples
#' KostkaNumbers(4)
KostkaNumbers <- function(n, alpha = NULL) {
  if(is.null(alpha)) {
    Knumbers <- SchurCoefficientsQ(n)
    stringParts <- paste0("(", gsub("(, 0| )", "", colnames(Knumbers)), ")")
  } else {
    alpha <- as.bigq(alpha)
    if(is.na(alpha)) {
      stop("Invalid `alpha`.")
    }
    Jcoeffs <- JackCoefficientsQ(n, alpha)
    JackPcoeffs <- c_bigq(lapply(rownames(Jcoeffs), function(lambda) {
      JackPcoefficient(fromString(lambda), alpha)
    }))
    Knumbers <- as.character(as.bigq(Jcoeffs)*JackPcoeffs)
    stringParts <- paste0("(", gsub("(, 0| )", "", colnames(Jcoeffs)), ")")
  }
  colnames(Knumbers) <- rownames(Knumbers) <- stringParts
  Knumbers
}
