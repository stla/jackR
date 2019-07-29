#' Kostka numbers
#'
#' The Kostka numbers for partitions of a given weight.
#'
#' @param n positive integer, the weight of the partitions
#'
#' @return A matrix of integers.
#' @export
#'
#' @examples
#' KostkaNumbers(4)
KostkaNumbers <- function(n){
  sc <- SchurCoefficientsQ(n)
  stringParts <- paste0("(", gsub("(, 0| )", "", colnames(sc)), ")")
  colnames(sc) <- rownames(sc) <- stringParts
  apply(sc, c(1L,2L), as.integer)
}
