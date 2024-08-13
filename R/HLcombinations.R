
## Symmetric polynomial as linear combination of Hall-Littlewood
##   P-polynomials or Q-polynomials.

#' @title Symmetric polynomial in terms of Hall-Littlewood 
#'   P-polynomials
#' @description Returns the expression of a symmetric polynomial as a 
#'   linear combination of some Hall-Littlewood P-polynomials.
#' 
#' @param Qspray a \code{qspray} polynomial or a 
#'   \code{symbolicQspray} polynomial, which must be symmetric
#' @param check Boolean, whether to check the symmetry
#' 
#' @return A list defining the linear combination. Each element of this list 
#'   is itself a list, with two elements: \code{coeff}, which is a 
#'   \code{ratioOfQsprays}, and the second element of this list is 
#'   \code{lambda}, an integer partition; then this list
#'   corresponds to the term \code{coeff * HallLittlewoodPol(n, lambda, "P")} 
#'   of the linear combination, where \code{n} is the number of variables in 
#'   the symmetric polynomial \code{Qspray}. 
#'   The output list defining the linear combination is named 
#'   by some character strings encoding the partitions \code{lambda}.
#' @export 
#' @importFrom qspray isQzero getConstantTerm isConstant
#' @importFrom symbolicQspray isQzero getConstantTerm isConstant
#' @importFrom ratioOfQsprays as.ratioOfQsprays
HLcombinationP <- function(Qspray, check = TRUE) {
  stopifnot(inherits(Qspray, "qspray") || inherits(Qspray, "symbolicQspray"))
  stopifnot(isBoolean(check))
  constantTerm <- getConstantTerm(Qspray)
  if(isConstant(Qspray)) {
    if(isQzero(Qspray)) {
      out <- list()
    } else {
      out <-
        list(
          list(
            "coeff" = as.ratioOfQsprays(constantTerm), 
            "lambda" = integer(0L)
          )
        )
      names(out) <- "[]"
    }
    return(out)
  }
  Qspray <- Qspray - constantTerm
  hlpCombo <- .HLcombinationP(Qspray, check = check, takeNumerators = FALSE)
  if(constantTerm != 0L) {
    hlpCombo <- c(
      list(
        "[]" = list(
          "coeff" = as.ratioOfQsprays(constantTerm), 
          "lambda" = integer(0L)
        )
      ),
      hlpCombo
    )
  }
  hlpCombo
}