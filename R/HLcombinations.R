
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


## whether a 'ratioOfQsprays' is zero
#' @importFrom qspray isQzero
#' @noRd
.isZeroROQ <- function(rOQ) {
  isQzero(rOQ@numerator)
}

## the `ms[lambda]` polynomial in the Hall-Littlewood Q-polynomials basis
#' @noRd
msPolynomialInHLQbasis <- function(lambda) {
  weight <- sum(lambda)
  msCombos <- msPolynomialsInSchurBasis(weight)
  lambdasAsStrings <- names(msCombos)
  lambdas <- lapply(lambdasAsStrings, fromPartitionAsString)
  lambdaAsString <- partitionAsString(lambda)
  msCombo <- msCombos[[lambdaAsString]]
  musAsStrings <- names(msCombo)
  hlqCombos <- lapply(musAsStrings, function(muAsString) {
    mu <- fromPartitionAsString(muAsString)
    r <- msCombo[muAsString] 
    lapply(lambdas, function(kappa) {
      r * KostkaFoulkesPolynomial(mu, kappa) / b(kappa)
    })
  })
  out <- Reduce(
    function(combo1, combo2) {
      mapply(
        `+`,
        combo1, combo2
      )
    },
    hlqCombos
  )
  names(out) <- lambdasAsStrings
  Filter(Negate(.isZeroROQ), out)
}


#' @title Symmetric polynomial in terms of Hall-Littlewood 
#'   Q-polynomials
#' @description Returns the expression of a symmetric polynomial as a 
#'   linear combination of some Hall-Littlewood Q-polynomials.
#' 
#' @param Qspray a \code{qspray} polynomial or a 
#'   \code{symbolicQspray} polynomial, which must be symmetric
#' @param check Boolean, whether to check the symmetry
#' 
#' @return A list defining the linear combination. Each element of this list 
#'   is itself a list, with two elements: \code{coeff}, which is a 
#'   \code{ratioOfQsprays}, and the second element of this list is 
#'   \code{lambda}, an integer partition; then this list
#'   corresponds to the term \code{coeff * HallLittlewoodPol(n, lambda, "Q")} 
#'   of the linear combination, where \code{n} is the number of variables in 
#'   the symmetric polynomial \code{Qspray}. 
#'   The output list defining the linear combination is named 
#'   by some character strings representing the partitions \code{lambda}.
#' @export 
#' @importFrom methods new
#' @importFrom qspray MSPcombination orderedQspray isQzero
#' @importFrom symbolicQspray Qzero isQzero
#' @importFrom ratioOfQsprays as.ratioOfQsprays
HLcombinationQ <- function(Qspray, check = TRUE) {
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
  fullMsCombo <- MSPcombination(Qspray, check = check)
  lambdas <- lapply(fullMsCombo, `[[`, "lambda")
  finalQspray <- Qzero()
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  for(lambda in lambdas) {
    hlqCombo <- msPolynomialInHLQbasis(lambda)
    kappas <- lapply(names(hlqCombo), fromPartitionAsString)
    msCombo <- fullMsCombo[[partitionAsString(lambda)]]
    sprays <- lapply(kappas, function(kappa) {
      new(
        "symbolicQspray",
        powers = list(kappa),
        coeffs = list(unitRatioOfQsprays)
      )
    })
    names(sprays) <- names(hlqCombo)
    spray <- Qzero()
    for(kappa in names(hlqCombo)) {
      coeff <- hlqCombo[[kappa]]
      #if(!.isZeroROQ(coeff)) { useless, there's Filter
      spray <- spray + coeff * sprays[[kappa]]
      #}
    }
    finalQspray <- finalQspray + msCombo[["coeff"]] * spray
  }
  finalQspray <- orderedQspray(finalQspray)
  powers <- finalQspray@powers
  coeffs <- finalQspray@coeffs
  hlqCombo <- mapply(
    function(lambda, coeff) {
      list("coeff" = coeff, "lambda" = lambda)
    },
    powers, coeffs,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  names(hlqCombo) <-
    vapply(powers, partitionAsString, character(1L), USE.NAMES = FALSE)
  if(constantTerm != 0L) {
    hlqCombo <- c(
      list(
        "[]" = list(
          "coeff" = as.ratioOfQsprays(constantTerm), 
          "lambda" = integer(0L)
        )
      ),
      hlqCombo
    )
  }
  hlqCombo
}
