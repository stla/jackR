
## whether a 'ratioOfQsprays' is zero
#' @importFrom qspray isQzero
#' @noRd
.isZeroROQ <- function(rOQ) {
  isQzero(rOQ@numerator)
}

## the `ms[lambda]` polynomial in the Hall-Littlewood Q-polynomials basis
#' @importFrom qspray isQzero
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

## the `Qspray` polynomial in the Hall-Littlewood Q-polynomials basis
#' @importFrom methods new
#' @importFrom qspray MSPcombination orderedQspray isQzero
#' @importFrom symbolicQspray Qzero
#' @importFrom ratioOfQsprays as.ratioOfQsprays
#' @noRd
HLcombinationQ <- function(Qspray, check = TRUE) {
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
  combo <- mapply(
    function(lambda, coeff) {
      list("coeff" = coeff, "lambda" = lambda)
    },
    powers, coeffs,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  names(combo) <-
    vapply(powers, partitionAsString, character(1L), USE.NAMES = FALSE)
  combo
}


test1 <- function(lambda, alpha, which) {
  n <- sum(lambda)
  p <- JackPol(n, lambda, alpha, which)
  combo <- HLcombinationQ(p)
  toAdd <- lapply(combo, function(lst) {
    mu <- lst$lambda
    coeff <- lst$coeff
    HallLittlewoodPol(n, mu, "Q") * coeff
  })
  obtained <- Reduce(`+`, toAdd)
  pp <- as(p, "symbolicQspray")
  obtained == pp
} 

test1(c(2, 2, 1), "3", "C")

lambda = c(2, 1)
alpha = "3"
which = "C"


test2 <- function(lambda, which) {
  n <- sum(lambda)
  p <- JackSymPol(n, lambda, which)
  combo <- HLcombinationQ(p)
  toAdd <- lapply(combo, function(lst) {
    mu <- lst$lambda
    coeff <- lst$coeff
    HallLittlewoodPol(n, mu, "Q") * coeff
  })
  obtained <- Reduce(`+`, toAdd)
  #pp <- as(p, "symbolicQspray")
  obtained == p
} 

test2(c(2, 2), "C")

