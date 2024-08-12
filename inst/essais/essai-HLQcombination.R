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
  hlpCombos <- lapply(musAsStrings, function(muAsString) {
    mu <- fromPartitionAsString(muAsString)
    r <- msCombo[muAsString] 
    lapply(lambdas, function(kappa) {
      rOQ <- r * KostkaFoulkesPolynomial(mu, kappa) / b(kappa)
      rOQ
    })
  })
  out <- Reduce(
    function(combo1, combo2) {
      mapply(
        `+`,
        combo1, combo2
      )
    },
    hlpCombos
  )
  #out <- lapply(out, function(rOQ) rOQ@numerator)
  names(out) <- lambdasAsStrings
  return(out)
  Filter(Negate(isQzero), out)
}

## the `Qspray` polynomial in the Hall-Littlewood Q-polynomials basis
#' @importFrom methods new
#' @importFrom qspray MSPcombination orderedQspray isQzero
#' @importFrom symbolicQspray Qzero
#' @importFrom ratioOfQsprays as.ratioOfQsprays
#' @noRd
HLQcombination <- function(Qspray) {
  fullMsCombo <- MSPcombination(Qspray, check = FALSE)
  lambdas <- lapply(fullMsCombo, `[[`, "lambda")
  finalQspray <- Qzero()
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  for(lambda in lambdas) {
    hlpCombo <- msPolynomialInHLQbasis(lambda)
    kappas <- lapply(names(hlpCombo), fromPartitionAsString)
    msCombo <- fullMsCombo[[partitionAsString(lambda)]]
    sprays <- lapply(kappas, function(kappa) {
      new(
        "symbolicQspray",
        powers = list(kappa),
        coeffs = list(unitRatioOfQsprays)
      )
    })
    names(sprays) <- names(hlpCombo)
    spray <- as.ratioOfQsprays(0L) #Qzero()
    for(kappa in names(hlpCombo)) {
      coeff <- hlpCombo[[kappa]]
      if(TRUE){#!isQzero(coeff)) {
        spray <- spray + coeff * sprays[[kappa]]
      }
    }
    finalQspray <- finalQspray + msCombo[["coeff"]]*spray
  }
  finalQspray <- orderedQspray(finalQspray)
  powers <- finalQspray@powers
  coeffs <- finalQspray@coeffs
  combo <- mapply(
    function(lambda, coeff) {
      qspray <- coeff#@numerator
      list("coeff" = qspray, "lambda" = lambda)
    },
    powers, coeffs,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  names(combo) <-
    vapply(powers, partitionAsString, character(1L), USE.NAMES = FALSE)
  combo
}



p <- JackPol(3, c(2, 1), alpha = "2")
co <- HLQcombination(p)

x <- qlone(1)
HallLittlewoodPol(3, c(2,1), "Q") * co[["[2, 1]"]]$coeff + HallLittlewoodPol(3, c(1,1,1), "Q") * co[["[1, 1, 1]"]]$coeff