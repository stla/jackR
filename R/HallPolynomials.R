#' @importFrom syt KostkaNumbersWithGivenLambda
#' @importFrom utils tail
#' @noRd
msPolynomialsInSchurBasis <- function(weight) {
  lambdas <- listOfPartitions(weight)
  nparts <- length(lambdas)
  lambdasAsStrings <-
    vapply(lambdas, partitionAsString, character(1L))
  KostkaMatrix <- matrix(0L, nrow = nparts, ncol = nparts)
  colnames(KostkaMatrix) <- lambdasAsStrings
  for(i in seq_len(nparts)) {
    kNumbers <- KostkaNumbersWithGivenLambda(lambdas[[i]], output = "vector")
    KostkaMatrix[i, names(kNumbers)] <- kNumbers
  }
  invKostkaMatrix <- backsolve(KostkaMatrix, diag(nparts))
  storage.mode(invKostkaMatrix) <- "integer"
  out <- lapply(seq_len(nparts), function(i) {
    coeffs <- tail(invKostkaMatrix[i, ], nparts - i + 1L)
    names(coeffs) <- tail(lambdasAsStrings, nparts - i + 1L)
    coeffs
  })
  names(out) <- lambdasAsStrings
  out
}

#' @importFrom qspray isQzero
#' @noRd
msPolynomialInHLPbasis <- function(lambda) {
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
      r * KostaFoulkesPolynomial(mu, kappa)
      # list(
      #   "partition" = kappa,
      #   "qspray" = r * KostaFoulkesPolynomial(mu, kappa)
      # )
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
  names(out) <- lambdasAsStrings
  Filter(Negate(isQzero), out)
}

#' @importFrom methods new
#' @importFrom qspray MSPcombination orderedQspray isQzero showQsprayOption<- showQsprayXYZ
#' @importFrom symbolicQspray Qzero
#' @importFrom ratioOfQsprays as.ratioOfQsprays
#' @noRd
HLPcombination <- function(Qspray) {
  fullMsCombo <- MSPcombination(Qspray, check = FALSE)
  lambdas <- lapply(fullMsCombo, `[[`, "lambda")
  finalQspray <- Qzero()
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  for(lambda in lambdas) {
    hlpCombo <- msPolynomialInHLPbasis(lambda)
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
    spray <- Qzero()
    for(kappa in names(hlpCombo)) {
      coeff <- hlpCombo[[kappa]]
      if(!isQzero(coeff)) {
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
      qspray <- coeff@numerator
      showQsprayOption(qspray, "showQspray") <- showQsprayXYZ("t")
      list("coeff" = qspray, "lambda" = lambda)
    },
    powers, coeffs,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  names(combo) <-
    vapply(powers, partitionAsString, character(1L), USE.NAMES = FALSE)
  combo
}

.substitute_invt <- function(qspray) {
  constantTerm <- getConstantTerm(qspray)
  if(isConstant(qspray)) {
    return(as.ratioOfQsprays(constantTerm))
  }
  qspray <- qspray - constantTerm
  powers <- qspray@powers
  coeffs <- qspray@coeffs
  t <- qlone(1L)
  rOQs <- mapply(
    function(coeff, power) {
      new(
        "ratioOfQsprays",
        numerator = as.qspray(coeff),
        denominator = t^power
      )
    },
    coeffs, powers,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  Reduce(`+`, rOQs) + constantTerm
}

HallPolynomials <- function(mu, nu) {
  stopifnot(isPartition(mu), isPartition(nu))
  n <- sum(mu) + sum(nu)
  Qspray <- HallLittlewoodPol(n, mu, "P") * HallLittlewoodPol(n, nu, "P")
  hlpCombo <- HLPcombination(Qspray)
  t <- qlone(1L)
  .n_mu_nu <- .n(mu) + .n(nu)
  lapply(hlpCombo, function(coeff_lambda) {
    lambda <- coeff_lambda[["lambda"]]
    qspray <- coeff_lambda[["coeff"]]
    rOQ <- t^(.n(lambda) - .n_mu_nu) * .substitute_invt(qspray)
    list(
      "lambda" = lambda,
      "polynomial" = rOQ@numerator
    )
  })
}
