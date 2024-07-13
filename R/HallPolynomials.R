#' @importFrom syt KostkaNumbersWithGivenLambda
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

# _msPolynomialInHLPbasis ::
#   Int -> Partition -> Map Partition (Spray Rational)
# _msPolynomialInHLPbasis n lambda =
#   DM.filter (not . isZeroSpray) (DM.unionsWith (^+^) hlpCombos)
#   where
#     weight = sum lambda
#     msCombos = msPolynomialsInSchurBasis n weight
#     lambdas = DM.keys msCombos
#     hlpCombo mu =
#       DM.filter (not . isZeroSpray) $
#         DM.fromDistinctAscList
#           (map (\kappa -> (kappa, _kostkaFoulkesPolynomial mu kappa)) lambdas)
#     msAssocs = DM.assocs (msCombos DM.! lambda)
#     hlpCombos =
#       map
#         (\(mu, r) ->
#           DM.map (\spray -> r *^ spray) (hlpCombo mu))
#         msAssocs

HLPcombination <- function(Qspray) {
  fullMsCombo <- MSPcombination(Qspray, check = FALSE)
  lambdas <- lapply(fullMsCombo, `[[`, "lambda")
  finalQspray <- Qzero()
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  for(lambda in lambdas) {
    hlpCombo <- msPolynomialInHLPbasis(lambda)
    kappas <- lapply(names(hlpCombo), fromPartitionAsString)
    msCombo <- fullMsCombo[[partitionAsString(lambda)]]
#    coeffs <- lapply(msCombo, `[[`, "coeff")
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

    # lambdas <- names(msCombo)
    # for(i in seq_along(lambdas)) {
    #   spray <- qzero()
    #   for(kappa in names(hlpCombo)) {
    #     coeff <- hlpCombo[[kappa]]
    #     if(!isQzero(coeff)) {
    #       spray <- spray + coeff * sprays[[kappa]]
    #     }
    #   }
    #   finalQspray <- finalQspray + coeffs[[i]]@numerator*spray
    # }
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
