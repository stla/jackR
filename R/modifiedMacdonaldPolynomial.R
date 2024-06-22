MacdonaldPolynomialJinMSPbasis <- function(lambda) {
  mus <- listOfDominatedPartitions(lambda)
  listsOfPairs <- lapply(mus, function(mu) {
    lapply(GelfandTsetlinPatterns(lambda, mu), function(pattern) {
      pairing(gtPatternDiagonals(pattern))
    })
  })
  allPairs <- unique(
    do.call(
      c,
      do.call(
        c,
        listsOfPairs
      )
    )
  )
  pairsMap <- lapply(allPairs, function(pair) {
    psiLambdaMu(pair[[1L]], pair[[2L]])
  })
  names(pairsMap) <- vapply(allPairs, toString, character(1L))
  c <- clambda(lambda)
  lapply(seq_along(mus), function(i) {
    mu <- mus[[i]]
    listOfPairs <- listsOfPairs[[i]]
    rOQ <- c * Reduce(`+`, lapply(listOfPairs, function(pairs) {
      makeRatioOfSprays(pairsMap, pairs)
    }))
    list(
      "mu" = mu,
      "coeff" = rOQ@numerator
    )
  })
}

# macdonaldJinPSbasis ::
#   (Eq a, AlgField.C a) => Partition -> Map Partition (Spray a)
# macdonaldJinPSbasis mu =
#   DM.filter (not . isZeroSpray)
#     (unionsWith (^+^) (DM.elems $ DM.mapWithKey combo_to_map macdonaldCombo))
#   where
#     macdonaldCombo = macdonaldJinMSPbasis mu
#     combo_to_map lambda spray =
#       DM.map
#         (\r -> fromRational r *^ spray)
#           (mspInPSbasis lambda)
MacdonaldPolynomialJinPSbasis <- function(mu) {
  macdonaldCombo <- MacdonaldPolynomialJinMSPbasis(mu)
  mapOfMaps <- lapply(macdonaldCombo, function(t1) {
    lambda <- t1[["mu"]]
    spray <- t1[["coeff"]]
    lapply(PSPcombination(MSFpoly(sum(lambda), lambda)), function(t2) {
      list(
        "lambda" = t2[["lambda"]],
        "lambdaAsString" = partitionAsString(t2[["lambda"]]),
        "coeff" = t2[["coeff"]] * spray
      )
    })
    # out <- Reduce(
    #   `+`,
    #   lapply(PSPcombination(MSFpoly(sum(lambda), lambda)), function(t2) {
    #     list("lambda" = t2[["lambda"]], "coeff" = t2[["coeff"]] * spray)
    #   })
    # )
    # if(isQzero(out)) {
    #   NULL
    # } else {
    #   list("lambda" = lambda, "coeff" = out)
    # }
  })
  cmapOfMaps <- do.call(c, mapOfMaps)
  f <- vapply(cmapOfMaps, `[[`, character(1L), "lambdaAsString")
  Filter(
    Negate(is.null),
    lapply(split(cmapOfMaps, f), function(l) {
      spray <- Reduce(`+`, lapply(l, `[[`, "coeff"))
      if(isQzero(spray)) {
        NULL
      } else {
        list(
          "lambda" = l[[1L]][["lambda"]],
          "coeff" = spray
        )
      }
    })
  )
}
