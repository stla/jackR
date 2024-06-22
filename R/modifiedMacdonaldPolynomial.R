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

# modifiedMacdonaldPolynomial n mu = jp

#     nmu = sum (zipWith (*) [1 .. ] (drop 1 mu))
#     jp = DM.foldlWithKey
#       (\spray lambda c ->
#           spray ^+^
#             _numerator (toROS (t' (nmu + sum lambda) ^*^ c) %/% den lambda)
#               *^ psPolynomial n lambda)
#       zeroSpray psCombo

#' @importFrom symbolicQspray showSymbolicQsprayOption<- Qone Qzero
#' @importFrom ratioOfQsprays showRatioOfQspraysXYZ
#' @importFrom methods as
#' @importFrom qspray qlone qone
modifiedMacdonaldPolynomial <- function(n, mu) {
  psCombo <- MacdonaldPolynomialJinPSbasis(mu)
  nmu <- sum(seq_len(length(mu) - 1L) * tail(mu, -1L))
  t <- qlone(2L)
  unitSpray <- qone()
  out <- Reduce(
    `+`,
    lapply(psCombo, function(term) {
      lambda <- term[["lambda"]]
      spray <- term[["coeff"]]
      den_lambda <- Reduce(
        `*`,
        lapply(lambda, function(k) {
          t^k - unitSpray
        })
      )
      rOS <- .toROS(t^(nmu + sum(lambda)) * spray) / den_lambda
      rOS@numerator * as(PSFpoly(n, lambda), "symbolicQspray")
    })
  )
  showSymbolicQsprayOption(out, "showRatioOfQsprays") <-
    showRatioOfQspraysXYZ(c("q", "t"))
  out
}
#   where
#     psCombo = macdonaldJinPSbasis mu
#     q' = lone' 1
#     t' = lone' 2
#     -- toROS spray =
#     --   evaluateAt [q', tm1] (HM.map constantRatioOfSprays spray)
#     num_and_den Empty = undefined
#     num_and_den (e :<| Empty) = (q' e, unitSpray)
#     num_and_den (e1 :<| (e2 :<| _)) = (q' e1, t' e2)
#     rOS_from_term powers coeff = coeff *^ RatioOfSprays spray1 spray2
#       where
#         (spray1, spray2) = num_and_den (exponents powers)
.rOS_from_term <- function(powers, coeff) {
  q <- qlone(1)
  t <- qlone(2)
  if(length(powers) == 1L) {
    coeff * as.ratioOfQsprays(q^powers)
  } else {
    coeff *
      new(
        "ratioOfQsprays",
        numerator = q^(powers[1L]),
        denominator = t^(powers[2L])
      )
  }
}
#     toROS spray =
#       HM.foldlWithKey'
#         (\ros powers coeff -> ros AlgAdd.+ rOS_from_term powers coeff)
#           zeroRatioOfSprays spray
.toROS <- function(spray) {
  Reduce(
    `+`,
    mapply(
      .rOS_from_term,
      spray@powers, spray@coeffs,
      USE.NAMES = FALSE, SIMPLIFY = FALSE
    )
  )
}
#     den lambda = productOfSprays [t' k ^-^ unitSpray | k <- lambda]
