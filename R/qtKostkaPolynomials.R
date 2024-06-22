# qtKostkaPolynomials ::
#   forall a. (Eq a, AlgField.C a)
#   => Partition
#   -> Map Partition (Spray a)
# qtKostkaPolynomials mu =  DM.map _numerator scs
#   where
#     psCombo = macdonaldJinPSbasis mu
#     t = lone' 2
# --    f = psCombo -- psCombination'' (macdonaldJpolynomial n mu)
#     den lambda = productOfSprays [unitSpray ^-^ t k | k <- lambda]
#     msCombo lambda =
#       msCombination (psPolynomial (length lambda) lambda)
#     ikn = inverseKostkaNumbers (sum mu)
#     coeffs lambda =
#       let combo = msCombo lambda in
#         DM.map
#           (\ikNumbers ->
#             DF.sum $ DM.intersectionWith (*) combo ikNumbers)
#           ikn
#     scs = DM.foldlWithKey
#       (\m lambda c ->
#         DM.unionWith (AlgAdd.+) m
#           (DM.map (\ikNumber -> (ikNumber .^ c) %//% den lambda) (coeffs lambda))
#       )
#       DM.empty psCombo


#' @importFrom qspray MSPcombination PSFpoly qlone qone showQsprayOption<- showQsprayXYZ
#' @importFrom RationalMatrix Qinverse
#' @importFrom partitions parts
#' @importFrom gmp c_bigq
qtKostkaPolynomials <- function(mu) {
  stopifnot(isPartition(mu))
  n <- sum(mu)
  if(n == 0L) {
    out <- list(
      list(
        "lambda" = integer(0L),
        "polynomial" = qone()
      )
    )
    names(out) <- partitionAsString(integer(0L))
    return(out)
  }
  psCombo <- MacdonaldPolynomialJinPSbasis(mu)
  iknMatrix <- Qinverse(KostkaNumbers(n))
  lambdas <- apply(parts(n), 2L, removeTrailingZeros, simplify = FALSE)
  lambdasAsStrings <-
    vapply(lambdas, partitionAsString, character(1L))
  rownames(iknMatrix) <- lambdasAsStrings
  coeffs <- function(lambda) {
    combo <- MSPcombination(PSFpoly(length(lambda), lambda), check = FALSE)
    out <- lapply(seq_along(lambdas), function(i) {
      list(
        "lambda" = lambdas[[i]],
        "coeff" = sum(c_bigq(lapply(names(combo), function(p) {
          combo[[p]][["coeff"]] * iknMatrix[p, i]
        })))
      )
    })
    names(out) <- lambdasAsStrings
    out
  }
  coeffsMap <- lapply(lambdas, coeffs)
  names(coeffsMap) <- lambdasAsStrings
  t <- qlone(2L)
  unitSpray <- qone()
  maps <- lapply(names(psCombo), function(p) {
    lambda <- psCombo[[p]][["lambda"]]
    spray <- psCombo[[p]][["coeff"]]
    coeffs_lambda <- coeffsMap[[p]]
    den_lambda <- Reduce(
      `*`,
      lapply(lambda, function(k) {
        unitSpray - t^k
      })
    )
    lapply(lambdasAsStrings, function(lambdaAsString) {
      coeff_lambda <- coeffs_lambda[[lambdaAsString]]
      list(
        "lambda" = coeff_lambda[["lambda"]],
        "lambdaAsString" = lambdaAsString,
        "coeff" =  coeff_lambda[["coeff"]] * spray / den_lambda
      )
    })
  })
  cmapOfMaps <- do.call(c, maps)
  f <- vapply(cmapOfMaps, `[[`, character(1L), "lambdaAsString")
  lapply(split(cmapOfMaps, f), function(l) {
    rOQ <- Reduce(`+`, lapply(l, `[[`, "coeff"))
    spray <- rOQ@numerator
    showQsprayOption(spray, "showQspray") <- showQsprayXYZ(c("q", "t"))
    list(
      "lambda" = l[[1L]][["lambda"]],
      "polynomial" = spray
    )
  })
}
