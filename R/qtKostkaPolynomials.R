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

qtSkewKostkaPolynomials <- function(lambda, mu) {
  stopifnot(isPartition(lambda))
  stopifnot(isPartition(mu))
  lambda <- as.integer(removeTrailingZeros(lambda))
  mu <- as.integer(removeTrailingZeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu || any(lambda[seq_len(ellMu)] < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  w <- sum(lambda) - sum(mu)
  if(w == 0L){
    out <- list(
      list(
        "nu" = integer(0L),
        "polynomial" = qone()
      )
    )
    names(out) <- partitionAsString(integer(0L))
    return(out)
  }
  lrCoeffs <- LRskew(lambda, mu)
  nus <- apply(parts(w), 2L, removeTrailingZeros, simplify = FALSE)
  out <- lapply(nus, function(nu) {
    qtKostkaPolys <- qtKostkaPolynomials(nu)
    pis <- intersect(names(lrCoeffs), names(qtKostkaPolys))
    poly <- Reduce(
      `+`,
      lapply(pis, function(pi) {
        lrCoeffs[[pi]][["coeff"]] * qtKostkaPolys[[pi]][["polynomial"]]
      })
    )
    showSymbolicQsprayOption(poly, "showRatioOfQsprays") <-
      showRatioOfQspraysXYZ(c("q", "t"))
    list(
      "nu" = nu,
      "polynomial" = poly
    )
  })
  names(out) <- vapply(nus, partitionAsString, character(1L))
  out
}
#   => Partition -- ^ outer partition of the skew partition
#   -> Partition -- ^ inner partition of the skew partition
#   -> Map Partition (Spray a)
# qtSkewKostkaPolynomials lambda mu
#   | not (isSkewPartition lambda mu) =
#       error "qtSkewKostkaPolynomials: invalid skew partition."
#   | lambda == mu =
#       DM.singleton [] unitSpray
#   | otherwise =
#       DM.fromList (map spray nus)
#   where
#     lrCoeffs = skewSchurLRCoefficients lambda mu
#     nus = partitions (sum lambda - sum mu)
#     spray nu =
#       let nu' = fromPartition nu in
#         (
#           nu',
#           foldl'
#             (^+^)
#               zeroSpray
#                 (DM.intersectionWith (.^) lrCoeffs (qtKostkaPolynomials nu'))
#         )
