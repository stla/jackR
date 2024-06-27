# _skewJackInMSPbasis ::
#   forall a. (AlgRing.C a)
#   => (([((Int, Int), Int)], [((Int, Int), Int)]) -> a)
#   -> (Partition -> Partition -> a)
#   -> Char
#   -> Partition
#   -> Partition
#   -> Map Partition a
.skewJackInMSPbasis <- function(func, ccoeff, which, lambda, mu) {
  nus <- listOfDominatedPartitions(
    lastSubpartition(sum(lambda) - sum(mu), lambda)
  )
  listsOfPatterns <- lapply(nus, function(nu) {
    skewGelfandTsetlinPatterns(lambda, mu, nu)
  })
  i_ <- which(lengths(listsOfPatterns) != 0L)
  listsOfPairs <- lapply(listsOfPatterns[i_], function(patterns) {
    lapply(patterns, function(pattern) {
      pairing(apply(pattern, 1L, removeTrailingZeros, simplify = FALSE))
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
  if(which == "Q") {
    funcLambdaMu <- phiLambdaMu
  } else {
    funcLambdaMu <- psiLambdaMu
  }
  pairsMap <- lapply(allPairs, function(pair) {
    funcLambdaMu(pair[[1L]], pair[[2L]])
  })
  names(pairsMap) <- vapply(allPairs, toString, character(1L))
  nus <- nus[i_]
  names(listsOfPairs) <- vapply(nus, partitionAsString, character(1L))
  makeAssocsFromPairs <- function(pairs) {
    pairsOfMatrices <- pairsMap[vapply(pairs, toString, character(1L))]
    matrix1 <- do.call(
      rbind,
      lapply(pairsOfMatrices, `[[`, 1L)
    )
    matrix2 <- do.call(
      rbind,
      lapply(pairsOfMatrices, `[[`, 2L)
    )
    simplifyTheTwoMatrices(matrix1, matrix2)
  }
  makeCoeffFromListOfPairs <- function(listOfPairs) {
    coeff <-
      Reduce(
        `+`,
        lapply(listOfPairs, function(pairs) {
          func(makeAssocsFromPairs(pairs))
        })
      )
    if(which == "J") {
      c <- func(.clambdamuMatrices(lambda, mu))
      c * coeff
    } else if(which == "C") {
      c <- func(.clambdamuMatrices(lambda, mu))
      cc <- ccoeff(lambda, mu)
      c * cc * coeff
    } else {
      coeff
    }
  }
  mapply(
    function(listOfPairs, nu) {
      list(
        "nu" = nu,
        "coeff" = makeCoeffFromListOfPairs(listOfPairs)
      )
    },
    listsOfPairs, nus,
    USE.NAMES = TRUE, SIMPLIFY = FALSE
  )
}
# _skewJackInMSPbasis func ccoeff which lambda mu =
#   DM.map makeCoeffFromListOfPairs mapOfPairs
#   where
#     nus =
#       dominatedPartitions
#         (toPartitionUnsafe (lastSubPartition (sum lambda - sum mu) lambda))
#     pairing lambdas = zip (drop1 lambdas) lambdas
#     mapOfPatterns = DM.filter (not . null)
#       (DM.fromList (map (\nu ->
#         let nu' = fromPartition nu in
#           (
#             nu'
#           , skewGelfandTsetlinPatterns lambda mu nu'
#           )
#         ) nus))
#     mapOfPairs = DM.map (map pairing) mapOfPatterns
#     listsOfPairs = DM.elems mapOfPairs
#     allPairs = nub $ concat (concat listsOfPairs)
#     funcLambdaMu = if which == 'Q' then phiLambdaMu else psiLambdaMu
#     pairsMap =
#       DM.fromList (zip allPairs (map funcLambdaMu allPairs))
#     makeAssocsFromPairs ::
#       [PartitionsPair] -> ([((Int, Int), Int)], [((Int, Int), Int)])
#     makeAssocsFromPairs pairs = assocsFromMaps num_map den_map
#       where
#         als =
#           both (S.fromList . concat)
#             (unzip (DM.elems $ DM.restrictKeys pairsMap (DS.fromList pairs)))
#         (num_map, den_map) = both alMapFromPairs als
#     makeCoeffFromListOfPairs :: [[PartitionsPair]] -> a
#     makeCoeffFromListOfPairs listOfPairs
#       | which == 'J' =
#           c AlgRing.* coeff
#       | which == 'C' =
#           ccoeff lambda mu AlgRing.* c AlgRing.* coeff
#       | otherwise =
#           coeff
#       where
#         c = func (clambdamuAssocs (S.fromList lambda) (S.fromList mu))
#         coeff = AlgAdd.sum (map (func . makeAssocsFromPairs) listOfPairs)
skewSymbolicJackInMSPbasis <- function(which, lambda, mu) {
  alpha <- qlone(1L)
  poly_from_alc <- function(alc) {
    (alc[1L] * alpha + alc[2L])^(alc[3L])
  }
  poly_from_alcs <- function(alcs) {
    Reduce(`*`, apply(alcs, 1L, poly_from_alc, simplify = FALSE))
  }
  rosFromMatrices <- function(alcsMatrices) {
    matrix1 <- alcsMatrices[[1L]]
    if(nrow(matrix1) >= 1L) {
      num <- poly_from_alcs(matrix1)
    } else {
      num <- qone()
    }
    matrix2 <- alcsMatrices[[2L]]
    if(nrow(matrix2) >= 1L) {
      den <- poly_from_alcs(matrix2)
    } else {
      den <- qone()
    }
    num / den
  }
  ccoeff <- function(.lambda, .mu) {
    symbolicJackCcoefficient(.lambda) / symbolicJackCcoefficient(.mu)
  }
  .skewJackInMSPbasis(rosFromMatrices, ccoeff, which, lambda, mu)
}

skewJackInMSPbasis <- function(alpha, which, lambda, mu) {
  coeff_from_alc <- function(alc) {
    (alc[1L] * alpha + alc[2L])^(alc[3L])
  }
  coeff_from_alcs <- function(alcs) {
    product(c_bigq(apply(alcs, 1L, coeff_from_alc, simplify = FALSE)))
  }
  ratioFromMatrices <- function(alcsMatrices) {
    matrix1 <- alcsMatrices[[1L]]
    if(nrow(matrix1) >= 1L) {
      num <- coeff_from_alcs(matrix1)
    } else {
      num <- 1L
    }
    matrix2 <- alcsMatrices[[2L]]
    if(nrow(matrix2) >= 1L) {
      den <- coeff_from_alcs(matrix2)
    } else {
      den <- 1L
    }
    num / den
  }
  ccoeff <- function(.lambda, .mu) {
    JackCcoefficient(.lambda, alpha) / JackCcoefficient(.mu, alpha)
  }
  .skewJackInMSPbasis(ratioFromMatrices, ccoeff, which, lambda, mu)
}
# skewJackInMSPbasis alpha =
#   _skewJackInMSPbasis ratioFromAssocs ccoeff
#   where
#     coeff ((a, l), c) =
#       (a .^ alpha AlgAdd.+ (_fromInt l)) AlgRing.^ (toInteger c)
#     ratioFromAssocs assocs = num AlgField./ den
#       where
#         (num, den) = both (AlgRing.product . (map coeff)) assocs
#     ccoeff lambda mu = jackCoeffC lambda alpha AlgField./ jackCoeffC mu alpha

#' Skew Jack polynomial
#' @description Computes a skew Jack polynomial.
#'
#' @param n positive integer, the number of variables
#' @param lambda outer integer partition of the skew partition
#' @param mu inner integer partition of the skew partition; it must be a
#'   subpartition of \code{lambda}
#' @param alpha the Jack parameter, an integer or a \code{bigq} number
#' @param which which Jack polynomial, \code{"J"}, \code{"P"}, \code{"Q"} or
#'   \code{"C"}
#'
#' @return A \code{qspray} polynomial.
#' @noRd
#' @importFrom gmp is.bigq as.bigq
#' @importFrom partitions parts
#' @importFrom qspray HallInnerProduct
#'
#' @examples
#' SkewJackPol(3, c(3,1), c(2), 2)
SkewJackPol <- function(n, lambda, mu, alpha, which = "J") {
  stopifnot(isPositiveInteger(n))
  stopifnot(isPartition(lambda), isPartition(mu))
  mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  if(any(lambda - mu < 0L)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  stopifnot(isInteger(alpha) || is.bigq(alpha))
  alpha <- as.bigq(alpha)
  Jlambda <- JackPol(n, lambda, alpha, which)
  Jmu     <- JackPol(n, mu, alpha, which)
  nus <- parts(sum(lambda) - sum(mu))
  terms <- apply(nus, 2L, function(nu) {
    if(length(lambda) < length(nu[nu>0L])) return(0L)
    Jnu <- JackPol(n, nu, alpha, which)
    coeff <- HallInnerProduct(Jlambda, Jmu * Jnu, alpha) /
      HallInnerProduct(Jnu, Jnu, alpha)
    coeff * Jnu
  }, simplify = FALSE)
  Reduce(`+`, terms)
}
