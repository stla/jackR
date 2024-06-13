lastSubpartition <- function(w, lambda) {
  if(length(lambda) == 0L) {
    integer(0L)
  } else {
    k <- lambda[1L]
    if(w <= k) {
      w
    } else {
      c(k, lastSubpartition(w - k, tail(lambda, -1L)))
    }
  }
}


#     nus =
#       filter ((<= n) . partitionWidth) $ dropEnd 1
#         (dominatedPartitions
#           (toPartitionUnsafe (lastSubPartition (sum lambda - sum mu) lambda)))
#     pairing lambdas = zip (drop1 lambdas) lambdas
#     listsOfPairs =
#       map (
#         map pairing
#           . (skewGelfandTsetlinPatterns lambda mu)
#           . fromPartition
#       ) nus
#' @importFrom DescTools Permn
#' @importFrom methods new
#' @importFrom syt skewGelfandTsetlinPatterns
#' @noRd
.SkewMacdonaldPolynomial <- function(f, n, lambda, mu) {
  nus <- Filter(
    function(nu) length(nu) <= n,
    listOfDominatedPartitions(
      lastSubpartition(sum(lambda) - sum(mu), lambda)
    )
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
  pairsMap <- lapply(allPairs, function(pair) {
    f(pair[[1L]], pair[[2L]])
  })
  names(pairsMap) <- vapply(allPairs, toString, character(1L))
  nus <- nus[i_]
  QSprays <- lapply(seq_along(nus), function(i) {
    nu <- nus[[i]]
    listOfPairs <- listsOfPairs[[i]]
    rOQ <- Reduce(`+`, lapply(listOfPairs, function(pairs) {
      makeRatioOfSprays(pairsMap, pairs)
    }))
    compos <- Permn(c(nu, rep(0L, n - length(nu))))
    powers <- apply(compos, 1L, removeTrailingZeros, simplify = FALSE)
    list(
      "powers" = powers,
      "coeffs" = rep(list(rOQ), length(powers))
    )
  })
  new(
    "symbolicQspray",
    powers = do.call(
      c,
      lapply(QSprays, `[[`, "powers")
    ),
    coeffs = do.call(
      c,
      lapply(QSprays, `[[`, "coeffs")
    )
  )
}

SkewMacdonaldPolynomial <- function(n, lambda, mu, which) {
  stopifnot(isPositiveInteger(n))
  stopifnot(isPartition(lambda))
  stopifnot(isPartition(mu))
  stopifnot(which %in% c("P", "Q"))
  lambda <- as.integer(removeTrailingZeros(lambda))
  mu <- as.integer(removeTrailingZeros(mu))
  if(which == "P") {
    .SkewMacdonaldPolynomial(psiLambdaMu, n, lambda, mu)
  } else {
    .SkewMacdonaldPolynomial(phiLambdaMu, n, lambda, mu)
  }
}
