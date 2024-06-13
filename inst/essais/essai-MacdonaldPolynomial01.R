# codedRatio ::
#   PartitionsPair -> PartitionsPair -> (Int, Int) -> ([(Int,Int)], [(Int,Int)])
# codedRatio (lambda, lambda') (mu, mu') (i, j)
#   | i <= ellMu && j <= mu_im1 =
#       ([(a+1, l), (a', l'+1)], [(a, l+1), (a'+1, l)])
#   | j <= lambda_im1 =
#       ([(a', l'+1)], [(a'+1, l')])
#   | otherwise =
#       ([], [])
#     where
#       ellMu = S.length mu
#       mu_im1 = mu `S.index` (i-1)
#       a = mu_im1 - j
#       l = mu' `S.index` (j-1) - i
#       lambda_im1 = lambda `S.index` (i-1)
#       a' = lambda_im1 - j
#       l' = lambda' `S.index` (j-1) - i

codedRatio <- function(
  lambda, lambdap, mu, mup, ij
) {
  i <- ij[1L]
  j <- ij[2L]
  ellMu <- length(mu)
  if(i <= ellMu && j <= (mu_i <- mu[i])) {
    a <- mu_i - j
    l <- mup[j] - i
    ap <- lambda[i] - j
    lp <- lambdap[j] - i
    list(
      rbind(
        c(a + 1L, l),
        c(ap, lp + 1L)
      ),
      rbind(
        c(a, l + 1L),
        c(ap + 1L, lp)
      )
    )
  } else if(j <= (lambda_i <- lambda[i])) {
    ap <- lambda_i - j
    lp <- lambdap[j] - i
    list(
      rbind(
        c(ap, lp + 1L)
      ),
      rbind(
        c(ap + 1L, lp)
      )
    )
  } else {
    list(
      matrix(NA_integer_, nrow = 0L, ncol = 2L),
      matrix(NA_integer_, nrow = 0L, ncol = 2L)
    )
  }
}

# psiLambdaMu (lambda, mu) =
#   both concat (unzip (map (swap . (codedRatio (lambda, lambda') (mu, mu'))) ss))
#   where
#     lambda' = _dualPartition' lambda
#     mu' = _dualPartition' mu
#     bools = S.zipWith (==) lambda mu >< S.replicate (S.length lambda - S.length mu) False
#     nonEmptyRows = S.elemIndicesL False bools
#     bools' = S.zipWith (==) lambda' mu'
#     emptyColumns = S.elemIndicesL True bools'
#     ss = [(i+1, j+1) | i <- nonEmptyRows, j <- emptyColumns]

#' @importFrom partitions conjugate
#' @noRd
psiLambdaMu <- function(lambda, mu) {
  lambdap <- conjugate(lambda)
  mup <- conjugate(mu)
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  nonEmptyRows <-
    c(which(lambda[seq_len(ellMu)] != mu), .rg(ellMu + 1L, ellLambda))
  emptyColumns <- which(lambdap[seq_along(mup)] == mup)
  ss <- do.call(
    rbind,
    lapply(nonEmptyRows, function(i) {
      columns <- intersect(emptyColumns, seq_len(lambda[i]))
      cbind(rep(i, length(columns)), columns)
    })
  )
  if(nrow(ss) >= 1L) {
    codedRatios <- apply(ss, 1L, function(ij) {
      codedRatio(lambda, lambdap, mu, mup, ij)
    }, simplify = FALSE)
    list(
      do.call(
        rbind,
        lapply(codedRatios, `[[`, 2L)
      ),
      do.call(
        rbind,
        lapply(codedRatios, `[[`, 1L)
      )
    )
  } else {
    list(
      matrix(NA_integer_, nrow = 0L, ncol = 2L),
      matrix(NA_integer_, nrow = 0L, ncol = 2L)
    )
  }
}


# gtPatternDiagonals' :: GT -> [Seq Int]
# gtPatternDiagonals' pattern = S.empty : [diagonal j | j <- [0 .. l]]
#   where
#     dropTrailingZeros = S.dropWhileR (== 0)
#     l = length pattern - 1
#     diagonal j =
#       dropTrailingZeros
#         (S.fromList
#           [pattern !! r !! c | (r, c) <- zip [l-j .. l] [0 .. j]])
gtPatternDiagonals <- function(pattern) {
  ell <- length(pattern)
  c(list(integer(0L)), lapply(seq_len(ell), function(j) {
    indices <- cbind(jack:::.rg(ell-j+1L, ell), jack:::.rg(1L, j))
    jack:::removeTrailingZeros(do.call(c, apply(indices, 1L, function(rc) {
      pattern[[rc[1L]]][rc[2L]]
    }, simplify = FALSE)))
  }))
}

simplifyTheTwoMatrices <- function(matrix1, matrix2) {
  pairs1 <- apply(matrix1, 1L, toString)
  pairs2 <- apply(matrix2, 1L, toString)
  allPairs <- union(pairs1, pairs2)
  table1 <- table(factor(pairs1, levels = allPairs))
  table2 <- table(factor(pairs2, levels = allPairs))
  diffs <- table1 - table2
  pairsNumerator <- Filter(function(count) count > 0L, diffs)
  pairsDenominator <- Filter(function(count) count > 0L, -diffs)
  rownames(matrix1) <- pairs1
  rownames(matrix2) <- pairs2
  colnames(matrix1) <- colnames(matrix2) <- c("i", "j")
  list(
    cbind(
      matrix1[names(pairsNumerator), , drop = FALSE],
      count = pairsNumerator,
      deparse.level = 0L
    ),
    cbind(
      matrix2[names(pairsDenominator), , drop = FALSE],
      count = pairsDenominator,
      deparse.level = 0L
    )
  )
}

matrix1 <- rbind(
  c(1, 2),
  c(1, 2),
  c(1, 2),
  c(1, 2),
  c(2, 2),
  c(3, 4)
)
matrix2 <- rbind(
  c(1, 2),
  c(1, 2),
  c(2, 2),
  c(2, 2),
  c(2, 2),
  c(4, 5)
)

pairing <- function(lambdas) {
  mapply(
    function(lambda1, lambda2) {
      list(lambda1, lambda2)
    },
    tail(lambdas, -1L), head(lambdas, -1L),
    USE.NAMES = FALSE, SIMPLIFY = FALSE
  )
}

#' @importFrom qspray qone qlone
#' @noRd
makeRatioOfSprays <- function(pairsMap, pairs) {
  pairsOfMatrices <- pairsMap[vapply(pairs, toString, character(1L))]
  matrix1 <- do.call(
    rbind,
    lapply(pairsOfMatrices, `[[`, 1L)
  )
  matrix2 <- do.call(
    rbind,
    lapply(pairsOfMatrices, `[[`, 2L)
  )
  simplifiedMatrices <- simplifyTheTwoMatrices(matrix1, matrix2)
  q <- qlone(1L)
  t <- qlone(2L)
  unitQSpray <- qone()
  num <- Reduce(
    `*`,
    apply(simplifiedMatrices[[1L]], 1L, function(alc) {
      (unitQSpray - q^(alc[1L]) * t^(alc[2L]))^alc[3L]
    }, simplify = FALSE)
  )
  den <- Reduce(
    `*`,
    apply(simplifiedMatrices[[2L]], 1L, function(alc) {
      (unitQSpray - q^(alc[1L]) * t^(alc[2L]))^alc[3L]
    }, simplify = FALSE)
  )
  num / den
}

# makeRatioOfSprays ::
#   (Eq a, AlgField.C a) =>
#   PairsMap -> [PartitionsPair] -> RatioOfSprays a
# makeRatioOfSprays pairsMap pairs = num %//% den
#   where
#     als = both concat (unzip (map ((DM.!) pairsMap) pairs))
#     (num_map, den_map) =
#       both (foldl' (\m k -> DM.insertWith (+) k 1 m) DM.empty) als
#     f k1 k2 = if k1 > k2 then Just (k1 - k2) else Nothing
#     assocs = both DM.assocs
#       (
#         DM.differenceWith f num_map den_map
#       , DM.differenceWith f den_map num_map
#       )
#     -- (num_als, den_als) = both concat (unzip (map ((DM.!) pairsMap) pairs))
#     -- als = (num_als \\ den_als, den_als \\ num_als)
#     -- assocs =
#     --   both (DM.assocs . (foldl' (\m k -> DM.insertWith (+) k 1 m) DM.empty)) als
#     q = lone' 1
#     t = lone' 2
#     poly ((a, l), c) = (unitSpray ^-^ q a ^*^ t l) ^**^ c
#     (num, den) = both (productOfSprays . (map poly)) assocs
#' @importFrom DescTools Permn
#' @importFrom methods new
#' @noRd
.MacdonaldPolynomial <- function(f, n, lambda) {
  mus <- Filter(
    function(mu) length(mu) <= n,
    apply(
      jack:::dominatedPartitions(lambda), 2L, jack:::removeTrailingZeros, simplify = FALSE
    )
  )
  listsOfPairs <- lapply(mus, function(mu) {
    lapply(syt::GelfandTsetlinPatterns(lambda, mu), function(pattern) {
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
    f(pair[[1L]], pair[[2L]])
  })
  names(pairsMap) <- vapply(allPairs, toString, character(1L))
  QSprays <- lapply(seq_along(mus), function(i) {
    mu <- mus[[i]]
    listOfPairs <- listsOfPairs[[i]]
    rOQ <- Reduce(`+`, lapply(listOfPairs, function(pairs) {
      makeRatioOfSprays(pairsMap, pairs)
    }))
    compos <- DescTools::Permn(c(mu, rep(0L, n - length(mu))))
    powers <- apply(compos, 1L, jack:::removeTrailingZeros, simplify = FALSE)
    list(
      "powers" = powers,
      "coeffs" = rep(list(rOQ), length(powers))
    )
  })
  new(
    "symbolicQSpray",
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
#     dropTrailingZeros = S.dropWhileR (== 0)
#     hashMaps =
#       map
#         (\mu ->
#           let mu' = fromPartition mu
#               mu'' = S.fromList mu'
#               mu''' = mu' ++ (replicate (n - S.length mu'') 0)
#               coeff = coeffs HM.! mu''
#               compos = permuteMultiset mu'''
#           in
#             HM.fromList
#               [let compo' = dropTrailingZeros (S.fromList compo) in
#                 (Powers compo' (S.length compo'), coeff) | compo <- compos]
#         ) mus

# _macdonaldPolynomial f n lambda = HM.unions hashMaps
#   where
#     lambda' = toPartitionUnsafe lambda
#     mus = filter (\mu -> partitionWidth mu <= n) (dominatedPartitions lambda')
#     pairing lambdas = zip (drop1 lambdas) lambdas
#     listsOfPairs =
#       map (
#         map (pairing . gtPatternDiagonals')
#           . (kostkaGelfandTsetlinPatterns lambda')
#       ) mus
#     allPairs = nub $ concat (concat listsOfPairs)
#     pairsMap = DM.fromList (zip allPairs (map f allPairs))
#     coeffs = HM.fromList
#       (zipWith
#         (\mu listOfPairs ->
#           (
#             S.fromList (fromPartition mu)
#           , AlgAdd.sum (map (makeRatioOfSprays pairsMap) listOfPairs)
#           )
#         ) mus listsOfPairs
#       )
#     dropTrailingZeros = S.dropWhileR (== 0)
#     hashMaps =
#       map
#         (\mu ->
#           let mu' = fromPartition mu
#               mu'' = S.fromList mu'
#               mu''' = mu' ++ (replicate (n - S.length mu'') 0)
#               coeff = coeffs HM.! mu''
#               compos = permuteMultiset mu'''
#           in
#             HM.fromList
#               [let compo' = dropTrailingZeros (S.fromList compo) in
#                 (Powers compo' (S.length compo'), coeff) | compo <- compos]
#         ) mus
#
