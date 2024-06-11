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

#   both concat (unzip (map (swap . (codedRatio (lambda, lambda') (mu, mu'))) ss))
