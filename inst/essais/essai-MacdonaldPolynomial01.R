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
  lambda, lambdap, mu, mup, i, j
) {
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

