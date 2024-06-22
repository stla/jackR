# clambdamu :: (Eq a, AlgField.C a) => Seq Int -> Seq Int -> RatioOfSprays a
# clambdamu lambda mu = num %//% den
#   where
#     lambda' = _dualPartition' lambda
#     mu' = _dualPartition' mu
.als <- function(lambda, lambdap) {
  do.call(
    rbind,
    apply(
      cbind(lambda, seq_along(lambda)),
      1L,
      function(mi) {
        m <- mi[1L]
        i <- mi[2L]
        t(apply(
          cbind(head(lambdap, m), seq_len(m)),
          1L,
          function(mpj) {
            mp <- mpj[1L]
            j <- mpj[2L]
            c(m - j, mp - i)
          }
        ))
      }, simplify = FALSE
    )
  )
}
#     rg = S.fromList [1 .. max (S.length lambda) (lambda `S.index` 0)]
#     als_lambda =
#        foldl'
#           (
#             \sq (m, i) ->
#               sq >< fmap (\(m', j) -> (m - j, m'- i)) (S.zip lambda' (S.take m rg))
#           )
#            S.empty (S.zip lambda rg)
#     als_mu =
#        foldl'
#           (
#             \sq (m, i) ->
#               sq >< fmap (\(m', j) -> (m - j, m'- i)) (S.zip mu' (S.take m rg))
#           )
#            S.empty (S.zip mu rg)
#     als =
#       (
#         als_lambda, als_mu
#       )
#     (num_map, den_map) =
#       both (foldl' (\i al -> DM.insertWith (+) al 1 i) DM.empty) als
#     f k1 k2 = if k1 > k2 then Just (k1 - k2) else Nothing
#     assocs = both DM.assocs
#       (
#         DM.differenceWith f num_map den_map
#       , DM.differenceWith f den_map num_map
#       )
.poly <- function(alc) {
  spray <- new(
    "qspray",
    powers = list(integer(0L), c(alc[1L], alc[2L] + 1L)),
    coeffs = c("1", "-1")
  )
  spray^(alc[3L])
}

#' @importFrom partitions conjugate
#' @importFrom qspray qone
#' @noRd
clambdamu <- function(lambda, mu) {
  lambdap <- conjugate(lambda)
  mup <- conjugate(mu)
  als_lambda <- .als(lambda, lambdap)
  als_mu <- .als(mu, mup)
  matrices <- simplifyTheTwoMatrices(als_lambda, als_mu)
  matrix1 <- matrices[[1L]]
  if(nrow(matrix1) >= 1L) {
    num <-
      Reduce(
        `*`,
        apply(matrix1, 1L, .poly, simplify = FALSE)
      )
  } else {
    num <- qone()
  }
  matrix2 <- matrices[[2L]]
  if(nrow(matrix2) >= 1L) {
    den <-
      Reduce(
        `*`,
        apply(matrix2, 1L, .poly, simplify = FALSE)
      )
  } else {
    den <- qone()
  }
  num / den
}

#     poly ((a, l), c) =
#       (HM.fromList
#         [
#           (Powers S.empty 0, AlgRing.one)
#         , (Powers (S.fromList [a, l+1]) 2, AlgAdd.negate AlgRing.one)
#         ]) ^**^ c -- (unitSpray ^-^ q a ^*^ t (l+1)) ^**^ c
#     (num, den) = both (productOfSprays . (map poly)) assocs
