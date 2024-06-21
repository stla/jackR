#     flambda p q r nu =
#       (S.take (p-1) nu |>
#         nu `S.index` (q-1) + p - q + r) ><
#           fmap (+1) (S.take (q-p) (S.drop (p-1) nu)) ><
#             S.drop q nu
.flambda <- function(pq, r, nu) {
  p <- pq[1L]
  q <- pq[2L]
  c(head(nu, p-1L), nu[q]+p-q+r, head(tail(nu, 1L-p), q-p) + 1L, tail(nu, -q))
}

#     ok q r nu =
#       let nu_qm1 = nu `S.index` (q-1) in
#         nu_qm1 - q + r > nu `S.index` 0 - 1
#           && nu_qm1 <= lambda `S.index` 0
.ok <- function(lambda, q, r, nu) {
  nu_q <- nu[q]
  nu_q - q + r >= nu[1L] && nu_q <= lambda[1L]
}
#     pairs' r nu =
#       [(p, q) | p <- [2 .. n], q <- [p .. n], ok' p q r nu]
#     ok' p q r nu =
#        let nu_qm1 = nu `S.index` (q-1) in
#         nu_qm1 - q + r > nu `S.index` (p-1) - p
#           && nu `S.index` (p-2) - p >= nu_qm1 - q + r
#           && nu_qm1 <= lambda `S.index` (p-1)
#           && all (uncurry (<))
#                   (S.zip (S.take (q-p) (S.drop (p-1) nu)) (S.drop p lambda))
.okp <- function(lambda, pq, r, nu) {
  p <- pq[1L]
  q <- pq[2L]
  nu_q <- nu[q]
  nu_q - q + r > nu[p] - p &&
    nu[p-1L] - p >= nu_q - q + r &&
      nu_q <= lambda[p] &&
        all(head(tail(nu, 1L-p), q-p) < tail(lambda, -p))
}
# sequencesOfRibbons :: Seq Int -> Seq Int -> Seq Int -> [Seq (Seq Int)]
# sequencesOfRibbons lambda mu rho =
#    foldr
#      (\r zs ->
#       [z |> lbda
#         | z <- zs
#         , lbda <- lambdas r (z `S.index` (S.length z - 1))
#         , and (S.zipWith (<=) lbda lambda)
#       ])
#         [S.singleton (mu >< (S.replicate (n - S.length mu) 0))]
#           rho
#    where
#     n = S.length lambda
#     lambdas r nu = [flambda p q r nu | (p, q) <- pairs r nu ++ pairs' r nu]
sequencesOfRibbons <- function(lambda, mu, rho) {
  f <- function(zs, r) {
    do.call(c, lapply(zs, function(z) {
      lapply(
          Filter(
            function(lbda) {
              all(lbda <= lambda)
            },
            .lambdas(lambda, r, z[[length(z)]])
          ),
          function(lbda) {
            c(z, list(lbda))
          }
      )
    }))
  }
  Reduce(f, rho, init = list(list(c(mu, rep(0L, length(lambda) - length(mu))))))
}
.lambdas <- function(lambda, r, nu) {
  apply(
    rbind(.pairs(lambda, r, nu), .pairsp(lambda, r, nu)),
    1L,
    function(pq) {
      .flambda(pq, r, nu)
    },
    simplify = FALSE
  )
}
#     flambda p q r nu =
#       (S.take (p-1) nu |>
#         nu `S.index` (q-1) + p - q + r) ><
#           fmap (+1) (S.take (q-p) (S.drop (p-1) nu)) ><
#             S.drop q nu
#     pairs r nu = [(1, q) | q <- [1 .. n], ok q r nu]
#     ok q r nu =
#       let nu_qm1 = nu `S.index` (q-1) in
#         nu_qm1 - q + r > nu `S.index` 0 - 1
.pairs <- function(lambda, r, nu) {
  cbind(
    1L,
    Filter(
      function(q) .ok(lambda, q, r, nu),
      seq_along(lambda)
    )
  )
}
#           && nu_qm1 <= lambda `S.index` 0
#     pairs' r nu =
#       [(p, q) | p <- [2 .. n], q <- [p .. n], ok' p q r nu]
.pairsp <- function(lambda, r, nu) {
  n <- length(lambda)
  Grid <- do.call(
    rbind,
    lapply(2L:n, function(p) cbind(p, p:n))
  )
  keep <- apply(Grid, 1L, function(pq) {
    .okp(lambda, pq, r, nu)
  })
  Grid[keep, , drop = FALSE]
}
#     ok' p q r nu =
#        let nu_qm1 = nu `S.index` (q-1) in
#         nu_qm1 - q + r > nu `S.index` (p-1) - p
#           && nu `S.index` (p-2) - p >= nu_qm1 - q + r
#           && nu_qm1 <= lambda `S.index` (p-1)
#           && all (uncurry (<))
#                   (S.zip (S.take (q-p) (S.drop (p-1) nu)) (S.drop p lambda))
