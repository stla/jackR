#     flambda p q r nu =
#       (S.take (p-1) nu |>
#         nu `S.index` (q-1) + p - q + r) ><
#           fmap (+1) (S.take (q-p) (S.drop (p-1) nu)) ><
#             S.drop q nu
.flambda <- function(pq, r, nu) {
  p <- pq[1L]
  q <- pq[2L]
  if(p == 1L) {
    tnu <- nu
  } else {
    tnu <- tail(nu, 1L-p)
  }
  c(head(nu, p-1L), nu[q]+p-q+r, head(tnu, q-p) + 1L, tail(nu, -q))
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
              all(lbda <= lambda[seq_along(lbda)])
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


# chi_lambda_mu_rho :: Seq Int -> Seq Int -> Seq Int -> Int
# chi_lambda_mu_rho lambda mu rho =
#   if S.null rho then 1 else 2 * nevens - length sequences
#   where
#     ribbonHeight :: Seq Int -> Seq Int -> Int
#     ribbonHeight kappa nu =
#       DF.sum
#         (S.zipWith (\k n -> fromEnum (k /= n)) kappa nu)
#           - 1
#       -- kappa and mu have same length so don't need to add S.length kappa - S.length mu
#     sequences = sequencesOfRibbons lambda mu rho
#     nevens =
#       sum $ map
#         (
#           \sq ->
#             (fromEnum . even . DF.sum) $
#               S.zipWith ribbonHeight (S.drop 1 sq) sq
#         )
#           sequences
chi_lambda_mu_rho <- function(lambda, mu, rho) {
  if(length(rho) == 0L) {
    1L
  } else {
    sequences <- sequencesOfRibbons(lambda, mu, rho)
    nevens <- sum(
      vapply(sequences, function(sq) {
        sum(
          vapply(seq_len(length(sq)-1L), function(i) {
            kappa <- sq[[i+1L]]
            nu <- sq[[i]]
            sum(kappa != nu) - 1L
          }, integer(1L))
        ) %% 2L == 0L
      }, logical(1L))
    )
    2L * nevens - length(sequences)
  }
}

zlambda <- function(lambda) {
  parts <- unique(lambda)
  mjs <- vapply(parts, function(j) {
    sum(lambda == j)
  }, integer(1L))
  prod(factorial(mjs) * parts^mjs)
}
# _tSkewSchurPolynomial ::
#   (Eq a, AlgField.C a)
#   => (Integer -> Integer -> a)
#   -> Int
#   -> Partition
#   -> Partition
#   -> SimpleParametricSpray a
.tSkewSchurPolynomial <- function(n, lambda, mu) {
  w <- sum(lambda) - sum(mu)
  rhos <- apply(parts(w), 2L, removeTrailingZeros, simplify = FALSE)
  unitSpray <- qone()
  t <- qlone(1L)
  mapOfSprays <- lapply(seq_len(w), function(r) {
    unitSpray - t^r
  })
  sprays <- lapply(rhos, function(rho) {
    c <- chi_lambda_mu_rho(lambda, mu, rho)
    if(c == 0L) {
      qzero()
    } else {
      psPoly <- PSFpoly(n, rho)
      coeffs <- lapply(psPoly@coeffs, function(coeff) {
        as.ratioOfQsprays(coeff * Reduce(`*`, mapOfSprays[rho]))
      })
      tPowerSumPol <- new(
        "symbolicQspray",
        powers = psPoly@powers,
        coeffs = coeffs
      )
      as.bigq(c, zlambda(rho)) * tPowerSumPol
    }
  })
  Reduce(`+`, sprays)
}
# _tSkewSchurPolynomial f n lambda mu = sumOfSprays sprays
#   where
#     w = sum lambda - sum mu
#     rhos = partitions w
#     t = lone' 1
#     mapOfSprays = IM.fromList (map (\r -> (r, unitSpray ^-^ t r)) [1 .. w])
#     tPowerSumPol rho =
#       HM.map
#         (flip (*^) (productOfSprays (map ((IM.!) mapOfSprays) rho)))
#           (psPolynomial n rho)
#     lambda' = S.fromList lambda
#     mu' = S.fromList mu
#     chi_lambda_mu_rhos =
#       [(rho', chi_lambda_mu_rho lambda' mu' (S.fromList rho'))
#         | rho <- rhos, let rho' = fromPartition rho]
#     sprays =
#       [
#         (f (toInteger c) (toInteger (zlambda rho)))
#          AlgMod.*> tPowerSumPol rho
#       | (rho, c) <- chi_lambda_mu_rhos, c /= 0
#       ]
