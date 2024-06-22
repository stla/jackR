#   , testCase "qt-Kostka polynomials" $ do
#     let
#       n = 4
#       mu = [2, 1, 1]
#       macJpoly = macdonaldJpolynomial' n mu
#       qtKostkaPolys = qtKostkaPolynomials' mu
#       expected =
#         sumOfSprays $ DM.elems $ DM.mapWithKey
#           (\lambda kp ->
#             kp *^ (HM.map (swapVariables (1, 2)) (tSchurPolynomial' n lambda)))
#         qtKostkaPolys
#     assertEqual "" macJpoly expected
test_that("qt-Kostka polynomials", {
  n <- 4
  mu <- c(2, 1, 1)
  macJpoly <- MacdonaldPol(n, mu, "J")
  qtKostkaPolys <- qtKostkaPolynomials(mu)
  t <- qlone(2)
  expected <- Reduce(
    `+`,
    lapply(qtKostkaPolys, function(lambda_kp) {
      lambda <- lambda_kp[["lambda"]]
      kp <- lambda_kp[["polynomial"]]
      tSchurPoly <- tSchurPol(n, lambda)
      kp * changeParameters(tSchurPoly, list(t))
    })
  )
  expect_true(macJpoly == expected)
})

test_that("Skew qt-Kostka polynomials", {
  lambda <- c(2, 1, 1)
  mu <- c(1, 1)
  qtSkewKostkaPolys <- qtSkewKostkaPolynomials(lambda, mu)
  q <- qlone(1)
  t <- qlone(2)
  expected <- lapply(qtSkewKostkaPolys, function(nu_poly) {
    changeParameters(
      changeParameters(
        nu_poly[["polynomial"]], list(qzero(), t)
      ),
      list(t, q)
    )
  })
  skewKFpolys <- lapply(qtSkewKostkaPolys, function(nu_poly) {
    SkewKostkaFoulkesPolynomial(lambda, mu, nu_poly[["nu"]])
  })
  check <- mapply(
    `==`,
    expected, skewKFpolys,
    SIMPLIFY = TRUE
  )
  expect_true(all(check))
})
#   , testCase "Skew qt-Kostka polynomials" $ do
#     let
#       lambda = [2, 1, 1]
#       mu = [1, 1]
#       qtSkewKostkaPolys = qtSkewKostkaPolynomials' lambda mu
#       expected =
#         map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing]))
#               (DM.elems qtSkewKostkaPolys)
#       skewKFpolys =
#         map (skewKostkaFoulkesPolynomial' lambda mu)
#               (DM.keys qtSkewKostkaPolys)
#     assertEqual "" skewKFpolys expected
