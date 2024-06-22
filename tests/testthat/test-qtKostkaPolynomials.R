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
