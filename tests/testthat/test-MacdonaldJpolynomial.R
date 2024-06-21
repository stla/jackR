#   , testCase "Macdonald J-polynomial" $ do
#     let
#       n = 3
#       macJpoly = macdonaldJpolynomial' n [2, 1]
#       q = qlone 1
#       t = qlone 2
#       mus = [[3], [2, 1], [1, 1, 1]]
#       qtKFpolys = [t, unitSpray ^+^ q ^*^ t, q]
#       tSchurPolys = map ((HM.map (flip changeVariables [t])) . (tSchurPolynomial' n)) mus
#       expected = sumOfSprays (zipWith (*^) qtKFpolys tSchurPolys)
#     assertEqual "" macJpoly expected
#
test_that("Macdonald J-polynomial", {
  n <- 3
  macJpoly <- MacdonaldPol(n, c(2, 1), "J")
  q <- qlone(1)
  t <- qlone(2)
  mus <- list(3, c(2, 1), c(1, 1, 1))
  qtKFpolys <- list(t, 1 + q*t, q)
  tSchurPolys <- lapply(1:3, function(i) {
    qtKFpolys[[i]] * changeParameters(tSchurPol(n, mus[[i]]), list(t))
  })
  expected <- Reduce(`+`, tSchurPolys)
  expect_true(macJpoly == expected)
})
#   , testCase "Skew Macdonald J-polynomial at q=0" $ do
#     let
#       n = 3
#       lambda = [3, 2]
#       mu = [2, 1]
#       skewHLpoly = skewHallLittlewoodPolynomial' n lambda mu 'Q'
#       expected = asSimpleParametricSpray $
#         HM.map ((swapVariables (1, 2)) . (substitute [Just 0, Nothing]))
#           (skewMacdonaldJpolynomial' n lambda mu)
#     assertEqual "" skewHLpoly expected
#
#   , testCase "Macdonald polynomial branching rule" $ do
#     let
#       nx = 2
#       ny = 2
#       lambda = [2, 2]
#       ys = [lone 3, lone 4]
#       macJpoly = macdonaldJpolynomial' (nx + ny) lambda
#       expected = asSimpleParametricSpray $
#         sumOfSprays
#           [
#             skewMacdonaldJpolynomial' nx lambda mu
#               ^*^ changeVariables (HM.map asRatioOfSprays $ macdonaldJpolynomial' ny mu) ys
#           | mu <- [[], [1], [2], [1,1], [2,1], [2,2]]
#           ]
#     assertEqual "" macJpoly expected
