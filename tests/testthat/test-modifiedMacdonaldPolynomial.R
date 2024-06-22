# testCase "Modified Macdonald polynomial mu (q, t) = Modified Macdonald polynomial mu' (t, q)" $ do
#     let
#       n = 4
#       mu = [2, 1, 1]
#       mu' = [3, 1]
#       macHpoly = modifiedMacdonaldPolynomial' n mu
#       macHpoly' = HM.map (swapVariables (1, 2)) (modifiedMacdonaldPolynomial' n mu')
#     assertEqual "" macHpoly macHpoly'
test_that("Modified Macdonald polynomial mu (q, t) = Modified Macdonald polynomial mu' (t, q)", {
  n <- 4
  mu <- c(2, 1, 1)
  mup <- c(3, 1)
  macHpoly <- modifiedMacdonaldPolynomial(n, mu)
  macHpolyp <- changeParameters(
    modifiedMacdonaldPolynomial(n, mup),
    list(qlone(2), qlone(1))
  )
  expect_true(macHpoly == macHpolyp)
})
