test_that("Factorial Schur polynomial with a=(0, ...) is Schur polynomial", {
  n <- 4
  lambda <- c(3, 3, 2, 2)
  a <- rep(0, n + lambda[1] -1)
  factorialSchurPoly <- factorialSchurPol(n, lambda, a)
  schurPoly <- SchurPol(n, lambda)
  expect_true(factorialSchurPoly == schurPoly)
})

#   testCase "Factorial Schur polynomial with y=[0 .. ] is Schur polynomial" $ do
#     let
#       n = 4
#       lambda = [3, 3, 2, 2]
#       y = replicate (n + lambda !! 0 - 1) 0
#       factorialSchurPoly = factorialSchurPol' n lambda y
#       schurPoly = schurPol' n lambda
#     assertEqual "" schurPoly factorialSchurPoly
#
