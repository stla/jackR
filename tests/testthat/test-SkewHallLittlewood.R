test_that("Skew Hall-Littlewood at t=1 is skew Schur", {
  n <- 3; lambda <- c(3, 2, 1); mu <- c(1, 1)
  skewHLpoly <- SkewHallLittlewoodPol(n, lambda, mu)
  skewSchurPoly <- SkewSchurPol(n, lambda, mu)
  expect_true(substituteParameters(skewHLpoly, 0) == skewSchurPoly)
})

test_that("Skew Hall-Littlewood for mu=[] is Hall-Littlewood", {
  n <- 3; lambda <- c(3, 2, 1); mu <- integer(0)
  skewHLpoly <- SkewHallLittlewoodPol(n, lambda, mu, "Q")
  HLpoly <- HallLittlewoodPol(n, lambda, "Q")
  expect_true(skewHLpoly == HLpoly)
})
