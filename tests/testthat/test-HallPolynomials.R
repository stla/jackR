test_that("A Hall polynomial (comparison with Sage)", {
  hallPolys <- HallPolynomials(c(3, 1, 1), c(2, 1))
  lambda <- c(4, 2, 1, 1)
  hallPoly <- hallPolys[[partitionAsString(lambda)]][["polynomial"]]
  t <- qlone(1)
  expected <- 2*t^3 + t^2 - t - 1
  expect_true(hallPoly == expected)
})

