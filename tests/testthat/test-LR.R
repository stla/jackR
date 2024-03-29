test_that("Littlewood-Richardson multiplication", {
  mu <- c(2, 1)
  nu <- c(3, 2, 1)
  LR <- LRmult(mu, nu, output = "list")
  LRcoeffs <- LR$coeff
  LRparts <- LR$lambda
  LRterms <- lapply(1:length(LRcoeffs), function(i) {
    LRcoeffs[i] * SchurPolCPP(3, LRparts[[i]])
  })
  smu_times_snu <- Reduce(`+`, LRterms)
  expect_true(smu_times_snu == SchurPolCPP(3, mu) * SchurPolCPP(3, nu))
})

test_that("Skew Schur polynomial", {
  sspol <- SkewSchurPol(3, lambda = c(3, 2, 1), mu = c(1, 1))
  expected <- SchurPolCPP(3, c(2, 1, 1)) + SchurPolCPP(3, c(2, 2)) +
    SchurPolCPP(3, c(3, 1))
  expect_true(sspol == expected)
})
