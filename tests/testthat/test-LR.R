test_that("Littlewood-Richardson", {
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
