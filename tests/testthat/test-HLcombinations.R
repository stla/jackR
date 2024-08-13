test_that("HLcombinationP", {
  lambda <- c(2, 2, 1)
  alpha <- "2/5"
  which <- "C"
  n <- sum(lambda)
  Jpoly <- JackPol(n, lambda, alpha, which)
  combo <- HLcombinationP(Jpoly)
  toAdd <- lapply(combo, function(lst) {
    mu <- lst[["lambda"]]
    coeff <- lst[["coeff"]]
    coeff * HallLittlewoodPol(n, mu, "P")
  })
  obtained <- Reduce(`+`, toAdd)
  Jpoly2 <- as(Jpoly, "symbolicQspray")
  expect_true(obtained == Jpoly2)
})
