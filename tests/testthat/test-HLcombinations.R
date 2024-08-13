test_that("HLcombinationP - Jack polynomial", {
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

test_that("HLcombinationP - Jack symbolic polynomial", {
  lambda <- c(2, 2, 1)
  which <- "Q"
  n <- sum(lambda)
  JsymPoly <- JackSymPol(n, lambda, which)
  combo <- HLcombinationP(JsymPoly)
  toAdd <- lapply(combo, function(lst) {
    mu <- lst[["lambda"]]
    coeff <- lst[["coeff"]]
    coeff * HallLittlewoodPol(n, mu, "P")
  })
  obtained <- Reduce(`+`, toAdd)
  expect_true(obtained == JsymPoly)
})


test_that("HLcombinationQ - Jack polynomial", {
  lambda <- c(3, 1, 1)
  alpha <- "22/5"
  which <- "C"
  n <- sum(lambda)
  Jpoly <- JackPol(n, lambda, alpha, which)
  combo <- HLcombinationQ(Jpoly)
  toAdd <- lapply(combo, function(lst) {
    mu <- lst[["lambda"]]
    coeff <- lst[["coeff"]]
    coeff * HallLittlewoodPol(n, mu, "Q")
  })
  obtained <- Reduce(`+`, toAdd)
  Jpoly2 <- as(Jpoly, "symbolicQspray")
  expect_true(obtained == Jpoly2)
})

test_that("HLcombinationQ - Jack symbolic polynomial", {
  lambda <- c(3, 1, 1)
  which <- "J"
  n <- sum(lambda)
  JsymPoly <- JackSymPol(n, lambda, which)
  combo <- HLcombinationQ(JsymPoly)
  toAdd <- lapply(combo, function(lst) {
    mu <- lst[["lambda"]]
    coeff <- lst[["coeff"]]
    coeff * HallLittlewoodPol(n, mu, "Q")
  })
  obtained <- Reduce(`+`, toAdd)
  expect_true(obtained == JsymPoly)
})
