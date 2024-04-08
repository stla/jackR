test_that("JackSymPol", {
  n <- 3
  lambda <- c(3, 1, 1)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda)
  JackPolynomial <- JackPol(n, lambda, alpha)
  x <- evalSymbolicQspray(symbolicJackPolynomial, a = alpha)
  expect_true(x == JackPolynomial)
})
