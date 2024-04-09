test_that("JackSymPol J", {
  n <- 3
  lambda <- c(3, 1, 1)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda)
  JackPolynomial <- JackPol(n, lambda, alpha)
  x <- evalSymbolicQspray(symbolicJackPolynomial, a = alpha)
  expect_true(x == JackPolynomial)
})

test_that("JackSymPol P", {
  n <- 3
  lambda <- c(3, 2, 1)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda, which = "P")
  JackPolynomial <- JackPol(n, lambda, alpha, which = "P")
  x <- evalSymbolicQspray(symbolicJackPolynomial, a = alpha)
  expect_true(x == JackPolynomial)
})

test_that("JackSymPol Q", {
  n <- 3
  lambda <- c(4, 2, 2)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda, which = "Q")
  JackPolynomial <- JackPol(n, lambda, alpha, which = "Q")
  x <- evalSymbolicQspray(symbolicJackPolynomial, a = alpha)
  expect_true(x == JackPolynomial)
})

test_that("JackSymPol is symmetric", {
  n <- 3
  lambda <- c(4, 2, 2)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda)
  expect_silent(MSPcombination(symbolicJackPolynomial, check = TRUE))
})
