test_that("JackSymPol J", {
  n <- 3
  lambda <- c(3, 1, 1)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda)
  JackPolynomial <- JackPol(n, lambda, alpha)
  x <- substituteParameters(symbolicJackPolynomial, alpha)
  expect_true(x == JackPolynomial)
})

test_that("JackSymPol P", {
  n <- 3
  lambda <- c(3, 2, 1)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda, which = "P")
  JackPolynomial <- JackPol(n, lambda, alpha, which = "P")
  x <- substituteParameters(symbolicJackPolynomial, alpha)
  expect_true(x == JackPolynomial)
})

test_that("JackSymPol Q", {
  n <- 3
  lambda <- c(4, 2, 2)
  alpha <- gmp::as.bigq("2/3")
  symbolicJackPolynomial <- JackSymPol(n, lambda, which = "Q")
  JackPolynomial <- JackPol(n, lambda, alpha, which = "Q")
  x <- substituteParameters(symbolicJackPolynomial, alpha)
  expect_true(x == JackPolynomial)
})

test_that("JackSymPol is symmetric", {
  n <- 3
  lambda <- c(4, 2, 2)
  symbolicJackPolynomial <- JackSymPol(n, lambda)
  expect_true(isSymmetricQspray(symbolicJackPolynomial))
})

test_that("JackSymPol J has polynomial coefficients only", {
  n <- 5
  lambda <- c(4, 4, 3, 2, 1)
  J <- JackSymPol(n, lambda)
  expect_true(hasPolynomialCoefficients(J))
})
