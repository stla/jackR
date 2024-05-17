test_that("Kostka numbers for alpha=1", {
  expect_identical(KostkaNumbers(4), KostkaNumbers(4, alpha = 1))
})

test_that("Kostka numbers are the coefficients of Jack P-polynomials", {
  n <- 4L
  alpha <- "3/2"
  Knumbers <- KostkaNumbers(4L, alpha = alpha)
  lambda <- c(3L, 1L)
  jp <- JackPol(n, lambda, alpha, which = "P")
  lambda <- paste0("(", paste0(lambda, collapse=","), ")")
  qspray <- qzero()
  for(mu in colnames(Knumbers)) {
    coeff <- Knumbers[lambda, mu]
    mu <- fromString(gsub("(\\(|\\))", "", mu))
    qspray <- qspray + coeff * MSFpoly(n, mu)
  }
  expect_true(qspray == jp)
})
