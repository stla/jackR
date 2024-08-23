test_that("SchurCombination", {
  combo <- SchurCombination(msPolynomial(4L, c(4L)))
  qspray <- comboToQspray(combo)
  expected <- new(
    "qspray",
    powers = list(c(1L,1L,1L,1L), c(2L,1L,1L), c(3L,1L), c(4L)),
    coeffs = c("-1", "1", "-1", "1")
  )
  expect_true(qspray == expected)
  #
  combo <- SchurCombination(psPolynomial(4, c(3, 1)) + esPolynomial(4, c(2, 1)))
  qspray <- comboToQspray(combo)
  expected <- new(
    "qspray",
    powers = list(c(1L,1L,1L), c(1L,1L,1L,1L), c(2L,1L), c(2L,2L), c(4L)),
    coeffs = c("1", "1", "1", "-1", "1")
  )
  expect_true(qspray == expected)
})
