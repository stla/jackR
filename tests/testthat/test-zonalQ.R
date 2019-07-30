test_that(
  "Quaternionic zonal polynomials sum to the trace - gmp", {
    x <- as.bigq(c(1L,2L,4L,7L), c(2L,3L,1L,2L))
    expected <- sum(x)^3
    obtained_DK <- ZonalQ(x, 3) + ZonalQ(x, c(2,1)) + ZonalQ(x, c(1,1,1))
    obtained_naive <- ZonalQ(x, 3, algorithm = "naive") +
      ZonalQ(x, c(2,1), algorithm = "naive") +
      ZonalQ(x, c(1,1,1), algorithm = "naive")
    expect_identical(obtained_DK, expected)
    expect_identical(obtained_naive, expected)
  }
)

test_that(
  "Quaternionic zonal polynomials sum to the trace - numeric", {
    x <- c(1,2,4,7) / c(2,3,1,2)
    expected <- sum(x)^3
    obtained_DK <- ZonalQ(x, 3) + ZonalQ(x, c(2,1)) + ZonalQ(x, c(1,1,1))
    obtained_naive <- ZonalQ(x, 3, algorithm = "naive") +
      ZonalQ(x, c(2,1), algorithm = "naive") +
      ZonalQ(x, c(1,1,1), algorithm = "naive")
    expect_equal(obtained_DK, expected)
    expect_equal(obtained_naive, expected)
  }
)
