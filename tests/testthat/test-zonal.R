test_that(
  "Zonal = 0 if l(lambda)>l(x)", {
    # numeric
    expect_equal(Zonal(c(1,2), c(3,2,1)), 0)
    expect_equal(Zonal(c(1,2), c(3,2,1), algorithm = "naive"), 0)
    # gmp
    x <- as.bigq(c(1L,2L))
    lambda <- c(3,2,1)
    expect_identical(Zonal(x, lambda), as.bigq(0L))
    expect_identical(Zonal(x, lambda, algorithm = "naive"), as.bigq(0L))
    # polynomial
    n <- 2
    lambda <- c(3,2,1)
    expect_true(ZonalPol(n, lambda) == as.qspray(0))
    expect_identical(ZonalPol(n, lambda, algorithm = "naive"),
                     as.qspray(0))
    expect_identical(ZonalPol(n, lambda, exact = FALSE, algorithm = "naive"),
                     mvp::constant(0))
    expect_identical(ZonalPol(n, lambda, algorithm = "naive", basis = "MSF"),
                     as.qspray(0))
    expect_identical(ZonalPol(n, lambda, exact = FALSE, algorithm = "naive",
                              basis = "MSF"),
                     mvp::constant(0))
  }
)


test_that(
  "Zonal polynomials sum to the trace - gmp", {
    x <- as.bigq(c(1L,2L,4L,7L), c(2L,3L,1L,2L))
    expected <- sum(x)^3
    obtained_DK <- Zonal(x, 3) + Zonal(x, c(2,1)) + Zonal(x, c(1,1,1))
    obtained_naive <- Zonal(x, 3, algorithm = "naive") +
      Zonal(x, c(2,1), algorithm = "naive") +
      Zonal(x, c(1,1,1), algorithm = "naive")
    expect_identical(obtained_DK, expected)
    expect_identical(obtained_naive, expected)
  }
)

test_that(
  "Zonal polynomials sum to the trace - numeric", {
    x <- c(1,2,4,7) / c(2,3,1,2)
    expected <- sum(x)^3
    obtained_DK <- Zonal(x, 3) + Zonal(x, c(2,1)) + Zonal(x, c(1,1,1))
    obtained_naive <- Zonal(x, 3, algorithm = "naive") +
      Zonal(x, c(2,1), algorithm = "naive") +
      Zonal(x, c(1,1,1), algorithm = "naive")
    expect_equal(obtained_DK, expected)
    expect_equal(obtained_naive, expected)
  }
)

test_that(
  "Zonal polynomials sum to the trace - complex", {
    x <- c(1i, 2+1i, 4-3i, 7-3i) / c(2,3,1,2)
    expected <- sum(x)^3
    obtained_DK <- Zonal(x, 3) + Zonal(x, c(2,1)) + Zonal(x, c(1,1,1))
    obtained_naive <- Zonal(x, 3, algorithm = "naive") +
      Zonal(x, c(2,1), algorithm = "naive") +
      Zonal(x, c(1,1,1), algorithm = "naive")
    expect_equal(obtained_DK, expected)
    expect_equal(obtained_naive, expected)
  }
)

test_that(
  "ZonalPol is correct", {
    lambda <- c(3,2)
    pol <- ZonalPol(4, lambda, algorithm = "naive")
    x <- as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L))
    polEval <- qspray::evalQspray(pol, x)
    expect_identical(polEval, Zonal(as.bigq(x), lambda))
  }
)

test_that(
  "Zonal polynomials sum to the trace - polynomial", {
    n <- 4
    expected <- (qlone(1) + qlone(2) + qlone(3) + qlone(4))^3
    obtained <- ZonalPol(n, 3) + ZonalPol(n, c(2,1)) + ZonalPol(n, c(1,1,1))
    expect_true(expected == obtained)
  }
)

test_that(
  "ZonalCPP is correct", {
    x <- as.bigq(c(6L, -7L, 8L, 9L), c(1L, 2L, 3L, 4L))
    lambda <- c(3, 2)
    res <- ZonalCPP(x, lambda)
    expect_identical(res, Zonal(x, lambda))
    #
    x <- c(6, -7, 8, 9) / c(1, 2, 3, 4)
    lambda <- c(3, 2)
    res <- ZonalCPP(x, lambda)
    expect_equal(res, Zonal(x, lambda))
  }
)
