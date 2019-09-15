test_that(
  "ZonalQ with lambda = (4)", {
    # gmp
    obtained <- ZonalQPol(4, c(4), algo="naive", basis="MSF")
    expected <- "M_(4) + 8/5 M_(3,1) + 9/5 M_(2,2) + 12/5 M_(2,1,1) + 16/5 M_(1,1,1,1)"
    expect_identical(obtained, expected)
    # numeric
    obtained <- ZonalQPol(4, c(4), algo="naive", basis="MSF", exact = FALSE)
    expected <- "M_(4) + 1.6 M_(3,1) + 1.8 M_(2,2) + 2.4 M_(2,1,1) + 3.2 M_(1,1,1,1)"
    expect_identical(obtained, expected)
  }
)

test_that(
  "An example of the quaternionic zonal polynomial", {
    poly <- ZonalQPol(3, c(2,1), algo = "naive", basis = "MSF")
    expect_identical(poly, "3/2 M_(2,1) + 18/5 M_(1,1,1)")
  }
)

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

test_that(
  "ZonalQPol is correct", {
    bigqMonomial <- function(vars, powers){
      do.call(prod, mapply(gmp::pow.bigq, vars, powers, SIMPLIFY = FALSE))
    }
    evalPol <- function(pol, x){
      vars <- paste0("x_", seq_along(x))
      polValue <- pol
      for(i in 1L:length(x)){
        polValue <- mvp::subsmvp(polValue, vars[i], mvp::mvp(x[i],1,1))
      }
      polValue$names <- lapply(polValue$names, as.bigq)
      terms <-
        mapply(bigqMonomial, polValue$names, polValue$power, SIMPLIFY = FALSE)
      value <-
        do.call(sum, mapply("*", polValue$coeffs, terms, SIMPLIFY = FALSE))
      value
    }
    #
    lambda <- c(3,2)
    pol <- ZonalQPol(4, lambda, algorithm = "naive")
    x <- as.character(as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L)))
    polEval <- evalPol(pol, x)
    expect_identical(polEval, ZonalQ(as.bigq(x), lambda))
  }
)

test_that(
  "ZonalQ polynomials sum to the trace - polynomial", {
    n <- 4
    expected <- (mvp("x_1",1,1)+mvp("x_2",1,1)+mvp("x_3",1,1)+mvp("x_4",1,1))^3
    obtained <- ZonalQPol(n, 3) + ZonalQPol(n, c(2,1)) + ZonalQPol(n, c(1,1,1))
    expect_identical(expected$names, obtained$names)
    expect_identical(expected$power, obtained$power)
    expect_equal(expected$coeffs, obtained$coeffs)
  }
)
