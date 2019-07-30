test_that(
  "Schur (3,2) - gmp", {
    x <- as.bigq(3L:5L, c(10L,2L,1L))
    expected <- x[1]^3*x[2]^2 + x[1]^3*x[3]^2 + x[1]^3*x[2]*x[3] +
      x[1]^2*x[2]^3 + x[1]^2*x[3]^3 + 2*x[1]^2*x[2]*x[3]^2 +
      2*x[1]^2*x[2]^2*x[3] + x[1]*x[2]*x[3]^3 + 2*x[1]*x[2]^2*x[3]^2 +
      x[1]*x[2]^3*x[3] + x[2]^2*x[3]^3 + x[2]^3*x[3]^2
    naive <- Schur(x, c(3,2), algorithm = "naive")
    DK <- Schur(x, c(3,2), algorithm = "DK")
    expect_identical(naive, expected)
    expect_identical(DK, expected)
  }
)

test_that(
  "Schur (3,2) - numeric", {
    x <- c(3L:5L) / c(10L,2L,1L)
    expected <- x[1]^3*x[2]^2 + x[1]^3*x[3]^2 + x[1]^3*x[2]*x[3] +
      x[1]^2*x[2]^3 + x[1]^2*x[3]^3 + 2*x[1]^2*x[2]*x[3]^2 +
      2*x[1]^2*x[2]^2*x[3] + x[1]*x[2]*x[3]^3 + 2*x[1]*x[2]^2*x[3]^2 +
      x[1]*x[2]^3*x[3] + x[2]^2*x[3]^3 + x[2]^3*x[3]^2
    naive <- Schur(x, c(3,2), algorithm = "naive")
    DK <- Schur(x, c(3,2), algorithm = "DK")
    expect_equal(naive, expected)
    expect_equal(DK, expected)
  }
)

test_that(
  "SchurPol is correct", {
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
    pol <- SchurPol(4, lambda)
    x <- as.character(as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L)))
    polEval <- evalPol(pol, x)
    expect_identical(polEval, Schur(as.bigq(x), lambda))
  }
)
