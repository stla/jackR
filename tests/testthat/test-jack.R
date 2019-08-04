test_that(
  "Jack (3,1) - gmp", {
    alpha <- as.bigq(5L,2L)
    x <- as.bigq(2L:5L)
    expected <- (2*alpha^2+4*alpha+2)*MSF(x,c(3,1)) +
      (6*alpha+10)*MSF(x,c(2,1,1)) + (4*alpha+4)*MSF(x,c(2,2)) +
      24*MSF(x,c(1,1,1,1))
    jack_naive <- Jack(x, c(3,1), alpha, algorithm = "naive")
    jack_DK <- Jack(x, c(3,1), alpha, algorithm = "DK")
    expect_identical(jack_naive, expected)
    expect_identical(jack_DK, expected)
  }
)

test_that(
  "Jack (3,1) - numeric", {
    alpha <- 5/2
    x <- 2L:5L
    expected <- (2*alpha^2+4*alpha+2)*MSF(x,c(3,1)) +
      (6*alpha+10)*MSF(x,c(2,1,1)) + (4*alpha+4)*MSF(x,c(2,2)) +
      24*MSF(x,c(1,1,1,1))
    jack_naive <- Jack(x, c(3,1), alpha, algorithm = "naive")
    jack_DK <- Jack(x, c(3,1), alpha, algorithm = "DK")
    expect_equal(jack_naive, expected)
    expect_equal(jack_DK, expected)
  }
)

test_that(
  "Jack (4) - gmp", {
    alpha <- as.bigq(5L,2L)
    x <- as.bigq(2L:5L)
    expected <- (1+alpha)*(1+2*alpha)*(1+3*alpha)*MSF(x,c(4)) +
      4*(1+alpha)*(1+2*alpha)*MSF(x,c(3,1)) +
      12*(1+alpha)*MSF(x,c(2,1,1)) +
      6*(1+alpha)^2*MSF(x,c(2,2)) +
      24*MSF(x,c(1,1,1,1))
    jack_naive <- Jack(x, c(4), alpha, algorithm = "naive")
    jack_DK <- Jack(x, c(4), alpha, algorithm = "DK")
    expect_identical(jack_naive, expected)
    expect_identical(jack_DK, expected)
  }
)

test_that(
  "Jack (4) - numeric", {
    alpha <- 5/2
    x <- 2L:5L
    expected <- (1+alpha)*(1+2*alpha)*(1+3*alpha)*MSF(x,c(4)) +
      4*(1+alpha)*(1+2*alpha)*MSF(x,c(3,1)) +
      12*(1+alpha)*MSF(x,c(2,1,1)) +
      6*(1+alpha)^2*MSF(x,c(2,2)) +
      24*MSF(x,c(1,1,1,1))
    jack_naive <- Jack(x, c(4), alpha, algorithm = "naive")
    jack_DK <- Jack(x, c(4), alpha, algorithm = "DK")
    expect_equal(jack_naive, expected)
    expect_equal(jack_DK, expected)
  }
)

test_that(
  "JackPol is correct", {
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
    alpha <- as.bigq(11L,3L)
    pol <- JackPol(4, lambda, alpha, algorithm = "naive")
    x <- as.character(as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L)))
    polEval <- evalPol(pol, x)
    expect_identical(polEval, Jack(as.bigq(x), lambda, alpha))
  }
)
