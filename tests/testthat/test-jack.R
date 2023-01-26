test_that("Jack at x = 0", {
  expect_equal(Jack(c(0,0), NULL, 2), 1)
  expect_equal(Jack(c(0,0,0), c(3,2), 2), 0)
})

test_that(
  "Jack = 0 if l(lambda)>l(x)", {
  # numeric
  expect_equal(Jack(c(1,2), c(3,2,1), 2), 0)
  expect_equal(Jack(c(1,2), c(3,2,1), 0), 0)
  expect_equal(Jack(c(1,2), c(3,2,1), 2, algorithm = "naive"), 0)
  expect_equal(Jack(c(1,2), c(3,2,1), 0, algorithm = "naive"), 0)
  # gmp
  x <- as.bigq(c(1L,2L))
  lambda <- c(3,2,1)
  alpha <- as.bigq(4L)
  expect_identical(Jack(x, lambda, alpha), as.bigq(0L))
  expect_identical(Jack(x, lambda, as.bigq(0L)), as.bigq(0L))
  expect_identical(Jack(x, lambda, alpha, algorithm = "naive"), as.bigq(0L))
  expect_identical(Jack(x, lambda, as.bigq(0L), algorithm = "naive"), as.bigq(0L))
  # polynomial
  n <- 2
  lambda <- c(3,2,1)
  expect_identical(JackPol(n, lambda, alpha = 4), mvp::constant(0))
  expect_identical(JackPol(n, lambda, alpha = 0), mvp::constant(0))
  expect_identical(JackPol(n, lambda, alpha = 4, algorithm = "naive"),
                   mvp::constant(0))
  expect_identical(JackPol(n, lambda, alpha = as.bigq(4L), algorithm = "naive"),
                   as.qspray(0))
  expect_identical(JackPol(n, lambda, alpha = 4, algorithm = "naive",
                           basis = "MSF"),
                   mvp::constant(0))
  expect_identical(JackPol(n, lambda, alpha = as.bigq(4L), algorithm = "naive",
                           basis = "MSF"),
                   as.qspray(0))
  }
)

test_that(
  "Jack - empty partition", {
    # numeric
    expect_equal(Jack(c(1,2), NULL, 2), 1)
    expect_equal(Jack(c(1,2), c(), 2), 1)
    expect_equal(Jack(c(1,2), c(0,0), 2), 1)
    expect_equal(Jack(c(1,2), NULL, 2, algorithm = "naive"), 1)
    expect_equal(Jack(c(1,2), c(), 2, algorithm = "naive"), 1)
    expect_equal(Jack(c(1,2), c(0,0), 2, algorithm = "naive"), 1)
    # gmp
    x <- as.bigq(1L, 2L)
    alpha <- as.bigq(3)
    expect_identical(Jack(x, NULL, alpha), as.bigq(1L))
    expect_identical(Jack(x, c(), alpha), as.bigq(1L))
    expect_identical(Jack(x, c(0,0), alpha), as.bigq(1L))
    expect_identical(Jack(x, NULL, alpha, algorithm = "naive"), as.bigq(1L))
    expect_identical(Jack(x, c(), alpha, algorithm = "naive"), as.bigq(1L))
    expect_identical(Jack(x, c(0,0), alpha, algorithm = "naive"), as.bigq(1L))
    # polynomial
    P <- JackPol(3, lambda = NULL, alpha = 4)
    expect_identical(P, mvp::constant(1))
    P <- JackPol(3, lambda = c(), alpha = 4)
    expect_identical(P, mvp::constant(1))
    P <- JackPol(3, lambda = c(0,0), alpha = 4)
    expect_identical(P, mvp::constant(1))
  }
)

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
  "Jack (3,1) - polynomial", {
    alpha <- 5/2
    m <- 4
    expected <- as_mvp_spray((2*alpha^2+4*alpha+2)*MSFspray(m,c(3,1)) +
      (6*alpha+10)*MSFspray(m,c(2,1,1)) + (4*alpha+4)*MSFspray(m,c(2,2)) +
      24*MSFspray(m,c(1,1,1,1)))
    obtained <- JackPol(m, c(3,1), alpha)
    expect_identical(expected$names, obtained$names)
    expect_identical(expected$power, obtained$power)
    expect_equal(expected$coeffs, obtained$coeffs)
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
    expect_true(jack_naive == expected)
    expect_true(jack_DK == expected)
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
    lambda <- c(3,2)
    alpha <- as.bigq(11L,3L)
    pol <- JackPol(4, lambda, alpha, algorithm = "naive")
    x <- as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L))
    polEval <- qspray::evalQspray(pol, x)
    expect_identical(polEval, Jack(as.bigq(x), lambda, alpha))
  }
)

test_that("JackPol gmp", {
  gmpol <-
    JackPol(3, lambda = c(3,1), alpha = gmp::as.bigq(2,3), algorithm = "DK")
  mvpol <-
    JackPol(3, lambda = c(3,1), alpha = 2/3, algorithm = "DK")
  gmvpol <- as_mvp_qspray(gmpol)
  expect_identical(gmvpol[["names"]], mvpol[["names"]])
  expect_identical(gmvpol[["power"]], mvpol[["power"]])
  expect_equal(gmvpol[["coeffs"]], mvpol[["coeffs"]])
})

test_that(
  "JackPolCPP is correct", {
    lambda <- c(3, 2)
    alpha <- as.bigq(11L, 3L)
    pol <- JackPolCPP(4, lambda, alpha)
    x <- as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L))
    polEval <- qspray::evalQspray(pol, x)
    expect_identical(polEval, Jack(as.bigq(x), lambda, alpha))
  }
)
