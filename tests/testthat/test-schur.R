test_that(
  # https://math.stackexchange.com/questions/3335885/expansion-of-sum-x-in-in-schur-polynomials
  "Schur expansion of (sum x_i)^n", {
    # numeric
    x <- c(3,4,5,6)
    e <- Schur(x, c(4)) + 3*Schur(x, c(3,1)) + 2*Schur(x, c(2,2)) +
      3*Schur(x, c(2,1,1)) + Schur(x, c(1,1,1,1))
    expect_equal(e, sum(x)^4)
    # gmp
    x <- as.bigq(c(3L,4L,5L,6L), c(4L,5L,6L,7L))
    e <- Schur(x, c(4)) + 3L*Schur(x, c(3,1)) + 2L*Schur(x, c(2,2)) +
      3L*Schur(x, c(2,1,1)) + Schur(x, c(1,1,1,1))
    expect_identical(e, sum(x)^4)
    # polynomial
    n <- 4
    P <- SchurPol(n, c(4)) + 3*SchurPol(n, c(3, 1)) + 2*SchurPol(n, c(2, 2)) +
      3*SchurPol(n, c(2, 1, 1)) + SchurPol(n, c(1, 1, 1, 1))
    Q <- (mvp("x_1", 1, 1) + mvp("x_2", 1, 1) + mvp("x_3", 1, 1) +
            mvp("x_4", 1, 1))^4
    expect_true(as_mvp_qspray(P) == Q)
  }
)

test_that(
  "Schur = 0 if l(lambda)>l(x)", {
    # numeric
    expect_equal(Schur(c(1,2), c(3,2,1)), 0)
    expect_equal(Schur(c(1,2), c(3,2,1), algorithm = "naive"), 0)
    # gmp
    x <- as.bigq(c(1L,2L))
    lambda <- c(3,2,1)
    expect_identical(Schur(x, lambda), as.bigq(0L))
    expect_identical(Schur(x, lambda, algorithm = "naive"), as.bigq(0L))
    # polynomial
    n <- 2
    lambda <- c(3,2,1)
    expect_true(SchurPol(n, lambda) == as.qspray(0))
    expect_identical(SchurPol(n, lambda, algorithm = "naive"),
                     as.qspray(0))
    expect_identical(SchurPol(n, lambda, exact = FALSE, algorithm = "naive"),
                     mvp::constant(0))
    expect_identical(SchurPol(n, lambda, algorithm = "naive",
                             basis = "MSF"),
                     as.qspray(0))
    expect_identical(SchurPol(n, lambda, exact = FALSE, algorithm = "naive",
                             basis = "MSF"),
                     mvp::constant(0))
  }
)


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
    lambda <- c(3,2)
    pol <- SchurPol(4, lambda, algorithm = "naive")
    x <- as.bigq(c(6L,-7L,8L,9L), c(1L,2L,3L,4L))
    polEval <- qspray::evalQspray(pol, x)
    expect_identical(polEval, Schur(as.bigq(x), lambda))
  }
)

test_that(
  "Pieri rule", {
    n <- 3
    P1 <- SchurPol(n, c(3, 2)) + 2 * SchurPol(n, c(2, 2, 1)) +
      SchurPol(n, c(3, 1, 1)) + 2 * SchurPol(n, c(2, 1, 1, 1)) +
      SchurPol(n, c(1, 1, 1, 1, 1))
    P2 <- ESFpoly(n, c(2, 2, 1))
    expect_true(P1 == P2)
  }
)
