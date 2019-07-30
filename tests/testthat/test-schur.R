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
