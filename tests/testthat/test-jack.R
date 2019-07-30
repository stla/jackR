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
