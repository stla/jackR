test_that(
  "Skew Jack polynomial for alpha=1 is collinear to skew Schur polynomial", {
    n <- 4L
    lambda <- c(4, 2)
    mu     <- c(2, 1)
    skjp <- SkewJackPol(n, lambda, mu, alpha = 1L)
    sksp <- SkewSchurPol(n, lambda, mu)
    expect_true(qspray::collinearQsprays(skjp, sksp))
  }
)
