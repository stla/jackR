test_that(
  "Skew Jack polynomial branching rule", {
    nx <- 2
    ny <- 2
    n <- nx + ny
    lambda <- c(2, 2)
    alpha <- 2L
    which <- "P"
    mus <- list(integer(0), 1, 2, c(1, 1), c(2, 1), c(2, 2))
    ys <- list(qlone(3), qlone(4))
    jp <- JackPol(n, lambda, alpha, which)
    expected <-
      Reduce(
        `+`,
        lapply(mus, function(mu) {
          JackPol(nx, mu, alpha, which) *
            changeVariables(SkewJackPol(ny, lambda, mu, alpha, which), ys)
        })
      )
    expect_true(jp == expected)
  }
)

test_that(
  "Skew symbolic Jack polynomial branching rule", {
    nx <- 2
    ny <- 2
    n <- nx + ny
    lambda <- c(2, 2)
    which <- "Q"
    mus <- list(integer(0), 1, 2, c(1, 1), c(2, 1), c(2, 2))
    ys <- list(Qlone(3), Qlone(4))
    jp <- JackSymPol(n, lambda, which)
    expected <-
      Reduce(
        `+`,
        lapply(mus, function(mu) {
          JackSymPol(nx, mu, which) *
            changeVariables(SkewJackSymPol(ny, lambda, mu, which), ys)
        })
      )
    expect_true(jp == expected)
  }
)
