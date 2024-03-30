test_that("Littlewood-Richardson multiplication", {
  mu <- c(2, 1)
  nu <- c(3, 2, 1)
  LR <- LRmult(mu, nu, output = "list")
  LRcoeffs <- LR$coeff
  LRparts <- LR$lambda
  LRterms <- lapply(1:length(LRcoeffs), function(i) {
    LRcoeffs[i] * SchurPolCPP(3, LRparts[[i]])
  })
  smu_times_snu <- Reduce(`+`, LRterms)
  expect_true(smu_times_snu == SchurPolCPP(3, mu) * SchurPolCPP(3, nu))
})

test_that("LR-rule and Standard Young Tableaux counting", {
  # https://math.stackexchange.com/q/3012744/38217
  skip_if_not_installed("syt")
  mu <- c(3, 2, 2, 1)
  nu <- c(4, 3, 2, 1)
  LR <- LRmult(mu, nu, output = "list")
  h <- function(p) as.integer(syt::count_sytx(p))
  counts <- vapply(LR$lambda, h, integer(1L))
  rhs <- sum(LR$coeff * counts)
  lhs <- h(mu) * h(nu) * choose(sum(mu) + sum(nu), sum(mu))
  expect_true(lhs == rhs)
})

test_that("Skew Schur polynomial", {
  sspol <- SkewSchurPol(3, lambda = c(3, 2, 1), mu = c(1, 1))
  expected <- SchurPolCPP(3, c(2, 1, 1)) + SchurPolCPP(3, c(2, 2)) +
    SchurPolCPP(3, c(3, 1))
  expect_true(sspol == expected)
})


test_that("Skew Schur polynomial and skew semistandard tableaux", {
  skip_if_not_installed("syt")
  lambda <- c(4, 3, 2, 1)
  mu <- c(2, 1)
  n <- 4
  sssktx <- syt::all_ssSkewTableaux(lambda, mu, n)
  wt <- function(ssyt) {
    ssyt <- unlist(ssyt)
    vapply(1:n, function(k) {
      length(which(ssyt == k))
    }, integer(1L))
  }
  qlones <- lapply(1:n, qspray::qlone)
  monomial <- function(ssyt) {
    powers <- wt(ssyt)
    Reduce(`*`, lapply(1:n, function(k) qlones[[k]]^powers[k]))
  }
  monomials <- lapply(sssktx, monomial)
  obtained <- Reduce(`+`, monomials)
  expected <- SkewSchurPol(n, lambda, mu)
  expect_true(obtained == expected)
})
