test_that("Some Green X-polynomials (comparison with Sage)", {
# Sym = SymmetricFunctions(FractionField(QQ['t']))
# HLP = Sym.hall_littlewood().P()
# p = Sym.power()
# HLP(p([2,2]))
# (t^6-t^5+t^4-2*t^3+t^2-t+1)*HLP[1, 1, 1, 1] + (t^3-t^2+t-1)*HLP[2, 1, 1] + (t^2-t+2)*HLP[2, 2] + (t-1)*HLP[3, 1] + HLP[4]
  GreenXpolys <- GreenXpolynomials(c(2, 2))
  t <- qlone(1)
  partitions <- c(
    "[1, 1, 1, 1]",
    "[2, 1, 1]",
    "[2, 2]",
    "[3, 1]",
    "[4]"
  )
  polynomials <- list(
    t^6-t^5+t^4-2*t^3+t^2-t+1, 
    t^3-t^2+t-1,
    t^2-t+2,
    t-1,
    qone()
  )
  names(polynomials) <- partitions
  checks <- vapply(partitions, function(part) {
    lst <- GreenXpolys[[part]]
    obtained <- lst[["polynomial"]]
    expected <- polynomials[[part]]
    obtained == expected
  }, logical(1L))
  expect_true(all(checks))
})


test_that("Green X-polynomial expansion", {
  
  mu  <- c(2, 1, 1, 1)
  rho <- c(2, 2, 1)

  lambdas <- listOfDominatingPartitions(mu)

  KFpolys <- lapply(lambdas, function(lambda) { 
    KostkaFoulkesPolynomial(lambda, mu)
  })

  chis <- lapply(lambdas, function(lambda) {
    chi_lambda_rho(lambda, rho)
  })

  toAdd <- mapply(
    function(chi, KFpoly) {
      chi * KFpoly
    },
    chis, KFpolys,
    USE.NAMES = FALSE, SIMPLIFY = FALSE
  )

  obtained <- Reduce(`+`, toAdd)

  muAsString <- partitionAsString(mu)
  expected <- GreenXpolynomials(rho = rho)[[muAsString]][["polynomial"]]

  expect_true(obtained == expected)
})


test_that("A Green Q-polynomial (rho = lambda = [2,2])", {
  # comparison with https://elad.zelingher.com/mathapps/gln/GreenPolynomials.html
  GreenQpolys <- GreenQpolynomials(c(2, 2))
  q <- qlone(1)
  obtained <- GreenQpolys[["[2, 2]"]][["polynomial"]]
  expected <- 2*q^2 - q + 1
  expect_true(obtained == expected)
})
