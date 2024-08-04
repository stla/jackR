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
