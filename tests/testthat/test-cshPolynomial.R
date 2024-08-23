test_that("CSH polynomials and Hall inner product", {
    lambda1 <- c(3, 1, 1, 1)
    lambda2 <- c(3, 2)
    n <- 4
    P <- 7 * msPolynomial(n, lambda1) - 9 * msPolynomial(n, lambda2)
    cshPol1 <- cshPolynomial(n, lambda1)
    cshPol2 <- cshPolynomial(n, lambda2)
    hp1 <- HallInnerProduct(P, cshPol1)
    hp2 <- HallInnerProduct(P, cshPol2)
    expect_equal(as.integer(hp1), 7L)
    expect_equal(as.integer(hp2), -9L)
})