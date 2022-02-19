test_that("Julia", {
  skip_if_not(JuliaConnectoR::juliaSetupOk(), "Julia setup is not OK")
  julia <- Jack_julia()
  n <- 3
  lambda <- c(3, 2)
  alpha <- "2/3"
  alphaq <- gmp::as.bigq(alpha)
  mvpol_julia <- julia$JackPol(n, lambda, alpha, poly = "mvp")
  gmpol_julia <- julia$JackPol(n, lambda, alpha, poly = "gmpoly")
  gmpol_r     <- JackPol(n, lambda, alphaq)
  expect_true(mvpEqual(mvpol_julia, gmpoly::gmpoly2mvp(gmpol_julia)))
  expect_true(mvpEqual(
    gmpoly::gmpoly2mvp(gmpol_r), gmpoly::gmpoly2mvp(gmpol_julia)
  ))
  JuliaConnectoR::stopJulia()
})
