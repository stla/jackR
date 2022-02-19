test_that("Julia", {
  skip_if_not(JuliaConnectoR::juliaSetupOk(), "Julia setup is not OK")
  julia <- Jack_julia()
  # numerical ####
  x <- c("1/2", "2/3", "5")
  xq <- gmp::as.bigq(x)
  lambda <- c(2, 1, 1)
  # jack
  alpha <- "2/3"
  alphaq <- gmp::as.bigq(alpha)
  expect_equal(
    julia$Jack(x, lambda, alpha), Jack(xq, lambda, alphaq)
  )
  # zonal
  expect_equal(
    julia$Zonal(x, lambda), Zonal(xq, lambda)
  )
  # zonalQ
  expect_equal(
    julia$ZonalQ(x, lambda), ZonalQ(xq, lambda)
  )
  # Schur
  expect_equal(
    julia$Schur(x, lambda), Schur(xq, lambda)
  )
  # polynomials ####
  n <- 3
  lambda <- c(3, 2)
  # jack
  alpha <- "2/3"
  alphaq <- gmp::as.bigq(alpha)
  mvpol_julia <- julia$JackPol(n, lambda, alpha, poly = "mvp")
  gmpol_julia <- julia$JackPol(n, lambda, alpha, poly = "gmpoly")
  gmpol_r     <- JackPol(n, lambda, alphaq)
  expect_true(mvpEqual(mvpol_julia, gmpoly::gmpoly2mvp(gmpol_julia)))
  expect_true(mvpEqual(
    gmpoly::gmpoly2mvp(gmpol_r), gmpoly::gmpoly2mvp(gmpol_julia)
  ))
  # zonal
  mvpol_julia <- julia$ZonalPol(n, lambda, poly = "mvp")
  gmpol_julia <- julia$ZonalPol(n, lambda, poly = "gmpoly")
  gmpol_r     <- ZonalPol(n, lambda)
  expect_true(mvpEqual(mvpol_julia, gmpoly::gmpoly2mvp(gmpol_julia)))
  expect_true(mvpEqual(
    gmpoly::gmpoly2mvp(gmpol_r), gmpoly::gmpoly2mvp(gmpol_julia)
  ))
  # zonalq
  mvpol_julia <- julia$ZonalQPol(n, lambda, poly = "mvp")
  gmpol_julia <- julia$ZonalQPol(n, lambda, poly = "gmpoly")
  gmpol_r     <- ZonalQPol(n, lambda)
  expect_true(mvpEqual(mvpol_julia, gmpoly::gmpoly2mvp(gmpol_julia)))
  expect_true(mvpEqual(
    gmpoly::gmpoly2mvp(gmpol_r), gmpoly::gmpoly2mvp(gmpol_julia)
  ))
  # schur
  mvpol_julia <- julia$SchurPol(n, lambda, poly = "mvp")
  gmpol_julia <- julia$SchurPol(n, lambda, poly = "gmpoly")
  gmpol_r     <- SchurPol(n, lambda)
  expect_true(mvpEqual(mvpol_julia, gmpoly::gmpoly2mvp(gmpol_julia)))
  expect_true(mvpEqual(
    gmpoly::gmpoly2mvp(gmpol_r), gmpoly::gmpoly2mvp(gmpol_julia)
  ))
  #
  JuliaConnectoR::stopJulia()
})
