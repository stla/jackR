test_that("Hall-Littlewood P (2,2,1)", {
  hlp <- HallLittlewoodP(5, c(2, 2, 1))
  t <- qlone(1)
  M221 <- as(MSFpoly(5, c(2, 2, 1)), "symbolicQspray")
  M2111 <- as(MSFpoly(5, c(2, 1, 1, 1)), "symbolicQspray")
  M11111 <- as(MSFpoly(5, c(1, 1, 1, 1, 1)), "symbolicQspray")
  expected <- M221 + (2-t-t^2)*M2111 + (5-4*t-4*t^2+t^3+t^4+t^5)*M11111
  expect_true(hlp == expected)
})
