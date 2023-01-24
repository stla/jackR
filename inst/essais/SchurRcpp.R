x <- jack:::SchurPolRcpp(n = 2L, lambda = 4L)
qspray::qsprayMaker(x[["exponents"]], x[["coeffs"]])
