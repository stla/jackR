x <- jack:::SchurPolRcpp(n = 3L, lambda = c(4L, 0L, 0L))
qspray::qsprayMaker(x[["exponents"]], x[["coeffs"]])
jack::SchurPol(3L, 4L)
