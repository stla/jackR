
x <- jack:::SchurPolRcpp(n = 3L, lambda = 4L)
qspray::qsprayMaker(x[["exponents"]], x[["coeffs"]])
jack::SchurPolR(3L, 4L)


library(microbenchmark)
library(jack)
library(qspray)

rcpp <- function(n, lambda) {
  x <- jack:::SchurPolRcpp(as.integer(n), as.integer(lambda))
  qsprayMaker(x[["exponents"]], x[["coeffs"]])
}

julia <- Jack_julia()

n <- 5
lambda <- c(4, 3, 3)
microbenchmark(
  Rcpp  = rcpp(n, lambda),
  Julia = julia$SchurPolR(n, lambda),
  times = 6L
)

JuliaConnectoR::stopJulia()
