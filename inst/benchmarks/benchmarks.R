library(jack)
library(microbenchmark)

microbenchmark(
  JackPol = JackPolCPP(n = n, lambda = lambda, alpha = alpha),
  JackSymPol = JackSymPol(n = n, lambda = lambda),
  JackPol2 = evalSymbolicQspray(JackSymPol(n = n, lambda = lambda), a = alpha),
  setup = {
    n <- 5
    lambda <- c(4, 2, 2, 1)
    alpha <- 2
  },
  times = 5
)
