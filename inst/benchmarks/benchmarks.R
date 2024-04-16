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

# Haskell is faster:
#
# benchmarking Jack/jackPol with the given alpha
# mean                 33.76 ms
# std dev              3.004 ms
#
# benchmarking Jack/jackSymbolicPol
# mean                 543.8 ms
# std dev              44.58 ms
#
# benchmarking Jack/jackSymbolicPol evaluated at alpha
# mean                 621.9 ms
# std dev              50.79 ms
