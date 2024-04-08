library(symbolicQspray)

n <- 2L
lambda <- c(3L, 1L)
Pl <- jack:::JackSymPolRcpp(n, lambda)
P <- symbolicQspray:::symbolicQspray_from_list(Pl)
P
jack::JackPol(n, lambda, gmp::as.bigq(2))
