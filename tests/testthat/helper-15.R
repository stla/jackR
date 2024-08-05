chi_lambda_rho <- function(lambda, rho) {
  chi_lambda_mu_rho(lambda, integer(0), rho)
}

listOfDominatingPartitions <- function(mu) {
  kappas <- listOfDominatedPartitions(partitions::conjugate(mu))
  lapply(kappas, partitions::conjugate)
}
