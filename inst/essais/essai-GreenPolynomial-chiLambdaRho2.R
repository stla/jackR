library(jack)

mu <- c(2, 1, 1, 1)
rho <- c(2, 2, 1)

chi_lambda_rho <- function(lambda, rho) {
  jack:::chi_lambda_mu_rho(lambda, integer(0), rho)
}

listOfDominatingPartitions <- function(mu) {
  kappas <- jack:::listOfDominatedPartitions(partitions::conjugate(mu))
  lapply(kappas, partitions::conjugate)
}

lambdas <- listOfDominatingPartitions(mu)

KFpolys <- lapply(lambdas, function(lambda) { 
  KostkaFoulkesPolynomial(lambda, mu)
})

chis <- lapply(lambdas, function(lambda) {
  chi_lambda_rho(lambda, rho)
})

toAdd <- mapply(
  function(chi, KFpoly) {
    chi * KFpoly
  },
  chis, KFpolys,
  USE.NAMES = FALSE, SIMPLIFY = FALSE
)

obtained <- Reduce(`+`, toAdd)

muAsString <- jack:::partitionAsString(mu)
expected <- GreenXpolynomials(rho = rho)[[muAsString]][["polynomial"]]

obtained == expected



