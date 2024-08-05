library(jack)

#lambda <- c(3, 1)
mu <- c(2, 1, 1)
rho <- c(2, 2)

#jack:::chi_lambda_mu_rho(lambda, integer(0), rho)
chi_lambda_rho <- function(lambda, rho) {
  jack:::chi_lambda_mu_rho(lambda, integer(0), rho)
}

Knumbers <- syt::KostkaNumbers(4)
lambdasAsStrings <- rownames(Knumbers)
lambdas <- lapply(lambdasAsStrings, jack:::fromPartitionAsString)

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

Reduce(`+`, toAdd)

GreenXpolynomials(rho = rho)

