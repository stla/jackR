JackCombinationToQspray <- function(combo, n, alpha, which) {
  lambdas <- lapply(combo, `[[`, "lambda")
  coeffs <- lapply(combo, `[[`, "coeff")
  Reduce(`+`, mapply(function(lambda, coeff) {
    coeff * JackPol(n, lambda, alpha, which)
  }, lambdas, coeffs, SIMPLIFY = FALSE))
}
