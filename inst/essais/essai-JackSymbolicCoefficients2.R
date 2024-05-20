library(gmp)
library(partitions)
library(ratioOfQsprays)

.eSymbolic <- function(lambda){
  jack:::.n(jack:::dualPartition(lambda))*qlone(1) - jack:::.n(lambda)
}


JackSymbolicCoefficients <- function(n, weight, which){
  #stopifnot(n > 0L, isPositiveInteger(n))
  zeroRatioOfQsprays <- as.ratioOfQsprays(0L)
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  if(weight == 1L) {
    return(list("1" = list("1" = unitRatioOfQsprays)))
  }
  allParts <- restrictedparts(weight, n)
  nParts <- ncol(allParts)
  stringParts <- apply(allParts, 2L, toString)
  row <- rep(list(zeroRatioOfQsprays), nParts)
  names(row) <- stringParts
  coefs <- rep(list(row), nParts)
  names(coefs) <- stringParts
  coefs[[1L]][[1L]] <- unitRatioOfQsprays
  for(m in seq_len(nParts-1L)){
    lambda0 <- allParts[, m]
    lambda <- lambda0[lambda0 != 0L]
    eSymbolic_lambda <- .eSymbolic(lambda)
    for(k in (m+1L):nParts){
      mu0 <- allParts[, k]
      mu <- mu0[mu0 != 0L]
      l <- length(mu)
      eSymbolic_mu <- .eSymbolic(mu)
      x <- zeroRatioOfQsprays
      for(i in 1L:(l-1L)){ # l is always >1
        for(j in (i+1L):l){
          for(t in seq_len(mu[j])){
            mucopy <- mu0
            mucopy[i] <- mucopy[i] + t
            mucopy[j] <- mucopy[j] - t
            muOrd <- sort(mucopy, decreasing = TRUE)
            if(isDominated(muOrd[muOrd != 0L], lambda)){
              x <- x + (mucopy[i] - mucopy[j]) /
                (eSymbolic_lambda - eSymbolic_mu) *
                coefs[[m]][[toString(muOrd)]]
            }
          }
        }
      }
      coefs[[m]][[toString(mu0)]] <- x
    }
    coefs[[m+1L]][[m+1L]] <- unitRatioOfQsprays
  }
  if(which != "P") {
    lambdas <- apply(allParts, 2L, function(part) {
      part[part != 0L]
    }, simplify = FALSE)
    invPcoeffs <- lapply(lambdas, jack:::symbolicJackPcoefficientInverse)
    names(invPcoeffs) <- names(coefs)
    if(which == "J") {
      coefs <- lapply(names(coefs), function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], rep(list(invPcoeffs[[lambda]]), nParts),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    } else if(which == "C") {
      Ccoefs <- lapply(lambdas, jack:::symbolicJackCcoefficient)
      factors <-
        mapply(`*`, Ccoefs, invPcoeffs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      coefs <- lapply(coefs, function(row) {
        mapply(
          `*`,
          row, factors,
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    } else {
      invQcoeffs <- lapply(lambdas, jack:::symbolicJackQcoefficientInverse)
      factors <-
        mapply(`/`, invPcoeffs, invQcoeffs, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      coefs <- lapply(coefs, function(row) {
        mapply(
          `*`,
          row, factors,
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    }
  }
  coefs
}

JackSymbolicCoefficients(3, 4, which = "J")
