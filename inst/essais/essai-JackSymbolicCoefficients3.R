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
  lambdas <- apply(allParts, 2L, function(part) {
    part[part != 0L]
  }, simplify = FALSE)
  stringParts <- vapply(lambdas, toString, character(1L))
  row <- rep(list(zeroRatioOfQsprays), nParts)
  names(row) <- stringParts
  coefs <- lapply(seq_len(nParts), function(i) {
    part <- list(unitRatioOfQsprays)
    names(part) <- stringParts[i]
    parts <- tail(stringParts, nParts - i)
    c(part, row[parts])
  })
  names(coefs) <- stringParts
  for(m in seq_len(nParts-1L)){
    lambda <- lambdas[[m]]
    eSymbolic_lambda <- .eSymbolic(lambda)
    for(k in (m+1L):nParts){
      mu <- lambdas[[k]]
      l <- length(mu)
      eSymbolic_mu <- .eSymbolic(mu)
      x <- zeroRatioOfQsprays
      for(i in 1L:(l-1L)){ # l is always >1
        for(j in (i+1L):l){
          for(t in seq_len(mu[j])){
            mucopy <- mu
            mucopy[i] <- mucopy[i] + t
            mucopy[j] <- mucopy[j] - t
            muOrd <- sort(mucopy[mucopy != 0L], decreasing = TRUE)
            if(isDominated(muOrd, lambda)){
              x <- x + (mucopy[i] - mucopy[j]) /
                (eSymbolic_lambda - eSymbolic_mu) *
                coefs[[m]][[toString(muOrd)]]
            }
          }
        }
      }
      coefs[[m]][[toString(mu)]] <- x
    }
  }
  if(which != "P") {
    names(lambdas) <- names(coefs)
    invPcoeffs <- lapply(lambdas, jack:::symbolicJackPcoefficientInverse)
    Names <- names(coefs)
    names(Names) <- Names
    if(which == "J") {
      coefs <- lapply(Names, function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], list(invPcoeffs[[lambda]]),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    } else if(which == "C") {
      Ccoefs <- lapply(lambdas, jack:::symbolicJackCcoefficient)
      factors <-
        mapply(`*`, Ccoefs, invPcoeffs, SIMPLIFY = FALSE, USE.NAMES = TRUE)
      coefs <- lapply(Names, function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], list(factors[[lambda]]),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    } else {
      invQcoeffs <- lapply(lambdas, jack:::symbolicJackQcoefficientInverse)
      factors <-
        mapply(`/`, invPcoeffs, invQcoeffs, SIMPLIFY = FALSE, USE.NAMES = TRUE)
      coefs <- lapply(Names, function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], list(factors[[lambda]]),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    }
  }
  coefs
}

JackSymbolicCoefficients(3, 4, which = "C")
