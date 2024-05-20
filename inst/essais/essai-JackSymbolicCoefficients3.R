library(gmp)
library(partitions)
library(ratioOfQsprays)

partitionAsString <- function(lambda) {
  paste0("[", toString(lambda), "]")
}

.eSymbolic <- function(lambda){
  jack:::.n(jack:::dualPartition(lambda))*qlone(1) - jack:::.n(lambda)
}


.symbolicKostkaNumbers <- function(n, weight, which){
  #stopifnot(n > 0L, isPositiveInteger(n))
  zeroRatioOfQsprays <- as.ratioOfQsprays(0L)
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  if(weight == 1L) {
    return(list("[1]" = list("[1]" = unitRatioOfQsprays)))
  }
  allParts <- restrictedparts(weight, n)
  nParts <- ncol(allParts)
  lambdas <- apply(allParts, 2L, function(part) {
    part[part != 0L]
  }, simplify = FALSE)
  stringParts <- vapply(lambdas, partitionAsString, character(1L))
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
                coefs[[m]][[partitionAsString(muOrd)]]
            }
          }
        }
      }
      coefs[[m]][[partitionAsString(mu)]] <- x
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

sKN <- .symbolicKostkaNumbers(4, 4, which = "J")

invTriMatrix <- function(L) {
  d <- length(L)
  f <- 1L / L[[d]][[1L]]
  if(d == 1L) {
    invL <- L
    invL[[1L]][[1L]] <- f
    return(invL)
  } else {
    B <- invTriMatrix(lapply(head(L, -1L), function(row) {
      head(row, -1L)
    }))
    newColumn <- lapply(seq_len(d-1L), function(i) {
      toAdd <- mapply(
        `*`,
        B[[i]], lapply(i:(d-1L), function(j) {
          L[[j]][[d-j+1L]]
        }),
        SIMPLIFY = FALSE, USE.NAMES = FALSE
      )
      -Reduce(`+`, toAdd) * f
    })
    B <- c(B, list(list()))
    newColumn <- c(newColumn, list(f))
    names(B) <- names(newColumn) <- names(L)
    Names <- lapply(L, names)
    mapply(
      function(row, x, nms) {
        out <- c(row, list(x))
        names(out) <- nms
        out
      },
      B, newColumn, Names,
      SIMPLIFY = FALSE, USE.NAMES = TRUE
    )
    # out <- lapply(seq_along(L), function(i) {
    #   row <- c(B[[i]], newColumn[i])
    #   names(row) <- names(L[[i]])
    #   row
    # })
    # names(out) <- names(L)
    # out
  }
}

invTriMatrix(sKN)

