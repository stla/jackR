library(gmp)
library(partitions)
library(ratioOfQsprays)

.eSymbolic <- function(lambda){
  jack:::.n(jack:::dualPartition(lambda))*qlone(1) - jack:::.n(lambda)
}


JackSymbolicCoefficients <- function(n, weight, which){
  #stopifnot(n > 0L, isPositiveInteger(n))
  if(weight == 1L) {
    out <- list(as.ratioOfQsprays(1L))
    dim(out) <- c(1L, 1L)
    rownames(out) <- colnames(out) <- "1"
    return(out)
  }
  allParts <- restrictedparts(weight, n)
  nParts <- ncol(allParts)
  stringParts <- apply(allParts, 2L, toString)
  coefs <- rep(list(as.ratioOfQsprays(0L)), nParts*nParts)
  dim(coefs) <- c(nParts, nParts)
  colnames(coefs) <- rownames(coefs) <- stringParts
  coefs[[1L, 1L]] <- as.ratioOfQsprays(1L)
  for(m in seq_len(nParts-1L)){
    lambda0 <- allParts[, m]
    lambda <- lambda0[lambda0 != 0L]
    eSymbolic_lambda <- .eSymbolic(lambda)
    for(k in (m+1L):nParts){
      mu0 <- allParts[, k]
      mu <- mu0[mu0 != 0L]
      l <- length(mu)
      eSymbolic_mu <- .eSymbolic(mu)
      x <- as.ratioOfQsprays(0L)
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
                coefs[[m, toString(muOrd)]]
            }
          }
        }
      }
      coefs[[m, toString(mu0)]] <- x
    }
    coefs[[m+1L,m+1L]] <- as.ratioOfQsprays(1L)
  }
  coefs <- lapply(coefs, showRatioOfQspraysXYZ("alpha"))

  dim(coefs) <- c(nParts, nParts)
  rownames(coefs) <- colnames(coefs) <- stringParts
  return(coefs)
  lastColumn <- coefs[indices[, nParts]]
  facto <- as.bigq(factorialZ(n))
  for(i in 1L:lastRow){
    f <- facto / lastColumn[i]
    coefs[indices[i,-nParts]] <- f * coefs[indices[i,-nParts]]
  }
  out <- cbind(matrix(as.character(coefs[indices[,-nParts]]),
                      nrow=lastRow, ncol = nParts-1L),
               rep(as.character(facto), lastRow))
  dimnames(out) <- dimnames(indices)
  out
}

JackSymbolicCoefficients(3, 4)
