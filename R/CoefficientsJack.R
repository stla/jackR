#' @importFrom multicool multinom
#' @importFrom gmp as.bigq asNumeric
NULL

#' Title
#'
#' @param n
#' @param alpha
#' @param until
#'
#' @return
#' @export
#' @importFrom gmp factorialZ
#'
#' @examples
JackCoefficientsQ <- function(n, alpha, until = NULL){
  allParts <- dominatedPartitions(n)
  nParts <- ncol(allParts)
  stringParts <- apply(allParts, 2L, toString)
  if(!is.null(until)){
    until <- toString(c(until, rep(0, n-length(until))))
    lastRow <- match(until, stringParts)
  }else{
    lastRow <- nParts
  }
  indices <- matrix(seq_len(lastRow*nParts),
                    nrow = lastRow, ncol = nParts)
  colnames(indices) <- stringParts
  rownames(indices) <- stringParts[1L:lastRow]
  coefs <- as.bigq(integer(lastRow*nParts))
  coefs[1L] <- as.bigq(1L)
  for(m in 1L:min(lastRow, nParts-1L)){
    kappa <- allParts[,m]
    for(k in (m+1L):nParts){
      lambda <- allParts[,k]
      btwn <- betweenPartitions(lambda, kappa)
      x <- as.bigq(0L)
      for(i in 1L:(n-1)){
        for(j in (i+1L):n){
          for(t in seq_len(lambda[j])){
            mu <- as.bigq(lambda)
            mu[i] <- mu[i] + as.bigq(t)
            mu[j] <- mu[j] - as.bigq(t)
            muOrd <- sort(asNumeric(mu), decreasing = TRUE)
            if(isDominated(muOrd, kappa)){
              x <- x + (mu[i]-mu[j]) /
                (.e(kappa, alpha)-.e(lambda, alpha)) *
                coefs[indices[m, toString(muOrd)]]
            }
          }
        }
      }
      coefs[indices[m, toString(lambda)]] <- x
    }
    if(m < lastRow){
      coefs[indices[m+1L,m+1L]] <- as.bigq(1L)
    }
  }
  # coefs <- matrix(coefs, nrow = lastRow, ncol = nParts)
  # lastColumn <- coefs[, nParts]
  # facto <- as.bigq(factorialZ(n))
  # for(i in 1L:lastRow){
  #   f <- facto / lastColumn[i]
  #   coefs[i,-nParts] <- f * coefs[i,-nParts]
  # }
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
