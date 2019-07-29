#' @importFrom gmp as.bigq
NULL

JackEvalNaive <- function(x, lambda, alpha){
  stopifnot(isPartition(lambda), alpha > 0)
  gmp <- is.bigq(x) || is.bigq(alpha)
  if(gmp){
    stopifnot(is.bigq(x), is.bigq(alpha))
  }
  mus <- dominatedPartitions(lambda)
  lambda <- mus[,1L] # to add trailing zeros
  coefs <- JackCoefficients(sum(lambda), until = lambda, alpha)
  if(gmp){
    out <- as.bigq(0L)
    for(i in 1L:ncol(mus)){
      out <- out + MSF(x, mus[,i]) *
        as.bigq(coefs[toString(lambda), toString(mus[,i])])
    }
  }else{
    out <- 0
    for(i in 1L:ncol(mus)){
      out <- out + MSF(x, mus[,i]) *
        coefs[toString(lambda), toString(mus[,i])]
    }
  }
  out
}
