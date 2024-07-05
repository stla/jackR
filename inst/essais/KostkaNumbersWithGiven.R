KostkaNumbersWithGivenLambda <- function(lambda) {
  lambda <- removeTrailingZeros(as.integer(lambda))
  mus <- rev(listOfDominatedPartitions(lambda))
  kNumbers <- seq_along(mus)
  musAsStrings <-
    vapply(mus, partitionAsString, character(1L), USE.NAMES = FALSE)
  names(mus) <- names(kNumbers) <- musAsStrings
  ellLambda <- length(lambda)
  if(ellLambda == 1L) {
    kNumbers <- 1L
    return(kNumbers)
  }
  # psums <- cumsum(tail(lambda, -2L))
  # ellPsums <- ellLambda - 2L
  # seq_Psums <- seq_len(ellPsums)
  # .isDominated <- function(mu){
  #   for(i in seq_Psums){
  #     if(sum(mu[1L:i]) > psums[i]){
  #       return(FALSE)
  #     }
  #   }
  #   TRUE
  # }
  lambdap <- conjugate(lambda)
  nlambda <- sum(seq_len(ellLambda - 1L) * tail(lambda, -1L))
  nlambdap <- sum(seq_len(lambda[1L] - 1L) * tail(lambdap, -1L))
  elambda <- as.bigq(nlambdap - nlambda)
  kNumbers[1L] <- 1L
  for(muAsString in tail(musAsStrings, -1L)) {
    mu <- mus[[muAsString]]
    mup <- conjugate(mu)
    ellMu <- length(mu)
    nmu <- sum(seq_len(ellMu - 1L) * tail(mu, -1L))
    nmup <- sum(seq_len(mu[1L] - 1L) * tail(mup, -1L))
    emu <- as.bigq(nmup - nmu)
    ee <- elambda - emu
    x <- as.bigq(0L)
    for(i in 1L:(ellMu-1L)) {
      for(j in (i+1L):ellMu) {
        dmuij <- mu[i] - mu[j]
        for(t in seq_len(mu[j])) {
          kappa <- mu
          kappa[i] <- mu[i] + t
          kappa[j] <- mu[j] - t
          muOrd <- removeTrailingZeros(sort(kappa, decreasing = TRUE))
          if(partitionAsString(muOrd) %in% musAsStrings){
            x <- x + (dmuij + 2*t) /
              ee *
              kNumbers[partitionAsString(muOrd)]
          }
        }
      }
    }
    kNumbers[muAsString] <- as.integer(x)
  }
  kNumbers
}


.isDominated <- function(mu, lambda){
  # assumptions:
  # sum(mu) == sum(lambda)
  # length(mu) == length(lambda) == sum(mu)
  ellMu <- match(0L, mu, nomatch = length(mu) + 1L) - 1L
  for(i in seq_len(ellMu)){
    if(sum(mu[1L:i]) > sum(lambda[1L:i])){
      return(FALSE)
    }
  }
  TRUE
}
