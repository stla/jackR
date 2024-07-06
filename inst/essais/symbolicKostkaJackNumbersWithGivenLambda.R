symbolicKostkaJackNumbersWithGivenLambda <- function(lambda) {
  stopifnot(isPartition(lambda))
  lambda <- removeTrailingZeros(as.integer(lambda))
  mus <- rev(listOfDominatedPartitions(lambda))
  nparts <- length(mus)
  musAsStrings <-
    vapply(mus, partitionAsString, character(1L), USE.NAMES = FALSE)
  kNumbers <- vector("list", nparts)
  names(kNumbers) <- musAsStrings
  kNumbers[[1L]] <- new(
    "ratioOfQsprays",
    numerator = qone(), denominator = qone()
  )
  if(nparts >= 2L) {
    alpha <- qlone(1L)
    names(mus) <- musAsStrings
    ellLambda <- length(lambda)
    lambdap <- conjugate(lambda)
    nlambda <- sum(seq_len(ellLambda - 1L) * tail(lambda, -1L))
    nlambdap <- sum(seq_len(lambda[1L] - 1L) * tail(lambdap, -1L))
    elambda <- alpha*nlambdap - nlambda
    for(muAsString in tail(musAsStrings, -1L)) {
      mu <- mus[[muAsString]]
      mup <- conjugate(mu)
      ellMu <- mup[1L]
      i_ <- seq_len(ellMu - 1L)
      nmu <- sum(i_ * tail(mu, -1L))
      nmup <- sum(seq_len(mu[1L] - 1L) * tail(mup, -1L))
      emu <- alpha*nmup - nmu
      ee <- elambda - emu
      x <- 0L
      for(i in i_) {
        mu_i <- mu[i]
        for(j in (i+1L):ellMu) {
          mu_j <- mu[j]
          dmuij <- mu_i - mu_j
          kappa <- mu
          for(t in seq_len(mu_j - 1L)) {
            kappa[i] <- mu_i + t
            kappa[j] <- mu_j - t
            kappaOrd <- sort(kappa, decreasing = TRUE)
            kappaOrdAsString <- partitionAsString(kappaOrd)
            if(kappaOrdAsString %in% musAsStrings){
              x <- x + kNumbers[[kappaOrdAsString]] * (dmuij + 2L*t)
            }
          }
          mu_i_plus_mu_j <- mu_i + mu_j
          kappa[i] <- mu_i_plus_mu_j
          kappaOrd <- sort(kappa[-j], decreasing = TRUE)
          kappaOrdAsString <- partitionAsString(kappaOrd)
          if(kappaOrdAsString %in% musAsStrings){
            x <- x + kNumbers[[kappaOrdAsString]] * mu_i_plus_mu_j
          }
        }
      }
      kNumbers[[muAsString]] <- x / ee
    }
  }
  kNumbers
}

symbolicKostkaJackNumbersWithGivenLambda(c(3,2))
symbolicKostkaJackNumbers(5)[["[3, 2]"]]

