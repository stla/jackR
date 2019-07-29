#' @importFrom partitions conjugate parts
NULL

dualPartition <- function(lambda){
  if(all(lambda == 0)) 0 else conjugate(lambda)
}

logHookLengths <- function(lambda, alpha){
  i <- rep(seq_along(lambda), times = lambda)
  j <- unlist(sapply(lambda, seq_len, simplify = FALSE))
  lambdaPrime <- dualPartition(lambda)
  upperHL <- lambdaPrime[j] - i + alpha*(lambda[i] - j + 1)
  lowerHL <- lambdaPrime[j] - i + 1 + alpha*(lambda[i] - j)
  log(c(upperHL,lowerHL))
}

.B <- function(nu, lambda, mu, alpha){
  if(all(nu == 0)) return(0)
  i <- rep(seq_along(nu), times = nu)
  j <- unlist(sapply(nu, seq_len, simplify = FALSE))
  nuPrime <- dualPartition(nu)
  lambdaPrime <- dualPartition(lambda); l1 <- length(lambdaPrime)
  muPrime <- dualPartition(mu); l2 <- length(muPrime)
  n <- length(i)
  out <- numeric(n)
  for(k in 1L:n){
    ii <- i[k]; jj <- j[k]
    lambdaPj <- if(jj <= l1) lambdaPrime[jj] else 0
    muPj <- if(jj <= l2) muPrime[jj] else 0
    out[k] <- if(lambdaPj == muPj){
      nuPrime[jj] - ii + alpha*(nu[ii] - jj + 1)
    }else{
      nuPrime[jj] - ii + 1 + alpha*(nu[ii] - jj)
    }
  }
  sum(log(out))
}

.beta <- function(lambda, mu, alpha){
  exp(.B(lambda, lambda, mu, alpha) - .B(mu, lambda, mu, alpha))
}

.N <- function(lambda, mu){
  n <- length(lambda)
  M <- sapply(1L:n, function(i) prod(tail(lambda+1,n-i)))
  sum(mu * M)
}

#####
isDominated <- function(mu, lambda){
  n <- sum(lambda)
  lambda <- c(lambda, rep(0,n-length(lambda)))
  mu <- mu[seq_len(match(0, mu, nomatch = n+1)-1L)]
  for(i in 1:length(mu)){
    if(sum(mu[1:i]) > sum(lambda[1:i])){
      return(FALSE)
    }
  }
  TRUE
}

dominatedPartitions <- function(lambda){
  allParts <- parts(sum(lambda))
  allParts[, apply(allParts, 2L, isDominated, lambda = lambda),
           drop = FALSE]
}

betweenPartitions <- function(mu, lambda){
  doms <- dominatedPartitions(lambda)
  n <- sum(mu)
  doms[,apply(doms, 2L, function(p){
    isDominated(mu,p) && !all(c(mu, rep(0,n-length(mu)))==p)
  }), drop = FALSE]
}

.rho <- function(lambda) sum(lambda*(lambda-seq_along(lambda)))
