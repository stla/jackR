#' @importFrom partitions conjugate parts
#' @importFrom gmp as.bigq is.bigq
NULL

isPartition <- function(lambda){
  all(floor(lambda) == lambda) && all(diff(lambda) <= 0)
}

dualPartition <- function(lambda){
  if(all(lambda == 0)) 0 else conjugate(lambda)
}

logHookLengths <- function(lambda, alpha){
  i <- rep(seq_along(lambda), times = lambda)
  j <- unlist(sapply(lambda, seq_len, simplify = FALSE))
  lambdaPrime <- dualPartition(lambda)
  upperHL <- lambdaPrime[j] - i + alpha*(lambda[i] - j + 1)
  lowerHL <- lambdaPrime[j] - i + 1 + alpha*(lambda[i] - j)
  log(c(upperHL, lowerHL))
}

hookLengths_gmp <- function(lambda, alpha){
  i <- rep(seq_along(lambda), times = lambda)
  j <- unlist(sapply(lambda, seq_len, simplify = FALSE))
  lambdaPrime <- as.bigq(dualPartition(lambda))
  lambda <- as.bigq(lambda)
  upperHL <- lambdaPrime[j] - i + alpha*(lambda[i] - j + 1L)
  lowerHL <- lambdaPrime[j] - i + 1L + alpha*(lambda[i] - j)
  c(upperHL, lowerHL)
}
.Blog <- function(nu, lambda, mu, alpha){
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
  if(all(lambda == mu)) return(1)
  exp(.Blog(lambda, lambda, mu, alpha) - .Blog(mu, lambda, mu, alpha))
}

.B_gmp <- function(nu, lambda, mu, alpha){
  if(all(nu == 0)) return(as.bigq(1L))
  i <- rep(seq_along(nu), times = nu)
  j <- unlist(sapply(nu, seq_len, simplify = FALSE))
  nuPrime <- as.bigq(dualPartition(nu))
  lambdaPrime <- dualPartition(lambda); l1 <- length(lambdaPrime)
  muPrime <- dualPartition(mu); l2 <- length(muPrime)
  n <- length(i)
  out <- as.bigq(integer(n))
  for(k in 1L:n){
    ii <- i[k]; jj <- j[k]
    lambdaPj <- if(jj <= l1) lambdaPrime[jj] else 0
    muPj <- if(jj <= l2) muPrime[jj] else 0
    out[k] <- if(lambdaPj == muPj){
      nuPrime[jj] - ii + alpha*(nu[ii] - jj + 1L)
    }else{
      nuPrime[jj] - ii + 1L + alpha*(nu[ii] - jj)
    }
  }
  sum(prod(out))
}

.beta_gmp <- function(lambda, mu, alpha){
  if(all(lambda == mu)) return(as.bigq(1L))
  .B_gmp(lambda, lambda, mu, alpha) / .B_gmp(mu, lambda, mu, alpha)
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

.rhoQ <- function(lambda) sum(lambda*(lambda-4*seq_along(lambda)))

#####
.n <- function(lambda){
  sum((seq_len(length(lambda))-1)*lambda)
}

.e <- function(lambda, alpha){
  if(is.bigq(alpha)){
    alpha * as.bigq(.n(dualPartition(lambda))) - .n(lambda)
  }else{
    alpha * .n(dualPartition(lambda)) - .n(lambda)
  }
}
