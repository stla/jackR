library(jack)


# skew Hall-Littlewood

horizontalStrip <- function(lambda, mu) {
  mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  x <- lambda[1] >= mu[1]
  i <- 1L
  while(x && i < length(lambda)) {
    x <- mu[i] >= lambda[i+1] && lambda[i+1] >= mu[i+1]
    i <- i+1
  }
  return(x)
  lambdaPrime <- conjugate(lambda)
  muPrime <- conjugate(mu)
  muPrime <- c(muPrime, rep(0L, length(lambdaPrime) - length(muPrime)))
  thetaPrime <- lambdaPrime - muPrime
  all(thetaPrime %in% c(0L, 1L))
}

columnStrictTableau <- function(tableau) {
  all(
    mapply(
      horizontalStrip,
      head(tableau, -1L), tail(tableau, -1L),
      SIMPLIFY = TRUE, USE.NAMES = FALSE
    )
  )
}

psi <- function(lambda, mu) {
  t <- qlone(1L)
  out <- qone()
  i_ <- seq_len(lambda[1L])
  mlambda <- vapply(i_, function(i) {
    sum(lambda == i)
  }, integer(1L))
  mmu <- vapply(i_, function(i) {
    sum(mu == i)
  }, integer(1L))
  for(i in i_) {
    mmu_i <- mmu[i]
    if(mmu_i == mlambda[i]+1) {
      out <- out * (1L - t^mmu_i)
    }
  }
  out
}

phi <- function(lambda, mu) {
  t <- qlone(1L)
  out <- qone()
  i_ <- seq_len(lambda[1L])
  mlambda <- vapply(i_, function(i) {
    sum(lambda == i)
  }, integer(1L))
  mmu <- vapply(i_, function(i) {
    sum(mu == i)
  }, integer(1L))
  for(i in i_) {
    mlambda_i <- mlambda[i]
    if(mmu[i]+1L == mlambda_i) {
      out <- out * (1L - t^mlambda_i)
    }
  }
  out
}

lambda <- c(3, 1)
mu <- c(2, 1)
jack:::b(lambda) / jack:::b(mu) * psi(lambda, mu)
phi(lambda, mu)

Combos <- function(a, b, n) {
  if(n == 0L) {
    return(matrix(NA_integer_, nrow = 1L, ncol = 0L))
  }
  if(n == 1L) {
    return(cbind(a:b))
  }
  do.call(rbind, lapply(a:b, function(i) {
    cbind(i, Combos(i, b, n-1L))
  }))
}

cartesianProduct <- function(diffs) {
  if(length(diffs) == 1L) {
    return(cbind(rev(c(0L, seq_len(diffs)))))
  }
  previous <- cartesianProduct(tail(diffs, -1L))
  do.call(rbind, lapply(rev(c(0L, seq_len(diffs[1L]))), function(i) {
    cbind(i, previous)
  }))
}

# assumes lambda is clean and length(mu)=length(lambda)
Paths <- function(n, lambda, mu) {
  # mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  diffs <- lambda - mu
  # Grid <- as.matrix(expand.grid(lapply(diffs, function(i) c(0L, seq_len(i)))))
  # o <- qspray:::lexorder(Grid)
  # Grid <- Grid[o, ]
  Grid <- cartesianProduct(diffs)
  kappas <- Filter(jack:::isDecreasing, apply(Grid, 1L, function(kappa) {
    kappa + mu
  }, simplify = FALSE))
  # combos <- arrangements::combinations(length(kappas), n-1L, replace = TRUE)
  combos <- Combos(1L, length(kappas), n - 1L)
  Filter(columnStrictTableau, apply(combos, 1L, function(combo) {
    c(list(lambda), lapply(combo, function(i) {
      kappas[[i]]
    }), list(mu))
  }, simplify = FALSE))
}

# diffs <- c(1,1,2)
# Grid <- as.matrix(expand.grid(lapply(diffs, function(i) c(0L, seq_len(i)))))
# o <- qspray:::lexorder(Grid)
# all(Grid[o, ] == cartesianProduct(diffs))


SkewHallLittlewoodP <- function(n, lambda, mu) {
  lambda <- as.integer(jack:::removeTrailingZeros(lambda))
  mu <- as.integer(jack:::removeTrailingZeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  mu <- c(mu, rep(0L, ellLambda - ellMu))
  if(any(lambda - mu < 0L)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  if(n == 0L){
    if(all(lambda == mu)) {
      return(Qone())
    } else {
      return(Qzero())
    }
  }
  paths <- Paths(n, lambda, mu)
  Pskew <- Qzero()
  lones <- lapply(1L:n, Qlone)
  for(j in 1:length(paths)) {
    nu <- rev(paths[[j]])
    Pskew <- Pskew + Reduce(`*`, lapply(1L+seq_len(length(nu)-1L), function(i) {
      psi(nu[[i]], nu[[i-1]]) * lones[[i-1]]^(sum(nu[[i]]-nu[[i-1]]))
    }))
  }
  Pskew
}

SkewHallLittlewoodQ <- function(n, lambda, mu) {
  lambda <- as.integer(jack:::removeTrailingZeros(lambda))
  mu <- as.integer(jack:::removeTrailingZeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  mu <- c(mu, rep(0L, ellLambda - ellMu))
  if(any(lambda - mu < 0L)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  if(n == 0L){
    if(all(lambda == mu)) {
      return(Qone())
    } else {
      return(Qzero())
    }
  }
  paths <- Paths(n, lambda, mu)
  Qskew <- Qzero()
  lones <- lapply(1L:n, Qlone)
  for(j in 1:length(paths)) {
    nu <- rev(paths[[j]])
    Qskew <- Qskew + Reduce(`*`, lapply(1L+seq_len(length(nu)-1L), function(i) {
      phi(nu[[i]], nu[[i-1]]) * lones[[i-1]]^(sum(nu[[i]]-nu[[i-1]]))
    }))
  }
  Qskew
}

lambda <- c(3, 2, 1)
mu <- c(2, 1)
n <- sum(lambda) - sum(mu)
Pskew <- SkewHallLittlewoodP(n, lambda, mu)
SkewSchurPol(n, lambda, mu)
substituteParameters(Pskew, 0)
Qskew <- SkewHallLittlewoodQ(n, lambda, mu)

jack:::b(lambda) / jack:::b(mu) * Pskew == Qskew

