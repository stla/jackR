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
  mlambda <- vapply(1:lambda[1], function(i) {
    sum(lambda == i)
  }, integer(1L))
  mmu <- vapply(1:lambda[1], function(i) {
    sum(mu == i)
  }, integer(1L))
  for(i in 1:lambda[1]) {
    if(mmu[i] == mlambda[i]+1) {
      out <- out * (1-t^mmu[i])
    }
  }
  out
}

phi <- function(lambda, mu) {
  t <- qlone(1L)
  out <- qone()
  mlambda <- vapply(1:lambda[1], function(i) {
    sum(lambda == i)
  }, integer(1L))
  mmu <- vapply(1:lambda[1], function(i) {
    sum(mu == i)
  }, integer(1L))
  for(i in 1:lambda[1]) {
    if(mmu[i]+1L == mlambda[i]) {
      out <- out * (1-t^mlambda[i])
    }
  }
  out
}

lambda <- c(3, 1)
mu <- c(2, 1)
jack:::b(lambda) / jack:::b(mu) * psi(lambda, mu)
phi(lambda, mu)

Combos <- function(a, b, n) {
  if(n == 1L) {
    return(cbind(a:b))
  }
  do.call(rbind, lapply(a:b, function(i) {
    cbind(i, Combos(i, b, n-1L))
  }))
}

Paths <- function(n, lambda, mu) {
  mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  diffs <- lambda - mu
  Grid <- as.matrix(expand.grid(lapply(diffs, function(i) c(0L, seq_len(i)))))
  o <- qspray:::lexorder(Grid)
  Grid <- Grid[o, ]
  kappas <- Filter(jack:::isDecreasing, apply(as.matrix(Grid), 1L, function(kappa) {
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

SkewHallLittlewoodP <- function(n, lambda, mu) {
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

