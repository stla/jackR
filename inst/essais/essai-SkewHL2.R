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

Paths <- function(n, lambda, mu) {
  mu0 <- mu
  mu <- c(mu, rep(0L, length(lambda) - length(mu)))
  diffs <- lambda - mu
  Grid <- as.matrix(expand.grid(lapply(diffs, function(i) c(0L, seq_len(i)))))
  o <- qspray:::lexorder(Grid)
  Grid <- Grid[o, ]
  kappas <- Filter(jack:::isDecreasing, apply(as.matrix(Grid), 1L, function(kappa) {
    kappa + mu
  }, simplify = FALSE))
  combos <- arrangements::combinations(length(kappas), n-1L, replace = TRUE)
  Filter(columnStrictTableau, apply(combos, 1L, function(combo) {
    Filter(function(nu) length(nu) > 0, c(list(lambda), lapply(combo, function(i) {
      kappas[[i]]
    }), list(mu)))
  }, simplify = FALSE))
}

lambda <- c(4, 3)
mu <- c(2, 2)
n <- sum(lambda) - sum(mu) + 1
paths <- Paths(n, lambda, mu)
Pskew <- Qzero()
lones <- lapply(1L:n, Qlone)
for(j in 1:length(paths)) {
  #  nu <- Filter(function(kappa) any(kappa)>0, paths[[i]])
  nu <- rev(paths[[j]])
  Pskew <- Pskew + Reduce(`*`, lapply(1L+seq_len(length(nu)-1L), function(i) {
    psi(nu[[i]], nu[[i-1]]) * lones[[i-1]]^(sum(nu[[i]]-nu[[i-1]]))
  }))
}
Pskew
SkewSchurPol(n, lambda, mu)
substituteParameters(Pskew, 0)


