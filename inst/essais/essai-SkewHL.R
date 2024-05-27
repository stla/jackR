library(jack)

# skew Schur ####
content <- function(tableau) {
  entries <- unlist(tableau)
  entries <- entries[!is.na(entries)]
  n <- max(entries)
  vapply(seq_len(n), function(i) {
    sum(entries == i)
  }, integer(1L))
}

lambda <- c(2, 2, 1)
mu <- c(2, 1)

n <- 4
sktx <- syt::all_ssSkewTableaux(lambda, mu, n)

p <- qzero()
lones <- lapply(1:n, qlone)
for(skt in sktx) {
  ct <- content(skt)
  q <- qone()
  for(i in seq_along(ct)) {
    q <- q * lones[[i]]^(ct[i])
  }
  p <- p + q
}
p
SkewSchurPol(n, lambda, mu)

# skew Hall-Littlewood
phi_theta <- function(skt) {
  lambda <- lengths(skt)
  mu <- jack:::removeTrailingZeros(vapply(skt, function(row) {
    sum(is.na(row))
  }, integer(1L)))
  lambdaPrime <- conjugate(lambda)
  muPrime <- conjugate(mu)
  l <- length(lambdaPrime)
  muPrime <- c(muPrime, rep(0L, l - length(muPrime)))
  thetaPrime <- lambdaPrime - muPrime
  i_ <- integer(0L)
  for(k in seq_len(length(thetaPrime) - 1L)) {
    if(thetaPrime[k] > thetaPrime[k+1]) {
      i_ <- c(i_, k)
    }
  }
  if(length(i_) > 0L) {
    m_ <- vapply(i_, function(i) {
      sum(lambda == i)
    }, integer(1L))
    t <- qlone(1L)
    Reduce(`*`, lapply(m_, function(m) {
      1L - t^m
    }))
  } else {
    qone()
  }
}

lambda <- c(2, 2, 1)
mu <- c(2, 2)
phi_theta(lambda, mu)

lambdais <- function(lambda, mu) {
  lambda0 <- mu
  mugrown <- c(mu, rep(0L, length(lambda) - length(mu)))
  diffs <- lambda - mugrown
  ws <- which(diffs > 0L)
  r <- length(ws)
  lambdas <- vector("list", r)
  lambdas[[1L]] <- mu
  for(i in 1L+seq_len(r-1L)) {
    kappa <- lambdas[[i-1L]]
    kappa[ws[i-1L]] <- lambda[ws[i-1L]]
    lambdas[[i]] <- kappa
  }
  c(lambdas, list(lambda))
}

lambdais(c(3,3,3,2), c(2,1,1))
