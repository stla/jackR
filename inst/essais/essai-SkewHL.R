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
phi_lambda_mu <- function(lambda, mu) {
  lambdaPrime <- conjugate(lambda)
  muPrime <- conjugate(mu)
  l <- length(lambdaPrime)
  muPrime <- c(muPrime, rep(0L, l - length(muPrime)))
  thetaPrime <- lambdaPrime - muPrime
  i_ <- which(head(thetaPrime, -1L) - tail(thetaPrime, -1L) == 1L)
  # i_ <- integer(0L)
  # for(k in seq_len(length(thetaPrime) - 1L)) {
  #   if(thetaPrime[k] > thetaPrime[k+1]) {
  #     i_ <- c(i_, k)
  #   }
  # }
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

psi_lambda_mu <- function(lambda, mu) {
  lambdaPrime <- conjugate(lambda)
  muPrime <- conjugate(mu)
  l <- length(lambdaPrime)
  muPrime <- c(muPrime, rep(0L, l - length(muPrime)))
  thetaPrime <- lambdaPrime - muPrime
  i_ <- which(head(thetaPrime, -1L) - tail(thetaPrime, -1L) == -1L)
  # i_ <- integer(0L)
  # for(k in seq_len(length(thetaPrime) - 1L)) {
  #   if(thetaPrime[k] > thetaPrime[k+1]) {
  #     i_ <- c(i_, k)
  #   }
  # }
  if(length(i_) > 0L) {
    m_ <- vapply(i_, function(i) {
      sum(mu == i)
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
phi_lambda_mu(lambda, mu)

lambdais <- function(lambda, mu) {
  mugrown <- c(mu, rep(0L, length(lambda) - length(mu)))
  diffs <- lambda - mugrown
  ws <- which(diffs > 0L)
  r <- length(ws)
  lambdas <- vector("list", r+1L)
  lambdas[[1L]] <- mu
  lambdas[[r+1L]] <- lambda
  for(i in seq_len(r-1L)) {
    kappa <- lambdas[[i]]
    kappa[ws[i]] <- lambda[ws[i]]
    lambdas[[i+1L]] <- kappa
  }
  lambdas
}

lambdais(c(3,3,3,2), c(2,1,1))

phiT <- function(skt) {
  lambda <- lengths(skt)
  mu <- jack:::removeTrailingZeros(vapply(skt, function(row) {
    sum(is.na(row))
  }, integer(1L)))
  lambdas <- lambdais(lambda, mu)
  r <- length(lambdas)
  Reduce(`*`, lapply(2L:r, function(i) {
    psi_lambda_mu(lambdas[[i]], lambdas[[i-1L]])
  }))

}


lambda <- c(4,3)
mu <- c(2,2)
n <- 7
sktx <- syt::all_ssSkewTableaux(lambda, mu, n)
length(sktx)
phiT(sktx[[3]])
lapply(sktx[1:36], phiT)

p <- Qzero()
lones <- lapply(1:n, Qlone)
for(skt in sktx) {
  ct <- content(skt)
  q <- Qone()
  for(i in seq_along(ct)) {
    q <- q * lones[[i]]^(ct[i])
  }
  p <- p + phiT(skt)*q
}
p
SkewSchurPol(n, lambda, mu)
substituteParameters(p, 0)
evalQspray(jack:::b(mu), 0) / evalQspray(jack:::b(lambda), 0)
jack:::b(lambda) / jack:::b(mu) * p
jack:::b(mu) / jack:::b(lambda) * p

