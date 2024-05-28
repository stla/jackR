R <- function(a, i, j) {
  aa <- a
  aa[i] <- a[i] + 1L
  aa[j] <- a[j] - 1L
  aa
}

library(jack)

q0 <- function(n, r, t) {
  if(r < 0) {
    qzero()
  } else if(r == 0) {
    qone()
  } else {
    (1-t)*substituteParameters(HallLittlewoodPol(n, r), t)
  }
}

qa <- function(n, alpha, t) {
  Reduce(`*`, lapply(alpha, function(a) q0(n, a, t)))
}

lambda <- c(1,1,1)
prod <- qone()
t <- -1
for(i in 1:2) {
  for(j in (i+1):3) {
    qq <- qa(3, lambda, t)
    p <- qa(3, R(lambda, i, j), t)
    prod <- prod * (qq - p) / (qq - t*p)
  }
}
prod
substituteParameters(HallLittlewoodPol(3, lambda, "Q"), t)

prod <- qone()
n <- 4
t <- 2
lambda <- c(n, 0, 0, 0)

summ <- qzero()
for(i in 0:(n-1)) {
  summ <- summ + (-t)^i*SchurPol(n, c(n-i, rep(1,i)))
}
summ
substituteParameters(HallLittlewoodPol(n, lambda, "P"), t)

n <- 4
lambda <- c(n, 0, 0, 0)
prod <- qone()
for(j in 2:n) {
  kappa <- R(lambda, j, 1)
  s <- SchurPol(n, lambda)
  ss <- SchurPol(n, kappa)
  prod <- (s - t*ss)
  lambda <- kappa
}
prod

n <- 4
lambda <- c(n, 0, 0, 0)
kappas <- vector("list", n)
kappas[[1]] <- lambda
for(j in 2:n) {
  kappas[[j]] <- R(kappas[[j-1]], j, 1)
}
kappas
p0 <- SchurPol(n, kappas[[1]])
p1 <- p0 - t * SchurPol(n, kappas[[2]])
# p1 = kappa1 - t*kappa2
p2 <- p1 - t * (SchurPol(n, kappas[[2]]) - t*SchurPol(n, kappas[[3]]))
# p2 = kappa1 - t*kappa2 - t*(kappa2 - t*kappa3)
p3 <- p2 - t * (SchurPol(n, kappas[[2]]) - t*SchurPol(n, kappas[[3]]) - t*(SchurPol(n, kappas[[3]]) - t*SchurPol(n, kappas[[4]])))
p3 <- p2 - t * (SchurPol(n, kappas[[4]]) - SchurPol(n, kappas[[3]]))# - t*(SchurPol(n, kappas[[3]]) - t*(SchurPol(n, kappas[[4]]) - t*SchurPol(n, kappas[[4]]))))



ss <- SchurPol(n, kappa)
(s - t*ss)


lambdas <- partitions::parts(3)
kappas <- jack:::Columns(lambdas)

mu1 <- kappas[[1]]
mu2 <- kappas[[1]]
t <- qlone(1)
prod <- qone()
for(i in 1:2) {
  for(j in (i+1):3) {
    x <- c(crossprod(R(kappas[[i]], i, j), R(kappas[[j]], i, j)))
    prod <- prod / (1- x*t)
  }
}
prod
