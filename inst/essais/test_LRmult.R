# https://math.stackexchange.com/q/3012744/38217
library(jack)
library(syt)

mu <- c(3, 2, 2, 1)
nu <- c(4, 3, 2, 1)
LR <- LRmult(mu, nu, output = "list")
counts <- vapply(LR$lambda, count_sytx, numeric(1L))
rhs <- sum(LR$coeff * counts)
lhs <- count_sytx(mu) * count_sytx(nu) * choose(sum(mu) + sum(nu), sum(mu))
lhs == rhs
