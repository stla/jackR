library(jack)
library(syt)
library(qspray)

lambda <- c(4, 3, 2, 1)
mu <- c(2, 1)
ssytx <- all_ssSkewTableaux(lambda, mu, 4)

# ! quand |lambda| > n, length(all_ssytx(lambda, n)) > count_ssytx(lambda, n)
# -> FIXED

x <- ssytx[[19]]

wt <- function(ssyt) {
  ssyt <- unlist(ssyt)
  vapply(1:4, function(k) {
    length(which(ssyt == k))
  }, integer(1L))
}

qlones <- lapply(1:4, qlone)

monomial <- function(ssyt) {
  powers <- wt(ssyt)
  Reduce(`*`, lapply(1:4, function(k) qlones[[k]]^powers[k]))
}

monomials <- lapply(ssytx, monomial)
Reduce(`+`, monomials)

SkewSchurPol(4, lambda, mu)


