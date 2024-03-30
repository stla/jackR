library(jack)
library(syt)
library(qspray)

lambda <- c(3, 1)
ssytx <- all_ssytx(lambda, 4)

# ! quand |lambda| > n, length(all_ssytx(lambda, n)) > count_ssytx(lambda, n)

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

SchurPol(4, lambda)


