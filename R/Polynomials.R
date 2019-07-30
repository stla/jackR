#' Jack polynomial
#'
#' Returns the Jack polynomial.
#'
#' @param m number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param alpha parameter of the Jack polynomial,
#' a positive \code{bigq} rational number
#' @param basis the polynomial basis, either \code{"canonical"} or
#' \code{"MSF"} (monomial symmetric functions)
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}).
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples JackPol(3, lambda = c(3,1), alpha = gmp::as.bigq(2,3))
#' JackPol(3, lambda = c(3,1), alpha= gmp::as.bigq(2,3), basis = "MSF")
JackPol <- function(m, lambda, alpha, basis = "canonical"){
  stopifnot(floor(m) == m, alpha >= 0, isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  lambda <- lambda[lambda>0]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- numeric(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  coefs <- JackCoefficientsQ(sum(lambda), alpha, until = lambda)
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    for(i in 1:ncol(mus)){
      l <- sum(mus[,i] > 0)
      if(l <= m){
        toAdd <- MSFpoly(m, mus[,i])
        if(coefs[toString(mus[,i])] != "1")
          toAdd <- toAdd * mvp(coefs[toString(mus[,i])], 1, 1)
        out <- out + toAdd
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

#' Zonal polynomial
#'
#' Returns the Zonal polynomial.
#'
#' @param m number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param basis the polynomial basis, either \code{"canonical"} or
#' \code{"MSF"} (monomial symmetric functions)
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}).
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples ZonalPol(3, lambda = c(3,1))
#' ZonalPol(3, lambda = c(3,1), basis = "MSF")
ZonalPol <- function(m, lambda, basis = "canonical"){
  stopifnot(floor(m) == m, isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  lambda <- lambda[lambda>0]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- numeric(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  coefs <- zonalCoefficientsQ(sum(lambda), until = lambda)
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    for(i in 1:ncol(mus)){
      l <- sum(mus[,i] > 0)
      if(l <= m){
        toAdd <- MSFpoly(m, mus[,i])
        if(coefs[toString(mus[,i])] != "1")
          toAdd <- toAdd * mvp(coefs[toString(mus[,i])], 1, 1)
        out <- out + toAdd
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

#' Schur polynomial
#'
#' Returns the Schur polynomial.
#'
#' @param m number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param basis the polynomial basis, either \code{"canonical"} or
#' \code{"MSF"} (monomial symmetric functions)
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}).
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples SchurPol(3, lambda = c(3,1))
#' SchurPol(3, lambda = c(3,1), basis = "MSF")
SchurPol <- function(m, lambda, basis = "canonical"){
  stopifnot(floor(m) == m, isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  lambda <- lambda[lambda>0]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- numeric(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  coefs <- SchurCoefficientsQ(sum(lambda), until = lambda)
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    for(i in 1:ncol(mus)){
      l <- sum(mus[,i] > 0)
      if(l <= m){
        toAdd <- MSFpoly(m, mus[,i])
        if(coefs[toString(mus[,i])] != "1")
          toAdd <- toAdd * mvp(coefs[toString(mus[,i])], 1, 1)
        out <- out + toAdd
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

