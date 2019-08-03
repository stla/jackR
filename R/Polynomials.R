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
    rowIdx <- match(toString(lambda00), names(coefs))
    coefs <- coefs[-seq_len(rowIdx-1L)]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}


JackPolDK <- function(n, lambda, alpha){
  jac <- function(m, k, mu, nu, beta){
    if(length(nu) == 0L || nu[1L]==0 || m == 0L) return(constant(1))
    if(length(nu) > m && nu[m+1L] > 0) return(constant(0))
    if(m == 1L) return(mvp("x_1", nu[1L], prod(alpha*seq_len(nu[1L]-1)+1)))
    if(k == 0L && !is.na(s <- S[[.N(lambda,nu),m]])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0L, nu, nu, 1) * beta *
      mvp(x[m], sum(mu)-sum(nu), 1)
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i && nu[i] > 0 || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        gamma <- beta * .betaratio(mu, nu, i, alpha)
        if(nu[i] > 1){
          s <- s + jac(m, i, mu, .nu, gamma)
        }else{
          s <- s + jac(m-1L, 0L, .nu, .nu, 1) * gamma *
            mvp(x[m], sum(mu)-sum(.nu), 1)
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[[.N(lambda,nu),m]] <- s
    return(s)
  }
  S <- as.list(rep(NA, .N(lambda,lambda) * n))
  dim(S) <- c(.N(lambda,lambda), n)
  x <- paste0("x_", 1L:n)
  jac(n, 0L, lambda, lambda, 1)
}

#' Zonal polynomial
#'
#' Returns the zonal polynomial.
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

ZonalPolDK <- function(n, lambda){
  jack <- JackPolDK(n, lambda, alpha= 2)
  jlambda <- sum(logHookLengths(lambda, alpha = 2))
  n <- sum(lambda)
  exp(n*log(2) + lfactorial(n) - jlambda) * jack
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

SchurPolDK <- function(n, lambda){
  sch <- function(m, k, nu){
    if(length(nu) == 0L || nu[1L]==0 || m == 0L){
      return(constant(1))
    }
    if(length(nu) > m && nu[m+1L] > 0) return(constant(0))
    if(m == 1L) return(mvp(x[1L], nu[1L], 1))
    if(!is.na(s <- S[[.N(lambda,nu),m]])) return(s)
    s <- sch(m-1L, 1L, nu)
    i <- k
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        if(nu[i] > 1){
          s <- s + mvp(x[m], 1, 1) * sch(m, i, .nu)
        }else{
          s <- s + mvp(x[m], 1, 1) * sch(m-1L, 1L, .nu)
        }
      }
      i <- i + 1L
    }
    if(k == 1L) S[[.N(lambda,nu),m]] <- s
    return(s)
  }
  x <- paste0("x_", 1L:n)
  S <- as.list(rep(NA, .N(lambda,lambda)*n))
  dim(S) <- c(.N(lambda,lambda), n)
  sch(n, 1L, lambda)
}

#' Quaternionic zonal polynomial
#'
#' Returns the quaternionic zonal polynomial.
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
#' @examples ZonalQPol(3, lambda = c(3,1))
#' ZonalQPol(3, lambda = c(3,1), basis = "MSF")
ZonalQPol <- function(m, lambda, basis = "canonical"){
  stopifnot(floor(m) == m, isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  lambda <- lambda[lambda>0]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- numeric(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  coefs <- zonalQCoefficientsQ(sum(lambda), until = lambda)
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

ZonalQPolDK <- function(n, lambda){
  jack <- JackPolDK(n, lambda, alpha= 1/2)
  jlambda <- sum(logHookLengths(lambda, alpha = 1/2))
  n <- sum(lambda)
  exp(-n*log(2) + lfactorial(n) - jlambda) * jack
}
