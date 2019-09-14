JackPolNaive <- function(n, lambda, alpha, basis = "canonical"){
  stopifnot(isPositiveInteger(n), alpha >= 0, isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  if(length(lambda) == 0L){
    if(basis == "canonical"){
      return(constant(1))
    }else{
      return("M_()")
    }
  }
  gmp <- is.bigq(alpha)
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > n) return(constant(0))
  lambda00 <- integer(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  if(gmp){
    coefs <- JackCoefficientsQ(sum(lambda), alpha, until = lambda)
  }else{
    coefs <- JackCoefficientsNum(sum(lambda), alpha, until = lambda)
  }
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    if(gmp){
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= n){
          toAdd <- MSFpoly(n, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * mvp(coefs[toString(mu)], 1, 1)
          out <- out + toAdd
        }
      }
    }else{
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= n){
          toAdd <- MSFpoly(n, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0L], collapse = ","), ")")
    })
    rowIdx <- match(toString(lambda00), names(coefs))
    coefs <- coefs[-seq_len(rowIdx-1L)]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

JackPolDK <- function(n, lambda, alpha){
  stopifnot(isPositiveInteger(n), alpha >= 0, isPartition(lambda))
  jac <- function(m, k, mu, nu, beta){
    if(length(nu) == 0L || nu[1L] == 0L || m == 0L) return(constant(1))
    if(length(nu) > m && nu[m+1L] > 0L) return(constant(0))
    if(m == 1L) return(mvp("x_1", nu[1L], prod(alpha*seq_len(nu[1L]-1L)+1)))
    if(k == 0L && !is.na(s <- S[[.N(lambda,nu),m]])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0L, nu, nu, 1) * beta *
      mvp(x[m], sum(mu)-sum(nu), 1)
    while(length(nu) >= i && nu[i] > 0L){
      if(length(nu) == i && nu[i] > 0L || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1L
        gamma <- beta * .betaratio(mu, nu, i, alpha)
        if(nu[i] > 1L){
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

#' Jack polynomial
#'
#' Returns the Jack polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param alpha parameter of the Jack polynomial, always a positive number
#' for \code{algorithm = "DK"}, a positive number or a positive \code{bigq}
#' rational number for \code{algorithm = "naive"}
#' @param algorithm the algorithm used, either \code{"DK"} or \code{"naive"}
#' @param basis the polynomial basis for \code{algorithm = "naive"},
#' either \code{"canonical"} or \code{"MSF"} (monomial symmetric functions);
#' for \code{algorithm = "DK"} the canonical basis is always used and
#' this parameter is ignored
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}) or a
#' character string if \code{basis = "MSF"}.
#' @importFrom mvp constant mvp
#' @importFrom gmp is.bigq
#' @export
#'
#' @examples JackPol(3, lambda = c(3,1), alpha = gmp::as.bigq(2,3),
#'                   algorithm = "naive")
#' JackPol(3, lambda = c(3,1), alpha = 2/3, algorithm = "DK")
#' JackPol(3, lambda = c(3,1), alpha= gmp::as.bigq(2,3),
#'         algorithm = "naive", basis = "MSF")
JackPol <- function(n, lambda, alpha, algorithm = "DK",
                    basis = "canonical"){
  algo <- match.arg(algorithm, c("DK", "naive"))
  stopifnot(
    is.numeric(alpha) || is.bigq(alpha),
    length(alpha) == 1L
  )
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    if(is.bigq(alpha)){
      stop("Algorithm `DK` is not implemented for rational `alpha`")
    }
    JackPolDK(n, lambda, alpha)
  }else{
    JackPolNaive(n, lambda, alpha, basis)
  }
}


ZonalPolNaive <- function(m, lambda, basis = "canonical", exact = TRUE){
  stopifnot(isPositiveInteger(m), isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  if(length(lambda) == 0L){
    if(basis == "canonical"){
      return(constant(1))
    }else{
      return("M_()")
    }
  }
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- numeric(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  if(exact){
    coefs <- zonalCoefficientsQ(sum(lambda), until = lambda)
  }else{
    coefs <- zonalCoefficientsNum(sum(lambda), until = lambda)
  }
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    if(exact){
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * mvp(coefs[toString(mu)], 1, 1)
          out <- out + toAdd
        }
      }
    }else{
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0L], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

ZonalPolDK <- function(m, lambda){
  jack <- JackPolDK(m, lambda, alpha= 2)
  jlambda <- sum(logHookLengths(lambda, alpha = 2))
  n <- sum(lambda)
  exp(n*log(2) + lfactorial(n) - jlambda) * jack
}

#' Zonal polynomial
#'
#' Returns the zonal polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param algorithm the algorithm used, either \code{"DK"} or \code{"naive"}
#' @param basis the polynomial basis for \code{algorithm = "naive"},
#' either \code{"canonical"} or \code{"MSF"} (monomial symmetric functions);
#' for \code{algorithm = "DK"} the canonical basis is always used and
#' this parameter is ignored
#' @param exact logical, whether to get rational coefficients when using
#' \code{algorithm = "naive"}; ignored if \code{algorithm = "DK"}
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}) or a
#' character string if \code{basis = "MSF"}.
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples ZonalPol(3, lambda = c(3,1), algorithm = "naive")
#' ZonalPol(3, lambda = c(3,1), algorithm = "DK")
#' ZonalPol(3, lambda = c(3,1), algorithm = "naive", basis = "MSF")
ZonalPol <- function(n, lambda, algorithm = "DK", basis = "canonical",
                     exact = TRUE){
  algo <- match.arg(algorithm, c("DK", "naive"))
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    ZonalPolDK(n, lambda)
  }else{
    ZonalPolNaive(n, lambda, basis, exact)
  }
}


SchurPolNaive <- function(m, lambda, basis = "canonical",
                          exact = TRUE){
  stopifnot(isPositiveInteger(m), isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  if(length(lambda) == 0L){
    if(basis == "canonical"){
      return(constant(1))
    }else{
      return("M_()")
    }
  }
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- integer(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  if(exact){
    coefs <- SchurCoefficientsQ(sum(lambda), until = lambda)
  }else{
    coefs <- SchurCoefficientsNum(sum(lambda), until = lambda)
  }
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    if(exact){
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * mvp(coefs[toString(mu)], 1, 1)
          out <- out + toAdd
        }
      }
    }else{
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0L], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

SchurPolDK <- function(n, lambda){
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  sch <- function(m, k, nu){
    if(length(nu) == 0L || nu[1L]==0L || m == 0L){
      return(constant(1))
    }
    if(length(nu) > m && nu[m+1L] > 0L) return(constant(0))
    if(m == 1L) return(mvp(x[1L], nu[1L], 1))
    if(!is.na(s <- S[[.N(lambda,nu),m]])) return(s)
    s <- sch(m-1L, 1L, nu)
    i <- k
    while(length(nu) >= i && nu[i] > 0L){
      if(length(nu) == i || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1L
        if(nu[i] > 1L){
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
  sch(n, 1L, as.integer(lambda))
}

#' Schur polynomial
#'
#' Returns the Schur polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param algorithm the algorithm used, either \code{"DK"} or \code{"naive"}
#' @param basis the polynomial basis for \code{algorithm = "naive"},
#' either \code{"canonical"} or \code{"MSF"} (monomial symmetric functions);
#' for \code{algorithm = "DK"} the canonical basis is always used and
#' this parameter is ignored
#' @param exact logical, whether to get rational coefficients when using
#' \code{algorithm = "naive"}; ignored if \code{algorithm = "DK"}
#'
#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}) or a
#' character string if \code{basis = "MSF"}.
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples SchurPol(3, lambda = c(3,1), algorithm = "naive")
#' SchurPol(3, lambda = c(3,1), algorithm = "DK")
#' SchurPol(3, lambda = c(3,1), algorithm = "naive", basis = "MSF")
SchurPol <- function(n, lambda, algorithm = "DK", basis = "canonical",
                     exact = TRUE){
  algo <- match.arg(algorithm, c("DK", "naive"))
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    SchurPolDK(n, lambda)
  }else{
    SchurPolNaive(n, lambda, basis, exact)
  }
}


ZonalQPolNaive <- function(m, lambda, basis = "canonical", exact = TRUE){
  stopifnot(isPositiveInteger(m), isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  if(length(lambda) == 0L){
    if(basis == "canonical"){
      return(constant(1))
    }else{
      return("M_()")
    }
  }
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(constant(0))
  lambda00 <- numeric(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  if(exact){
    coefs <- zonalQCoefficientsQ(sum(lambda), until = lambda)
  }else{
    coefs <- zonalQCoefficientsNum(sum(lambda), until = lambda)
  }
  coefs <- coefs[toString(lambda00),]
  if(basis == "canonical"){
    out <- constant(0)
    if(exact){
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * mvp(coefs[toString(mu)], 1, 1)
          out <- out + toAdd
        }
      }
    }else{
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }
    out
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0L], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

ZonalQPolDK <- function(m, lambda){
  jack <- JackPolDK(m, lambda, alpha= 1/2)
  jlambda <- sum(logHookLengths(lambda, alpha = 1/2))
  n <- sum(lambda)
  exp(-n*log(2) + lfactorial(n) - jlambda) * jack
}

#' Quaternionic zonal polynomial
#'
#' Returns the quaternionic (or symplectic) zonal polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param algorithm the algorithm used, either \code{"DK"} or \code{"naive"}
#' @param basis the polynomial basis for \code{algorithm = "naive"},
#' either \code{"canonical"} or \code{"MSF"} (monomial symmetric functions);
#' for \code{algorithm = "DK"} the canonical basis is always used and
#' this parameter is ignored
#' @param exact logical, whether to get rational coefficients when using
#' \code{algorithm = "naive"}; ignored if \code{algorithm = "DK"}

#' @return A polynomial (\code{mvp} object; see \link[mvp]{mvp-package}) or a
#' character string if \code{basis = "MSF"}.
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples ZonalQPol(3, lambda = c(3,1), algorithm = "naive")
#' ZonalQPol(3, lambda = c(3,1), algorithm = "DK")
#' ZonalQPol(3, lambda = c(3,1), algorithm = "naive", basis = "MSF")
ZonalQPol <- function(n, lambda, algorithm = "DK", basis = "canonical",
                     exact = TRUE){
  algo <- match.arg(algorithm, c("DK", "naive"))
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    ZonalQPolDK(n, lambda)
  }else{
    ZonalQPolNaive(n, lambda, basis, exact)
  }
}
