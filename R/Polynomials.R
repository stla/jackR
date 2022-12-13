#' @importFrom qspray as.qspray
#' @importFrom spray zero one lone
NULL

#' @importFrom mvp mvp
as_mvp_spray <- function(s) {
  powers <- s[["index"]]
  m <- nrow(powers)
  n <- ncol(powers)
  vars <- replicate(m, paste0("x_", 1L:n), simplify = FALSE)
  powers <- lapply(1L:m, function(i) powers[i, ])
  mvp(vars, powers, s[["value"]])
}

JackPolNaive <- function(n, lambda, alpha, basis = "canonical"){
  stopifnot(isPositiveInteger(n), alpha >= 0, isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  gmp <- is.bigq(alpha)
  if(length(lambda) == 0L) {
    if(basis == "canonical") {
      return(if(gmp) as.qspray(1) else one(n))
    } else {
      return("M_()")
    }
  }
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > n) return(if(gmp) as.qspray(0) else zero(n))
  lambda00 <- integer(sum(lambda))
  lambda00[seq_along(lambda)] <- lambda
  mus <- dominatedPartitions(lambda)
  if(gmp) {
    coefs <- JackCoefficientsQ(sum(lambda), alpha, until = lambda)
  } else {
    coefs <- JackCoefficientsNum(sum(lambda), alpha, until = lambda)
  }
  coefs <- coefs[toString(lambda00), ]
  if(basis == "canonical") {
    if(gmp) {
      out <- as.qspray(0)
      for(i in 1L:ncol(mus)){
        mu <- mus[, i]
        l <- sum(mu > 0L)
        if(l <= n) {
          toAdd <- MSFpoly(n, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    } else {
      out <- zero(n)
      for(i in 1L:ncol(mus)) {
        mu <- mus[, i]
        l <- sum(mu > 0L)
        if(l <= n) {
          toAdd <- MSFspray(n, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
      as_mvp_spray(out)
    }
    out
  } else {
    vars <- apply(mus, 2L, function(mu) {
      paste0("M_(", paste0(mu[mu > 0L], collapse = ","), ")")
    })
    rowIdx <- match(toString(lambda00), names(coefs))
    coefs <- coefs[-seq_len(rowIdx-1L)]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

JackPolDK <- function(n, lambda, alpha) {
  stopifnot(isPositiveInteger(n), alpha >= 0, isPartition(lambda))
  jac <- function(m, k, mu, nu, beta) {
    if(length(nu) == 0L || nu[1L] == 0L || m == 0L) return(one(n))
    if(length(nu) > m && nu[m+1L] > 0L) return(zero(n))
    if(m == 1L) return(prod(alpha*seq_len(nu[1L]-1L)+1) * lone(1, n)^nu[1L])
    if(k == 0L && !is.na(s <- S[[.N(lambda,nu),m]])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0L, nu, nu, 1) * beta * lone(m, n)^(sum(mu)-sum(nu))
    while(length(nu) >= i && nu[i] > 0L){
      if(length(nu) == i && nu[i] > 0L || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1L
        gamma <- beta * .betaratio(mu, nu, i, alpha)
        if(nu[i] > 1L){
          s <- s + jac(m, i, mu, .nu, gamma)
        }else{
          s <- s + jac(m-1L, 0L, .nu, .nu, 1) * gamma *
            lone(m, n)^(sum(mu)-sum(.nu))
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[[.N(lambda,nu),m]] <- s
    return(s)
  }
  S <- as.list(rep(NA, .N(lambda,lambda) * n))
  dim(S) <- c(.N(lambda,lambda), n)
  as_mvp_spray(jac(n, 0L, lambda, lambda, 1))
}

#' @importFrom qspray as.qspray qsprayMaker qlone
#' @importFrom gmp as.bigq
#' @noRd
JackPolDK_gmp <- function(n, lambda, alpha) {
  stopifnot(isPositiveInteger(n), alpha >= 0, isPartition(lambda))
  jac <- function(m, k, mu, nu, beta) {
    if(length(nu) == 0L || nu[1L] == 0L || m == 0L) {
      return(as.qspray(1))
    }
    if(length(nu) > m && nu[m+1L] > 0L) return(as.qspray(0))
    if(m == 1L) {
      return(
        prod(alpha * seq_len(nu[1L]-1L) + 1L) * qlone(1)^nu[1L]
      )
    }
    if(k == 0L && inherits(s <- S[[.N(lambda, nu), m]], "qspray")) return(s)
    i <- max(1L, k)
    s <- jac(m-1L, 0L, nu, nu, oneq) * as.qspray(beta) *
      qsprayMaker(
        coeffs = "1",
        powers = list(c(rep(0L, m-1L), sum(mu)-sum(nu)))
      )
    while(length(nu) >= i && nu[i] > 0L) {
      if(length(nu) == i && nu[i] > 0L || nu[i] > nu[i+1L]) {
        .nu <- nu; .nu[i] <- nu[i]-1L
        gamma <- beta * .betaratio(mu, nu, i, alpha)
        if(nu[i] > 1L) {
          s <- s + jac(m, i, mu, .nu, gamma)
        } else {
          s <- s + jac(m-1L, 0L, .nu, .nu, oneq) * as.qspray(gamma) *
            qsprayMaker(
              coeffs = "1",
              powers = list(c(rep(0L, m-1L), sum(mu)-sum(.nu)))
            )
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[[.N(lambda, nu), m]] <- s
    return(s)
  }
  Nlambdalambda <- .N(lambda, lambda)
  S <- as.list(rep(NA, Nlambdalambda * n))
  dim(S) <- c(Nlambdalambda, n)
  oneq <- as.bigq(1L)
  jac(n, 0L, lambda, lambda, oneq)
}

#' Jack polynomial
#'
#' Returns the Jack polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#' integers
#' @param alpha parameter of the Jack polynomial, a positive number, possibly
#'   a \code{\link[gmp]{bigq}} rational number
#' @param algorithm the algorithm used, either \code{"DK"} or \code{"naive"}
#' @param basis the polynomial basis for \code{algorithm = "naive"},
#' either \code{"canonical"} or \code{"MSF"} (monomial symmetric functions);
#' for \code{algorithm = "DK"} the canonical basis is always used and
#' this parameter is ignored
#'
#' @return A \code{mvp} multivariate polynomial (see \link[mvp]{mvp-package}),
#'  or a \code{\link[gmpoly]{gmpoly}} multivariate polynomial if \code{alpha}
#'  is a \code{bigq} rational number and \code{algorithm = "DK"}, or a
#'  character string if \code{basis = "MSF"}.
#' @importFrom gmp is.bigq
#' @export
#'
#' @examples JackPol(3, lambda = c(3,1), alpha = gmp::as.bigq(2,3),
#'                   algorithm = "naive")
#' JackPol(3, lambda = c(3,1), alpha = 2/3, algorithm = "DK")
#' JackPol(3, lambda = c(3,1), alpha = gmp::as.bigq(2,3), algorithm = "DK")
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
      JackPolDK_gmp(n, lambda, alpha)
    }else{
      JackPolDK(n, lambda, alpha)
    }
  }else{
    JackPolNaive(n, lambda, alpha, basis)
  }
}


ZonalPolNaive <- function(m, lambda, basis = "canonical", exact = TRUE){
  stopifnot(isPositiveInteger(m), isPartition(lambda))
  basis <- match.arg(basis, c("canonical", "MSF"))
  if(length(lambda) == 0L){
    if(basis == "canonical"){
      return(if(exact) as.qspray(1) else one(m))
    }else{
      return("M_()")
    }
  }
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(if(exact) as.qspray(0) else zero(m))
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
    if(exact){
      out <- as.qspray(0)
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
      out
    }else{
      out <- zero(m)
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFspray(m, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }
    as_mvp_spray(out)
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
  jack <- JackPolDK(m, lambda, alpha = 2)
  jlambda <- sum(logHookLengths(lambda, alpha = 2))
  n <- sum(lambda)
  exp(n*log(2) + lfactorial(n) - jlambda) * jack
}

#' @importFrom gmp as.bigq factorialZ
#' @importFrom gmpoly gmpolyConstant
#' @noRd
ZonalPolDK_gmp <- function(m, lambda){
  twoq <- as.bigq(2)
  jack <- JackPolDK_gmp(m, lambda, alpha = twoq)
  jlambda <- prod(hookLengths_gmp(lambda, alpha = twoq))
  n <- sum(lambda)
  (twoq^n * factorialZ(n) / jlambda) * jack
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
#' @param exact logical, whether to get rational coefficients
#'
#' @return A \code{mvp} multivariate polynomial (see \link[mvp]{mvp-package}),
#'  or a \code{\link[gmpoly]{gmpoly}} multivariate polynomial if
#'  \code{exact = TRUE} and \code{algorithm = "DK"}, or a
#'  character string if \code{basis = "MSF"}.
#'
#' @importFrom mvp constant mvp
#' @export
#'
#' @examples ZonalPol(3, lambda = c(3,1), algorithm = "naive")
#' ZonalPol(3, lambda = c(3,1), algorithm = "DK")
#' ZonalPol(3, lambda = c(3,1), algorithm = "DK", exact = FALSE)
#' ZonalPol(3, lambda = c(3,1), algorithm = "naive", basis = "MSF")
ZonalPol <- function(n, lambda, algorithm = "DK", basis = "canonical",
                     exact = TRUE){
  algo <- match.arg(algorithm, c("DK", "naive"))
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    if(exact){
      ZonalPolDK_gmp(n, lambda)
    }else{
      ZonalPolDK(n, lambda)
    }
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
      return(if(exact) as.qspray(1) else one(m))
    }else{
      return("M_()")
    }
  }
  lambda <- lambda[lambda > 0L]
  if(length(lambda) > m) return(if(exact) as.qspray(0) else zero(m))
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
    if(exact){
      out <- as.qspray(0)
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFpoly(m, mu)
          if(coefs[toString(mu)] != "1")
            toAdd <- toAdd * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }else{
      out <- zero(m)
      for(i in 1L:ncol(mus)){
        mu <- mus[,i]
        l <- sum(mu > 0L)
        if(l <= m){
          toAdd <- MSFspray(m, mu) * coefs[toString(mu)]
          out <- out + toAdd
        }
      }
    }
    as_mvp_spray(out)
  }else{
    vars <- apply(mus, 2L, function(mu){
      paste0("M_(", paste0(mu[mu>0L], collapse = ","), ")")
    })
    coefs <- coefs[coefs != "0"]
    coefs <- ifelse(coefs == "1", "", paste0(coefs, " "))
    paste0(coefs, vars, collapse = " + ")
  }
}

#' @importFrom mvp constant mvp
#' @noRd
SchurPolDK <- function(n, lambda){
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  sch <- function(m, k, nu){
    if(length(nu) == 0L || nu[1L] == 0L || m == 0L){
      return(one(n))
    }
    if(length(nu) > m && nu[m+1L] > 0L) return(zero(n))
    if(m == 1L) return(lone(1, n)^nu[1L])
    if(!is.na(s <- S[[.N(lambda, nu), m]])) return(s)
    s <- sch(m-1L, 1L, nu)
    i <- k
    while(length(nu) >= i && nu[i] > 0L){
      if(length(nu) == i || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1L
        if(nu[i] > 1L){
          s <- s + lone(m, n) * sch(m, i, .nu)
        }else{
          s <- s + lone(m, n) * sch(m-1L, 1L, .nu)
        }
      }
      i <- i + 1L
    }
    if(k == 1L) S[[.N(lambda, nu), m]] <- s
    return(s)
  }
  Nlambdalambda <- .N(lambda,lambda)
  S <- as.list(rep(NA, Nlambdalambda*n))
  dim(S) <- c(Nlambdalambda, n)
  sch(n, 1L, as.integer(lambda))
}

#' @importFrom gmpoly gmpolyConstant gmpoly gmpolyGrow
SchurPolDK_gmp <- function(n, lambda){
  stopifnot(isPositiveInteger(n), isPartition(lambda))
  sch <- function(m, k, nu){
    if(length(nu) == 0L || nu[1L] == 0L || m == 0L){
      return(as.qspray(1))
    }
    if(length(nu) > m && nu[m+1L] > 0L) return(as.qspray(0))
    if(m == 1L) return(qlone(1)^nu[1L])
    if(inherits(s <- S[[.N(lambda, nu), m]], "qspray")) return(s)
    s <- sch(m-1L, 1L, nu)
    i <- k
    while(length(nu) >= i && nu[i] > 0L){
      if(length(nu) == i || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1L
        if(nu[i] > 1L){
          s <- s + x[[m]] * sch(m, i, .nu)
        }else{
          s <- s + x[[m]] * sch(m-1L, 1L, .nu)
        }
      }
      i <- i + 1L
    }
    if(k == 1L) S[[.N(lambda, nu), m]] <- s
    return(s)
  }
  Nlambdalambda <- .N(lambda,lambda)
  S <- as.list(rep(NA, Nlambdalambda*n))
  dim(S) <- c(Nlambdalambda, n)
  oneq <- as.bigq(1L)
  x <- lapply(1L:n, function(m){
    qsprayMaker(
      coeffs = "1",
      powers = list(c(rep(0L, m-1L), 1L))
    )
  })
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
#' @param exact logical, whether to use exact arithmetic
#'
#' @return A \code{mvp} multivariate polynomial (see \link[mvp]{mvp-package}),
#'  or a \code{\link[gmpoly]{gmpoly}} multivariate polynomial if
#'  \code{exact = TRUE} and \code{algorithm = "DK"}, or a
#'  character string if \code{basis = "MSF"}.
#'
#' @export
#'
#' @examples SchurPol(3, lambda = c(3,1), algorithm = "naive")
#' SchurPol(3, lambda = c(3,1), algorithm = "DK")
#' SchurPol(3, lambda = c(3,1), algorithm = "DK", exact = FALSE)
#' SchurPol(3, lambda = c(3,1), algorithm = "naive", basis = "MSF")
SchurPol <- function(n, lambda, algorithm = "DK", basis = "canonical",
                     exact = TRUE){
  algo <- match.arg(algorithm, c("DK", "naive"))
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    if(exact){
      SchurPolDK_gmp(n, lambda)
    }else{
      SchurPolDK(n, lambda)
    }
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
  jack <- JackPolDK(m, lambda, alpha = 1/2)
  jlambda <- sum(logHookLengths(lambda, alpha = 1/2))
  n <- sum(lambda)
  exp(-n*log(2) + lfactorial(n) - jlambda) * jack
}

#' @importFrom gmp as.bigq factorialZ
#' @noRd
ZonalQPolDK_gmp <- function(m, lambda){
  onehalfq <- as.bigq(1L, 2L)
  jack <- JackPolDK_gmp(m, lambda, alpha = onehalfq)
  jlambda <- prod(hookLengths_gmp(lambda, alpha = onehalfq))
  n <- sum(lambda)
  as.qspray(onehalfq^n * factorialZ(n) / jlambda) * jack
}

#' Quaternionic zonal polynomial
#'
#' Returns the quaternionic (or symplectic) zonal polynomial.
#'
#' @param n number of variables, a positive integer
#' @param lambda an integer partition, given as a vector of decreasing
#'   integers
#' @param algorithm the algorithm used, either \code{"DK"} or \code{"naive"}
#' @param basis the polynomial basis for \code{algorithm = "naive"},
#' either \code{"canonical"} or \code{"MSF"} (monomial symmetric functions);
#' for \code{algorithm = "DK"} the canonical basis is always used and
#' this parameter is ignored
#' @param exact logical, whether to get rational coefficients
#'
#' @return A \code{mvp} multivariate polynomial (see \link[mvp]{mvp-package}),
#'  or a \code{\link[gmpoly]{gmpoly}} multivariate polynomial if
#'  \code{exact = TRUE} and \code{algorithm = "DK"}, or a
#'  character string if \code{basis = "MSF"}.
#'
#' @export
#'
#' @examples ZonalQPol(3, lambda = c(3,1), algorithm = "naive")
#' ZonalQPol(3, lambda = c(3,1), algorithm = "DK")
#' ZonalQPol(3, lambda = c(3,1), algorithm = "DK", exact = FALSE)
#' ZonalQPol(3, lambda = c(3,1), algorithm = "naive", basis = "MSF")
ZonalQPol <- function(n, lambda, algorithm = "DK", basis = "canonical",
                     exact = TRUE){
  algo <- match.arg(algorithm, c("DK", "naive"))
  lambda <- as.integer(lambda)
  if(algo == "DK"){
    if(exact){
      ZonalQPolDK_gmp(n, lambda)
    }else{
      ZonalQPolDK(n, lambda)
    }
  }else{
    ZonalQPolNaive(n, lambda, basis, exact)
  }
}
