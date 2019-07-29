#' Jack polynomial
#'
#' Evaluates the Jack polynomials.
#'
#' @param x numeric vector
#' @param lambda integer partition, a vector of decreasing integers
#' @param alpha positive number
#'
#' @return A number.
#' @export
Jack <- function(x, lambda, alpha){
  stopifnot(alpha > 0)
  if(any(diff(lambda) > 0) || any(floor(lambda) != lambda)){
    stop("`lambda` is not a valid integer partition")
  }
  jac <- function(m, k, mu, nu){
    if(all(nu==0) || length(nu) == 0 || nu[1L]==0 || m == 0) return(1)
    if(length(nu) > m && nu[m+1] > 0) return(0)
    if(m == 1) return(x[1]^nu[1] * prod(alpha*seq_len(nu[1L]-1)+1))
    if(k == 0 && !is.na(s <- S[.N(lambda,nu),m])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0, nu, nu) * .beta(mu,nu,alpha) *
      x[m]^(sum(mu)-sum(nu))
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i && nu[i] > 0 || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        if(nu[i] > 1){
          s <- s + jac(m, i, mu, .nu)
        }else{
          s <- s + jac(m-1L, 0L, .nu, .nu) * .beta(mu,.nu,alpha) *
            x[m]^(sum(mu)-sum(.nu))
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[.N(lambda,nu),m] <- s
    return(s)
  }
  S <- matrix(NA_real_, nrow = .N(lambda,lambda), ncol = length(x))
  jac(length(x), 0L, lambda, lambda)
}

#' Evaluation of Jack polynomial with rational parameters
#'
#' Evaluates the Jack polynomials with rational parameters.
#'
#' @param x \link[gmp]{bigq} vector
#' @param lambda integer partition, a vector of decreasing integers
#' @param alpha positive \code{bigq} rational
#'
#' @return A \code{bigq} rational number.
#' @export
#' @importFrom gmp as.bigq is.bigq
#'
#' @examples x <- c(1, 2, 5/2)
#' lambda <- c(3, 1)
#' Jack(x, lambda, alpha = 2)
#' library(gmp)
#' x <- c(as.bigq(1), as.bigq(2), as.bigq(5,2))
#' JackQ(x, lambda, alpha = as.bigq(2))
JackQ <- function(x, lambda, alpha){
  stopifnot(alpha > 0, is.bigq(alpha), is.bigq(x))
  if(any(diff(lambda) > 0) || any(floor(lambda) != lambda)){
    stop("`lambda` is not a valid integer partition")
  }
  jac <- function(m, k, mu, nu){
    if(all(nu==0) || length(nu) == 0 || nu[1L]==0 || m == 0){
      return(as.bigq(1L))
    }
    if(length(nu) > m && nu[m+1] > 0) return(as.bigq(0L))
    if(m == 1) return(x[1L]^nu[1L] * prod(alpha*seq_len(nu[1L]-1)+1L))
    if(k == 0 && !is.na(s <- S[.N(lambda,nu),m])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0, nu, nu) * .beta_gmp(mu,nu,alpha) *
      x[m]^(sum(mu)-sum(nu))
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i && nu[i] > 0 || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        if(nu[i] > 1){
          s <- s + jac(m, i, mu, .nu)
        }else{
          s <- s + jac(m-1L, 0L, .nu, .nu) * .beta_gmp(mu,.nu,alpha) *
            x[m]^(sum(mu)-sum(.nu))
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[.N(lambda,nu),m] <- s
    return(s)
  }
  S <- as.bigq(matrix(NA_real_, nrow = .N(lambda,lambda), ncol = length(x)))
  jac(length(x), 0L, lambda, lambda)
}
