#' @importFrom gmp as.bigq is.bigq factorialZ
NULL

# JackEvalNum <- function(x, lambda, alpha){
#   jac <- function(m, k, mu, nu){
#     if(all(nu==0) || length(nu) == 0 || nu[1L]==0 || m == 0) return(1)
#     if(length(nu) > m && nu[m+1] > 0) return(0)
#     if(m == 1) return(x[1]^nu[1] * prod(alpha*seq_len(nu[1L]-1)+1))
#     if(k == 0 && !is.na(s <- S[.N(lambda,nu),m])) return(s)
#     i <- max(1L,k)
#     s <- jac(m-1L, 0, nu, nu) * .beta(mu,nu,alpha) *
#       x[m]^(sum(mu)-sum(nu))
#     while(length(nu) >= i && nu[i] > 0){
#       if(length(nu) == i && nu[i] > 0 || nu[i] > nu[i+1L]){
#         .nu <- nu; .nu[i] <- nu[i]-1
#         if(nu[i] > 1){
#           s <- s + jac(m, i, mu, .nu)
#         }else{
#           s <- s + jac(m-1L, 0L, .nu, .nu) * .beta(mu,.nu,alpha) *
#             x[m]^(sum(mu)-sum(.nu))
#         }
#       }
#       i <- i + 1L
#     }
#     if(k == 0L) S[.N(lambda,nu),m] <- s
#     return(s)
#   }
#   S <- matrix(NA_real_, nrow = .N(lambda,lambda), ncol = length(x))
#   jac(length(x), 0L, lambda, lambda)
# }

JackEvalNum <- function(x, lambda, alpha){
  jac <- function(m, k, mu, nu, beta){
    if(length(nu) == 0 || nu[1L]==0 || m == 0L) return(1)
    if(length(nu) > m && nu[m+1L] > 0) return(0)
    if(m == 1L) return(x[1L]^nu[1L] * prod(alpha*seq_len(nu[1L]-1)+1))
    if(k == 0L && !is.na(s <- S[.N(lambda,nu),m])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0L, nu, nu, 1) * beta * x[m]^(sum(mu)-sum(nu))
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i && nu[i] > 0 || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        gamma <- beta * .betaratio(mu, nu, i, alpha)
        if(nu[i] > 1){
          s <- s + jac(m, i, mu, .nu, gamma)
        }else{
          s <- s + jac(m-1L, 0L, .nu, .nu, 1) * gamma *
            x[m]^(sum(mu)-sum(.nu))
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[.N(lambda,nu),m] <- s
    return(s)
  }
  S <- matrix(NA_real_, nrow = .N(lambda,lambda), ncol = length(x))
  jac(length(x), 0L, lambda, lambda, 1)
}

JackEvalQ <- function(x, lambda, alpha){
  jac <- function(m, k, mu, nu, beta){
    if(length(nu) == 0 || nu[1L]==0 || m == 0){
      return(as.bigq(1L))
    }
    if(length(nu) > m && nu[m+1L] > 0) return(as.bigq(0L))
    if(m == 1L) return(x[1L]^nu[1L] * prod(alpha*seq_len(nu[1L]-1)+1L))
    if(k == 0L && !is.na(s <- S[.N(lambda,nu),m])) return(s)
    i <- max(1L,k)
    s <- jac(m-1L, 0L, nu, nu, as.bigq(1L)) * beta * x[m]^(sum(mu)-sum(nu))
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i|| nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        gamma <- beta * .betaratio(mu, nu, i, alpha)
        if(nu[i] > 1){
          s <- s + jac(m, i, mu, .nu, gamma)
        }else{
          s <- s + jac(m-1L, 0L, .nu, .nu, as.bigq(1L)) * gamma *
            x[m]^(sum(mu)-sum(.nu))
        }
      }
      i <- i + 1L
    }
    if(k == 0L) S[.N(lambda,nu),m] <- s
    return(s)
  }
  S <- as.bigq(matrix(NA_integer_, nrow = .N(lambda,lambda), ncol = length(x)))
  jac(length(x), 0L, lambda, lambda, as.bigq(1L))
}

# JackEvalQ <- function(x, lambda, alpha){ # TODO: use .betaratio
#   jac <- function(m, k, mu, nu){
#     if(length(nu) == 0 || nu[1L]==0 || m == 0){
#       return(as.bigq(1L))
#     }
#     if(length(nu) > m && nu[m+1L] > 0) return(as.bigq(0L))
#     if(m == 1L) return(x[1L]^nu[1L] * prod(alpha*seq_len(nu[1L]-1)+1L))
#     if(k == 0L && !is.na(s <- S[.N(lambda,nu),m])) return(s)
#     i <- max(1L,k)
#     s <- jac(m-1L, 0L, nu, nu) * .beta_gmp(mu,nu,alpha) *
#       x[m]^(sum(mu)-sum(nu))
#     while(length(nu) >= i && nu[i] > 0){
#       if(length(nu) == i|| nu[i] > nu[i+1L]){
#         .nu <- nu; .nu[i] <- nu[i]-1
#         if(nu[i] > 1){
#           s <- s + jac(m, i, mu, .nu)
#         }else{
#           s <- s + jac(m-1L, 0L, .nu, .nu) * .beta_gmp(mu,.nu,alpha) *
#             x[m]^(sum(mu)-sum(.nu))
#         }
#       }
#       i <- i + 1L
#     }
#     if(k == 0L) S[.N(lambda,nu),m] <- s
#     return(s)
#   }
#   S <- as.bigq(matrix(NA_integer_, nrow = .N(lambda,lambda), ncol = length(x)))
#   jac(length(x), 0L, lambda, lambda)
# }

JackEval <- function(x, lambda, alpha){
  stopifnot(isPartition(lambda), alpha >= 0)
  gmp <- is.bigq(x) || is.bigq(alpha)
  if(gmp){
    stopifnot(is.bigq(x), is.bigq(alpha))
    JackEvalQ(x, lambda, alpha)
  }else{
    JackEvalNum(x, lambda, alpha)
  }
}

ZonalEvalNum <- function(x, lambda){
  jack <- JackEvalNum(x, lambda, alpha= 2)
  jlambda <- sum(logHookLengths(lambda, alpha = 2))
  n <- sum(lambda)
  exp(n*log(2) + lfactorial(n) - jlambda) * jack
}

ZonalEvalQ <- function(x, lambda){
  jack <- JackEvalQ(x, lambda, alpha= as.bigq(2L))
  jlambda <- prod(hookLengths_gmp(lambda, alpha = as.bigq(2L)))
  n <- sum(lambda)
  as.bigq(2L)^n * as.bigq(factorialZ(n)) / jlambda * jack
}

ZonalEval <- function(x, lambda){
  stopifnot(isPartition(lambda))
  gmp <- is.bigq(x)
  if(gmp){
    ZonalEvalQ(x, lambda)
  }else{
    ZonalEvalNum(x, lambda)
  }
}

SchurEvalNum <- function(x, lambda){
  # i <- rep(seq_along(lambda), times = lambda)
  # j <- unlist(sapply(lambda, seq_len, simplify = FALSE))
  # lambdaPrime <- dualPartition(lambda)
  # hookslengths <- lambdaPrime[j] - i + lambda[i] - j + 1
  # JackEvalNum(x, lambda, alpha = 1) / prod(hookslengths)
  sch <- function(m, k, nu){
    if(length(nu) == 0L || nu[1L]==0 || m == 0L){
      return(1)
    }
    if(length(nu) > m && nu[m+1L] > 0) return(0)
    if(m == 1L) return(x[1L]^nu[1L])
    if(!is.na(s <- S[.N(lambda,nu),m])) return(s)
    s <- sch(m-1L, 1L, nu)
    i <- k
    # if(length(nu) > 1L && nu[1L] > nu[2L]){
    #   .nu <- nu; .nu[1L] <- nu[1L] - 1
    #   s <- s + 1 * sch(m, 1L, .nu)
    #   i <- 2L
    # }
    while(length(nu) >= i && nu[i] > 0){
      if(length(nu) == i || nu[i] > nu[i+1L]){
        .nu <- nu; .nu[i] <- nu[i]-1
        if(nu[i] > 1){
          s <- s + x[m] * sch(m, i, .nu)
        }else{
          s <- s + x[m] * sch(m-1L, 1L, .nu)
        }
      }
      i <- i + 1L
    }
    if(k == 1L) S[.N(lambda,nu),m] <- s
    return(s)
  }
  S <- matrix(NA_real_, nrow = .N(lambda,lambda), ncol = length(x))
  sch(length(x), 1L, lambda)
}

SchurEvalQ <- function(x, lambda){ # TODO: implement as SchurEvalNum
  jack <- JackEvalQ(x, lambda, alpha = as.bigq(1L))
  i <- rep(seq_along(lambda), times = lambda)
  j <- unlist(sapply(lambda, seq_len, simplify = FALSE))
  lambdaPrime <- as.bigq(dualPartition(lambda))
  lambda <- as.bigq(lambda)
  hookslengths <- lambdaPrime[j] - i + lambda[i] - j + 1L
  jack / prod(hookslengths)
}

SchurEval <- function(x, lambda){
  stopifnot(isPartition(lambda))
  gmp <- is.bigq(x)
  if(gmp){
    SchurEvalQ(x, lambda)
  }else{
    SchurEvalNum(x, lambda)
  }
}

ZonalQEvalNum <- function(x, lambda){
  jack <- JackEvalNum(x, lambda, alpha = 1/2)
  jlambda <- sum(logHookLengths(lambda, alpha = 1/2))
  n <- sum(lambda)
  exp(-n*log(2) + lfactorial(n) - jlambda) * jack
}

ZonalQEvalQ <- function(x, lambda){
  jack <- JackEvalQ(x, lambda, alpha = as.bigq(1L,2L))
  jlambda <- prod(hookLengths_gmp(lambda, alpha = as.bigq(1L,2L)))
  n <- sum(lambda)
  as.bigq(1L,2L)^n * as.bigq(factorialZ(n)) / jlambda * jack
}

ZonalQEval <- function(x, lambda){
  stopifnot(isPartition(lambda))
  gmp <- is.bigq(x)
  if(gmp){
    ZonalQEvalQ(x, lambda)
  }else{
    ZonalQEvalNum(x, lambda)
  }
}
