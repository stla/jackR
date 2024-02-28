insertWith <- function(f, mp, key, value) {
  if(key %in% names(mp)) {
    mp[[key]] <- f(value, mp[[key]])
  } else {
    mp[[key]] <- value
  }
  mp
}

#' @importFrom utils tail
drop <- function(n, x) {
  tail(x, max(length(x) - n, 0L))
}

#' @title Littlewood-Richardson rule
#' @description Expression of the product of two Schur polynomials as a linear
#'   combination of Schur polynomials.
#'
#' @param mu,nu integer partitions, given as vectors of decreasing integers
#' @param output the type of the output, \code{"dataframe"} or \code{"list"}
#'
#' @return This computes the expression of the two Schur polynomials
#'   associated to \code{mu} and \code{nu} as a linear combination of Schur
#'   polynomials. If \code{output="dataframe"}, the output is a dataframe with
#'   two columns: the column \code{coeff} gives the coefficients of this
#'   linear combination, and the column \code{lambda} gives the partitions
#'   defining the Schur polynomials of this linear combination as character
#'   strings, e.g. the partition \code{c(4, 3, 1)} is given by \code{"4, 3, 1"}.
#'   If \code{output="list"}, the output is a list with two fields: the field
#'   \code{coeff} is the vector made of the coefficients of the linear
#'   combination, and the field \code{lambda} is the list of partitions
#'   defining the Schur polynomials of the linear combination given as
#'   integer vectors.
#' @export
#'
#' @examples
#' library(jack)
#' mu <- c(2, 1)
#' nu <- c(3, 2, 1)
#' LR <- LRmult(mu, nu, output = "list")
#' LRcoeffs <- LR$coeff
#' LRparts <- LR$lambda
#' LRterms <- lapply(1:length(LRcoeffs), function(i) {
#'   LRcoeffs[i] * SchurPolCPP(3, LRparts[[i]])
#' })
#' smu_times_snu <- Reduce(`+`, LRterms)
#' smu_times_snu == SchurPolCPP(3, mu) * SchurPolCPP(3, nu)
LRmult <- function(mu, nu, output = "dataframe") {
  stopifnot(isPartition(mu), isPartition(nu))
  output <- match.arg(output, c("dataframe", "list"))
  add <- function(old, lambda) {
    insertWith(`+`, old, toString(lambda), 1L)
  }
  v <- Reduce(add, addMu(mu, nu), init = integer(0L))
  if(output == "dataframe") {
    data.frame("coeff" = v, "lambda" = names(v))
  } else {
    partitions <- lapply(strsplit(names(v), ","), as.integer)
    list("coeff" = unname(v), "lambda" = partitions)
  }
}

addMu <- function(mu, part) {
  go <- function(ubs, mms, dds, part) {
    if(length(mms) == 0L) {
      list(part)
    } else {
      d <- dds[1L]
      ds <- dds[-1L]
      ms <- mms[-1L]
      L <- addRowOf(ubs, part)
      LL <- lapply(L, function(x) {
        ubsprime <- x[[1L]]
        partprime <- x[[2L]]
        go(drop(d, ubsprime), ms, ds, partprime)
      })
      do.call(c, LL)
    }
  }
  ubs0 <- headOrZero(part) + seq_len(headOrZero(mu))
  dmu <- diffSeq(mu)
  go(ubs0, mu, dmu, part)
}

addRowOf <- function(pcols, part) {
  go <- function(lb, ububs, p, ncols) {
    if(length(ububs) == 0L) {
      list(list(rev(ncols), p))
    } else {
      ub <- ububs[1L]
      ubs <- ububs[-1L]
      LL <- lapply(newBoxes(lb + 1L, ub, p), function(ij) {
        col <- ij[2L]
        go(col, ubs, addBox(ij, p), c(col, ncols))
      })
      do.call(c, LL)
    }
  }
  go(0L, pcols, part, integer(0L))
}

newBoxes <- function(lb, ub, part) {
  go <- function(i, jjs, lp) {
    if(length(jjs) == 0L) {
      if(lb <= 1L && 1L <= ub && lp > 0L) {
        list(c(i, 1L))
      } else {
        list()
      }
    } else {
      j <- jjs[1L]
      js <- jjs[-1L]
      j1 <- j + 1L
      if(j1 < lb) {
        list()
      } else if(j1 <= ub && lp > j) {
        c(list(c(i, j1)), go(i + 1L, js, j))
      } else {
        go(i + 1L, js, j)
      }
    }
  }
  rev(go(1L, part, headOrZero(part) + 1L))
}

addBox <- function(kl, part) {
  k <- kl[1L]
  go <- function(i, pps) {
    if(length(pps) == 0L) {
      if(i == k) {
        1L
      } else {
        stop("addBox: shouldn't happen")
      }
    } else {
      p <- pps[1L]
      ps <- pps[-1L]
      if(i == k) {
        c(p + 1L, ps)
      } else {
        c(p, go(i + 1L, ps))
      }
    }
  }
  go(1L, part)
}

headOrZero <- function(xs) {
  if(length(xs) == 0L) {
    0L
  } else {
    xs[1L]
  }
}

diffSeq <- function(x) {
  diff(-c(x, 0L))
}
