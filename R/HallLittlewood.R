SSYTXwithGivenShapeAndContent <- function(lambda, mu) { # mu partition
  if(sum(lambda) != sum(mu)) {
    return(list())
  }
  if(!isDominated(mu, lambda)) {
    return(list())
  }
  if(all(lambda == 1L)) {
    if(all(mu == 1L)) {
      return(list(as.list(seq_along(lambda))))
    } else {
      return(list())
    }
  }
  l <- length(mu)
  mu[l] <- mu[l] - 1L
  mu <- removeTrailingZeros(mu)
  out <- list()
  for(i in seq_along(lambda)) {
    kappa <- lambda
    kappa[i] <- kappa[i] - 1L
    if(isDecreasing(kappa)) {
      kappa <- removeTrailingZeros(kappa)
      ssytx <- SSYTXwithGivenShapeAndContent(kappa, mu)
      ssytx <- lapply(ssytx, function(ssyt) {
        copy <- ssyt
        if(i > length(ssyt)) {
          row <- integer(0L)
        } else {
          row <- ssyt[[i]]
        }
        copy[[i]] <- c(row, l)
        copy
      })
      out <- c(out, unique(ssytx))
    }
  }
  out
}

#' @importFrom utils tail
#' @noRd
charge <- function(w) {
  l <- length(w)
  if(l == 0L) {
    return(0L)
  }
  n <- max(w)
  if(n == 1L) {
    return(0L)
  }
  pos <- match(1L, w)
  index <- 0L
  indices <- integer(n-1L)
  positions <- integer(n)
  positions[1L] <- pos
  for(r in 1L + seq_len(n-1L)) {
    v <- tail(w, l - pos)
    pos <- match(r, v) + pos
    if(is.na(pos)) {
      pos <- match(r, w)
      index <- index + 1L
    }
    indices[r-1L] <- index
    positions[r] <- pos
  }
  sum(indices) + charge(w[-positions])
}

#' @title Kostka-Foulkes polynomial
#' @description Kostka-Foulkes polynomial for two given partitions.
#'
#' @param lambda,mu integer partitions
#'
#' @return The Kostka-Foulkes polynomial associated to \code{lambda} and
#'   \code{mu}.
#' @export
#' @importFrom qspray qlone qzero showQsprayOption<- showQsprayXYZ
#'
#' @examples
KostaFoulkesPolynomial <- function(lambda, mu) {
  if(isDominated(mu, lambda)) {
    Tx <- SSYTXwithGivenShapeAndContent(lambda, mu)
    ws <- lapply(Tx, function(t) unlist(lapply(t, rev)))
    charges <- vapply(ws, charge, integer(1L))
    t <- qlone(1L)
    out <- Reduce(`+`, lapply(charges, function(e) t^e))
  } else {
    out <- qzero()
  }
  showQsprayOption(out, "showQspray") <- showQsprayXYZ("t")
  out
}

#' @importFrom partitions parts
#' @importFrom symbolicQspray Qzero
HallLittlewoodP <- function(n, lambda) {
  weight <- sum(lambda)
  lambdas <- lapply(Columns(parts(weight)), removeTrailingZeros)
  lambdaStrings <- vapply(lambdas, partitionAsString, character(1L))
  names(lambdas) <- lambdaStrings
  i <- match(partitionAsString(lambda), lambdaStrings)
  lambdas <- lambdas[i:length(lambdas)]
  kfs <- lapply(lambdas, function(kappa) {
    dom <- lapply(Columns(dominatedPartitions(kappa)), removeTrailingZeros)
    names(dom) <- vapply(dom, partitionAsString, character(1L))
    lapply(dom, function(mu) {
      KostaFoulkesPolynomial(kappa, mu)
    })
  })
  coeffs <- invTriMatrix(kfs)
  coeffs <- coeffs[[lambdaStrings[i]]]
  hlp <- Qzero()
  for(mu in names(coeffs)) {
    c <- coeffs[[mu]]
    mu <- fromPartitionAsString(mu)
    hlp <- hlp + c * as(SchurPol(n, mu), "symbolicQspray")
  }
  hlp
}
