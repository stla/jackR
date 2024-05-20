#' @title Kostka numbers
#'
#' @description The Kostka numbers for partitions of a given weight.
#'
#' @param n positive integer, the weight of the partitions
#' @param alpha the Jack parameter, a \code{bigq} number or an object coercible
#'   to a \code{bigq} number; setting \code{alpha=NULL} is equivalent to set
#'   \code{alpha=1}
#'
#' @return A matrix of character strings representing integers or fractions.
#' @export
#' @importFrom gmp as.bigq c_bigq
#'
#' @examples
#' KostkaNumbers(4)
KostkaNumbers <- function(n, alpha = NULL) {
  if(is.null(alpha)) {
    Knumbers <- SchurCoefficientsQ(n)
    stringParts <- paste0("(", gsub("(, 0| )", "", colnames(Knumbers)), ")")
  } else {
    alpha <- as.bigq(alpha)
    if(is.na(alpha)) {
      stop("Invalid `alpha`.")
    }
    Jcoeffs <- JackCoefficientsQ(n, alpha)
    JackPcoeffs <- c_bigq(lapply(rownames(Jcoeffs), function(lambda) {
      JackPcoefficient(fromString(lambda), alpha)
    }))
    Knumbers <- as.character(as.bigq(Jcoeffs)*JackPcoeffs)
    stringParts <- paste0("(", gsub("(, 0| )", "", colnames(Jcoeffs)), ")")
  }
  colnames(Knumbers) <- rownames(Knumbers) <- stringParts
  Knumbers
}

#' @importFrom ratioOfQsprays as.ratioOfQsprays showRatioOfQspraysOption
#' @importFrom partitions restrictedparts
#' @importFrom utils tail
#' @noRd
.symbolicKostkaNumbers <- function(n, weight, which){
  stopifnot(n > 0L, isPositiveInteger(n))
  stopifnot(isPositiveInteger(weight), weight > 0L)
  zeroRatioOfQsprays <- as.ratioOfQsprays(0L)
  unitRatioOfQsprays <- as.ratioOfQsprays(1L)
  if(weight == 1L) {
    return(list("[1]" = list("[1]" = unitRatioOfQsprays)))
  }
  allParts <- restrictedparts(weight, n)
  nParts <- ncol(allParts)
  lambdas <- apply(allParts, 2L, function(part) {
    part[part != 0L]
  }, simplify = FALSE)
  stringParts <- vapply(lambdas, partitionAsString, character(1L))
  row <- rep(list(zeroRatioOfQsprays), nParts)
  names(row) <- stringParts
  coefs <- lapply(seq_len(nParts), function(i) {
    part <- list(unitRatioOfQsprays)
    names(part) <- stringParts[i]
    parts <- tail(stringParts, nParts - i)
    c(part, row[parts])
  })
  names(coefs) <- stringParts
  for(m in seq_len(nParts-1L)){
    lambda <- lambdas[[m]]
    eSymbolic_lambda <- .eSymbolic(lambda)
    for(k in (m+1L):nParts){
      mu <- lambdas[[k]]
      l <- length(mu)
      eSymbolic_mu <- .eSymbolic(mu)
      x <- zeroRatioOfQsprays
      for(i in 1L:(l-1L)){ # l is always >1
        for(j in (i+1L):l){
          for(t in seq_len(mu[j])){
            mucopy <- mu
            mucopy[i] <- mucopy[i] + t
            mucopy[j] <- mucopy[j] - t
            muOrd <- sort(mucopy[mucopy != 0L], decreasing = TRUE)
            if(isDominated(muOrd, lambda)){
              x <- x + (mucopy[i] - mucopy[j]) /
                (eSymbolic_lambda - eSymbolic_mu) *
                coefs[[m]][[partitionAsString(muOrd)]]
            }
          }
        }
      }
      showRatioOfQspraysOption(x, "showRatioOfQsprays") <-
        showRatioOfQspraysXYZ("alpha")
      coefs[[m]][[partitionAsString(mu)]] <- x
    }
  }
  if(which != "P") {
    names(lambdas) <- names(coefs)
    invPcoeffs <- lapply(lambdas, symbolicJackPcoefficientInverse)
    Names <- names(coefs)
    names(Names) <- Names
    if(which == "J") {
      coefs <- lapply(Names, function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], list(invPcoeffs[[lambda]]),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    } else if(which == "C") {
      Ccoefs <- lapply(lambdas, symbolicJackCcoefficient)
      factors <-
        mapply(`*`, Ccoefs, invPcoeffs, SIMPLIFY = FALSE, USE.NAMES = TRUE)
      coefs <- lapply(Names, function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], list(factors[[lambda]]),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    } else {
      invQcoeffs <- lapply(lambdas, symbolicJackQcoefficientInverse)
      factors <-
        mapply(`/`, invPcoeffs, invQcoeffs, SIMPLIFY = FALSE, USE.NAMES = TRUE)
      coefs <- lapply(Names, function(lambda) {
        mapply(
          `*`,
          coefs[[lambda]], list(factors[[lambda]]),
          SIMPLIFY = FALSE, USE.NAMES = TRUE
        )
      })
    }
  }
  coefs
}

#' @title Symbolic Kostka numbers
#'
#' @description Kostka numbers with symbolic parameter for partitions of a
#'   given weight.
#'
#' @param n positive integer, the weight of the partitions
#'
#' @return A named list of named lists of \code{ratioOfQsprays} objects.
#'   Denoting the Kostka numbers by \eqn{K_{\lambda,\mu}(\alpha)}, the names
#'   of the outer list correspond to the partitions \eqn{\lambda}, and the
#'   names of the inner lists correspond to the partitions \eqn{\mu}.
#'   This list contains only the non-zero Kostka numbers.
#' @export
#'
#' @examples
#' symbolicKostkaNumbers(3)
symbolicKostkaNumbers <- function(n) {
  .symbolicKostkaNumbers(n, n, which = "P")
}
