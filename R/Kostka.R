#' @title Kostka-Jack numbers
#'
#' @description The Kostka-Jack numbers for partitions of a given weight.
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
#' @details The Kostka-Jack numbers generalize the Kostka numbers: these ones
#'   are obtained when \code{alpha=1}.
#'
#' @examples
#' KostkaJackNumbers(4)
KostkaJackNumbers <- function(n, alpha = NULL) {
  stopifnot(isPositiveInteger(n))
  if(n == 0L) {
    Knumbers <- as.matrix("1")
    colnames(Knumbers) <- rownames(Knumbers) <- "()"
    return(Knumbers)
  }
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

#' @importFrom ratioOfQsprays as.ratioOfQsprays showRatioOfQspraysXYZ showRatioOfQspraysOption<-
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

#' @title Symbolic Kostka-Jack numbers
#'
#' @description Kostka-Jack numbers with symbolic Jack parameter for partitions
#'   of a given weight.
#'
#' @param n positive integer, the weight of the partitions
#'
#' @return A named list of named lists of \code{ratioOfQsprays} objects.
#'   Denoting the Kostka numbers by \eqn{K_{\lambda,\mu}(\alpha)}, the names
#'   of the outer list correspond to the partitions \eqn{\lambda}, and the
#'   names of the inner lists correspond to the partitions \eqn{\mu}.
#' @export
#'
#' @examples
#' symbolicKostkaJackNumbers(3)
symbolicKostkaJackNumbers <- function(n) {
  if(n == 0L) {
    list("[]" = list("[]" = as.ratioOfQsprays(1L)))
  } else {
    .symbolicKostkaNumbers(n, n, which = "P")
  }
}

#' @title Skew Kostka-Jack numbers with given Jack parameter
#' @description Skew Kostka-Jack numbers associated to a given skew partition
#'   and a given Jack parameter.
#'
#' @param lambda,mu integer partitions defining the skew partition:
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition
#'   (so \code{mu} must be a subpartition of \code{lambda})
#' @param alpha the Jack parameter, a \code{bigq} number or an object coercible
#'   to a \code{bigq} number; setting \code{alpha=NULL} is equivalent to set
#'   \code{alpha=1}
#' @param output the format of the output, either \code{"vector"} or
#'   \code{"list"}
#'
#' @return If \code{output="vector"}, the function returns a named vector.
#'   This vector is made of the non-zero skew Kostka numbers
#'   \eqn{K_{\lambda/\mu,\nu}(\alpha)} given as character strings and its names
#'   encode the partitions \eqn{\nu}.
#'   If \code{ouput="list"}, the function returns a list. Each element of this
#'   list is a named list with two elements: an integer partition \eqn{\nu}
#'   in the field named \code{"nu"}, and the corresponding skew Kostka number
#'   \eqn{K_{\lambda/\mu,\nu}(\alpha)} in the field named \code{"value"}. Only
#'   the non-null skew Kostka numbers are provided by this list.
#' @export
#' @importFrom gmp as.bigq c_bigq
#' @importFrom utils head
#'
#' @examples
#' skewKostkaJackNumbers(c(4,2,2), c(2,2))
skewKostkaJackNumbers <- function(lambda, mu, alpha = NULL, output = "vector") {
  stopifnot(isPartition(lambda), isPartition(mu))
  output <- match.arg(output, c("vector", "list"))
  lambda <- as.integer(removeTrailingZeros(lambda))
  mu <- as.integer(removeTrailingZeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu || any(head(lambda, ellMu) < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  if(is.null(alpha)) {
    alpha <- "1"
  }
  alpha <- as.bigq(alpha)
  if(is.na(alpha)) {
    stop("Invalid `alpha`.")
  }
  listOfKnumbers <- skewJackInMSPbasis(alpha, "P", lambda, mu)
  if(output == "vector") {
    values <- lapply(listOfKnumbers, `[[`, "coeff")
    kNumbers <- as.character(c_bigq(values))
    names(kNumbers) <- names(listOfKnumbers)
  } else {
    kNumbers <- lapply(
      listOfKnumbers,
      function(lst) {
        list("nu" = lst[["nu"]], "value" = lst[["coeff"]])
      }
    )
  }
  kNumbers
}

#' @title Skew Kostka-Jack numbers with symbolic Jack parameter
#' @description Skew Kostka-Jack numbers associated to a given skew partition
#'   with a symbolic Jack parameter.
#'
#' @param lambda,mu integer partitions defining the skew partition:
#'   \code{lambda} is the outer partition and \code{mu} is the inner partition
#'   (so \code{mu} must be a subpartition of \code{lambda})
#'
#' @return The function returns a list. Each element of this
#'   list is a named list with two elements: an integer partition \eqn{\nu}
#'   in the field named \code{"nu"}, and the corresponding skew Kostka number
#'   \eqn{K_{\lambda/\mu,\nu}(\alpha)} in the field named \code{"value"}, a
#'   \code{ratioOfQsprays} object.
#' @export
#' @importFrom utils head
#' @importFrom ratioOfQsprays showRatioOfQspraysXYZ showRatioOfQspraysOption<-
#'
#' @examples
#' symbolicSkewKostkaJackNumbers(c(4,2,2), c(2,2))
symbolicSkewKostkaJackNumbers <- function(lambda, mu) {
  stopifnot(isPartition(lambda), isPartition(mu))
  lambda <- as.integer(removeTrailingZeros(lambda))
  mu <- as.integer(removeTrailingZeros(mu))
  ellLambda <- length(lambda)
  ellMu <- length(mu)
  if(ellLambda < ellMu || any(head(lambda, ellMu) < mu)) {
    stop("The partition `mu` is not a subpartition of the partition `lambda`.")
  }
  listOfKnumbers <- skewSymbolicJackInMSPbasis("P", lambda, mu)
  lapply(
    listOfKnumbers,
    function(lst) {
      rOS <- lst[["coeff"]]
      showRatioOfQspraysOption(rOS, "showRatioOfQsprays") <-
        showRatioOfQspraysXYZ("alpha")
      list("nu" = lst[["nu"]], "value" = rOS)
    }
  )
}
