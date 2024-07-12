.invKostkaMatrix <- function(weight, alpha, which) {
  # if(alpha == 1L) {
  #   Pcoeffs <- SchurCoefficientsQ(weight)
  #   dimNames <- gsub(", 0", "", colnames(Pcoeffs), fixed = TRUE)
  #   kappas <- lapply(dimNames, fromString)
  #   dimNames <- paste0("[", dimNames, "]")
  #   KostkaMatrix <- as.bigq(Pcoeffs)
  #   if(which != "P") {
  #     JackPcoeffs <- c_bigq(lapply(kappas, function(lambda) {
  #       JackPcoefficient(lambda, 1L)
  #     }))
  #     if(which == "Q") {
  #       JackQcoeffs <- c_bigq(lapply(kappas, function(lambda) {
  #         JackQcoefficient(lambda, 1L)
  #       }))
  #       factors <- JackQcoeffs / JackPcoeffs
  #     } else if(which == "C") {
  #       JackCcoeffs <- c_bigq(lapply(kappas, function(lambda) {
  #         JackCcoefficient(lambda, 1L)
  #       }))
  #       factors <- JackCcoeffs / JackPcoeffs
  #     } else {
  #       factors <- 1L / JackPcoeffs
  #     }
  #     KostkaMatrix <- KostkaMatrix * factors
  #   }
  # } else {
  #   Jcoeffs <- JackCoefficientsQ(weight, alpha)
  #   dimNames <- gsub(", 0", "", colnames(Jcoeffs), fixed = TRUE)
  #   kappas <- lapply(dimNames, fromString)
  #   dimNames <- paste0("[", dimNames, "]")
  #   KostkaMatrix <- as.bigq(Jcoeffs)
  #   if(which != "J") {
  #     if(which == "Q") {
  #       factors <- c_bigq(lapply(kappas, function(lambda) {
  #         JackQcoefficient(lambda, alpha)
  #       }))
  #     } else if(which == "C") {
  #       factors <- c_bigq(lapply(kappas, function(lambda) {
  #         JackCcoefficient(lambda, alpha)
  #       }))
  #     } else {
  #       factors <- c_bigq(lapply(kappas, function(lambda) {
  #         JackPcoefficient(lambda, alpha)
  #       }))
  #     }
  #     KostkaMatrix <- KostkaMatrix * factors
  #   }
  # }
  Pcoeffs <- KostkaJackNumbers(weight, alpha)
  dimNames <- colnames(Pcoeffs)
  kappas <- lapply(dimNames, fromPartitionAsString)
  KostkaMatrix <- as.bigq(Pcoeffs)
  if(which != "P") {
    JackPcoeffs <- c_bigq(lapply(kappas, function(lambda) {
      JackPcoefficient(lambda, alpha)
    }))
    if(which == "Q") {
      JackQcoeffs <- c_bigq(lapply(kappas, function(lambda) {
        JackQcoefficient(lambda, alpha)
      }))
      factors <- JackQcoeffs / JackPcoeffs
    } else if(which == "C") {
      JackCcoeffs <- c_bigq(lapply(kappas, function(lambda) {
        JackCcoefficient(lambda, alpha)
      }))
      factors <- JackCcoeffs / JackPcoeffs
    } else {
      factors <- 1L / JackPcoeffs
    }
    KostkaMatrix <- KostkaMatrix * factors
  }
  invKostkaMatrix <- Qinverse(KostkaMatrix)
  rownames(invKostkaMatrix) <- dimNames
  list(
    "matrix" = invKostkaMatrix,
    "kappas" = kappas
  )
}

#' @title Symmetric polynomial in terms of Jack polynomials
#' @description Expression of a symmetric polynomial as a linear combination
#'   of Jack polynomials.
#'
#' @param qspray a \code{qspray} object defining a symmetric polynomial
#' @param alpha Jack parameter, must be coercible to a \code{bigq} number
#' @param which which Jack polynomials, \code{"J"}, \code{"P"}, \code{"Q"} or
#'   \code{"C"}
#' @param check Boolean, whether to check the symmetry
#'
#' @return A list defining the combination. Each element of this list is a
#'   list with two elements: \code{coeff}, a \code{bigq} number, and
#'   \code{lambda}, an integer partition; then this list corresponds to the
#'   term \code{coeff * JackPol(n, lambda, alpha, which)}, where \code{n} is
#'   the number of variables in the symmetric polynomial.
#' @export
#' @importFrom methods new
#' @importFrom gmp as.bigq c_bigq
#' @importFrom qspray MSPcombination qzero orderedQspray isConstant isQzero getConstantTerm
#' @importFrom RationalMatrix Qinverse
JackCombination <- function(qspray, alpha, which = "J", check = TRUE) {
  if(isConstant(qspray)) {
    if(isQzero(qspray)) {
      out <- list()
    } else {
      out <-
        list(list("coeff" = getConstantTerm(qspray), "lambda" = integer(0L)))
      names(out) <- "[]"
    }
    return(out)
  }
  constantTerm <- getConstantTerm(qspray)
  alpha <- as.bigq(alpha)
  if(is.na(alpha)) {
    stop("Invalid `alpha`.")
  }
  which <- match.arg(which, c("J", "P", "Q", "C"))
  fullMsCombo <- MSPcombination(qspray - constantTerm, check = check)
  lambdas <- lapply(fullMsCombo, `[[`, "lambda")
  weights <- unique(vapply(lambdas, sum, integer(1L)))
  # invKostkaMatrices <- lapply(weights, function(weight) {
  #   .invKostkaMatrix(weight, alpha, which)
  # })
  finalQspray <- qzero()
  for(weight in weights) {
    invKostkaMatrix <- .invKostkaMatrix(weight, alpha, which)
    kappas <- invKostkaMatrix[["kappas"]]
    invKostkaMatrix <- invKostkaMatrix[["matrix"]]
    msCombo <- Filter(function(t) {sum(t[["lambda"]]) == weight}, fullMsCombo)
    coeffs <- c_bigq(lapply(msCombo, `[[`, "coeff"))
    sprays <- lapply(kappas, function(kappa) {
      new("qspray", powers = list(kappa), coeffs = "1")
    })
    lambdas <- names(msCombo)
    range <- seq_along(kappas)
    for(i in seq_along(lambdas)) {
      invKostkaNumbers <- invKostkaMatrix[lambdas[i], ]
      spray <- qzero()
      for(j in range) {
        coeff <- invKostkaNumbers[j]
        if(coeff != "0") {
          spray <- spray + coeff * sprays[[j]]
        }
      }
      finalQspray <- finalQspray + coeffs[i]*spray
    }
  }
  finalQspray <- orderedQspray(finalQspray)
  powers <- finalQspray@powers
  coeffs <- as.bigq(finalQspray@coeffs)
  combo <- mapply(
    function(lambda, coeff) {
      list("coeff" = coeff, "lambda" = lambda)
    },
    powers, coeffs,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  names(combo) <-
    vapply(powers, partitionAsString, character(1L), USE.NAMES = FALSE)
  if(constantTerm != 0L) {
    combo <-
      c(combo, list("[]" = list("coeff" = constantTerm, "lambda" = integer(0L))))
  }
  combo
}
