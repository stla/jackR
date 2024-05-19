library(jack)

qspray <- ESFpoly(4, c(3, 1)) + PSFpoly(4, c(2, 2))
alpha <- "2"
which <- "P"

alpha <- gmp::as.bigq(alpha)
if(is.na(alpha)) {
  stop("Invalid `alpha`.")
}
msCombo <- MSPcombination(qspray, check = FALSE)
lambdas <- lapply(msCombo, `[[`, "lambda")
coeffs <- gmp::c_bigq(lapply(msCombo, `[[`, "coeff"))
totalDegree <- max(vapply(lambdas, sum, integer(1L)))
if(alpha == 1L) {
  Pcoeffs <- jack:::SchurCoefficientsQ(totalDegree)
  dimNames <- gsub(", 0", "", colnames(Pcoeffs), fixed = TRUE)
  kappas <- lapply(dimNames, jack:::fromString)
  dimNames <- paste0("[", dimNames, "]")
  KostkaMatrix <- gmp::as.bigq(Pcoeffs)
  if(which != "P") {
    JackPcoeffs <- gmp::c_bigq(lapply(kappas, function(lambda) {
      jack:::JackPcoefficient(lambda, 1L)
    }))
    if(which == "Q") {
      JackQcoeffs <- gmp::c_bigq(lapply(kappas, function(lambda) {
        jack:::JackQcoefficient(lambda, 1L)
      }))
      factors <- JackQcoeffs / JackPcoeffs
    } else if(which == "C") {
      JackCcoeffs <- gmp::c_bigq(lapply(kappas, function(lambda) {
        jack:::JackCcoefficient(lambda, 1L)
      }))
      factors <- JackCcoeffs / JackPcoeffs
    } else {
      factors <- 1L / JackPcoeffs
    }
    KostkaMatrix <- KostkaMatrix * factors
  }
} else {
  Jcoeffs <- jack:::JackCoefficientsQ(totalDegree, alpha)
  dimNames <- gsub(", 0", "", colnames(Jcoeffs), fixed = TRUE)
  kappas <- lapply(dimNames, jack:::fromString)
  dimNames <- paste0("[", dimNames, "]")
  KostkaMatrix <- gmp::as.bigq(Jcoeffs)
  if(which != "J") {
    if(which == "Q") {
      factors <- gmp::c_bigq(lapply(kappas, function(lambda) {
        jack:::JackQcoefficient(lambda, alpha)
      }))
    } else if(which == "C") {
      factors <- gmp::c_bigq(lapply(kappas, function(lambda) {
        jack:::JackCcoefficient(lambda, alpha)
      }))
    } else {
      factors <- gmp::c_bigq(lapply(kappas, function(lambda) {
        jack:::JackPcoefficient(lambda, alpha)
      }))
    }
    KostkaMatrix <- KostkaMatrix * factors
  }
}
invKostkaMatrix <- RationalMatrix::Qinverse(KostkaMatrix)
rownames(invKostkaMatrix) <- colnames(invKostkaMatrix) <- dimNames
sprays <- lapply(kappas, function(kappa) {
  new("qspray", powers = list(kappa), coeffs = "1")
})
finalQspray <- qzero()
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
finalQspray <- orderedQspray(finalQspray)
powers <- finalQspray@powers
coeffs <- gmp::as.bigq(finalQspray@coeffs)
combo <- mapply(
  function(lambda, coeff) {
    list("coeff" = coeff, "lambda" = lambda)
  },
  powers, coeffs,
  SIMPLIFY = FALSE, USE.NAMES = FALSE
)
names(combo) <- paste0(
  "[",
  vapply(powers, toString, character(1L)),
  "]"
)
combo

comboToQspray <- function(combo) {
  lambdas <- lapply(combo, `[[`, "lambda")
  coeffs <- lapply(combo, `[[`, "coeff")
  Reduce(`+`, mapply(function(lambda, coeff) {
    coeff * JackPol(4, lambda, alpha, which)
  }, lambdas, coeffs, SIMPLIFY = FALSE))
}


comboToQspray(combo)
