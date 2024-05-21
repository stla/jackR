f <- function(qspray, alpha, which = "J", check = TRUE) {
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
  finalQspray <- Qzero()
  for(weight in weights) {
    invKostkaMatrix <- .invKostkaMatrix(weight, alpha, which)
    kappas <- invKostkaMatrix[["kappas"]]
    invKostkaMatrix <- invKostkaMatrix[["matrix"]]
    msCombo <- Filter(function(t) {sum(t[["lambda"]]) == weight}, fullMsCombo)
    coeffs <- lapply(msCombo, `[[`, "coeff")
    sprays <- lapply(kappas, function(kappa) {
      new("symbolicQspray", powers = list(kappa), coeffs = list(as.ratioOfQsprays(1L)))
    })
    lambdas <- names(msCombo)
    range <- seq_along(kappas)
    for(i in seq_along(lambdas)) {
      invKostkaNumbers <- invKostkaMatrix[lambdas[i], ]
      spray <- Qzero()
      for(j in range) {
        coeff <- invKostkaNumbers[j]
        if(coeff != "0") {
          spray <- spray + coeff * sprays[[j]]
        }
      }
      finalQspray <- finalQspray + coeffs[[i]]*spray
    }
  }
  finalQspray <- orderedQspray(finalQspray)
  powers <- finalQspray@powers
  coeffs <- finalQspray@coeffs
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
  if(constantTerm != 0L) {
    combo <-
      c(combo, list("[]" = list("coeff" = constantTerm, "lambda" = integer(0L))))
  }
  combo
}

a <- qlone(1)
jp <- a*JackSymPol(3, c(2,1)) + 4*JackSymPol(3, c(1,1,1))
combo <- f(jp, alpha = 2)

JackSymbolicCombinationToQspray <- function(combo, n, which) {
  lambdas <- lapply(combo, `[[`, "lambda")
  coeffs <- lapply(combo, `[[`, "coeff")
  Reduce(`+`, mapply(function(lambda, coeff) {
    coeff * JackSymPol(n, lambda, which)
  }, lambdas, coeffs, SIMPLIFY = FALSE))
}

jp2 <- JackSymbolicCombinationToQspray(combo, 3, "J")

substituteParameters(jp, 2)
substituteParameters(jp2, 2)

