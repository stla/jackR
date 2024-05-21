f <- function(qspray, which = "J", check = TRUE) {
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
  which <- match.arg(which, c("J", "P", "Q", "C"))
  fullMsCombo <- MSPcombination(qspray - constantTerm, check = check)
  lambdas <- lapply(fullMsCombo, `[[`, "lambda")
  weights <- unique(vapply(lambdas, sum, integer(1L)))
  n <- numberOfVariables(qspray)
  finalQspray <- Qzero()
  for(weight in weights) {
    invKostkaMatrix <- jack:::msPolynomialsInJackSymbolicBasis(which, n, weight)
    kappas <- lapply(names(invKostkaMatrix), jack:::fromPartitionAsString)
    msCombo <- Filter(function(t) {sum(t[["lambda"]]) == weight}, fullMsCombo)
    coeffs <- lapply(msCombo, `[[`, "coeff")
    sprays <- lapply(kappas, function(kappa) {
      new(
        "symbolicQspray",
        powers = list(kappa),
        coeffs = list(as.ratioOfQsprays(1L))
      )
    })
    names(sprays) <- names(invKostkaMatrix)
    lambdas <- names(msCombo)
    for(i in seq_along(lambdas)) {
      invKostkaNumbers <- invKostkaMatrix[[lambdas[i]]]
      spray <- Qzero()
      for(kappa in names(invKostkaNumbers)) {
        coeff <- invKostkaNumbers[[kappa]]
        if(coeff != 0L) {
          spray <- spray + coeff * sprays[[kappa]]
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
      c(
        combo,
        list("[]" = list(
          "coeff" = as.ratioOfQsprays(constantTerm),
          "lambda" = integer(0L)
        )
        )
      )
  }
  combo
}

a <- qlone(2)
jp <- a*JackSymPol(3, c(2,1)) + 4*JackSymPol(3, c(1,1,1))
f(jp)
