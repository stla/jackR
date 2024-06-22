# qtKostkaPolynomials ::
#   forall a. (Eq a, AlgField.C a)
#   => Partition
#   -> Map Partition (Spray a)
# qtKostkaPolynomials mu =  DM.map _numerator scs
#   where
#     psCombo = macdonaldJinPSbasis mu
#     t = lone' 2
# --    f = psCombo -- psCombination'' (macdonaldJpolynomial n mu)
#     den lambda = productOfSprays [unitSpray ^-^ t k | k <- lambda]
#     msCombo lambda =
#       msCombination (psPolynomial (length lambda) lambda)
#     ikn = inverseKostkaNumbers (sum mu)
#     coeffs lambda =
#       let combo = msCombo lambda in
#         DM.map
#           (\ikNumbers ->
#             DF.sum $ DM.intersectionWith (*) combo ikNumbers)
#           ikn
#     scs = DM.foldlWithKey
#       (\m lambda c ->
#         DM.unionWith (AlgAdd.+) m
#           (DM.map (\ikNumber -> (ikNumber .^ c) %//% den lambda) (coeffs lambda))
#       )
#       DM.empty psCombo
.den2 <- function(lambda) {
  t <- qlone(2)
  unitSpray <- qone()
  Reduce(
    `*`,
    lapply(lambda, function(k) {
      unitSpray - t^k
    })
  )
}

qtKostkaPolynomials <- function(mu) {
  psCombo <- MacdonaldPolynomialJinPSbasis(mu)
  n <- sum(mu)
  iknMatrix <- Qinverse(KostkaNumbers(n))
  lambdas <- apply(parts(n), 2L, removeTrailingZeros, simplify = FALSE)
  colnames(iknMatrix) <- rownames(iknMatrix) <- vapply(lambdas, partitionAsString, character(1L))
  coeffs <- function(lambda) {
    combo <- MSPcombination(PSFpoly(length(lambda), lambda), check = FALSE)
    out <- lapply(seq_along(lambdas), function(i) {
      list(
        "lambda" = lambdas[[i]],
        "coeff" = sum(c_bigq(lapply(names(combo), function(p) {
          combo[[p]][["coeff"]] * iknMatrix[p, i]
        })))
      )
    })
    names(out) <- vapply(lambdas, partitionAsString, character(1L))
    out
  }
  maps <- lapply(names(psCombo), function(p) {
    lambda <- psCombo[[p]][["lambda"]]
    spray <- psCombo[[p]][["coeff"]]
    coeffs_lambda <- coeffs(lambda)
    lapply(names(coeffs_lambda), function(x) {
      list(
        "lambda" = coeffs_lambda[[x]][["lambda"]],
        "lambdaAsString" = x,
        "coeff" =  coeffs_lambda[[x]][["coeff"]] * spray / .den2(lambda)
      )
    })
  })
  cmapOfMaps <- do.call(c, maps)
  f <- vapply(cmapOfMaps, `[[`, character(1L), "lambdaAsString")
  lapply(split(cmapOfMaps, f), function(l) {
    spray <- Reduce(`+`, lapply(l, `[[`, "coeff"))
    list(
      "lambda" = l[[1L]][["lambda"]],
      "coeff" = spray
    )
  })


  # ikn <- lapply(seq_along(lambdas), function(i) {
  #   # ikNumbers <- iknMatrix[, i]
  #   # names(ikNumbers) <- vapply(lambdas, partitionAsString, character(1L))
  #   lambda <- lambdas[[i]]
  #   combo <- MSPcombination(PSFpoly(length(lambda), lambda), check = FALSE)
  #   # out <- lapply(names(combo), function(p) {
  #   #   list(
  #   #     "lambda" = combo[[p]][["lambda"]],
  #   #     "coeff" = sum(c_bigq(lapply(names(combo), function(p) {
  #   #       combo[[p]][["coeff"]] * iknMatrix[i,p]
  #   #     })))
  #   #   )
  #   # })
  #   # names(out) <- names(combo)#vapply(lambdas, partitionAsString, character(1L))
  #
  #   out <- lapply(seq_along(lambdas), function(j) {
  #     # ikNumbers <- iknMatrix[, j]
  #     # names(ikNumbers) <- vapply(lambdas, partitionAsString, character(1L))
  #     kappa <- lambdas[[j]]
  #     list(
  #       "lambda" = kappa,
  #       "coeff" = sum(c_bigq(lapply(names(combo), function(p) {
  #         combo[[p]][["coeff"]] * iknMatrix[i,p]
  #       })))
  #     )
  #   })
  #   names(out) <- vapply(lambdas, partitionAsString, character(1L))
  #   out
  #   # lapply(names(combo), function(p) {
  #   #   list(
  #   #     "lambda" = combo[[p]][["lambda"]],
  #   #   )
  #   # })
  #   # list(
  #   #   "lambda" = lambda,
  #   #   "coeff" = sum(c_bigq(lapply(names(combo), function(p) {
  #   #     combo[[p]][["coeff"]] * iknMatrix[p,i]
  #   #   })))
  #   # )
  # })
  # names(ikn) <- vapply(lambdas, partitionAsString, character(1L))
  # maps <- lapply(psCombo, function(term) {
  #   lambda <- term[["lambda"]]
  #   p <- partitionAsString(lambda)
  #   spray <- psCombo[[p]][["coeff"]]
  #   # lambda <- psCombo[[p]][["lambda"]]
  #   # list(
  #   #   "lambdaAsString" = partitionAsString(ikn[[p]][["lambda"]]),
  #   #   "lambda" = ikn[[p]][["lambda"]],
  #   #   "coeff" = ikn[[p]][["coeff"]] * spray / .den(lambda)
  #   # )
  #   # coeff <- Reduce(`+`, lapply(ikn[[p]], function(t) {
  #   #   t[["coeff"]] * spray / .den(lambda)
  #   #   # list(
  #   #   #   "lambda" = t[["lambda"]],
  #   #   #   "coeff" = t[["coeff"]] * spray / .den(lambda)
  #   #   # )
  #   # }))
  #   #    names(coeff) <- vapply(ikn[[p]], function(t) partitionAsString(t[["lambda"]]), character(1L))
  #   # list(
  #   #   "lambdaAsString" = p,
  #   #   "lambda" = lambda,
  #   #   "coeff" = ikn[[p]][["coeff"]] * spray / .den(lambda)
  #   #   # lapply(ikn, function(t) {
  #   #   #   t[["coeff"]] * spray / .den(lambda)
  #   #   # })
  #   #
  #   #   # "coeff" = list(
  #   #   #   "lambda" = ikn[[p]][["lambda"]],
  #   #   #   "coeff" = ikn[[p]][["coeff"]] * spray / .den(lambda)
  #   #   # )
  #   # )
  #   list(
  #     "lambdaAsString" = p,
  #     "coeff" = lapply(ikn, function(l) {
  #       lapply(l, function(t) {
  #         list("lambda" = t[["lambda"]], "coeff" = t[["coeff"]] * spray / .den(lambda))
  #         #        function(t) t[["coeff"]] * spray / .den(lambda)
  #       })
  #     })
  #   )
  #   # lapply(ikn, function(l) {
  #   #   list(
  #   #     "lambdaAsString" = p,#partitionAsString(t[["lambda"]]),
  #   #     "lambda" = lambda,
  #   #     "coeff" = lapply(l, function(t) t[["coeff"]] * spray / .den(lambda))
  #   #   )
  #   #
  #   #   # list(
  #   #   #   "lambdaAsString" = p,#partitionAsString(t[["lambda"]]),
  #   #   #   "lambda" = t[["lambda"]],
  #   #   #   "coeff" = t[["coeff"]] * spray / .den(lambda)
  #   #   # )
  #   # })
  # })
  # cmapOfMaps <- maps#do.call(c, maps)
  # f <- vapply(cmapOfMaps, `[[`, character(1L), "lambdaAsString")
  # lapply(maps, function(l) {
  #   rOS <- lapply(do.call(c, l[["coeff"]]), `[[`, "coeff")
  #   list(
  #     "lambda" = NULL,#l[["lambdaAsString"]],
  #     "coeff" = rOS
  #   )
  # })
}
