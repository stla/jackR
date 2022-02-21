rationalMonomial <- function(variables, powers){
  factors <- gsub(
    "\\^1( ?)$", "\\1",
    paste0(variables, "^", powers, c(rep(" ", length(powers)-1L), ""))
  )
  paste0(factors, collapse = "")
}

#' @importFrom Ryacas yac_str
rationalPolynomial <- function(powers, coeffs, stars = FALSE){
  nterms <- length(coeffs)
  variables <- vector(mode = "list", length = nterms)
  for(i in 1L:nterms){
    pwrs <- powers[[i]]
    zs <- which(pwrs != 0L)
    powers[[i]] <- pwrs[zs]
    variables[[i]] <- paste0("x", zs)
  }
  monomials <- mapply(rationalMonomial, variables, powers, USE.NAMES = FALSE)
  # reversed
  terms <-
    paste0(gsub(" ", "*", monomials, fixed = TRUE), " * z Where z==", coeffs)
  yacpol <- yac_str(
    paste0(vapply(terms, yac_str, character(1L)), collapse = " + ")
  )
  # if(grepl("^\\(", yacpol)){
  #   yacpol <- sub(")", "", sub("(", "", yacpol, fixed = TRUE), fixed = TRUE)
  # }
  #
  spaces <- rep(" ", length(coeffs))
  ones <- coeffs == "1"
  minusones <- coeffs == "-1"
  spaces[ones] <- ""
  coeffs[ones] <- ""
  spaces[minusones] <- ""
  coeffs[minusones] <- "-"
  if(stars){
    out <- gsub("(\\d) x", "\\1 * x",
         gsub("+  -", "-  ",
              paste0(
                paste0(coeffs, spaces, monomials), collapse = "  +  "
              ),
              fixed = TRUE
         )
    )
    attr(out, "yacas") <- yacpol
  }else{
    out <- gsub("+ -", "- ",
         paste0(
           paste0(coeffs, spaces, monomials), collapse = " + "
         ),
         fixed = TRUE
    )
  }
  out
}


#' @title Print an \code{exactmvp} object
#' @description Print an \code{exactmvp} object.
#'
#' @param x object of class \code{exactmvp}; the functions returned by
#'   \code{\link{Jack_julia}} can return such objects
#' @param ... arguments passed to \code{\link[mvp]{print.mvp}}
#'
#' @return Nothing.
#' @export
#'
#' @importFrom mvp print.mvp
print.exactmvp <- function(x, ...){
  print.mvp(x, ...)
  cat("\nExact expression:\n")
  cat(attr(x, "exact"))
  cat("\n")
}

#' @title Exact multivariate polynomial as function
#' @description Coerces an exact multivariate polynomial into a function.
#'
#' @param x object of class \code{exactmvp}; the functions returned by
#'   \code{\link{Jack_julia}} can return such objects
#' @param ... ignored
#'
#' @return A function having the same variables as the polynomial.
#' @export
#'
#' @importFrom Ryacas yac_str
#'
#' @examples # library(jack)
#' \donttest{if(JuliaConnectoR::juliaSetupOk()){
#'   julia <- Jack_julia()
#'   ( pol <- julia$JackPol(m = 2, lambda = c(3, 1), alpha = "3/2") )
#'   f <- as.function(pol)
#'   f(2, "3/7")
#'   # the evaluation is performed by (R)yacas and complex numbers are
#'   # allowed; the imaginary unit is denoted by `I`
#'   f("2 + 2*I", "1/4")
#'   JuliaConnectoR::stopJulia()
#' }}
as.function.exactmvp <- function(x, ...){
  expr <- sprintf("Simplify(%s)", attr(x, "exact"))
  nvars <- attr(x, "nvars")
  vars <- paste0("x", seq_len(nvars))
  values <- paste0(paste0(vars, "==%s"), collapse = " And ")
  yacas <- paste0(expr, " Where ", values)
  f <- function(){
    yac_str(
      do.call(function(...) sprintf(yacas, ...), lapply(vars, function(xi){
        eval(parse(text = xi))
      }))
    )
  }
  formals(f) <- sapply(vars, function(xi){
    `names<-`(alist(y=), xi)
  }, USE.NAMES = FALSE)
  f
}

rationalize <- function(x){
  if(grepl("^ *\\d+ */ *\\d+ *$", x)){
    return(as.list(as.integer(strsplit(x, "/")[[1L]])))
  }
  if(grepl("^ *\\d+ *$", x)){
    return(list(as.integer(x), 1L))
  }
  stop(
    sprintf("The string '%s' cannot be identified to a fraction.", x)
  )
}

asIntegerList <- function(lambda){
  if(length(lambda) == 0L){
    integer(0L)
  }else{
    unname(as.list(as.integer(lambda)))
  }
}

#' @title Evaluation with Julia
#' @description Evaluate the Jack polynomials with Julia. This is highly faster.
#'
#' @return A list of functions having the same names as the R functions of this
#'   package (\code{Jack}, \code{JackPol}, \code{Schur}, etc). The
#'   \code{XXXPol} functions have an argument \code{poly}, whose possible
#'   value is \code{"mvp"} (default) or \code{"gmpoly"}, and this is the
#'   class of the polynomial returned by these functions. See the examples
#'   and the \href{https://github.com/stla/jackR#readme}{README} file.
#'
#' @importFrom JuliaConnectoR juliaSetupOk juliaCall juliaImport juliaGet juliaEval
#' @importFrom mvp mvp print.mvp
#' @importFrom gmpoly gmpoly
#' @importFrom gmp as.bigq
#'
#' @export
#'
#' @seealso \code{\link{as.function.exactmvp}}
#'
#' @note See \code{\link[JuliaConnectoR]{JuliaConnectoR-package}} for
#'   information about setting up Julia. If you want to directly use Julia,
#'   you can use \href{https://github.com/stla/JackPolynomials.jl}{my package}.
#'
#' @examples library(jack)
#' \donttest{if(JuliaConnectoR::juliaSetupOk()){
#'   julia <- Jack_julia()
#'   # numerical evaluation ####
#'   julia$Jack(x = c(2, 2/3), lambda = c(3, 1), alpha = 3/2)
#'   # to pass rational numbers, use strings:
#'   julia$Jack(x = c("2", "2/3"), lambda = c(3, 1), alpha = "3/2")
#'   # symbolic polynomials ####
#'   # for `JackPol`, you can pass a rational `alpha` as a string:
#'   ( pol <- julia$JackPol(m = 2, lambda = c(3, 1), alpha = "3/2") )
#'   class(pol)
#'   # you _must_ give `alpha` as a string if you choose `poly = "gmpoly"`
#'   julia$JackPol(m = 2, lambda = c(3, 1), alpha = "3/2", poly = "gmpoly")
#'   JuliaConnectoR::stopJulia()
#' }}
Jack_julia <- function(){
  if(!juliaSetupOk()){
    stop("Julia setup is not OK.")
  }
  toEval <- paste0(
    'if isnothing(Base.find_package("DynamicPolynomials")) ',
    'using Pkg; Pkg.add("DynamicPolynomials") ',
    'end'
  )
  juliaEval(toEval)
  module <- system.file("julia", "JackPolynomials.jl", package = "jack")
  . <- juliaCall("include", module)
  JackPolynomials <- juliaImport(".JackPolynomials", all = FALSE)
  Jack <- function(x, lambda, alpha = 2){
    rational <- FALSE
    if(is.character(alpha)){
      if(!is.character(x)){
        stop(
          "If you want to use a rational `alpha`, you have to use rational ",
          "numbers for the components of `x` as well."
        )
      }
      alpha <- rationalize(alpha)
      x <- lapply(x, rationalize)
      rational <- TRUE
    }
    if(is.character(x)){
      if(!is.character(alpha)){
        stop(
          "If you want to use rational numbers in `x`, you have to use a ",
          "rational number for `alpha` as well."
        )
      }
      alpha <- rationalize(alpha)
      x <- lapply(x, rationalize)
      rational <- TRUE
    }
    result <- JackPolynomials$Jack(
      unname(as.list(x)), asIntegerList(lambda), unname(alpha)
    )
    if(rational){
      result <- juliaGet(result)
      result <- as.bigq(result[["num"]], result[["den"]])
    }
    result
  }
  JackPol <- function(m, lambda, alpha, poly = "mvp"){
    poly <- match.arg(poly, c("mvp", "gmpoly"))
    rational <- FALSE
    if(is.character(alpha)){
      if(!grepl("^\\d+/\\d+$", alpha)){
        stop(
          "If you want to supply a rational `alpha`, ",
          "it must be a character string of the form `p/q`."
        )
      }
      alpha <- as.integer(strsplit(alpha, "/")[[1L]])
      rational <- TRUE
    }else{
      if(poly == "gmpoly"){
        stop(
          "If you want a `gmpoly` polynomial, you have ",
          "to supply a rational `alpha`; ",
          "it must be a character string of the form `p/q`."
        )
      }
    }
    J <- juliaGet(JackPolynomials$JackPolynomial(
      unname(as.integer(m)), asIntegerList(lambda), unname(alpha),
      TRUE
    ))
    if(poly == "mvp"){
      coefficients <- J[["coefficients"]]
      vars <- paste0("x_", seq_len(m))
      vars <- rep(list(vars), length(coefficients))
      powers <- J[["powers"]]
      poly <- mvp(vars, powers, coefficients)
      if(rational){
        coeffs <- vapply(J[["qcoefficients"]], function(f){
          den <- f[["den"]]
          if(den == 1L){
            as.character(f[["num"]])
          }else{
            paste0(f[["num"]], "/", den)
          }
        }, character(1L))
        attr(poly, "exact") <-
          rationalPolynomial(powers, coeffs, stars = TRUE)
        attr(poly, "nvars") <- m
        class(poly) <- c("exactmvp", class(poly))
      }
    }else{ # gmpoly
      powers <- do.call(rbind, J[["powers"]])
      coeffs <- vapply(J[["qcoefficients"]], function(f){
        paste0(f[["num"]], "/", f[["den"]])
      }, character(1L))
      poly <- gmpoly(coeffs = as.bigq(coeffs), powers = powers)
    }
    poly
  }
  Zonal <- function(x, lambda){
    if(rational <- is.character(x)){
      x <- lapply(x, rationalize)
    }
    result <- JackPolynomials$Zonal(
      unname(as.list(x)), asIntegerList(lambda)
    )
    if(rational){
      result <- juliaGet(result)
      result <- as.bigq(result[["num"]], result[["den"]])
    }
    result
  }
  ZonalPol <- function(m, lambda, poly = "mvp"){
    poly <- match.arg(poly, c("mvp", "gmpoly"))
    J <- juliaGet(JackPolynomials$ZonalPolynomial(
      unname(as.integer(m)), asIntegerList(lambda)
    ))
    if(poly == "mvp"){
      coefficients <- J[["coefficients"]]
      vars <- paste0("x_", seq_len(m))
      vars <- rep(list(vars), length(coefficients))
      powers <- J[["powers"]]
      poly <- mvp(vars, powers, coefficients)
      coeffs <- vapply(J[["qcoefficients"]], function(f){
        den <- f[["den"]]
        if(den == 1L){
          as.character(f[["num"]])
        }else{
          paste0(f[["num"]], "/", den)
        }
      }, character(1L))
      attr(poly, "exact") <-
        rationalPolynomial(powers, coeffs, stars = TRUE)
      attr(poly, "nvars") <- m
      class(poly) <- c("exactmvp", class(poly))
    }else{ # gmpoly
      powers <- do.call(rbind, J[["powers"]])
      coeffs <- vapply(J[["qcoefficients"]], function(f){
        paste0(f[["num"]], "/", f[["den"]])
      }, character(1L))
      poly <- gmpoly(coeffs = as.bigq(coeffs), powers = powers)
    }
    poly
  }
  ZonalQ <- function(x, lambda){
    if(rational <- is.character(x)){
      x <- lapply(x, rationalize)
    }
    result <- JackPolynomials$ZonalQ(
      unname(as.list(x)), asIntegerList(lambda)
    )
    if(rational){
      result <- juliaGet(result)
      result <- as.bigq(result[["num"]], result[["den"]])
    }
    result
  }
  ZonalQPol <- function(m, lambda, poly = "mvp"){
    poly <- match.arg(poly, c("mvp", "gmpoly"))
    J <- juliaGet(JackPolynomials$ZonalQPolynomial(
      unname(as.integer(m)), asIntegerList(lambda)
    ))
    if(poly == "mvp"){
      coefficients <- J[["coefficients"]]
      vars <- paste0("x_", seq_len(m))
      vars <- rep(list(vars), length(coefficients))
      powers <- J[["powers"]]
      poly <- mvp(vars, powers, coefficients)
      coeffs <- vapply(J[["qcoefficients"]], function(f){
        den <- f[["den"]]
        if(den == 1L){
          as.character(f[["num"]])
        }else{
          paste0(f[["num"]], "/", den)
        }
      }, character(1L))
      attr(poly, "exact") <-
        rationalPolynomial(powers, coeffs, stars = TRUE)
      attr(poly, "nvars") <- m
      class(poly) <- c("exactmvp", class(poly))
    }else{ # gmpoly
      powers <- do.call(rbind, J[["powers"]])
      coeffs <- vapply(J[["qcoefficients"]], function(f){
        paste0(f[["num"]], "/", f[["den"]])
      }, character(1L))
      poly <- gmpoly(coeffs = as.bigq(coeffs), powers = powers)
    }
    poly
  }
  Schur <- function(x, lambda){
    if(rational <- is.character(x)){
      x <- lapply(x, rationalize)
    }
    result <- JackPolynomials$Schur(
      unname(as.list(x)), asIntegerList(lambda)
    )
    if(rational){
      result <- juliaGet(result)
      result <- as.bigq(result[["num"]], result[["den"]])
    }
    result
  }
  SchurPol <- function(m, lambda, poly = "mvp"){
    poly <- match.arg(poly, c("mvp", "gmpoly"))
    J <- juliaGet(JackPolynomials$SchurPolynomial(
      unname(as.integer(m)), asIntegerList(lambda)
    ))
    if(poly == "mvp"){
      coefficients <- J[["coefficients"]]
      vars <- paste0("x_", seq_len(m))
      vars <- rep(list(vars), length(coefficients))
      powers <- J[["powers"]]
      poly <- mvp(vars, powers, coefficients)
      coeffs <- as.character(coefficients)
      attr(poly, "exact") <-
        rationalPolynomial(powers, coeffs, stars = TRUE)
      attr(poly, "nvars") <- m
      class(poly) <- c("exactmvp", class(poly))
    }else{ # gmpoly
      powers <- do.call(rbind, J[["powers"]])
      coeffs <- unlist(J[["coefficients"]])
      poly <- gmpoly(coeffs = as.bigq(coeffs), powers = powers)
    }
    poly
  }
  list(
    Jack = Jack,
    JackPol = JackPol,
    Zonal = Zonal,
    ZonalPol = ZonalPol,
    ZonalQ = ZonalQ,
    ZonalQPol = ZonalQPol,
    Schur = Schur,
    SchurPol = SchurPol
  )
}

#' @title Pretty exact expression
#' @description Pretty form of the exact expression of a polynomial.
#'
#' @param poly an \code{exactmvp} object, that is, a polynomial with an exact
#'   expression
#' @param asCharacter Boolean, whether to return a character string; if
#'   \code{FALSE}, the pretty form is printed
#'
#' @return A character string if \code{asCharacter=TRUE}, otherwise it is also
#'   returned but invisibly, and it is printed in the console.
#' @export
#'
#' @importFrom Ryacas yac_str
#'
#' @examples library(jack)
#' \donttest{if(JuliaConnectoR::juliaSetupOk()){
#'   julia <- Jack_julia()
#'   ( pol <- julia$ZonalPol(m = 2, lambda = c(3, 1)) )
#'   prettyForm(pol)
#'   JuliaConnectoR::stopJulia()
#' }}
prettyForm <- function(poly, asCharacter = FALSE){
  if(!inherits(poly, "exactmvp")){
    stop("The 'prettyForm' function is not applicable to this object.")
  }
  p <- yac_str(sprintf("PrettyForm(%s)", attr(attr(poly, "exact"), "yacas")))
  if(asCharacter){
    p
  }else{
    cat(p)
    invisible(p)
  }
}

#' @title Exact expression to LaTeX
#' @description LaTeX form of the exact expression of a polynomial.
#'
#' @param poly an \code{exactmvp} object, that is, a polynomial with an exact
#'   expression
#' @param asCharacter Boolean, whether to return a character string; if
#'   \code{FALSE}, the LaTeX code is printed
#'
#' @return A character string if \code{asCharacter=TRUE}, otherwise it is also
#'   returned but invisibly, and it is printed in the console.
#' @export
#'
#' @importFrom Ryacas yac_str
#'
#' @examples library(jack)
#' \donttest{if(JuliaConnectoR::juliaSetupOk()){
#'   julia <- Jack_julia()
#'   ( pol <- julia$ZonalQPol(m = 2, lambda = c(3, 2)) )
#'   toLaTeX(pol)
#'   JuliaConnectoR::stopJulia()
#' }}
toLaTeX <- function(poly, asCharacter = FALSE){
  if(!inherits(poly, "exactmvp")){
    stop("The 'toLaTeX' function is not applicable to this object.")
  }
  p <- yac_str(sprintf("TexForm(%s)", attr(attr(poly, "exact"), "yacas")))
  p <- gsub(" ^", "^", p, fixed = TRUE)
  if(asCharacter){
    p
  }else{
    cat(p)
    invisible(p)
  }
}
