rationalMonomial <- function(variables, powers){
  paste0(
    paste0(variables, "^", powers, c(rep(" ", length(powers)-1L), "")),
    collapse = ""
  )
}

rationalPolynomial <- function(variables, powers, coeffs, stars = FALSE){
  monomials <- mapply(rationalMonomial, variables, powers, USE.NAMES = FALSE)
  spaces <- rep(" ", length(coeffs))
  ones <- coeffs == "1"
  minusones <- coeffs == "-1"
  spaces[ones] <- ""
  coeffs[ones] <- ""
  spaces[minusones] <- ""
  coeffs[minusones] <- "-"
  if(stars){
    gsub("(\\d) x", "\\1 * x",
         gsub("+  -", "-  ",
              paste0(
                paste0(coeffs, spaces, gsub("^1", "", monomials, fixed = TRUE)),
                collapse = "  +  "
              ),
              fixed = TRUE
         )
    )
  }else{
    gsub("+ -", " - ",
         paste0(
           paste0(coeffs, spaces, gsub("^1", "", monomials, fixed = TRUE)),
           collapse = "+ "
         ),
         fixed = TRUE
    )
  }
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
#'   JuliaConnectoR::stopJulia()
#' }}
as.function.exactmvp <- function(x, ...){
  expr <- attr(x, "exact")
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

#' @title Evaluation with Julia
#' @description Evaluate the Jack polynomials with Julia. This is highly faster.
#'
#' @return A list of functions having the same names as the R functions of this
#'   package (\code{Jack}, \code{JackPol}, \code{Schur}, etc).
#'
#' @importFrom JuliaConnectoR juliaSetupOk juliaCall juliaImport juliaGet juliaEval
#' @importFrom mvp mvp print.mvp
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
#'   # for `JackPol`, you can pass a rational `alpha` as a string:
#'   ( pol <- julia$JackPol(m = 2, lambda = c(3, 1), alpha = "3/2") )
#'   class(pol)
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
    JackPolynomials$Jack(
      unname(as.list(x)), unname(as.list(as.integer(lambda))), unname(alpha)
    )
  }
  JackPol <- function(m, lambda, alpha){
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
    }
    J <- juliaGet(JackPolynomials$JackPolynomial(
      unname(as.integer(m)), unname(as.list(as.integer(lambda))), unname(alpha),
      TRUE
    ))
    coefficients <- J[["coefficients"]]
    vars <- paste0("x", seq_len(m))
    vars <- rep(list(vars), length(coefficients))
    poly <- mvp(vars, J[["powers"]], coefficients)
    if(rational){
      variables <- poly[["names"]]
      powers <- poly[["power"]]
      qcoefficients <- J[["qcoefficients"]]
      coeffs <- vapply(qcoefficients, function(f){
        den <- f[["den"]]
        if(den == 1L){
          as.character(f[["num"]])
        }else{
          paste0(f[["num"]], "/", den)
        }
      }, character(1L))
      attr(poly, "exact") <-
        rationalPolynomial(variables, powers, coeffs, stars = TRUE)
      attr(poly, "nvars") <- m
      class(poly) <- c("exactmvp", class(poly))
    }
    poly
  }
  Zonal <- function(x, lambda){
    JackPolynomials$Zonal(
      unname(as.list(x)), unname(as.list(as.integer(lambda)))
    )
  }
  ZonalPol <- function(m, lambda){
    J <- juliaGet(JackPolynomials$ZonalPolynomial(
      unname(as.integer(m)), unname(as.list(as.integer(lambda)))
    ))
    coefficients <- J[["coefficients"]]
    vars <- paste0("x", seq_len(m))
    vars <- rep(list(vars), length(coefficients))
    poly <- mvp(vars, J[["powers"]], coefficients)
    variables <- poly[["names"]]
    powers <- poly[["power"]]
    qcoefficients <- J[["qcoefficients"]]
    coeffs <- vapply(qcoefficients, function(f){
      den <- f[["den"]]
      if(den == 1L){
        as.character(f[["num"]])
      }else{
        paste0(f[["num"]], "/", den)
      }
    }, character(1L))
    attr(poly, "exact") <- rationalPolynomial(variables, powers, coeffs)
    attr(poly, "nvars") <- m
    class(poly) <- c("exactmvp", class(poly))
    poly
  }
  ZonalQ <- function(x, lambda){
    JackPolynomials$ZonalQ(
      unname(as.list(x)), unname(as.list(as.integer(lambda)))
    )
  }
  ZonalQPol <- function(m, lambda){
    J <- juliaGet(JackPolynomials$ZonalQPolynomial(
      unname(as.integer(m)), unname(as.list(as.integer(lambda)))
    ))
    coefficients <- J[["coefficients"]]
    vars <- paste0("x", seq_len(m))
    vars <- rep(list(vars), length(coefficients))
    poly <- mvp(vars, J[["powers"]], coefficients)
    variables <- poly[["names"]]
    powers <- poly[["power"]]
    qcoefficients <- J[["qcoefficients"]]
    coeffs <- vapply(qcoefficients, function(f){
      den <- f[["den"]]
      if(den == 1L){
        as.character(f[["num"]])
      }else{
        paste0(f[["num"]], "/", den)
      }
    }, character(1L))
    attr(poly, "exact") <- rationalPolynomial(variables, powers, coeffs)
    attr(poly, "nvars") <- m
    class(poly) <- c("exactmvp", class(poly))
    poly
  }
  Schur <- function(x, lambda){
    JackPolynomials$Schur(
      unname(as.list(x)), unname(as.list(as.integer(lambda)))
    )
  }
  SchurPol <- function(m, lambda){
    J <- juliaGet(JackPolynomials$SchurPolynomial(
      unname(as.integer(m)), unname(as.list(as.integer(lambda)))
    ))
    coefficients <- J[["coefficients"]]
    vars <- paste0("x", seq_len(m))
    vars <- rep(list(vars), length(coefficients))
    poly <- mvp(vars, J[["powers"]], coefficients)
    variables <- poly[["names"]]
    powers <- poly[["power"]]
    coeffs <- as.character(coefficients)
    attr(poly, "exact") <- rationalPolynomial(variables, powers, coeffs)
    attr(poly, "nvars") <- m
    class(poly) <- c("exactmvp", class(poly))
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
