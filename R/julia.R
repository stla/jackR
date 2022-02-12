#' @title Evaluation with Julia
#' @description Evaluate the Jack polynomials with Julia. This is highly faster.
#'
#' @return A list of functions having the same names as the R functions of this
#'   package (\code{Jack}, \code{JackPol}, \code{Schur}, etc).
#'
#' @importFrom JuliaConnectoR juliaSetupOk juliaCall juliaImport juliaGet
#' @import mvp mvp
#' @export
#'
#' @note See \code{\link[JuliaConnectoR]{JuliaConnectoR-package}} for
#'   information about setting up Julia. If you want to directly use Julia,
#'   you can use \href{https://github.com/stla/JackPolynomials.jl}{my package}.
#'
#' @examples library(jack)
#' \donttest{if(JuliaConnectoR::juliaSetupOk()){
#'   JuliaConnectoR::stopJulia()
#' }}
Jack_julia <- function(){
  if(!juliaSetupOk()){
    stop("Julia setup is not OK.")
  }
  module <- system.file("julia", "JackPolynomials.jl", package = "jack")
  . <- juliaCall("include", module)
  JackPolynomials <- juliaImport(".JackPolynomials", all = FALSE)
  Jack <- function(x, lambda, alpha = 2){
    JackPolynomials$Jack(
      unname(as.list(x)), unname(as.list(as.integer(lambda))), unname(alpha)
    )
  }
  JackPol <- function(m, lambda, alpha = 2){
    J <- juliaGet(JackPolynomials$JackPolynomial(
      unname(as.integer(m)), unname(as.list(as.integer(lambda))), unname(alpha)
    ))
    coefficients <- J[["coefficients"]]
    vars <- paste0("x", seq_len(m))
    vars <- rep(list(vars), length(coefficients))
    mvp(vars, J[["powers"]], coefficients)
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
    mvp(vars, J[["powers"]], coefficients)
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
    mvp(vars, J[["powers"]], coefficients)
  }
  list(
    Jack = Jack,
    JackPol = JackPol,
    Zonal = Zonal,
    ZonalPol = ZonalPol,
    Schur = Schur,
    SchurPol = SchurPol
  )
}
