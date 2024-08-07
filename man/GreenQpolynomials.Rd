% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GreenPolynomials.R
\name{GreenQpolynomials}
\alias{GreenQpolynomials}
\title{Green Q-polynomials}
\usage{
GreenQpolynomials(rho)
}
\arguments{
\item{rho}{an integer partition}
}
\value{
The Green Q-polynomials, usually denoted by \eqn{Q_\rho^\lambda}, 
  for all integer partitions \eqn{\lambda} of the same weight as the 
  given integer partition \eqn{\rho}. 
  They are returned in a list. Each element of this list is itself a list, 
  with two elements. The first one, called \code{lambda}, represents the 
  partition \eqn{\lambda}. The second one, called \code{polynomial}, 
  represents the Green Q-polynomial \eqn{Q_\rho^\lambda}. This is a 
  univariate \code{qspray} polynomial whose variable is denoted by 
  \code{q}.
  The names of the list encode the partitions \eqn{\lambda}.
}
\description{
Computes the Green Q-polynomials for a given integer partition.
}
\note{
The Green Q-polynomials are the "true" Green polynomials. 
  The Green X-polynomials (\code{\link{GreenXpolynomials}}) are a 
  variant of the Green Q-polynomials.
}
\examples{
GreenQpolynomials(c(2, 1)) 
}
\seealso{
\code{\link{GreenXpolynomials}}.
}
