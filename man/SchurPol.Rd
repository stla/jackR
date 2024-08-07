% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CPP.R
\name{SchurPol}
\alias{SchurPol}
\title{Schur polynomial - C++ implementation}
\usage{
SchurPol(n, lambda)
}
\arguments{
\item{n}{number of variables, a positive integer}

\item{lambda}{an integer partition, given as a vector of decreasing
integers}
}
\value{
A \code{qspray} multivariate polynomial.
}
\description{
Returns a Schur polynomial. The Schur polynomials are the
  Jack \eqn{P}-polynomials with Jack parameter \eqn{\alpha=1}.
}
\examples{
( schur <- SchurPol(3, lambda = c(3, 1)) )
schur == JackPol(3, lambda = c(3, 1), alpha = "1", which = "P")
}
