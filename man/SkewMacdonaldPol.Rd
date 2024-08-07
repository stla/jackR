% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SkewMacdonaldPolynomial.R
\name{SkewMacdonaldPol}
\alias{SkewMacdonaldPol}
\title{Skew Macdonald polynomial}
\usage{
SkewMacdonaldPol(n, lambda, mu, which = "P")
}
\arguments{
\item{n}{number of variables, a positive integer}

\item{lambda, mu}{integer partitions defining the skew partition:
\code{lambda} is the outer partition and \code{mu} is the inner partition
(so \code{mu} must be a subpartition of \code{lambda})}

\item{which}{which skew Macdonald polynomial, \code{"P"}, \code{"Q"}
or \code{"J"}}
}
\value{
A \code{symbolicQspray} multivariate polynomial, the skew
  Macdonald polynomial associated to the skew partition defined by
  \code{lambda} and \code{mu}. It has two parameters usually
  denoted by \eqn{q} and \eqn{t}. Substituting \eqn{q} with \eqn{0}
  yields the skew Hall-Littlewood polynomials.
}
\description{
Returns the skew Macdonald polynomial associated to
  the given skew partition.
}
