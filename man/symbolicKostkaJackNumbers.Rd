% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Kostka.R
\name{symbolicKostkaJackNumbers}
\alias{symbolicKostkaJackNumbers}
\title{Symbolic Kostka-Jack numbers}
\usage{
symbolicKostkaJackNumbers(n)
}
\arguments{
\item{n}{positive integer, the weight of the partitions}
}
\value{
A named list of named lists of \code{ratioOfQsprays} objects.
  Denoting the Kostka numbers by \eqn{K_{\lambda,\mu}(\alpha)}, the names
  of the outer list correspond to the partitions \eqn{\lambda}, and the
  names of the inner lists correspond to the partitions \eqn{\mu}.
}
\description{
Kostka-Jack numbers with symbolic Jack parameter for partitions
  of a given weight.
}
\examples{
symbolicKostkaJackNumbers(3)
}
