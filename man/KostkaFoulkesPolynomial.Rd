% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HallLittlewood.R
\name{KostkaFoulkesPolynomial}
\alias{KostkaFoulkesPolynomial}
\title{Kostka-Foulkes polynomial}
\usage{
KostkaFoulkesPolynomial(lambda, mu)
}
\arguments{
\item{lambda, mu}{integer partitions; in order for the Kostka-Foulkes
polynomial to be non-zero, a necessary condition is that \code{lambda}
and \code{mu} have the same weight; more precisely, \code{mu} must
be dominated by \code{lambda}}
}
\value{
The Kostka-Foulkes polynomial associated to \code{lambda} and
  \code{mu}. This is a univariate \code{qspray} polynomial whose value
  at \code{1} is the Kostka number associated to \code{lambda} and
  \code{mu}.
}
\description{
Kostka-Foulkes polynomial for two given partitions.
}
