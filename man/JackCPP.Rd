% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CPP.R
\name{JackCPP}
\alias{JackCPP}
\title{Evaluation of Jack polynomial - C++ implementation}
\usage{
JackCPP(x, lambda, alpha)
}
\arguments{
\item{x}{variables, a vector of \code{bigq} numbers, or a vector that can
be coerced as such (e.g. \code{c("2", "5/3")})}

\item{lambda}{an integer partition, given as a vector of decreasing
integers}

\item{alpha}{positive rational number, given as a string such as
\code{"2/3"} or as a \code{bigq} number}
}
\value{
A \code{bigq} number.
}
\description{
Evaluates the Jack polynomial.
}
\examples{
JackCPP(c("1", "3/2", "-2/3"), lambda = c(3, 1), alpha = "1/4")
}
