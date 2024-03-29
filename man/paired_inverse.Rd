% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unbalanced_EM.R
\name{paired_inverse}
\alias{paired_inverse}
\title{Compute (A^(-1) + n E^(-1))^(-1) for a vector of ns}
\usage{
paired_inverse(Sigma_E, Sigma_A, ns, E_type = c("cov", "prec"))
}
\arguments{
\item{ns}{A vector of doubles}

\item{E_type}{If this equals \code{"prec"}, the first argument is instead \code{E^(-1)}}

\item{A, E}{Symmetric nonnegative-definite square matrices}
}
\value{
A list of inverted matrices indexed by the vector \code{ns}
}
\description{
Compute (A^(-1) + n E^(-1))^(-1) for a vector of ns
}
