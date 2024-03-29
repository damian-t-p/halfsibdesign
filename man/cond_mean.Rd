% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cond_norm.R
\name{cond_mean}
\alias{cond_mean}
\title{Compute the conditional mean of two-way MANOVA random effects}
\usage{
cond_mean(init_covs, cond_cov, data, prior_mean = rep(0, data$dims$q))
}
\arguments{
\item{init_covs}{A list of prior covariances. Must have an entry \verb{$ind}.}

\item{cond_cov}{A function that returns conditional covariance matrices as
created by \code{halfsibdesign::cond_cov}}

\item{data}{An object inheriting \code{halfsibdata}}

\item{prior_mean}{A vector of the prior global mean}
}
\value{
A list with entries \code{sire} and \code{dam} whose rows are the posterior
means of \code{alpha[i]} and \code{beta[ij]} respectively.
}
\description{
Under the model \code{y[ijk] == mu + a[i] + beta[ij] + epsilon[ijk]}, where
each \code{alpha[i]}, \code{beta[ij]} and \code{epsilon[ijk]} are independent mean-0,
\code{q}-dimensional normal random vectors with with covariance matrices
\code{Sigma[A]}, \code{Sigma[B]} and \code{Sigma[E]} respectively, compute the means
of \verb{(alpha[i], beta[i1], ..., beta[iJ])} conditional on the observed data
for each \code{i}.
}
