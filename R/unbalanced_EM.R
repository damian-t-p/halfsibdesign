# y_data should be a list indexed by i, whose elements are O_i * p
# matrices, where each row is an observation y_ij

EM_oneway <- function(y_data, Sigma_E_init, Sigma_A_init) {
}

n_observed <- function(y_data) {
  sapply(y_data, nrow)
}

#' Compute (A^(-1) + n E^(-1))^(-1) for a vector of ns
#'
#' @param A,E Symmetric nonnegative-definite square matrices
#' @param ns A vector of doubles
#'
#' @return A list of inverted matrices indexed by the vector `ns`
#'
#' @export
paired_inverse <- function(Sigma_E, Sigma_A, ns) {
  U <- chol(Sigma_E)
  W <- U %*% solve(Sigma_A, t(U))

  W_eigen <- eigen(W, symmetric = TRUE)

  UtP <- t(U) %*% W_eigen$vectors

  inv_mats <- list()
  for(n in ns) {
    inv_mats[[paste(n)]] <- UtP %*% diag(1/(W_eigen$values + n)) %*% t(UtP)
  }

  return(inv_mats)
}

#' Compute conditional distribution of sires effect
#'
#' Given the model
#' \deqn{
#' y_{ij} | \alpha_i \sim N(\alpha_i, \Sigma_E)
#' \alpha_i \sim N(0, \Sigma_A),
#' }{
#' y[ij] | alpha[i] ~ N(alpha[i], Sigma[E])
#' alpha[i] ~ N(0, Sigma[A]),
#' }
#' computes the parameters of the Gaussian conditional distribution
#' \deqn{
#' \alpha_i | y_i.
#' }{
#' alpha[i] | y_i.
#' }
#'
#' @param y_data Observed data that inherits "`fullsibdata`"
#' @param Sigma_E Between-individuals covariance matrix
#' @param Sigma_A Between-sires covariance matrix
#'
#' @return A list of multivariate Gaussian parameters for the conditional
#' distribution of the sires random effect. Each entry is a list with entries
#' `mean` and `cov`, which are the mean vector and covariance matrix respectively.
#'
#' @export
alpha_cond_params <- function(y_data, Sigma_E, Sigma_A) {
  observed_n   <- unique(n_observed(y_data))
  covariances <- paired_inverse(Sigma_E, Sigma_A, observed_n)
  
  params <- list()

  for(i in 1:length(y_data)) {

    n   <- nrow(y_data[[i]])
    cov <- covariances[[paste(n)]]

    params[[i]] <- list(
      mean = cov %*% solve(Sigma_E, colSums(y_data[[i]])),
      cov  = cov
    )
  }

  return(params)
}
