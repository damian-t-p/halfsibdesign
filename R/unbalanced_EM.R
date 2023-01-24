# y_data should be a list indexed by i, whose elements are O_i * p
# matrices, where each row is an observation y_ij

EM_oneway <- function(y_data, Sigma_E_init, Sigma_A_init,
                      max_iter = 1000,
                      err.tol  = 1e-6) {

  ## Sigma_E <- Sigma_E_init
  ## Sigma_A <- Sigma_A_init

  prev_covs <- list(Sigma_E_init, Sigma_A_init)
  
  for(iter in 1:max_iter) {

    Sigma_E <- prev_covs[[1]]
    Sigma_A <- prev_covs[[2]]

    # M step
    cond_params   <- alpha_cond_params(y_data, Sigma_E, Sigma_A)
    balanced_data <- balance_data(y_data, cond_params)
    ss_mats_base  <- ss_oneway(balanced_data)

    J         <- y_data$n_ind
    I         <- y_data$n_sires
    n_missing <- sapply(y_data$tables, \(X) {J - nrow(X)})

    # compute sum-of-squares matrices
    A_E <- ss_mats_base$A_E + sum(n_missing) * (1 - 1/J) * Sigma_E
    A_A <- ss_mats_base$A_A + sum(n_missing) * (1 - 1/I)/J * Sigma_E

    for(i in 1:length(y_data$tables)) {
      A_E <- A_E + n_missing[i] * (1 - n_missing[i]/J) * cond_params[[i]]$cov
      A_A <- A_A + n_missing[i]^2 * (1 - 1/I)/J * cond_params[[i]]$cov
    }

    # E_step
    curr_covs <- stepreml_1way(A_E, I*(J-1), A_A, I - 1)$primal

    # check for convergence
    if(iter > 1) {
      err <- mat_err(prev_covs, curr_covs, list(I*(J-1), I-1))
      if(err < err.tol) {break}
    } else {
      err <- NA
    }

    prev_covs <- curr_covs
  }

  list(
    Sigma_E = curr_covs[[1]],
    Sigma_A = curr_covs[[2]]
  )
  
}

n_observed <- function(y_data) {
  sapply(y_data$tables, nrow)
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
  covariances  <- paired_inverse(Sigma_E, Sigma_A, observed_n)
  
  params <- list()

  for(i in 1:length(y_data$tables)) {

    n   <- nrow(y_data$tables[[i]])
    cov <- covariances[[paste(n)]]

    params[[i]] <- list(
      mean = cov %*% solve(Sigma_E, colSums(y_data$tables[[i]])),
      cov  = cov
    )
  }

  return(params)
}
