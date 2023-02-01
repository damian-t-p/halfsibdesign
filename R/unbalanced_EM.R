# y_data should be a list indexed by i, whose elements are O_i * p
# matrices, where each row is an observation y_ij

#'
#' @export
EM_oneway <- function(y_data,
                      Sigma_E_init, Sigma_A_init,
                      method   = c("REML", "ML"),
                      max_iter = 1000,
                      err.tol  = 1e-6) {

  method = match.arg(method)
  
  Sigma_E <- Sigma_E_init
  Sigma_A <- Sigma_A_init

  if(method == "ML") {
    mu <- rowMeans(sapply(y_data$tables, colMeans))
  } else {
    mu <- rep(0, nrow(Sigma_E_init))
  }
  
  for(iter in 1:max_iter) {

    # M step
    cond_params   <- alpha_cond_params(y_data, Sigma_E, Sigma_A, mu)
    balanced_data <- balance_data(y_data, cond_params)
    ss_mats_base  <- ss_oneway(balanced_data)

    if(method == "ML") {
      mu <- rowMeans(sapply(balanced_data$tables, colMeans))
    } else {
      mu <- rep(0, nrow(Sigma_E_init))
    }
    
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
    degf_E <- I * (J - 1)
    degf_A <- I - (method == "REML")
    curr_primal <- stepreml_1way(A_E, degf_E, A_A, degf_A)$primal
    
    Sigma_E = curr_primal[[1]]
    Sigma_A = (curr_primal[[2]] - curr_primal[[1]])/J

    # check for convergence
    if(iter > 1) {
      err <- mat_err(prev_primal, curr_primal, list(I*(J-1), I-1))
      if(err < err.tol) {break}
    } else {
      err <- NA
    }

    prev_primal <- curr_primal
  }

  list(
    Sigma_E = curr_primal[[1]],
    Sigma_A = (curr_primal[[2]] - curr_primal[[1]])/J
  )
  
}

n_observed <- function(y_data) {
  sapply(y_data$tables, nrow)
}

#' Use svd to compute a square root of the non-negative definite symmetric matrix A
svdsqrt <- function(A) {
  decomp <- svd(A, nu=0)

  V <- decomp$v
  d <- decomp$d

  stopifnot(all(d > -1e-6))
  
  return(t(V) * sqrt(abs(d)))
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
  ## U <- Re(expm::sqrtm(Sigma_A))
  U <- svdsqrt(Sigma_A)
  W <- U %*% solve(Sigma_E, t(U))

  W_eigen <- eigen(W, symmetric = TRUE)

  UtP <- t(U) %*% W_eigen$vectors

  inv_mats <- list()
  for(n in ns) {
    inv_mats[[paste(n)]] <- UtP %*% diag(1/(W_eigen$values * n + 1)) %*% t(UtP)
  }

  return(inv_mats)
}

## paired_inverse <- function(Sigma_E, Sigma_A, ns) {
##   U <- chol(Sigma_E)
##   W <- U %*% solve(Sigma_A, t(U))

##   W_eigen <- eigen(W, symmetric = TRUE)

##   UtP <- t(U) %*% W_eigen$vectors

##   inv_mats <- list()
##   for(n in ns) {
##     inv_mats[[paste(n)]] <- UtP %*% diag(1/(W_eigen$values + n)) %*% t(UtP)
##   }

##   return(inv_mats)
## }

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
#' @param mu Global mean
#'
#' @return A list of multivariate Gaussian parameters for the conditional
#' distribution of the sires random effect. Each entry is a list with entries
#' `mean` and `cov`, which are the mean vector and covariance matrix respectively.
#'
#' @export
alpha_cond_params <- function(y_data, Sigma_E, Sigma_A, mu = rep(0, nrow(Sigma_E))) {
  observed_n   <- unique(n_observed(y_data))
  covariances  <- paired_inverse(Sigma_E, Sigma_A, observed_n)
  
  params <- list()

  for(i in 1:length(y_data$tables)) {

    n   <- nrow(y_data$tables[[i]])
    cov <- covariances[[paste(n)]]

    params[[i]] <- list(
      mean = mu + cov %*% solve(Sigma_E, colSums(y_data$tables[[i]]) - n * mu),
      cov  = cov
    )
  }

  return(params)
}
