#' @export
EM_fit.halfsibdata <- function(data,
                               prior_covs,
                               method   = c("REML", "ML"),
                               max_iter = 1000,
                               err.tol  = 1e-6) {

  method <- match.arg(method)

  if(method == "ML") {
    mu <- colMeans(data$dam_sums) / data$dims$K
  } else {
    mu <- rep(0, data$dims$q)
  }

  I <- data$dims$I
  J <- data$dims$J
  K <- data$dims$K
  
  n_missing <- K - data$n.observed$inds
  
  for(iter in 1:max_iter) {
    ccov  <- cond_cov(prior_covs, data)
    cmean <- cond_mean(ccov, data, prior_mean = mu)

    balanced_data <- balance(data, cmean, globmean = mu)

    # Naive sum-of-squares M-matrices
    ss_base <- ss_mats.halfsibdata(balanced_data)

    M_E <- ss_base$m_ind + sum(n_missing) * (1 - 1/K) * prior_covs$ind * 1/(I * J * (K-1))
    M_B <- ss_base$m_dam + sum(n_missing) * (1 - 1/J) * prior_covs$ind * 1/(I * (J-1) * K)
    M_A <- ss_base$m_sire + sum(n_missing) * (1 - 1/I) * prior_covs$ind * 1/((I-1) * J * K)

    sire_names <- names(data$sires)
    dam_names  <- split(names(data$sires), data$sires)

    ccomp <- function(sire, dam1, dam2) {
      ccov(sire, "group", "group") + ccov(sire, "group", dam2) +
        ccov(sire, dam1, "group") + ccov(sire, dam1, dam2)
    }
    
    for(sire in sire_names) {
      for(dam in dam_names[[sire]]) {
        M_E <- M_E +
          n_missing[[dam]] * (1 - n_missing[[dam]]/K) * ccomp(sire, dam, dam) *
          1/(I * J * (K-1))
        
        M_B <- M_B +
          n_missing[[dam]]^2 * ccomp(sire, dam, dam) *
          1/K * 1 /(I * (J-1)) 

        for(dam2 in dam_names[[sire]]) {
          M_B <- M_B -
            n_missing[[dam]] * n_missing[[dam2]] * ccomp(sire, dam, dam2) *
            1/(J * K) * 1/(J-1)

          M_A <- M_A +
            n_missing[[dam]] * n_missing[[dam2]] * ccomp(sire, dam, dam2) *
            (1 - 1/I) * 1/(J * K) * 1/(I - 1)
        }
      }

      # E step
      curr_primal <- stepreml_2way_mat(M_E, K, M_B, J, M_A, I + (method == "ML"), log_crit = "never")

      prior_covs <- list(
        ind = curr_primal$S1,
        dam = curr_primal$S2,
        sire = curr_primal$S3
      )
      
      if(iter > 1) {
        err <- mat_err(prev_primal, curr_primal, list(I * J * (K-1), I * J, I - (method == "REML")))
        if(err < err.tol) {break}
      } else {
        err <- NA
      }

      prev_primal <- curr_primal
    }

    return(prior_covs)
  }
  
}
