#' @export
EM_fit <- function(data, ...) {
  UseMethod("EM_fit", data)
}

#' @export
EM_fit.halfsibdata <- function(data,
                               prior_covs = list(
                                 ind  = diag(1, nrow = data$dims$q),
                                 dam  = diag(1, nrow = data$dims$q),
                                 sire = diag(1, nrow = data$dims$q)
                               ),
                               method     = c("ML", "REML", "ML_nofix"),
                               flat_sire  = (method == "ML_nofix"),
                               max_iter   = 1000,
                               err.tol    = 1e-6) {

  method <- match.arg(method)  

  if(!is.dam_balanced(data)) {
    data <- include_unobs_dams(data)
  }
  
  I <- data$dims$I
  J <- data$dims$J
  K <- data$dims$K
  
  n_missing <- K - data$n.observed$inds

  if(method == "ML") {
    # Might have 0 observed individuals per dam, in which case the mean should be 0.
    mu <- colMeans(data$dam_sums / pmax(data$n.observed$inds, 1))
  } else {
    mu <- rep(0, data$dims$q)
  }
  
  for(iter in 1:max_iter) {

    ccov  <- cond_cov(prior_covs, data, flat_sire = flat_sire)
    if(method == "REML") {
      ccov_raw  <- cond_cov_counts(prior_covs, data)
      ccov_reml <- cond_cov_reml(prior_covs, ccov, ccov_raw, data)
      cmean     <- cond_mean_reml(prior_covs, ccov_reml, data)
    } else {
      cmean <- cond_mean(prior_covs, ccov, data, prior_mean = mu)
    }

    balanced_data <- balance(data, cmean, globmean = mu)

    if(method == "ML") {
      mu <- colMeans(balanced_data$dam_sums) / data$dims$K
    }
    
    # Naive sum-of-squares M-matrices
    ss_base <- ss_mats.halfsibdata(balanced_data)

    M_E <- ss_base$m_ind + sum(n_missing) * (1 - 1/K) * prior_covs$ind * 1/(I * J * (K-1))
    M_B <- ss_base$m_dam + sum(n_missing) * (1 - 1/J) * prior_covs$ind * 1/(I * (J-1) * K)
    M_A <- ss_base$m_sire + sum(n_missing) * (1 - 1/I) * prior_covs$ind * 1/((I-1) * J * K)

    sire_names <- names(data$sires)
    dam_names  <- split(names(data$sires), data$sires)

    if(method == "REML") {
      ccomp_reml <- function(sire1, sire2, dam1, dam2) {
        ccov_reml("group", "group", NA, NA) +
          ccov_reml("group", sire2, NA, "group") + ccov_reml(sire1, "group", "group", NA) +
          ccov_reml("group", sire2, NA, dam2) + ccov_reml(sire2, "group", dam2, NA) +
          ccov_reml(sire1, sire2, "group", "group") +
          ccov_reml(sire1, sire2, "group", dam2) + ccov_reml(sire1, sire2, dam1, "group") +
          ccov_reml(sire1, sire2, dam1, dam2)
      }

      ccomp <- function(sire, dam1, dam2) {
        ccov_reml(sire, sire, dam1, dam2)
      }
    } else {
      ccomp <- function(sire, dam1, dam2) {
        ccov(sire, "group", "group") + ccov(sire, "group", dam2) +
          ccov(sire, dam1, "group") + ccov(sire, dam1, dam2)
      }
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
            1/(J * K) * 1/(I * (J-1))

          M_A <- M_A +
            n_missing[[dam]] * n_missing[[dam2]] * ccomp(sire, dam, dam2) *
            (1 - (method != "REML")/I) * 1/(J * K) * 1/(I-1)
        }
      }
    }

    if(method == "REML") {
      for(sire in sire_names) {
        for(sire2 in sire_names) {
          for(dam in dam_names[[sire]]) {
            for(dam2 in dam_names[[sire2]]) {
              M_A <- M_A -
                n_missing[[dam]] * n_missing[[dam2]] * ccomp_reml(sire, sire2, dam, dam2) *
                (1 - 1/I) * 1/(J * K) * 1/(I-1)
            }
          }
        }
      }
    }

    if (I == 1) {
      M_A <- matrix(0, nrow = data$dims$q, ncol = data$dims$q)
    }
    
    # E step
    if (method %in% c("ML", "ML_nofix")) {
      E_method <- "ML"
    } else {
      E_method <- method
    }
    curr_primal <- stepreml_2way_mat(M_E, K, M_B, J, M_A, I,
                                     log_crit = "never",
                                     method = E_method)

    prior_covs <- list(
      ind = curr_primal$S1,
      dam = curr_primal$S2,
      sire = curr_primal$S3
    )

    ## print("-----------------------------------")
    ## print(prior_covs)
    
    if(iter > 1) {
      err <- mat_err(prev_primal, curr_primal, list(I * J * (K-1), I * (J-1), I - (E_method == "ML")))
      if(err < err.tol) {break}
    } else {
      err <- NA
    }

    prev_primal <- curr_primal
  }

  out_covs <- rlang::list2(
    !!data$level_names$ind_name  := prior_covs$ind,
    !!data$level_names$dam_name  := prior_covs$dam,
    !!data$level_names$sire_name := prior_covs$sire
  )

  if(isTRUE(flat_sire)) {
    out_covs[[3]] <- NULL
  }
  
  lapply(
    out_covs,
    \(A) {dimnames(A) <- dimnames(M_E); A}
  )
}
