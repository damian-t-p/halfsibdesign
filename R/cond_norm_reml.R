#' Conjugate each block of cond_cov by a matrix
block_conjugate <- function(cond_cov, conj_mat, one_side = TRUE) {

  # We will be re-assigning variables within the cond_cov environment, so create a cloned
  # copy of the closure to leave the original unaffected
  new_env                 <- rlang::env_clone(rlang::fn_env(cond_cov))
  rlang::fn_env(cond_cov) <- new_env

  if (isTRUE(one_side)) {
    conj_by <- function(M) conj_mat %*% M
  } else {
    conj_by <- function(M) conj_mat %*% M %*% conj_mat
  }

  assign("W_blocks",
         rapply(new_env$W_blocks, conj_by, how = "replace"),
         envir = new_env)

  assign("top_blocks",
         rapply(new_env$top_blocks, conj_by, how = "replace"),
         envir = new_env)

  assign("main_blocks",
         rapply(new_env$main_blocks, conj_by, how = "replace"),
         envir = new_env)

  assign("D_inv",
         rapply(new_env$D_inv, conj_by, how = "replace"),
         envir = new_env)

  return(cond_cov)
}

cond_cov_reml_conj <- function(init_covs, cond_cov, data) {
  
  Omega_E       <- solve(init_covs$ind)
  cond_cov_conj <- block_conjugate(cond_cov, Omega_E)
  
  n_obs_total <- sum(data$n.observed$inds)
  W_inv       <- n_obs_total * Omega_E

  dam_names <- split(names(data$sires), data$sires)
  
  for(sire in names(dam_names)) {
    n_obs_sire <- sum(data$n.observed$inds[dam_names[[sire]]])

    W_inv <- W_inv - n_obs_sire^2 * cond_cov_conj(sire, "group", "group")

    for(dam in dam_names[[sire]]) {

      n_obs_dam <- data$n.observed$inds[dam]
      
      W_inv <- W_inv - n_obs_dam * n_obs_sire * cond_cov_conj(sire, "group", dam)
      W_inv <- W_inv - n_obs_dam * n_obs_sire * cond_cov_conj(sire, dam, "group")

      for(dam2 in dam_names[[sire]]) {

        n_obs_dam2 <- data$n.observed$inds[dam2]
        
        W_inv <- W_inv <- n_obs_dam * n_obs_dam2 * cond_cov_conj(sire, dam, dam2)
      }
    }
  }

  W <- solve(W_inv)

  return(W)
  
}

cond_cov_reml <- function(init_covs, cond_cov, data) {
  
  Sigma_E     <- init_covs$ind
  
  n_obs_total <- sum(data$n.observed$inds)
  W_inv       <- n_obs_total * Sigma_E

  dam_names <- split(names(data$sires), data$sires)
  
  for(sire in names(dam_names)) {
    n_obs_sire <- sum(data$n.observed$inds[dam_names[[sire]]])

    W_inv <- W_inv - n_obs_sire^2 * cond_cov(sire, "group", "group")

    for(dam in dam_names[[sire]]) {

      n_obs_dam <- data$n.observed$inds[dam]
      
      W_inv <- W_inv - n_obs_dam * n_obs_sire * cond_cov(sire, "group", dam)
      W_inv <- W_inv - n_obs_dam * n_obs_sire * cond_cov(sire, dam, "group")

      for(dam2 in dam_names[[sire]]) {

        n_obs_dam2 <- data$n.observed$inds[dam2]
        
        W_inv <- W_inv - n_obs_dam * n_obs_dam2 * cond_cov(sire, dam, dam2)
      }
    }
  }

  W <- Sigma_E %*% solve(W_inv) %*% Sigma_E

  return(W)
  
}
