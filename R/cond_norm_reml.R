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

cond_cov_counts_reml <- function(init_covs, cond_cov, ccov_raw, data) {
  
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

  W_orig <- solve(W_inv)
  W_right <- W_orig %*% Sigma_E
  W <- Sigma_E %*% W_orig %*% Sigma_E


  # List of numbers of individuals observed per dam indexed by sire
  dam_counts  <- split(data$n.observed$inds[names(data$sires)], data$sires)

  # All observation count vectors that are distince up to reordering
  unique_count_vecs <- dam_counts %>%
    lapply(\(v) unname(sort(v))) %>%
    unique()

  # All observation counts per dat
  unique_counts <- unique(unname(data$n.observed$inds))

  # All observation counts per sire
  sire_counts <- sapply(dam_counts, sum)
  
  # Indexing format
  #  $sort(n[i1], ..., n[iJ])
  #    $n[ij]
  DC_sire <- list()
  
  # Indexing format
  #  $sort(n[i1], ..., n[iJ])
  #    $n[ij]
  DC_dam <- list()

  for(ns in unique_count_vecs) {

    DC_sire[[toString(ns)]] <- sum(ns) * ccov_raw(toString(ns), "group", "group")
    DC_dam_curr <- list()

    for(n in ns) {
      DC_sire[[toString(ns)]] <- DC_sire[[toString(ns)]] +
        n * ccov_raw(toString(ns), "group", paste(n))
      
      DC_dam_curr[[paste(n)]] <- n * ccov_raw(toString(ns), paste(n), "group")
      
      for(m in ns) {
        DC_dam_curr[[paste(n)]] <- DC_dam_curr[[paste(n)]] +
          m * ccov_raw(toString(ns), paste(n), paste(m))
      }

      DC_dam[[toString(ns)]] <- DC_dam_curr
    }
    
  }

  function(sire1_ns, sire2_ns, dam1_n, dam2_n, equal_sire_idx = FALSE, equal_dam_idx = FALSE) {

    if(sire1_ns == "group") {
      if(sire2_ns == "group") {
        return(W_orig)
      } else {
        if(dam2_n == "group") {
          return(-W_right %*% t(DC_sire[[sire2_ns]]))
        } else {
          return(-W_right %*% t(DC_dam[[sire2_ns]][[dam2_n]]))
        }
      }
    } else if(sire2_ns == "group") {
      if(dam1_n == "group") {
        return(-DC_sire[[sire1_ns]] %*% t(W_right))
      } else {
        return(-DC_dam[[sire1_ns]][[dam1_n]] %*% t(W_right))
      }
    }
    
    if(dam1_n == "group") {

      if(dam2_n == "group") {

        out <- DC_sire[[sire1_ns]] %*% W %*% t(DC_sire[[sire2_ns]])
        if(isTRUE(equal_sire_idx)) {
          out <- out + ccov_raw(sire1_ns, "group", "group")
        }
        
        return(out)

      } else {

        return(DC_sire[[sire1_ns]] %*% W %*% t(DC_dam[[sire2_ns]][[dam2_n]]))
        
      }

    } else {

      if(dam2_n == "group") {

        return(DC_dam[[sire1_ns]][[dam1_n]] %*% W %*% t(DC_sire[[sire2_ns]]))
        
      } else {

        out <- DC_dam[[sire1_ns]][[dam1_n]] %*% W %*% t(DC_dam[[sire2_ns]][[dam2_n]])
        if(isTRUE(equal_sire_idx) & isTRUE(equal_dam_idx)) {
          out <- out + ccov_raw(sire1_ns, dam1_n, dam1_n, equal_idx = TRUE)
        }

        return(out)

      }
      
    }
  }
  
}

cond_cov_reml <- function(init_covs, cond_cov, ccov_raw, data) {

  dam_counts  <- split(data$n.observed$inds[names(data$sires)], data$sires)
  
  ccov_counts <- cond_cov_counts_reml(init_covs, cond_cov, ccov_raw, data)
  
  function(sire1, sire2, dam1, dam2) {

    if(sire1 == "group") {
      sire1_ns <- "group"
      dam1_n   <- NA
    } else {
      sire1_ns_num <- dam_counts[[sire1]]
      sire1_ns     <- toString(sort(sire1_ns_num))
      
      if(dam1 == "group") {
        dam1_n <- "group"
      } else {
        dam1_n <- paste(sire1_ns_num[dam1])
      }
    }

    if(sire2 == "group") {
      sire2_ns <- "group"
      dam2_n   <- NA
    } else {
      sire2_ns_num <- dam_counts[[sire2]]
      sire2_ns     <- toString(sort(sire2_ns_num))
      
      if(dam2 == "group") {
        dam2_n <- "group"
      } else {
        dam2_n <- paste(sire2_ns_num[dam2])
      }
    }

    ccov_counts(
      sire1_ns, sire2_ns, dam1_n, dam2_n,
      equal_sire_idx = (sire1 == sire2),
      equal_dam_idx  = (dam1 == dam2)
    )

  }
  
}


cond_mean_reml <- function(init_covs, cond_cov, data) {

  Omega_E <- solve(init_covs$ind)

  dam_skew   <- data$dam_sums %*% Omega_E
  sire_skew  <- rowsum(dam_skew, data$sires)
  grand_skew <- colSums(sire_skew)

  # list indexed by sires with vectors of dam names
  dam_idxs <- split(names(data$sires), data$sires)
  
  # Empty matrices to be populated by outputs
  sire_means <- matrix(
    0,
    nrow     = data$dims$I,
    ncol     = data$dims$q,
    dimnames = list(names(dam_idxs))
  )
  
  dam_means <- matrix(
    0,
    nrow     = nrow(data$dam_sums),
    ncol     = data$dims$q,
    dimnames = list(rownames(data$dam_sums))
  )

  for(sire in names(dam_idxs)) {
    
    sire_means[sire, ] <- sire_means[sire, ] +
      cond_cov(sire, "group", "group", NA) %*% grand_skew

    for(dam in dam_idxs[[sire]]) {
      dam_means[dam, ] <- dam_means[dam, ] +
        cond_cov(sire, "group", dam, NA) %*% dam_skew[dam, ]
    }
    
    for(sire2 in names(dam_idxs)) {
      sire_means[sire, ] <- sire_means[sire, ] +
        cond_cov(sire, sire2, "group", "group") %*% sire_skew[sire2, ]

      for(dam in dam_idxs[[sire]]) {
        dam_means[dam, ] <- dam_means[dam, ] +
          cond_cov(sire, sire2, dam, "group") %*% sire_skew[sire2, ]
      }

      for(dam2 in dam_idxs[[sire2]]) {
        sire_means[sire, ] <- sire_means[sire, ] +
          cond_cov(sire, sire2, "group", dam2) %*% dam_skew[dam2, ]

        for(dam in dam_idxs[[sire]]) {
          dam_means[dam, ] <- dam_means[dam, ] +
            cond_cov(sire, sire2, dam, dam2) %*% dam_skew[dam2, ]
        }
        
      }
    }
     
  }

  list(
    sire = sire_means,
    dam  = dam_means
  )
}
