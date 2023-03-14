cov_idx <- function(data, method = "REML") {
  I <- data$dims$I
  J <- data$dims$J
  K <- data$dims$K
  q <- data$dims$q

  sorted_sires <- 1:I
  names(sorted_sires) <- sort(names(data$n.observed$dams))

  sorted_dams <- rep(1:J, times = I)
  names(sorted_dams)  <- names(data$sires[order(data$sires, names(data$sires))])

  function(sire, dam) {
    if(sire == "group") {
      if(method == "ML") {
        stop("Group random effect not estimated")
      } else {
        init_idx <- 1
      }
    } else if(dam == "group") {
      init_idx <- (sorted_sires[sire] - 1) * (J + 1) + 1 + (method == "REML")
    } else {
      init_idx <- (sorted_sires[sire] - 1) * (J + 1) + 1 + (method == "REML") + sorted_dams[dam]
    }

    (init_idx - 1) * q + 1:q
      
  }
}

full_prec <- function(init_covs, data, method = "REML") {

  I <- data$dims$I
  J <- data$dims$J
  K <- data$dims$K
  q <- data$dims$q

  Omega_E <- solve(init_covs$ind)
  Omega_A <- solve(init_covs$sire)
  Omega_B <- solve(init_covs$dam)
  
  prec <- matrix(
    0,
    nrow = q * (I * (J + 1) + (method == "REML")),
    ncol = q * (I * (J + 1) + (method == "REML"))
  )

  dam_names <- split(names(data$sires), data$sires)

  idx <- cov_idx(data, method = method)
  if(method == "REML") {
    prec[idx("group"), idx("group")] <- sum(data$n.observed$inds) * Omega_E
  }
  
  for(sire in names(dam_names)) {
    n_obs_sire <- sum(data$n.observed$inds[dam_names[[sire]]])

    if(method == "REML") {
      prec[idx("group"), idx(sire, "group")] <- n_obs_sire * Omega_E
      prec[idx(sire, "group"), idx("group")] <- n_obs_sire * Omega_E
    }
    
    prec[idx(sire, "group"), idx(sire, "group")] <- Omega_A + n_obs_sire * Omega_E
    
    for(dam in dam_names[[sire]]) {
      n_obs_dam <- data$n.observed$inds[dam] 

      if(method == "REML") {
        prec[idx("group"), idx(sire, dam)] <- n_obs_dam * Omega_E
        prec[idx(sire, dam), idx("group")] <- n_obs_dam * Omega_E
      }
      
      prec[idx(sire, "group"), idx(sire, dam)] <- n_obs_dam * Omega_E
      prec[idx(sire, dam), idx(sire, "group")] <- n_obs_dam * Omega_E

      prec[idx(sire, dam), idx(sire, dam)] <- Omega_B + n_obs_dam * Omega_E
    }
    
  }

  return(prec)
  
}


cond_cov_counts_reml <- function(init_covs, ccov_raw, ccov = cond_cov(ccov_raw), data) {
  
  Sigma_E     <- init_covs$ind
  
  n_obs_total <- sum(data$n.observed$inds)
  W_inv       <- n_obs_total * Sigma_E

  dam_names <- split(names(data$sires), data$sires)
  
  for(sire in names(dam_names)) {
    n_obs_sire <- sum(data$n.observed$inds[dam_names[[sire]]])

    W_inv <- W_inv - n_obs_sire^2 * ccov(sire, "group", "group")

    for(dam in dam_names[[sire]]) {

      n_obs_dam <- data$n.observed$inds[dam]
      
      W_inv <- W_inv - n_obs_dam * n_obs_sire * ccov(sire, "group", dam)
      W_inv <- W_inv - n_obs_dam * n_obs_sire * ccov(sire, dam, "group")

      for(dam2 in dam_names[[sire]]) {

        n_obs_dam2 <- data$n.observed$inds[dam2]
        
        W_inv <- W_inv - n_obs_dam * n_obs_dam2 * ccov(sire, dam, dam2)
      }
    }
  }

  W <- solve(W_inv)

  # List of numbers of individuals observed per dam indexed by sire
  dam_counts  <- split(data$n.observed$inds[names(data$sires)], data$sires)

  # All observation count vectors that are distince up to reordering
  unique_count_vecs <- dam_counts %>%
    lapply(\(v) unname(sort.int(v))) %>%
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

    ns_str <- to_str(ns)
    
    DC_sire[[ns_str]] <- sum(ns) * ccov_raw(ns_str, "group", "group")
    DC_dam_curr <- list()

    for(n in ns) {
      DC_sire[[ns_str]] <- DC_sire[[ns_str]] +
        n * ccov_raw(ns_str, "group", paste(n))
      
      DC_dam_curr[[paste(n)]] <- sum(ns) * ccov_raw(ns_str, paste(n), "group")
      
      for(m in ns) {
        DC_dam_curr[[paste(n)]] <- DC_dam_curr[[paste(n)]] +
          m * ccov_raw(ns_str, paste(n), paste(m))
      }

      # Add the repeated term in the DC multiplication
      DC_dam_curr[[paste(n)]] <- DC_dam_curr[[paste(n)]] -
        n * ccov_raw(ns_str, paste(n), paste(n)) +
        n * ccov_raw(ns_str, paste(n), paste(n), equal_idx = TRUE)
    }

    DC_dam[[ns_str]] <- DC_dam_curr
    
  }  

  function(sire1_ns, sire2_ns, dam1_n, dam2_n, equal_sire_idx = FALSE, equal_dam_idx = FALSE) {

    if(sire1_ns == "group") {
      if(sire2_ns == "group") {
        return(Sigma_E %*% W %*% Sigma_E)
      } else {
        if(dam2_n == "group") {
          return(-Sigma_E %*% W %*% t(DC_sire[[sire2_ns]]))
        } else {
          return(-Sigma_E %*% W %*% t(DC_dam[[sire2_ns]][[dam2_n]]))
        }
      }
    } else if(sire2_ns == "group") {
      if(dam1_n == "group") {
        return(-DC_sire[[sire1_ns]] %*% t(W) %*% Sigma_E)
      } else {
        return(-DC_dam[[sire1_ns]][[dam1_n]] %*% t(W) %*% Sigma_E)
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

        out <- DC_sire[[sire1_ns]] %*% W %*% t(DC_dam[[sire2_ns]][[dam2_n]])
        if(isTRUE(equal_sire_idx)) {
          out <- out + ccov_raw(sire1_ns, "group", dam2_n)
        }
        return(out)
      
      }

    } else {

      if(dam2_n == "group") {
        
        out <- DC_dam[[sire1_ns]][[dam1_n]] %*% W %*% t(DC_sire[[sire2_ns]])
        if(isTRUE(equal_sire_idx)) {
          out <- out + ccov_raw(sire1_ns, dam1_n, "group")
        }
        return(out)
        
      } else {

        out <- DC_dam[[sire1_ns]][[dam1_n]] %*% W %*% t(DC_dam[[sire2_ns]][[dam2_n]])
        if(isTRUE(equal_sire_idx)) {
          out <- out + ccov_raw(sire1_ns, dam1_n, dam2_n, equal_idx = equal_dam_idx)
        }
        return(out)


      }
      
    }
  }
  
}

cond_cov_reml <- function(init_covs, ccov_raw, ccov, data) {

  dam_counts  <- split(data$n.observed$inds[names(data$sires)], data$sires)
  ccov_counts <- cond_cov_counts_reml(init_covs, ccov_raw, ccov, data)
  
  function(sire1, sire2, dam1, dam2) {

    if(sire1 == "group") {
      sire1_ns <- "group"
      dam1_n   <- NA
    } else {
      sire1_ns_num <- dam_counts[[sire1]]
      sire1_ns     <- to_str(sort.int(sire1_ns_num))
      
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
      sire2_ns     <- to_str(sort.int(sire2_ns_num))
      
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

cond_mean_reml <- function(init_covs, ccov_reml, data) {

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

  glob_mean <- ccov_reml("group", "group", NA, NA) %*% grand_skew

  for(sire2 in names(dam_idxs)) {
    glob_mean <- glob_mean + ccov_reml("group", sire2, NA, "group") %*% sire_skew[sire2, ]

    for(dam2 in dam_idxs[[sire2]]) {
      glob_mean <- glob_mean + ccov_reml("group", sire2, NA, dam2) %*% dam_skew[dam2, ]
    }
  }

  for(sire in names(dam_idxs)) {
    
    sire_means[sire, ] <- sire_means[sire, ] +
      ccov_reml(sire, "group", "group", NA) %*% grand_skew

    for(dam in dam_idxs[[sire]]) {
      dam_means[dam, ] <- dam_means[dam, ] +
        ccov_reml(sire, "group", dam, NA) %*% grand_skew
    }
    
    for(sire2 in names(dam_idxs)) {
      sire_means[sire, ] <- sire_means[sire, ] +
        ccov_reml(sire, sire2, "group", "group") %*% sire_skew[sire2, ]

      for(dam in dam_idxs[[sire]]) {
        dam_means[dam, ] <- dam_means[dam, ] +
          ccov_reml(sire, sire2, dam, "group") %*% sire_skew[sire2, ]
      }

      for(dam2 in dam_idxs[[sire2]]) {
        sire_means[sire, ] <- sire_means[sire, ] +
          ccov_reml(sire, sire2, "group", dam2) %*% dam_skew[dam2, ]

        for(dam in dam_idxs[[sire]]) {
          dam_means[dam, ] <- dam_means[dam, ] +
            ccov_reml(sire, sire2, dam, dam2) %*% dam_skew[dam2, ]
        }
        
      }
    }
     
  }

  list(
    grand = c(glob_mean),
    sire  = sire_means,
    dam   = dam_means
  )
}
