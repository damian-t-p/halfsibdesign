cond_cov_mats <- function(init_covs, data, flat_sire = FALSE) {

  Omega_E <- solve(init_covs$ind)
  Sigma_B <- init_covs$dam
  Sigma_A <- init_covs$sire

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

  ## Compute matrices required for block inversion
  W     <- make_W(init_covs, data, flat_sire)
  D_inv <- paired_inverse(Omega_E, Sigma_B, unique_counts, E_type = "prec")

  CD_inv <- list()
  for(n in unique_counts) {
    CD_inv[[paste(n)]] <- n * Omega_E %*% D_inv[[paste(n)]]
  }

  ## Build tables for the blocks of the conditional covariance using the block inversion formula

  # The blocks (1, 1)
  # Indexing format:
  #   $sort(n[i1], ..., n[iJ])
  #     - a list for each combinations of observation numbers present among the sires
  W_blocks <- list()
  
  # The blocks (1, 2), ..., (1, J + 1) running along the top of the covariance matrix for each i
  # Indexing format:
  #   $sort(n[i1], ..., n[iJ])
  #     - a list for each combinations of observation numbers present among the sires
  #   $sort(n[i1], ..., n[iJ])$n[ij]
  #     - a q*q block for each n[ij] present in (n[i1], ..., n[iJ])
  #     - conditional covariance between alpha[i] and beta[ij]
  top_blocks  <- list()

  # The block in the main pars of the matrix, lying between indices (2, 2) and (J+1, J+1)
  # Indexing format:
  #   $sort(n[i1], ..., n[iJ])
  #     - a list for each combinations of observation numbers present among the sires
  #   $sort(n[i1], ..., n[iJ])$n[ij]
  #     - a list for each n[ij] present in (n[i1], ..., n[iJ])
  #   $sort(n[i1], ..., n[iJ])$n[ij]$n[ij']
  #     - a q*q block for each n[ij'] >= n[ij] present in (n[i1], ..., n[iJ])
  #     - conditional covariance between beta[ij] and beta[ij']
  #     - excludes the additional D[ij]^(-1) present in the full covariance matrix, as this
  #       is different from the covariance between to dam effects with the same observation
  #       numbers. This value must be manually added in the closure function.
  main_blocks <- list()

  
  for(ns in unique_count_vecs) {
    top_blocks_curr <- list()
    main_blocks_l1  <- list()

    for(n in ns) {
      top_blocks_curr[[paste(n)]] <- -W(ns) %*% CD_inv[[paste(n)]]

      main_blocks_l2 <- list()
      for(m in ns) {
        if(m >= n) {
          main_blocks_l2[[paste(m)]] <- t(CD_inv[[paste(n)]]) %*% W(ns) %*% CD_inv[[paste(m)]]
        }
      }

      main_blocks_l1[[paste(n)]] <- main_blocks_l2
      
    }

    top_blocks[[to_str(ns)]]  <- top_blocks_curr
    main_blocks[[to_str(ns)]] <- main_blocks_l1
    W_blocks[[to_str(ns)]]    <- W(ns)
  }

  list(
    top_blocks  = top_blocks,
    main_blocks = main_blocks,
    W_blocks    = W_blocks,
    D_inv       = D_inv
  )
  
}

cond_cov_counts <- function(init_covs, data, flat_sire = FALSE) {

  mats <- cond_cov_mats(init_covs, data, flat_sire)

  top_blocks  <- mats$top_blocks
  main_blocks <- mats$main_blocks
  W_blocks    <- mats$W_blocks
  D_inv       <- mats$D_inv

  dam_counts  <- split(data$n.observed$inds[names(data$sires)], data$sires)
  
  # The closure function looks up the relevant entry in `top_blocks` and `main_blocks`,
  # additionally adding D[ij]^(-1) as required along the block diagonal.
  #
  # i, j, k should all be sire and dam labels, not integers
  function(ns, n, m, equal_idx = FALSE) {
    
    if(n == "group") {
      # correlation with the sire effect along first row
      
      if(m == "group") {
        return(W_blocks[[ns]])
      } else {
        return(top_blocks[[ns]][[m]])
      }
    } else if (m == "group") {
      # correlation with the sire effect along the first column
      
      return(t(top_blocks[[ns]][[n]]))
      
    } else {
      # correlations between dam effects
      
      if(as.numeric(n) <= as.numeric(m)) {
        out_block <- main_blocks[[ns]][[n]][[m]]
      } else {
        out_block <- t(main_blocks[[ns]][[m]][[n]])
      }

      # Add additional blocks along main block diagonal
      if(equal_idx) {
        out_block <- out_block + D_inv[[n]]
      }

      return(out_block)
      
    }
  }
}

#' Compute the conditional covariances of two-way MANOVA random effects
#'
#' Under the model `y[ijk] == mu + a[i] + beta[ij] + epsilon[ijk]`, where
#' each `alpha[i]`, `beta[ij]` and `epsilon[ijk]` are independent mean-0,
#' `q`-dimensional normal random vectors with with covariance matrices
#' `Sigma[A]`, `Sigma[B]` and `Sigma[E]` respectively, compute the covariance
#' matrices of `(alpha[i], beta[i1], ..., beta[iJ])` conditional on the
#' observed data for each `i`.
#' 
#' @param init_covs A list of prior covariance matrices with entries named
#' `ind`, `dam` and `sire` indicating the between-individuals, between-dams
#' and between-sires covariances,
#' @param data An object inheriting `halfsibdata`
#' @param flat_sire A logical indicating whether a flat prior should be used
#' for the sire effect.
#'
#' @return A closure that takes the following arguments:
#'  - `i`: The name of sire
#'  - `j`, `k`: The name of dam, or `"group"`
#'
#' The closure function returns the posterior covariance between the random
#' effects `beta[ij]` and `beta[ik]`. If `j == "group"`, returns the posterior
#' covariance of `alpha[i]` and `beta[ik]`. and similarly if `j == "group"`.
#' 
#' @export
cond_cov <- function(ccov_counts, data, ...) {

  dam_counts  <- split(data$n.observed$inds[names(data$sires)], data$sires)
  
  function(sire, dam1, dam2) {

    sire_ns_num <- dam_counts[[sire]]
    sire_ns     <- to_str(sort.int(sire_ns_num))
    
    if(dam1 == "group") {
      dam1_n <- "group"
    } else {
      dam1_n <- paste(sire_ns_num[dam1])
    }

    if(dam2 == "group") {
      dam2_n <- "group"
    } else {
      dam2_n <- paste(sire_ns_num[dam2])
    }

    ccov_counts(
      sire_ns, dam1_n, dam2_n,
      equal_idx  = (dam1 == dam2),
      ...
    )

  }
}

#' Compute the conditional mean of two-way MANOVA random effects
#' 
#' Under the model `y[ijk] == mu + a[i] + beta[ij] + epsilon[ijk]`, where
#' each `alpha[i]`, `beta[ij]` and `epsilon[ijk]` are independent mean-0,
#' `q`-dimensional normal random vectors with with covariance matrices
#' `Sigma[A]`, `Sigma[B]` and `Sigma[E]` respectively, compute the means
#' of `(alpha[i], beta[i1], ..., beta[iJ])` conditional on the observed data
#' for each `i`.
#'
#' @param init_covs A list of prior covariances. Must have an entry `$ind`.
#' @param cond_cov A function that returns conditional covariance matrices as
#' created by `halfsibdesign::cond_cov`
#' @param data An object inheriting `halfsibdata`
#' @param prior_mean A vector of the prior global mean
#'
#' @return A list with entries `sire` and `dam` whose rows are the posterior
#' means of `alpha[i]` and `beta[ij]` respectively.
#' 
#' @export
cond_mean <- function(init_covs, ccov, data, prior_mean = rep(0, data$dims$q)) {

  Omega_E <- solve(init_covs$ind)
  
  # centre dam sums by prior mean
  dam_sums  <- data$dam_sums - data$n.observed$inds %o% prior_mean

  # Skew sire and dam sums by precision
  dam_skew  <- dam_sums %*% Omega_E
  sire_skew <- rowsum(dam_skew, data$sires)
  
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

    # Indices of dams whose sire is `sire`
    dams <- dam_idxs[[sire]]
    
    sire_means[sire, ] <- sire_means[sire, ] + ccov(sire, "group", "group") %*% sire_skew[sire, ]
    for(dam in dams) {
      sire_means[sire, ] <- sire_means[sire, ] + ccov(sire, "group", dam) %*% dam_skew[dam, ]
    }
    
    for(dam in dams) {

      dam_means[dam, ] <- dam_means[dam, ] + ccov(sire, dam, "group") %*% sire_skew[sire, ]
      for(dam2 in dams) {
        dam_means[dam, ] <- dam_means[dam, ] + ccov(sire, dam, dam2) %*% dam_skew[dam2, ]
      }
        
    }
    
  }

  list(
    sire = sire_means,
    dam  = dam_means
  )
}

#' Compute (Sigma[A] + sum((Sigma[B] + Sigma[E]/n[j])^(-1)))^(-1) for
#' values of n[j] given in `data`
#' 
#' @export
make_W <- function(init_covs, data, flat_sire = FALSE) {

  Sigma_E <- init_covs$ind
  Sigma_B <- init_covs$dam
  Sigma_A <- init_covs$sire

  q <- data$dims$q
  
  # Get the unique (up to re-ordering) values of (|O_i|)_i present in the data
  ob_counts <- split(data$n.observed$inds[names(data$sires)], data$sires)

  unique_count_vecs <- ob_counts %>%
    lapply(\(v) unname(sort.int(v))) %>%
    unique()

  unique_counts <- unique(unname(data$n.observed$inds))
  
  # Compute (Sigma[B] + Sigma[E] / n[j])^(-1) for each n[j]
  # Sigma[E] = U' U
  decomp <- eigen(Sigma_E, symmetric = TRUE)
  V <- decomp$vectors
  d <- decomp$values

  stopifnot(all(d > -1e-6))
  U     <- t(V) * sqrt(abs(d))
  U_inv <- t(t(V) * 1/sqrt(abs(d)))
  
  # inv(U') Sigma[B] inv(U) = Q D Q'
  eigs <- eigen(t(U_inv) %*% Sigma_B %*% U_inv, symmetric = TRUE)
  d    <- eigs$values
  Q    <- eigs$vectors

  ## U_invQ_inv <- t(Q) %*% U
  U_invQ <- U_inv %*% Q
  
  D <- list()
  for(n in unique_counts) {
    D[[paste(n)]] <- diag(1 / (d + 1/n))
  }

  W_list <- list()
  for(n_vec in unique_count_vecs) {
    D_curr <- Reduce(`+`, D[paste(n_vec)])
    
    if(isTRUE(flat_sire)) {
      W_list[[to_str(n_vec)]] <- solve(U_invQ %*% D_curr %*% t(U_invQ))
    } else {
      W_list[[to_str(n_vec)]] <- solve(diag(q) + Sigma_A %*% U_invQ %*% D_curr %*% t(U_invQ), Sigma_A)
    }
  }

  function(n_vec) {
    W_list[[to_str(sort.int(n_vec))]]
  }
}

