#' @export
cond_cov <- function(init_covs, data) {

  Omega_E <- solve(init_covs$ind)
  Sigma_B <- init_covs$dam
  Sigma_A <- init_covs$sire

  dam_counts  <- split(data$n.observed$inds, data$sires)

  unique_count_vecs <- dam_counts %>%
    lapply(\(v) unname(sort(v))) %>%
    unique()

  unique_counts <- unique(unname(data$n.observed$inds))
  
  sire_counts <- sapply(dam_counts, sum)
  
  W     <- make_W(init_covs, data)
  D_inv <- paired_inverse(Omega_E, Sigma_B, unique_counts, E_type = "prec")

  CD_inv <- list()
  for(n in unique_counts) {
    CD_inv[[paste(n)]] <- Omega_E %*% D_inv[[paste(n)]] / n
  }

  top_blocks  <- list()
  main_blocks <- list()
  for(ns in unique_count_vecs) {
    top_blocks_curr <- list()
    main_blocks_l1  <- list()

    for(n in ns) {
      top_blocks_curr[[paste(n)]] <- -W(ns) %*% CD_inv[[paste(n)]]

      main_blocks_l2 <- list()
      for(m in ns) {
        if(m >= n) {
          main_blocks_l2[[paste(m)]] <- t(CD_inv[[paste(m)]] %*% W(ns) %*% CD_inv[[paste(n)]])

          ## if(m == n) {
          ##   main_blocks_l2[[paste(n)]] <- main_blocks_l2[[paste(n)]] + D_inv[[paste(n)]]
          ## }
        }
      }

      main_blocks_l1[[paste(n)]] <- main_blocks_l2
      
    }

    top_blocks[[toString(ns)]] <- top_blocks_curr
    main_blocks[[toString(ns)]] <- main_blocks_l1
  }

  # CREATE THE CLOSURE

  function(i, j, k) {
    ns <- dam_counts[[i]]

    ns_sort <- sort(ns)
    
    if(j == "group") {
      if(k == "group") {
        return(W(ns))
      } else {
        return(top_blocks[[toString(ns_sort)]][[paste(ns[k])]])
      }
    } else if (k == "group") {
      return(t(top_blocks[[toString(ns_sort)]][[paste(ns[j])]]))
    } else {
      n <- ns[j]
      m <- ns[k]
      
      if(n <= m) {
        out_block <- main_blocks[[toString(ns_sort)]][[paste(n)]][[paste(m)]]
      } else {
        out_block <- t(main_blocks[[toString(ns_sort)]][[paste(m)]][[paste(n)]])
      }

      if(j == k) {
        out_block <- out_block + D_inv[[paste(n)]]
      }

      return(out_block)
      
    }

  }
  
}

cond_mean <- function(cond_cov, data) {

  sire_sums <- rowsum(data$dam_sums, data$sires)
  
  # List indexed by sires with vectors of dam names
  dam_idxs <- split(names(data$sires), data$sires)

  sire_means <- matrix(
    nrow     = data$dims$I,
    ncol     = data$dims$q,
    dimnames = list(names(dam_idxs))
  )
  
  dam_means <- matrix(
    nrow     = nrow(data$dam_sums),
    ncol     = data$dims$q,
    dimnames = list(rownames(data$dam_sums))
  )
  
  for(sire in names(dam_idxs)) {

    dams <- dam_idxs[[sire]]

    sire_means[sire, ] <- cond_cov(sire, "group", "group") %*% sire_sums[sire, ]

    for(dam in dams) {
      sire_means[sire, ] <- sire_means[sire, ] + cond_cov(sire, "group", dam) %*% data$dam_sums[dam, ]
    }

    for(dam in dams) {

      dam_means[dam, ] <- cond_cov(sire, dam, "group") %*% sire_sums[sire, ]
      
      for(dam_mult in dams) {
        dam_means[dam, ] <- dam_means[dam, ] + cond_cov(sire, dam, dam_mult) %*% data$dam_sums[dam_mult, ]
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
make_W <- function(init_covs, data) {

  Sigma_E <- init_covs$ind
  Sigma_B <- init_covs$dam
  Sigma_A <- init_covs$sire
  
  # Get the unique (up to re-ordering) values of (|O_i|)_i present in the data
  ob_counts <- split(data$n.observed$inds, data$sires)

  unique_count_vecs <- ob_counts %>%
    lapply(\(v) unname(sort(v))) %>%
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

  U_invQ_inv <- t(Q) %*% U

  
  D <- list()
  for(n in unique_counts) {
    D[[paste(n)]] <- diag(1 / (d + 1/n))
  }

  A <- U_invQ_inv %*% solve(Sigma_A, t(U_invQ_inv))

  W_list <- list()
  for(n_vec in unique_count_vecs) {
    D_curr <- Reduce(`+`, D[paste(n_vec)])

    W_list[[toString(n_vec)]] <- t(U_invQ_inv) %*% solve(A + D_curr, U_invQ_inv)
  }

  function(n_vec) {
    W_list[[toString(sort(n_vec))]]
  }
}

