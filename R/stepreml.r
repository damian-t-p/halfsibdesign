#' Compute 2-way REML estimates with a stepwise method
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' REML estimates for the covariances using a stepwise algorithm.
#'
#' The stepwise algorithm is described in \[1\]. Its convergence is proved in
#' that paper. Although the speed of convergence is not established therein, in
#' practice it is linear.
#'
#' @name stepreml
#'
#' @return A list with entries `S1`, `S2` and `S3`, which are the estimates of
#' the sires, dams and individual covariance matrices.
#'
#' @references \[1\] J. A. Calvin and R. L. Dykstra. "Maximum Likelihood Estimation
#' of a Set of Covariance Matrices Under Lowner Order Restrictions with Applications
#' to Balanced Multivariate Variance Components Models". In: *The Annals of Staistics*
#' 19.2 (1991).
#'
#' @seealso \code{\link{stepreml_2way}} for generic version.
NULL
#> NULL

stepreml_1way <- function(A1, n1, A2, n2) {

  # sample covariances
  S1 <- A1/n1
  S2 <- A2/n2

  # S1 = U' U
  decomp <- eigen(S1, symmetric = TRUE)
  V <- decomp$vectors
  d <- decomp$values

  stopifnot(all(d > -1e-6))
  U     <- t(V) * sqrt(abs(d))
  U_inv <- t(t(V) * 1/sqrt(abs(d)))

  # inv(U') S2 inv(U) = Q D Q'
  eigs <- eigen(t(U_inv) %*% S2 %*% U_inv, symmetric = TRUE)
  d <- eigs$values
  Q <- eigs$vectors

  # B S1 B' = I
  # B S2 B' = D
  B_inv <- t(U) %*% Q

  c_hat = 1 - n2/(n1 + n2) * pmax(1 - d, 0)
  d_hat = d + n1/(n1 + n2) * pmax(1 - d, 0)

  # primal solutions
  Sig1_hat <- B_inv %*% diag(c_hat, nrow = length(c_hat)) %*% t(B_inv)
  Sig2_hat <- B_inv %*% diag(d_hat, nrow = length(d_hat)) %*% t(B_inv)

  # dual solutions
  Psi1_hat <- (n1/2) * Sig1_hat - A1/2
  Psi2_hat <- (n2/2) * Sig2_hat - A2/2

  list(
    primal = list(Sig1_hat, Sig2_hat),
    dual = list(Psi1_hat, Psi2_hat)
  )
}

#' @rdname stepreml
#' @param M1,M2,M3 Sum-of-squares matrices within individual, dam and sire groups.
#' @param I1,I2,I3 Number of individuals per dam, dams per sire and sires.
#' @param max_iter Number of iterations allowed before termination.
#' @param err.tol Algorithm terminates when the differences between subsequent
#' matrices in the matric described in \[1\] fall below this threshold.
#' @param verbose If this is `TRUE`, the output will contain a field
#' `conv_details` containing information about the algorithm's convergence.
#' @param log_crit When should the REML criterion be computed. If `"output"`, it will
#' be computed for the final fitted value. If `"always"`, will be computed at each
#' iteration step.
#' 
#' @export
stepreml_2way_mat <- function(M1, I1, M2, I2, M3, I3,
                              max_iter = 50,
                              err.tol  = 1e-6,
                              verbose  = FALSE,
                              log_crit = c("output", "never", "always")) {

  log_crit <- match.arg(log_crit)

  n1 <- (I1-1)*I2*I3
  n2 <- (I2-1)*I3
  n3 <- I3-1

  A1 <- M1 * n1
  A2 <- M2 * n2
  A3 <- M3 * n3

  p <- nrow(A1)

  zero <- matrix(0, nrow = p, ncol = p)

  A_curr <- list(A1, A2, A3)
  n_list <- list(n1, n2, n3)

  primal_curr <- list(zero, zero, zero)
  primal_prev <- list(zero, zero, zero)

  dual_1_prev <- list(zero, zero, zero)
  dual_2_prev <- list(zero, zero, zero)

  conv_df <- tibble::tibble(
    iter       = integer(),
    est_diff   = numeric(),
    reml_score = numeric()
  )

  for(i in 1:max_iter) {

    # S1 < S2
    A_curr[[1]] <- A_curr[[1]] + 2*dual_1_prev[[1]] - 2*dual_2_prev[[1]]
    A_curr[[2]] <- A_curr[[2]] + 2*dual_1_prev[[2]] - 2*dual_2_prev[[2]]
    A_curr[[3]] <- A_curr[[3]] + 2*dual_1_prev[[3]] - 2*dual_2_prev[[3]]


    constr1_sol <- stepreml_1way(A_curr[[1]], n1 , A_curr[[2]], n2)


    dual_2_prev <- dual_1_prev

    dual_1_prev[[1]] <- constr1_sol$dual[[1]]
    dual_1_prev[[2]] <- constr1_sol$dual[[2]]
    dual_1_prev[[3]] <- zero

    primal_curr[[1]] <- constr1_sol$primal[[1]]
    primal_curr[[2]] <- constr1_sol$primal[[2]]

    # S2 < S3
    A_curr[[1]] <- A_curr[[1]] + 2*dual_1_prev[[1]] - 2*dual_2_prev[[1]]
    A_curr[[2]] <- A_curr[[2]] + 2*dual_1_prev[[2]] - 2*dual_2_prev[[2]]
    A_curr[[3]] <- A_curr[[3]] + 2*dual_1_prev[[3]] - 2*dual_2_prev[[3]]

    constr2_sol <- stepreml_1way(A_curr[[2]], n2, A_curr[[3]], n3)

    dual_2_prev <- dual_1_prev

    dual_1_prev[[1]] <- zero
    dual_1_prev[[2]] <- constr2_sol$dual[[1]]
    dual_1_prev[[3]] <- constr2_sol$dual[[2]]

    primal_curr[[2]] <- constr2_sol$primal[[1]]
    primal_curr[[3]] <- constr2_sol$primal[[2]]

    # update REML
    S1 <- primal_curr[[1]]
    S2 <- (primal_curr[[2]] - primal_curr[[1]])/I1
    S3 <- (primal_curr[[3]] - primal_curr[[2]])/(I1*I2)

    # check for convergence
    if(i > 1) {
      err <- mat_err(primal_prev, primal_curr, n_list)
      if(err < err.tol) {break}
    } else {
      err <- NA
    }

    if(log_crit == "always") {
      reml_score <- reml_crit(M1, S1, I1, M2, S2, I2, M3, S3, I3)
    } else {
      reml_score <- NA
    }

    if(verbose == TRUE) {
      conv_df <- conv_df %>%
        tibble::add_row(
          iter       = i,
          est_diff   = err,
          reml_score = reml_score
        )
    }

    primal_prev <- primal_curr
  }


  out <- list(
    S1 = primal_curr[[1]],
    S2 = (primal_curr[[2]] - primal_curr[[1]])/I1,
    S3 = (primal_curr[[3]] - primal_curr[[2]])/(I1*I2)
  )

  if(log_crit %in% c("always", "output")) {
    out$reml_crit <- reml_crit(M1, out$S1, I1, M2, out$S2, I2, M3, out$S3, I3)
  }

  if(verbose == TRUE) {
    out$conv_details <- conv_df
  }

  return(out)
}

#' @rdname stepreml
#' @export
#' @param df Data frame for 2-way design
#' @param ... Other arguments.
stepreml_2way_df <- df_cov_estimate(stepreml_2way_mat)

mat_err <- function(Sig_curr_list, Sig_prev_list, n_list) {
  runsum <- 0

  for(i in length(n_list)) {
    D <- Sig_curr_list[[i]] - Sig_prev_list[[i]]
    runsum <- runsum + n_list[[i]] * sum(D^2)
  }

  return(sqrt(runsum))
}
