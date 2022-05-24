#' Compute 1- or 2-way REML estimates with a stepwise method
#'
#' Given either a tuple of sum-of-squares matrices or a dataframe, computes
#' stepwise REML estimates for the covariances
#' @name stepreml
#' @seealso \code{\link{stepreml_2way}} for generic version.
NULL
#> NULL

#' @rdname stepreml
#' @param A1,A2 Sum-of-squares matrices
#' @param n1,n2 number of lines
#'
#' @export
stepreml_1way <- function(A1, n1, A2, n2) {

  # sample covariances
  S1 <- A1/n1
  S2 <- A2/n2

  # S1 = U' U
  U <- chol(S1)
  U_inv <- solve(U)

  # inv(U') S2 inv(U) = Q D Q'
  eigs <- eigen(t(U_inv) %*% S2 %*% U_inv, symmetric = TRUE)
  d <- eigs$values
  Q <- eigs$vectors

  # B S1 B' = I
  # B S2 B' = D
  B <- t(Q) %*% t(U_inv)
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
#' @export
stepreml_2way_mat <- function(M1, I1, M2, I2, M3, I3, max_iter=50, err.tol=1e-6, verbose = FALSE) {

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
    iter = integer(),
    est_diff = numeric(),
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
    reml_score <- reml_crit(M1, S1, I1, M2, S2, I2, M3, S3, I3)

    # check for convergence
    if(i > 1) {
      err <- mat_err(primal_prev, primal_curr, n_list)
      if(err < err.tol) {break}
    } else {
      err <- NA
    }

    conv_df <- conv_df %>%
      tibble::add_row(
        iter = i,
        est_diff = err,
        reml_score = reml_score)

    primal_prev <- primal_curr
  }


  out <- list(
    S1 = primal_curr[[1]],
    S2 = (primal_curr[[2]] - primal_curr[[1]])/I1,
    S3 = (primal_curr[[3]] - primal_curr[[2]])/(I1*I2)
  )

  out$reml_crit <- reml_crit(M1, out$S1, I1, M2, out$S2, I2, M3, out$S3, I3)

  if(verbose == TRUE) {
    out$conv_details <- conv_df
  }

  return(out)
}

#' @rdname stepreml
#' @export
#' @param df Data frame for 2-way design
stepreml_2way_df <- df_cov_estimate(stepreml_2way_mat)

mat_err <- function(Sig_curr_list, Sig_prev_list, n_list) {
  runsum <- 0

  for(i in length(n_list)) {
    D <- Sig_curr_list[[i]] - Sig_prev_list[[i]]
    runsum <- runsum + n_list[[i]] * sum(D^2)
  }

  return(sqrt(runsum))
}
